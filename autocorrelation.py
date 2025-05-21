import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
from csv import writer
from scipy.optimize import curve_fit
import hl_errors as err 


def Get_corr(filename):
    """read in corrs from file and average"""
    with open(filename) as file:
        data = pd.read_csv(file,  header=None)
        data = data.to_numpy()
        file.close()

    corr = np.mean(data, axis=0)

    return corr, data

def Func_exp(x, a, b):
    return a*np.exp(-x/b)

Ns=16
Nt=128
Nt_2 = 64
part='ps'
N_max = 200
ntau = 16
t = np.arange(Nt)
corr_dir = f'/home/rachel/PhD/Bmeson/hl_corrs'
corr, data = Get_corr(corr_dir+f'/hl_{part}_{Ns}x{Nt}_new.csv')

#print('len(data[:,0])', data[:,0], 'data.shape', data.shape)
auto_corr = np.correlate(data[:,ntau]-corr[ntau], data[:,ntau]-corr[ntau], 'full')
#print(len(auto_corr))
N = len(auto_corr)
N_2 = N/2
N1 = int(N_2 - (1/2))
int_autotime = 0
for i in range(300):
    if auto_corr[i] > 0:
        int_autotime += auto_corr[i]

norm_int_autotime = int_autotime/auto_corr[0]
#print(norm_int_autotime, auto_corr[0])

#autocovariance function 
time_series = data[:,ntau]
nts = len(time_series)
time_series_centered = time_series - np.average(time_series)
autocov = np.empty(N_max)

for j in range(N_max):
    autocov[j] = np.dot(time_series_centered[:nts - j], time_series_centered[j:])
autocov /= nts

fig = plt.figure(figsize=(10, 6))
plt.plot(autocov)
plt.axhline(0, color="gray", linewidth=1)
plt.xlabel("lag time $j$")
plt.ylabel(r"$\hat{K}^{XX}_j$")

#exponential fit
x = np.arange(N_max)
j_log = np.logspace(0, 3, 100)
popt, pcov = curve_fit(Func_exp, x[2:200], auto_corr[2:200])
#print(popt[1])
plt.figure()
plt.plot(x, autocov)
#plt.plot(j, auto_corr[N1:], 'g*')
plt.plot(j_log, Func_exp(j_log, popt[0], popt[1]), 'y--')
plt.xscale('log')

#acf
acf = autocov / autocov[0]
tau_int_v = np.zeros(N_max)
for j in range(N_max):
    tau_int_v[j] = 0.5 + np.sum(acf[1:j+1])
j_max = 0 
while j_max < 10*tau_int_v[j_max]:
    j_max += 1
tau_int = np.sqrt(2*tau_int_v[j_max])
print(tau_int)
plt.figure()
plt.plot(x[1:], tau_int_v[1:])
plt.xscale("log")
plt.xlabel(r"sum length $j_\mathrm{max}$")
plt.ylabel(r"$\hat{\tau}_{X, \mathrm{int}}$")
plt.legend()


plt.show()
