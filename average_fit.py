import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
from csv import writer 

def Func_1exp(t,a,b):
    return a*np.exp(-b*t)
def Func_2exp(t,a,b,c,d):
    return a*np.exp(-b*t)+c*np.exp(-d*t)
def Func(t,a,b,c):
    return a*np.exp(-b*t)+c

file_dir=f'/home/rachel/PhD/Bmeson/data_files'
Ns=16
Nt=128
part='v'
x_1,y_1 = 15, Nt-60
x_2,y_2 = 5, Nt-53
x_3,y_3 = 3, 69

with open(file_dir+f'/corr_expfits_{part}_{Ns}x{Nt}.csv', 'r') as file1, open(file_dir+f'/corr_effmass_funcs_{part}_{Ns}x{Nt}.csv', 'r') as file2:
    param1 = pd.read_csv(file1, header=None)
    param1 = param1.to_numpy()[-1,1:]
    param1 = np.array(param1, dtype=float)
    y1 = pd.read_csv(file2, header=None)
    y1 = y1.to_numpy()[0,:]
    file1.close()
    file2.close()

with open(file_dir+f'/corr_2expfits_{part}_{Ns}x{Nt}.csv', 'r') as file1, open(file_dir+f'/corr_effmass_funcs_{part}_{Ns}x{Nt}.csv', 'r') as file2:
    param2 = pd.read_csv(file1, header=None)
    param2 = param2.to_numpy()[-1,1:]
    param2 = np.array(param2, dtype=float)
    y2 = pd.read_csv(file2, header=None)
    y2 = y2.to_numpy()[0,:]
    file1.close()
    file2.close()

with open(file_dir+f'/effmass_contexpfits_{part}_{Ns}x{Nt}.csv', 'r') as file1, open(file_dir+f'/corr_effmass_funcs_{part}_{Ns}x{Nt}.csv', 'r') as file2:
    param3 = pd.read_csv(file1, header=None)
    param3 = param3.to_numpy()[-1,1:]
    param3 = np.array(param3, dtype=float)
    y3 = pd.read_csv(file2, header=None)
    y3 = y3.to_numpy()[1,:]
    file1.close()
    file2.close()
print('corr2 fit parameters',param2)

x_data=np.arange(Nt)
res1 = 0
y1 = y1[x_1:y_1]
for i in range(len(y1)):
    res1 += ((Func_1exp(x_data[x_1:y_1], param1[0], param1[1])[i]-y1[i])**2)
AIC1 = np.exp(4-2*np.log(res1/len(y1)))

res2 = 0
y2 = y2[x_2:y_2]
for i in range(len(y2)):
    res2 += ((Func_2exp(x_data[x_2:y_2], param2[0], param2[1], param2[2], param2[3])[i]-y2[i])**2)
AIC2 = np.exp(8-2*np.log(res2/len(y2)))

res3 = 0
y3 = y3[x_3:y_3]
for i in range(len(y3)):
    res3 += ((Func(x_data[x_3:y_3], param3[0], param3[1], param3[2])[i]-y3[i])**2)
AIC3 = np.exp(6-2*np.log(res3/len(y3)))

AIC = np.array([AIC1,AIC2,AIC3])
AICmin = np.min(AIC)
print('AIC values',AIC1, AIC2, AIC3)

prob1 = np.exp((AICmin-AIC1)/2)
prob2 = np.exp((AICmin-AIC2)/2)
prob3 = np.exp((AICmin-AIC3)/2)
print('probabilities',prob1, prob2, prob3)

#mass = (prob1*param1[1]+prob2*param2[3]+prob3*param3[2])/(prob1+prob2+prob3)
#mass_err1 = np.sqrt((prob1*param1[3]**2+prob2*param2[7]**2+prob3*param3[5]**2)/(prob1+prob2+prob3))
#mass_err2 = np.sqrt(((mass-param1[1])**2 + (mass-param2[3])**2 + (mass-param3[2])**2)/(3))

mass = (prob1*param1[1]+prob2*param2[1]+prob3*param3[2])/(prob1+prob2+prob3)
mass_err1 = np.sqrt((prob1*param1[3]**2+prob2*param2[5]**2+prob3*param3[5]**2)/(prob1+prob2+prob3))
mass_err2 = np.sqrt(((mass-param1[1])**2 + (mass-param2[1])**2 + (mass-param3[2])**2)/(3))

mass_err=np.sqrt(mass_err1**2+mass_err2**2)
print('mass',mass, mass_err)

mass_save = []
mass_save.append([Ns, Nt, mass, mass_err])
with open(file_dir+f'/mass_{part}.csv', 'a') as file:
    write = writer(file)
    write.writerow(['Ns','Nt','m', 'm error'])
    write.writerow(mass_save[0])
    file.close()

