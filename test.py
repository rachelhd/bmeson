import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
from csv import writer
from scipy.optimize import curve_fit
from scipy.stats import chisquare
import hl_errors as err
import AIC_fitting as aic
###################################################
#average 
#eff mass
#fit corr - exp & exp+exp
#fit eff mass - constant(lowtemp) & constant+exp
####################################################


def Get_corr(filename):
    """read in corrs from file and average"""
    with open(filename) as file:
        data = pd.read_csv(file,  header=None)
        data = data.to_numpy()
        file.close()

    corr = np.mean(data, axis=0)

    return corr[1:], data[:,1:]

def Effective_mass(y):
    m = np.zeros(len(y))
    for i in range(len(y)-1):
        m[i] = np.log(y[i]/y[i+1])
    
    return m

def Func_exp(t, a,b):
    return a*np.exp(-b*t)

def Func_con_exp(t, a,b,c):
    return a*np.exp(-b*t)+c

def AIC_fitting(Func, x_1,x_2,y_1,y_2, x_data, y_data, y_error, data, Nt):

    AIC = []
    for i in range(x_1,x_2):
        for j in range(y_1,y_2):
            try:
                covmat = np.cov(data[i:j,i:j])
                param, pocv = curve_fit(Func, x_data[i:j], y_data[i:j], sigma=covmat)
                if np.all(np.isinf(pocv) == False):
                    perr = np.sqrt(np.diag(pocv))
                    res1 = 0
                    y1 = y_data[i:j]
                    yer = y_error[i:j]
                    n = len(y1)
                    rrs = 0

                    for k in range(n):
                        #rrs += ((Func(x_data[i:j], *param)[k]-y1[k])**2)/yer[k]
                        rrs += (y1[k] - Func(x_data[i:j], *param)[k])**2
                    AICval = 2*len(param) + n*np.log(rrs/n)
                    AIC.append([i,j,AICval, *param, *perr])
            except:
                continue

    AIC = np.array(AIC)
    AICmin = np.min(AIC[:,2])
    for i in range(len(AIC[:,0])):
        if AIC[i,2] == AICmin:
            minrow = i
    summass = 0 
    probsum = 0
    sumerr = 0
    for i in range(len(AIC[:,0])):
        if 0 <= AIC[i,2]-AICmin < 5:
            #print(AICmin-AIC[i,2])
            probsum += np.exp((AICmin - AIC[i,2])/2)
            if len(AIC[0,:]) == 7: #one exponential fit
                summass += np.exp((AICmin - AIC[i,2])/2)*AIC[i,4]
                sumerr += (AIC[i,6])**2
            if len(AIC[0,:]) == 11: #two expontential fit
                summass += np.exp((AICmin - AIC[i,2])/2)*np.min([AIC[i,4], AIC[i,6]])
                sumerr += np.min(AIC[i,8],AIC[i,10])**2
            if len(AIC[0,:]) == 9: #exponential + c
                summass += np.exp((AICmin - AIC[i,2])/2)*AIC[i,5]
                sumerr += (AIC[i,8])**2
    mass = summass/probsum
    masserr = np.sqrt(sumerr)
    for i in range(len(AIC[:,0])):
      if AIC[i,2] == AICmin:
            fit_params = AIC[i,:]
    
    if len(AIC[0,:]) == 11:
        mass_err = err.Eff_mass_curve_fit_error(Func, int(fit_params[0]), int(fit_params[1]), 
                                                Effective_mass, y_error, data, 
                                                len(data[:,0]), Nt)

    else:
        mass_err = err.Corr_curve_fit_error(Func, int(fit_params[0]), int(fit_params[1]), 
                                            y_error, data, len(data[:,0]), Nt)

    return mass, masserr, fit_params


Ns=16
Nt=128
Nt_2 = 64
part='v'
ss='s0'
lq='strange'
tau_int = 1.14
t = np.arange(Nt)
corr_dir = f'/home/rachel/PhD/Bmeson/hl_corrs/{lq}'
file_dir = f'/home/rachel/PhD/Bmeson/data_files/{lq}'
plot_dir = f'/home/rachel/PhD/Bmeson/plots/{lq}'
corr, data = Get_corr(corr_dir+f'/hl_{part}_{Ns}x128_new.csv')
effmass = Effective_mass(np.abs(corr))
#corrcov = np.cov(data)
number_of_files = len(data[:,0])

corr_err = err.Corr_error(number_of_files, data)*tau_int
effmass_err = err.Mass_error_effective(Effective_mass, number_of_files, data)*tau_int



corrcov=np.cov(data[1:Nt_2,1:Nt_2], rowvar=False)
popt, pcov = curve_fit(Func_con_exp, t[1:Nt_2], effmass[1:Nt_2], sigma=corrcov)
perr=np.sqrt(np.diag(pcov))
print('cov mat', popt, perr)
popt, pcov = curve_fit(Func_con_exp, t[1:Nt_2], effmass[1:Nt_2], sigma=effmass_err[1:Nt_2])
perr=np.sqrt(np.diag(pcov))
print('normal err', popt, perr)


mass, masserr, fits = AIC_fitting(Func_con_exp, 1,int(Nt_2)-10,Nt_2//2+10,Nt_2, \
                                t, effmass, effmass_err, data, Nt)



print('AIC',mass, masserr)


plt.figure()
plt.errorbar(5.63*1000/Nt, 5.63*1000*mass+4126, yerr=5.63*1000*masserr, fmt='g*')
plt.axhline(y=5366, color='b', linestyle='--', label=r'$B^0_S(exp)$')
plt.axhline(y=5415, color='r', linestyle='dashed', label=r'$B^*_S(exp)$')

plt.show()












