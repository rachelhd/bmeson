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

def Func_exp_exp(t, a,b, c,d):
    return a*np.exp(-b*t)+c*np.exp(-d*t)

def Func_con(t, a):
    return a

def Func_con_exp(t, a,b,c):
    return a*np.exp(-b*t)+c

def Curvefit_corr(func,Nt,Ns,x_1,x_2, x,y,a, x_data,y_data,y_error, Label):
    param, pcov = curve_fit(func, x_data[x_1:x_2], y_data[x_1:x_2], sigma=y_error[x_1:x_2])
    perr = np.sqrt(np.diag(pcov))

    plt.figure()
    plt.errorbar(x_data[x:y], y_data[x:y]*a, yerr=y_error[x:y]*a, fmt='b*', label=Label)
    if len(param) == 2:
        plt.plot(x_data[x_1:x_2], func(x_data[x_1:x_2], param[0], param[1]), 'r--')
    else:
        plt.plot(x_data[x_1:x_2], func(x_data[x_1:x_2], param[0], param[1], param[2], param[3]), 'r--')
    plt.legend()
    plt.yscale('log')
    plt.xlabel(r'$\tau/a_{\tau}$')
    plt.ylabel(r'$log(G(\tau)/G(\tau+1))$')
    #plt.savefig(f'{y_data}_{Nt}x{Ns}.png')
    return param

def Curvefit_effmass(func,Nt,Ns, x_1, x_2, x, y, a, x_data, y_data, y_error, Label):
    param, pcov = curve_fit(func, x_data[x_1:x_2], y_data[x_1:x_2], sigma=y_error[x_1:x_2])
    perr = np.sqrt(np.diag(pcov))

    plt.figure()
    plt.errorbar(x_data[x:y], y_data[x:y]*a, yerr=y_error[x:y]*a, fmt='b*', label=Label)
    if len(param) == 3:
        plt.plot(x_data[x_1:x_2], func(x_data[x_1:x_2], param[0], param[1], param[2]), 'r--')
    else: 
        plt.plot(x_data[x_1:x_2], func(x_data[x_1:x_2], param[0]), 'r--')
    plt.legend()
    plt.xlabel(r'$\tau/a_{\tau}$')
    plt.ylabel(r'$log(G(\tau)/G(\tau+1))$')
    plt.yscale('log')
    #plt.savefig(f'{y_data}_{Nt}x{Ns}.png')

    return param

def AIC_fitting(Func, x_1,x_2,y_1,y_2, x_data, y_data, y_error, data, Nt):

    AIC = []
    for i in range(x_1,x_2):
        for j in range(y_1,y_2):
            try:    
                param, pocv = curve_fit(Func, x_data[i:j], y_data[i:j], sigma=y_error[i:j]) 
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
                    AIC.append([i,j,AICval, *param])
            except:
           #     print("no", i,j)
                continue

    AIC = np.array(AIC)
#    print("len(AIC)", len(AIC[0,:]))
    AICmin = np.min(AIC[:,2])
    for i in range(len(AIC[:,0])):
        if AIC[i,2] == AICmin:
            minrow = i
    summass = 0 
    probsum = 0
    for i in range(len(AIC[:,0])):
        #if AICmin-AIC[i,2] == 0:
         #   probsum += 1
         #   if len(AIC[0,:]) == 5: #one exponential fit
         #       summass += AIC[i,4]
         #   if len(AIC[0,:]) == 11: #two expontential fit
         #       summass += np.min([AIC[i,4], AIC[i,6]])
         #   if len(AIC[0,:]) == 6: #exponential + c
         #       summass += AIC[i,5]

        if 0 <= AIC[i,2]-AICmin < 5:
            #print(AICmin-AIC[i,2])
            probsum += np.exp((AICmin - AIC[i,2])/2)
            if len(AIC[0,:]) == 5: #one exponential fit
                summass += np.exp((AICmin - AIC[i,2])/2)*AIC[i,4]
            if len(AIC[0,:]) == 11: #two expontential fit
                summass += np.exp((AICmin - AIC[i,2])/2)*np.min([AIC[i,4], AIC[i,6]])
            if len(AIC[0,:]) == 6: #exponential + c
                summass += np.exp((AICmin - AIC[i,2])/2)*AIC[i,5]

    mass = summass/probsum
    
    for i in range(len(AIC[:,0])):
 #       if 0 <= AIC[i,2]-AICmin < 5:
            
#            plt.figure(figsize=(16,9))
#            plt.errorbar(t, effmass, yerr=effmass_err[1:], fmt='b^')
#            plt.plot(t[int(AIC[i,0]):int(AIC[i,1])], Func(t[int(AIC[i,0]):int(AIC[i,1])], *AIC[i,3:]), 'r--')
#            plt.ylabel(r'$M(\tau)=log(\frac{G(\tau)}{G(\tau+1)})$')
#            plt.xlabel(r'$frac{\tau}{a_{\tau}}$')
#            plt.title(f'{Ns}x{Nt} effective mass')
#            plt.yscale('log')
            
 #           print(i)
        if AIC[i,2] == AICmin:
            fit_params = AIC[i,:]
    
    if len(AIC[0,:]) == 11:
        mass_err = err.Eff_mass_curve_fit_error(Func, int(fit_params[0]), int(fit_params[1]), 
                                                Effective_mass, y_error, data, 
                                                len(data[:,0]), Nt)

    else:
        mass_err = err.Corr_curve_fit_error(Func, int(fit_params[0]), int(fit_params[1]), 
                                            y_error, data, len(data[:,0]), Nt)

    return mass, mass_err*tau_int, fit_params

############################################################################################################


Ns=16
Nt=128
Nt_2 = 64 
part='v'
ss='s0'
lq='light'
tau_int = 1
t = np.arange(Nt)
corr_dir = f'/home/rachel/PhD/Bmeson/hl_corrs/{lq}'
file_dir = f'/home/rachel/PhD/Bmeson/data_files/{lq}'
plot_dir = f'/home/rachel/PhD/Bmeson/plots/{lq}'
corr, data = Get_corr(corr_dir+f'/hl_{part}_{Ns}x{Nt}_s0.csv')
#print(corr)
#corr = corr[:Nt]
#data = data[:,:Nt]
effmass = Effective_mass(np.abs(corr))
corrfuncs=[]
#corrfuncs.append(corr[1:])
#corrfuncs.append(effmass)

number_of_files = len(data[:,0])

#plt.figure()
#plt.plot(t, corr[1:], 'b.')
plt.figure(figsize=(16,9))
corr_err = err.Corr_error(number_of_files, data)*tau_int
effmass_err = err.Mass_error_effective(Effective_mass, number_of_files, data)*tau_int
nt = np.array([Nt])
nt_corr = np.hstack((nt, corr))
corrfuncs.append(nt_corr)
nt_correrr = np.hstack((nt, corr_err))
corrfuncs.append(nt_correrr)
nt_efmass = np.hstack((nt, effmass))
corrfuncs.append(nt_efmass)
nt_efmasserr = np.hstack((nt, effmass_err))
corrfuncs.append(nt_efmasserr)

#plot corr & mass
plt.errorbar(t, corr, yerr=corr_err, fmt='b^')
plt.ylabel(r'$G(\tau)$')
plt.xlabel(r'$\frac{\tau}{a_{\tau}}$')
#plt.title(f'{Ns}x{Nt} correlator')
plt.yscale('log')
#plt.savefig(plot_dir+f'/corr_{part}_{Ns}x{Nt}_{ss}.pdf')

plt.figure(figsize=(16,9))
plt.errorbar(t, effmass, yerr=effmass_err, fmt='b^')
plt.ylabel(r'$M(\tau)=log(\frac{G(\tau)}{G(\tau+1)})$')
plt.xlabel(r'$frac{\tau}{a_{\tau}}$')
plt.title(f'{Ns}x{Nt} effective mass')
plt.yscale('log')
#plt.savefig(plot_dir+f'/effmass_{part}_{Ns}x{Nt}_{ss}.pdf')

AIC_corr1exp, AIC_corr1experr, AIC_corr1exp_fitparams = \
                AIC_fitting(Func_exp, 1, int(Nt_2/2)-3, int(Nt_2/2)+3, Nt_2, t, corr, corr_err, data, Nt)
#AIC_corr2exp, AIC_corr2experr, AIC_corr2exp_fitparams = \
 #               AIC_fitting(Func_exp_exp, 1, int(Nt_2/2), int(Nt_2/2), Nt_2+10, t, corr, corr_err, data, Nt)
AIC_effmass, AIC_effmasserr, AIC_effmass_fitparams = \
                AIC_fitting(Func_con_exp, 1, int(Nt_2/2)-3, int(Nt_2/2)+3,Nt_2, t, effmass, effmass_err, data, Nt)
print("corr1", AIC_corr1exp, AIC_corr1experr, AIC_corr1exp_fitparams, \
      "effmass", AIC_effmass, AIC_effmasserr, AIC_effmass_fitparams)
   #   "corr2", AIC_corr2exp, AIC_corr2experr, AIC_corr2exp_fitparams)

##average both mass values
aicmin = np.min([AIC_corr1exp_fitparams[2], AIC_effmass_fitparams[2]])
probc1 = np.exp((aicmin-AIC_corr1exp_fitparams[2])/2)
probef = np.exp((aicmin-AIC_effmass_fitparams[2])/2)
mass = (AIC_corr1exp*probc1  + AIC_effmass*probef)/(probc1+probef)
masserr = np.sqrt(AIC_corr1experr[1]**2 +  AIC_effmasserr[2]**2)
print(mass, masserr, AIC_effmass_fitparams[4])
"""
plt.figure(figsize=(16,9))
#plot corr & mass
plt.errorbar(t, corr[1:], yerr=corr_err[1:], fmt='b^')
plt.plot(t[int(AIC_corr1exp_fitparams[0]):int(AIC_corr1exp_fitparams[1])], \
        Func_exp(t[int(AIC_corr1exp_fitparams[0]):int(AIC_corr1exp_fitparams[1])], *AIC_corr1exp_fitparams[3:]), 'r--')
plt.ylabel(r'$G(\tau)$')
plt.xlabel(r'$\frac{\tau}{a_{\tau}}$')
#plt.title(f'{Ns}x{Nt} correlator')
plt.yscale('log')
#plt.savefig(plot_dir+f'/corr_{part}_{Ns}x{Nt}_{ss}.pdf')
"""
plt.figure(figsize=(16,9))
plt.errorbar(t, effmass, yerr=effmass_err, fmt='b^')
plt.plot(t[int(AIC_effmass_fitparams[0]):int(AIC_effmass_fitparams[1])], \
        Func_con_exp(t[int(AIC_effmass_fitparams[0]):int(AIC_effmass_fitparams[1])], *AIC_effmass_fitparams[3:]), 'r--')
plt.ylabel(r'$M(\tau)=log(\frac{G(\tau)}{G(\tau+1)})$')
plt.xlabel(r'$frac{\tau}{a_{\tau}}$')
plt.title(f'{Ns}x{Nt} effective mass')
plt.yscale('log')
#plt.savefig(plot_dir+f'/effmassfit_{part}_{Ns}x{Nt}_{ss}.pdf')

plt.figure(figsize=(16,9))
plt.errorbar(t, corr, yerr=corr_err, fmt='b^')
plt.plot(t[int(AIC_corr1exp_fitparams[0]):int(AIC_corr1exp_fitparams[1])], \
        Func_exp(t[int(AIC_corr1exp_fitparams[0]):int(AIC_corr1exp_fitparams[1])], *AIC_corr1exp_fitparams[3:]), 'r--')
plt.ylabel(r'$G(\tau)$')
plt.xlabel(r'$frac{\tau}{a_{\tau}}$')
plt.title(f'{Ns}x{Nt} corr')
plt.yscale('log')
#plt.savefig(plot_dir+f'corrfit_{part}_{Ns}x{Nt}_{ss}.pdf')

#save fits to files
corr_ef = []
corr_ef.append([Nt, AIC_corr1exp, AIC_corr1experr[1]])
efmass_cef = []
efmass_cef.append([Nt, AIC_effmass, AIC_effmasserr[2]])
mass_f = []
mass_f.append([Nt, mass, masserr])

with open(file_dir+f'/corr_funcs_{part}_{Ns}_{ss}.csv', 'a') as file:
    write = writer(file)
    #write.writerow(corrfuncs[0])
    #write.writerow(corrfuncs[1])
    file.close()

with open(file_dir+f'/effmass_funcs_{part}_{ss}.csv', 'a') as file:
    write = writer(file)
    #write.writerow(corrfuncs[2])
    #write.writerow(corrfuncs[3])
    file.close()

with open(file_dir+f'/corr_fits_{part}_T0_{ss}_half.csv', 'a') as file:
    write = writer(file)
    write.writerow(['Nt','m','m error'])
    write.writerow(corr_ef[0])
    file.close()

with open(file_dir+f'/effmass_fits_{part}_T0_{ss}_half.csv', 'a') as file:
    write = writer(file)
    write.writerow(['Nt','m','m error'])
    write.writerow(efmass_cef[0])
    file.close()

with open(file_dir+f'/massave_{part}_T0_{ss}_half.csv', 'a') as file:
    write = writer(file)
    write.writerow(['Nt','m','m error'])
    write.writerow(mass_f[0])
    file.close()



plt.show()
