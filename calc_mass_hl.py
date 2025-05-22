import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
from csv import writer
from scipy.optimize import curve_fit
from scipy.stats import chisquare
import hl_errors as err
import AIC_fitting as aic
import os 

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

def Func_exp(t, a,b): return a*np.exp(-b*t)
def Func_exp_exp(t, a,b, c,d): return a*np.exp(-b*t)+c*np.exp(-d*t)
def Func_con(t, a): return a
def Func_con_exp(t, a,b,c): return a*np.exp(-b*t)+c

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
    return param

def AIC_fitting(Func, x_1,x_2,y_1,y_2, x_data, y_data, y_error, data, Nt):
    AIC = []
    for i in range(x_1,x_2):
        for j in range(y_1,y_2):
            try:    
                param, pocv = curve_fit(Func, x_data[i:j], y_data[i:j], sigma=y_error[i:j]) 
                if np.all(np.isinf(pocv) == False):
                    perr = np.sqrt(np.diag(pocv))
                    n = len(y_data[i:j])
                    rrs = np.sum((y_data[i:j] - Func(x_data[i:j], *param))**2)
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
            probsum += np.exp((AICmin - AIC[i,2])/2)
            if len(AIC[0,:]) == 5:
                summass += np.exp((AICmin - AIC[i,2])/2)*AIC[i,4]
            if len(AIC[0,:]) == 11:
                summass += np.exp((AICmin - AIC[i,2])/2)*np.min([AIC[i,4], AIC[i,6]])
            if len(AIC[0,:]) == 6:
                summass += np.exp((AICmin - AIC[i,2])/2)*AIC[i,5]
                sumerr += (AIC[i,8])**2
    mass = summass/probsum
    masserr = np.sqrt(sumerr)
    for i in range(len(AIC[:,0])):
        if AIC[i,2] == AICmin:
            fit_params = AIC[i,:]
    if len(AIC[0,:]) == 11:
        mass_err = err.Eff_mass_curve_fit_error(Func, int(fit_params[0]), int(fit_params[1]), Effective_mass, y_error, data, len(data[:,0]), Nt)
    else:
        mass_err = err.Corr_curve_fit_error(Func, int(fit_params[0]), int(fit_params[1]), y_error, data, len(data[:,0]), Nt)
    return mass, masserr, fit_params

def autocorrelation_analysis(data, corr, N_max=200, ntau=7):
    # autocorrelation time from autocorrelation.py
    # autocovariance function
    time_series = data[:,ntau]
    nts = len(time_series)
    time_series_centered = time_series - np.average(time_series)
    autocov = np.empty(N_max)
    for j in range(N_max):
        autocov[j] = np.dot(time_series_centered[:nts - j], time_series_centered[j:])
    autocov /= nts

    # acf
    acf = autocov / autocov[0]
    tau_int_v = np.zeros(N_max)
    for j in range(N_max):
        tau_int_v[j] = 0.5 + np.sum(acf[1:j+1])
    j_max = 0
    while j_max < len(tau_int_v) - 1 and j_max < 10*tau_int_v[j_max]:
        j_max += 1
    tau_int = np.sqrt(2*tau_int_v[j_max])
    return tau_int

#Define the lattice sizes 
lattice_sizes = [
    (16, 128),
    (32, 48),
    (32, 32),
    (24, 40),
    (24, 36),
    (24, 28),
    (24, 24),
]

parts=['ps','v']
ss='s0'
lq='strange'
#tau_int = 1
corr_dir = f'/home/rachel/PhD/Bmeson/hl_corrs/{lq}'
file_dir = f'/home/rachel/PhD/Bmeson/data_files/{lq}'
plot_dir = f'/home/rachel/PhD/Bmeson/plots/{lq}/'

mass_results = {}
for Ns, Nt in lattice_sizes:
    for part in parts:
        Nt_2 = Nt // 2
        t = np.arange(Nt)
        corr_path = f"{corr_dir}/hl_{part}_{Ns}x{Nt}_s0.csv"
        try:
            corr, data = Get_corr(corr_path)
        except Exception as e:
            print(f"Could not read {corr_path}: {e}")
            continue
        effmass = Effective_mass(np.abs(corr))
        number_of_files = len(data[:,0])

        #Autocorrelation analysis (ntau=7)
        tau_int = autocorrelation_analysis(
            data, np.hstack(([0],corr)), N_max=200, ntau=7
        )

        #errors
        corr_err = err.Corr_error(number_of_files, data) * tau_int
        effmass_err = err.Mass_error_effective(Effective_mass, number_of_files, data) * tau_int

        # Fitting
        try:
            #AIC_corr1exp, AIC_corr1experr, AIC_corr1exp_fitparams = \
             #   aic.AIC_fitting(Func_exp, 1, int(Nt_2/2)-3, int(Nt_2/2)+3, Nt_2, t, corr, corr_err, data, Nt)
            AIC_effmass, AIC_effmasserr, AIC_effmass_fitparams = \
                AIC_fitting(Func_con_exp, 1, Nt_2//2-3, Nt_2//2+3, Nt_2, t, effmass, effmass_err, data, Nt)
        except Exception as e:
            print(f"Fit failed for {part} {Ns}x{Nt}: {e}")
            continue

        mass = AIC_effmass
        masserr = AIC_effmasserr
        print(f"Lattice {Ns}x{Nt}, part={part}: mass={mass}, masserr={masserr}")
        if Ns not in mass_results:
            mass_results[Ns] = []
        mass_results[Ns].append([part, Nt, mass, masserr])
        # Save fits to files with lattice size and part in the filename
        #with open(f"{file_dir}/corr_fits_{part}_{ss}_half_{Ns}x{Nt}.csv", 'a') as file:
        #    write = writer(file)
        #    write.writerow(['Nt','m','m error'])
        #    write.writerow([Nt, AIC_corr1exp, AIC_corr1experr[1]])
for Ns in mass_results:
    mass_path = f"{file_dir}/mass_{ss}_half_Ns{Ns}.csv"
    with open(mass_path, 'w', newline='') as file:
        write = writer(file)
        write.writerow(['part','Nt','m','m error'])
        #write.writerow(mass_results[Ns])
        for row in mass_results[Ns]:
            write.writerow(row)


        #with open(f"{file_dir}/massave_{part}_{ss}_half_{Ns}x{Nt}.csv", 'a') as file:
        #    write = writer(file)
        #    write.writerow(['Nt','m','m error'])
        #    write.writerow([Nt, mass, masserr])
        #with open(f"{file_dir}/autocorr_tauint_{part}_{Ns}x{Nt}.csv", 'a') as file:
        #    write = writer(file)
        #    write.writerow(['ntau','tau_int'])
        #    write.writerow([16, tau_int])

plt.show()
