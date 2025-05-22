import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
from csv import writer
from scipy.optimize import curve_fit
import hl_errors as err
import AIC_fitting as aic



def Get_data(filename):
    """read in corrs from file and average"""
    with open(filename) as file:
        data = pd.read_csv(file,  header=None)
        data = data.to_numpy()
        file.close()

    return data




file_dir = "/home/rachel/PhD/Bmeson/data_files"

pscorr16 = Get_data(file_dir+f"/corr_funcs_ps_16_s0.csv")
psefm = Get_data(file_dir+f"/effmass_funcs_ps_s0.csv")
vcorr16 = Get_data(file_dir+f"/corr_funcs_v_16_s0.csv")
vefm = Get_data(file_dir+f"/effmass_funcs_v_s0.csv")
pscorr24 = Get_data(file_dir+f"/corr_funcs_ps_24_s0.csv")
vcorr24 = Get_data(file_dir+f"/corr_funcs_v_24_s0.csv")
pscorr32 = Get_data(file_dir+f"/corr_funcs_ps_32_s0.csv")
vcorr32 = Get_data(file_dir+f"/corr_funcs_v_32_s0.csv")
"""
#corr plot
plt.figure(figsize=(16,9))
plt.errorbar(np.arange(pscorr16[0,0]), pscorr16[0,1:], yerr=pscorr16[1,1:], fmt='b.', label=r'44 MeV')
plt.errorbar(np.arange(pscorr32[0,0]), pscorr32[0,1:49], yerr=pscorr32[1,1:49], fmt='r*', label=r'117 MeV')
plt.errorbar(np.arange(pscorr24[0,0]), pscorr24[0,1:41], yerr=pscorr24[1,1:41], fmt='y<', label=r'140 MeV')
plt.errorbar(np.arange(pscorr32[2,0]), pscorr32[2,1:33], yerr=pscorr32[3,1:33], fmt='g^', label=r'175 MeV')
plt.errorbar(np.arange(pscorr24[2,0]), pscorr24[2,1:25], yerr=pscorr24[3,1:25], fmt='mx', label=r'235 MeV')
plt.legend(fontsize=15)
plt.xlabel(r'$\tau/a_\tau$', fontsize=20)
plt.ylabel(r'$G(\tau)$', fontsize=20)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.yscale('log')
plt.xlim(-5,64)
#plt.savefig(file_dir+f"../plots/.pdf")

plt.figure(figsize=(16,9))
plt.errorbar(np.arange(vcorr16[0,0]), vcorr16[0,1:], yerr=vcorr16[0,1:], fmt='b.', label=r'44 MeV')
plt.errorbar(np.arange(vcorr32[0,0]), vcorr32[0,1:49], yerr=vcorr32[1,1:49], fmt='r*', label=r'117 MeV')
plt.errorbar(np.arange(vcorr24[0,0]), vcorr24[0,1:41], yerr=vcorr24[1,1:41], fmt='g^', label=r'140 MeV')
plt.errorbar(np.arange(vcorr32[2,0]), vcorr32[2,1:33], yerr=vcorr32[3,1:33], fmt='y<', label=r'175 MeV')
plt.errorbar(np.arange(vcorr24[2,0]), vcorr24[2,1:25], yerr=vcorr24[3,1:25], fmt='mx', label=r'235 MeV')
plt.legend(fontsize=15)
plt.xlabel(r'$\tau/a_\tau$',fontsize=20)
plt.ylabel(r'$G(\tau)$',fontsize=20)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.yscale('log')
plt.xlim(-5,64)
#plt.savefig(file_dir+f"../plots/.pdf")


#effmass plot
plt.figure(figsize=(16,9))
plt.errorbar(np.arange(psefm[0,0]), psefm[0,1:], yerr=psefm[1,1:], fmt='b.', label=r'44 MeV')
plt.errorbar(np.arange(psefm[2,0])[:-2], psefm[2,1:47], yerr=psefm[3,1:47], fmt='r*', label=r'117 MeV')
plt.errorbar(np.arange(psefm[6,0])[:-2], psefm[6,1:39], yerr=psefm[7,1:39], fmt='y<', label=r'140 MeV')
plt.errorbar(np.arange(psefm[4,0])[:-2], psefm[4,1:31], yerr=psefm[5,1:31], fmt='g^', label=r'175 MeV')
plt.errorbar(np.arange(psefm[8,0])[:-2], psefm[8,1:23], yerr=psefm[9,1:23], fmt='mx', label=r'235 MeV')
plt.legend(fontsize=15)
plt.xlabel(r'$\tau/a_\tau$', fontsize=20)
plt.ylabel(r'$M(\tau)$', fontsize=20)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
#plt.yscale('log')
plt.xlim(-5,64)
#plt.savefig(file_dir+f"../plots/.pdf")

plt.figure(figsize=(16,9))
plt.errorbar(np.arange(vefm[0,0]),vefm[0,1:], yerr=vefm[1,1:], fmt='b.', label=r'44 MeV')
plt.errorbar(np.arange(vefm[2,0])[:-2], vefm[2,1:47], yerr=vefm[3,1:47], fmt='r*', label=r'117 MeV')
plt.errorbar(np.arange(vefm[6,0])[:-2], vefm[6,1:39], yerr=vefm[7,1:39], fmt='y<', label=r'140 MeV')
plt.errorbar(np.arange(vefm[4,0])[:-2], vefm[4,1:31], yerr=vefm[5,1:31], fmt='g^', label=r'175 MeV')
plt.errorbar(np.arange(vefm[8,0])[:-2], vefm[8,1:23], yerr=vefm[9,1:23], fmt='mx', label=r'235 MeV')
plt.legend(fontsize=15)
plt.xlabel(r'$\tau/a_\tau$', fontsize=20)
plt.ylabel(r'$M(\tau)$', fontsize=20)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.yscale('log')
plt.xlim(-5,64)
#plt.savefig(file_dir+f"../plots/.pdf")
"""

with open(file_dir+"/massave_ps_16_s0.csv") as file:
    psm16 = pd.read_csv(file)
    psm16 = psm16.to_numpy()
    file.close()

with open(file_dir+"/massave_v_16_s0.csv") as file:
    vm16 = pd.read_csv(file)
    vm16 = vm16.to_numpy()
    file.close()

with open(file_dir+"/massave_ps_24_s0.csv") as file:
    psm24 = pd.read_csv(file)
    psm24 = psm24.to_numpy()
    file.close()

with open(file_dir+"/massave_v_24_s0.csv") as file:
    vm24 = pd.read_csv(file)
    vm24 = vm24.to_numpy()
    file.close()

with open(file_dir+"/massave_ps_32_s0.csv") as file:
    psm32 = pd.read_csv(file)
    psm32 = psm32.to_numpy()
    file.close()

with open(file_dir+"/massave_v_32_s0.csv") as file:
    vm32 = pd.read_csv(file)
    vm32 = vm32.to_numpy()
    file.close()

with open(file_dir+"/massave_ps_T0_s0.csv") as file:
    psmT0 = pd.read_csv(file)
    psmT0 = psmT0.to_numpy()
    file.close()

with open(file_dir+"/massave_v_T0_s0.csv") as file:
    vmT0 = pd.read_csv(file)
    vmT0 = vmT0.to_numpy()
    file.close()

#plot mass with Temp
plt.figure(figsize=(16,9))
plt.errorbar(5.63*1000/psm16[:,0], 5.63*1000*psm16[:,1]+4126, yerr=5.63*1000*psm16[:,2], fmt='bX',markersize=12)
plt.errorbar(5.63*1000/vm16[:,0], 5.63*1000*vm16[:,1]+4126, yerr=5.63*1000*vm16[:,2], fmt='r*',markersize=12)
plt.errorbar(5.63*1000/psm24[:,0], 5.63*1000*psm24[:,1]+4126, yerr=5.63*1000*psm24[:,2], fmt='bX',markersize=12, label=r'$B_s^0 Ns=24$')
plt.errorbar(5.63*1000/vm24[:,0], 5.63*1000*vm24[:,1]+4126, yerr=5.63*1000*vm24[:,2], fmt='r*',markersize=12, label=r'$B_s^* Ns=24$')
plt.errorbar(5.63*1000/psm32[:,0], 5.63*1000*psm32[:,1]+4126, yerr=5.63*1000*psm32[:,2], fmt='b+',markersize=12, label=r'$B_s^0 Ns=32$')
plt.errorbar(5.63*1000/vm32[:,0], 5.63*1000*vm32[:,1]+4126, yerr=5.63*1000*vm32[:,2], fmt='r.',markersize=12, label=r'$B_s^* Ns=32$')
plt.errorbar(5.63*1000/psmT0[:,0], 5.63*1000*psmT0[:,1]+4126, yerr=5.63*1000*psmT0[:,2], fmt='g^', markersize=12,label=r'$B_s^0 T_0$')
plt.errorbar(5.63*1000/vmT0[:,0], 5.63*1000*vmT0[:,1]+4126, yerr=5.63*1000*vmT0[:,2], fmt='m>',markersize=12, label=r'$B_s^* T_0$')
plt.legend(fontsize=15)
plt.xlabel(r'T MeV', fontsize=20)
plt.ylabel(r'M(T) MeV', fontsize=20)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.axhline(y=5366, color='b', linestyle='--', label=r'$B^0_S(exp)$')
plt.axhline(y=5415, color='r', linestyle='dashed', label=r'$B^*_S(exp)$')
plt.ylim(4600, 6200)
#plt.xlim(25,210)

with open(file_dir+"/massave_ps_16_s0_half.csv") as file:
    psm16 = pd.read_csv(file)
    psm16 = psm16.to_numpy()
    file.close()

with open(file_dir+"/massave_v_16_s0_half.csv") as file:
    vm16 = pd.read_csv(file)
    vm16 = vm16.to_numpy()
    file.close()

with open(file_dir+"/massave_ps_24_s0_half.csv") as file:
    psm24 = pd.read_csv(file)
    psm24 = psm24.to_numpy()
    file.close()

with open(file_dir+"/massave_v_24_s0_half.csv") as file:
    vm24 = pd.read_csv(file)
    vm24 = vm24.to_numpy()
    file.close()

with open(file_dir+"/massave_ps_32_s0_half.csv") as file:
    psm32 = pd.read_csv(file)
    psm32 = psm32.to_numpy()
    file.close()

with open(file_dir+"/massave_v_32_s0_half.csv") as file:
    vm32 = pd.read_csv(file)
    vm32 = vm32.to_numpy()
    file.close()
with open(file_dir+"/massave_ps_T0_s0_half.csv") as file:
    psmT0 = pd.read_csv(file)
    psmT0 = psmT0.to_numpy()
    file.close()

with open(file_dir+"/massave_v_T0_s0_half.csv") as file:
    vmT0 = pd.read_csv(file)
    vmT0 = vmT0.to_numpy()
    file.close()

#plot mass with Temp
plt.figure(figsize=(16,9))
plt.errorbar(5.63*1000/psm16[:,0], 5.63*1000*psm16[:,1]+4126, yerr=5.63*1000*psm16[:,2], fmt='bX',markersize=12)
plt.errorbar(5.63*1000/vm16[:,0], 5.63*1000*vm16[:,1]+4126, yerr=5.63*1000*vm16[:,2], fmt='r*',markersize=12)
plt.errorbar(5.63*1000/psm24[:,0], 5.63*1000*psm24[:,1]+4126, yerr=5.63*1000*psm24[:,2], fmt='bX',markersize=12, label=r'$B_s^0 Ns=24$')
plt.errorbar(5.63*1000/vm24[:,0], 5.63*1000*vm24[:,1]+4126, yerr=5.63*1000*vm24[:,2], fmt='r*',markersize=12, label=r'$B_s^* Ns=24$')
plt.errorbar(5.63*1000/psm32[:,0], 5.63*1000*psm32[:,1]+4126, yerr=5.63*1000*psm32[:,2], fmt='b+',markersize=12, label=r'$B_s^0 Ns=32$')
plt.errorbar(5.63*1000/vm32[:,0], 5.63*1000*vm32[:,1]+4126, yerr=5.63*1000*vm32[:,2], fmt='r.',markersize=12, label=r'$B_s^* Ns=32$')
plt.errorbar(5.63*1000/psmT0[:,0], 5.63*1000*psmT0[:,1]+4126, yerr=5.63*1000*psmT0[:,2], fmt='g^', markersize=12,label=r'$B_s^0 T_0$')
plt.errorbar(5.63*1000/vmT0[:,0], 5.63*1000*vmT0[:,1]+4126, yerr=5.63*1000*vmT0[:,2], fmt='m>',markersize=12, label=r'$B_s^* T_0$')
plt.legend(fontsize=15)
plt.xlabel(r'T MeV', fontsize=20)
plt.ylabel(r'M(T) MeV', fontsize=20)
plt.xticks(fontsize=15) 
plt.yticks(fontsize=15)
plt.axhline(y=5366, color='b', linestyle='--', label=r'$B^0_S(exp)$')
plt.axhline(y=5415, color='r', linestyle='dashed', label=r'$B^*_S(exp)$')
plt.ylim(4600, 6200)
#plt.xlim(25,210)
plt.title('half')


#############eff mass & corr vs temp

with open(file_dir+"/effmass_fits_ps_16_s0.csv") as file:
    psm16 = pd.read_csv(file)
    psm16 = psm16.to_numpy()
    file.close()

with open(file_dir+"/effmass_fits_v_16_s0.csv") as file:
    vm16 = pd.read_csv(file)
    vm16 = vm16.to_numpy()
    file.close()

with open(file_dir+"/effmass_fits_ps_24_s0.csv") as file:
    psm24 = pd.read_csv(file)
    psm24 = psm24.to_numpy()
    file.close()

with open(file_dir+"/effmass_fits_v_24_s0.csv") as file:
    vm24 = pd.read_csv(file)
    vm24 = vm24.to_numpy()
    file.close()

with open(file_dir+"/effmass_fits_ps_32_s0.csv") as file:
    psm32 = pd.read_csv(file)
    psm32 = psm32.to_numpy()
    file.close()

with open(file_dir+"/effmass_fits_v_32_s0.csv") as file:
    vm32 = pd.read_csv(file)
    vm32 = vm32.to_numpy()
    file.close()

with open(file_dir+"/effmass_fits_ps_T0_s0.csv") as file:
    psmT0 = pd.read_csv(file)
    psmT0 = psmT0.to_numpy()
    file.close()

with open(file_dir+"/effmass_fits_v_T0_s0.csv") as file:
    vmT0 = pd.read_csv(file)
    vmT0 = vmT0.to_numpy()
    file.close()

#plot corr with Temp
plt.figure(figsize=(16,9))
plt.errorbar(5.63*1000/psm16[:,0], 5.63*1000*psm16[:,1]+4126, yerr=5.63*1000*psm16[:,2], fmt='bX',markersize=12)
plt.errorbar(5.63*1000/vm16[:,0], 5.63*1000*vm16[:,1]+4126, yerr=5.63*1000*vm16[:,2], fmt='r*',markersize=12)
plt.errorbar(5.63*1000/psm24[:,0], 5.63*1000*psm24[:,1]+4126, yerr=5.63*1000*psm24[:,2], fmt='bX',markersize=12, label=r'$B_s^0 Ns=24$')
plt.errorbar(5.63*1000/vm24[:,0], 5.63*1000*vm24[:,1]+4126, yerr=5.63*1000*vm24[:,2], fmt='r*',markersize=12, label=r'$B_s^* Ns=24$')
plt.errorbar(5.63*1000/psm32[:,0], 5.63*1000*psm32[:,1]+4126, yerr=5.63*1000*psm32[:,2], fmt='b+',markersize=12, label=r'$B_s^0 Ns=32$')
plt.errorbar(5.63*1000/vm32[:,0], 5.63*1000*vm32[:,1]+4126, yerr=5.63*1000*vm32[:,2], fmt='r.',markersize=12, label=r'$B_s^* Ns=32$')
plt.errorbar(5.63*1000/psmT0[:,0], 5.63*1000*psmT0[:,1]+4126, yerr=5.63*1000*psmT0[:,2], fmt='g^', markersize=12,label=r'$B_s^0 T_0$')
plt.errorbar(5.63*1000/vmT0[:,0], 5.63*1000*vmT0[:,1]+4126, yerr=5.63*1000*vmT0[:,2], fmt='m>',markersize=12, label=r'$B_s^* T_0$')
plt.legend(fontsize=15)
plt.xlabel(r'T MeV', fontsize=20)
plt.ylabel(r'M(T) MeV', fontsize=20)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.axhline(y=5366, color='b', linestyle='--', label=r'$B^0_S(exp)$')
plt.axhline(y=5415, color='r', linestyle='dashed', label=r'$B^*_S(exp)$')
plt.ylim(4600, 6200)
#plt.xlim(25,210)
plt.title('effmass')

with open(file_dir+"/effmass_fits_ps_16_s0_half.csv") as file:
    psm16 = pd.read_csv(file)
    psm16 = psm16.to_numpy()
    file.close()

with open(file_dir+"/effmass_fits_v_16_s0_half.csv") as file:
    vm16 = pd.read_csv(file)
    vm16 = vm16.to_numpy()
    file.close()

with open(file_dir+"/effmass_fits_ps_24_s0_half.csv") as file:
    psm24 = pd.read_csv(file)
    psm24 = psm24.to_numpy()
    file.close()

with open(file_dir+"/effmass_fits_v_24_s0_half.csv") as file:
    vm24 = pd.read_csv(file)
    vm24 = vm24.to_numpy()
    file.close()

with open(file_dir+"/effmass_fits_ps_32_s0_half.csv") as file:
    psm32 = pd.read_csv(file)
    psm32 = psm32.to_numpy()
    file.close()

with open(file_dir+"/effmass_fits_v_32_s0_half.csv") as file:
    vm32 = pd.read_csv(file)
    vm32 = vm32.to_numpy()
    file.close()
with open(file_dir+"/effmass_fits_ps_T0_s0_half.csv") as file:
    psmT0 = pd.read_csv(file)
    psmT0 = psmT0.to_numpy()
    file.close()

with open(file_dir+"/effmass_fits_v_T0_s0_half.csv") as file:
    vmT0 = pd.read_csv(file)
    vmT0 = vmT0.to_numpy()
    file.close()

#plot mass with Temp
plt.figure(figsize=(16,9))
plt.errorbar(5.63*1000/psm16[:,0], 5.63*1000*psm16[:,1]+4126, yerr=5.63*1000*psm16[:,2], fmt='bX',markersize=12)
plt.errorbar(5.63*1000/vm16[:,0], 5.63*1000*vm16[:,1]+4126, yerr=5.63*1000*vm16[:,2], fmt='r*',markersize=12)
plt.errorbar(5.63*1000/psm24[:,0], 5.63*1000*psm24[:,1]+4126, yerr=5.63*1000*psm24[:,2], fmt='bX',markersize=12, label=r'$B_s^0 Ns=24$')
plt.errorbar(5.63*1000/vm24[:,0], 5.63*1000*vm24[:,1]+4126, yerr=5.63*1000*vm24[:,2], fmt='r*',markersize=12, label=r'$B_s^* Ns=24$')
plt.errorbar(5.63*1000/psm32[:,0], 5.63*1000*psm32[:,1]+4126, yerr=5.63*1000*psm32[:,2], fmt='b+',markersize=12, label=r'$B_s^0 Ns=32$')
plt.errorbar(5.63*1000/vm32[:,0], 5.63*1000*vm32[:,1]+4126, yerr=5.63*1000*vm32[:,2], fmt='r.',markersize=12, label=r'$B_s^* Ns=32$')
plt.errorbar(5.63*1000/psmT0[:,0], 5.63*1000*psmT0[:,1]+4126, yerr=5.63*1000*psmT0[:,2], fmt='g^', markersize=12,label=r'$B_s^0 T_0$')
plt.errorbar(5.63*1000/vmT0[:,0], 5.63*1000*vmT0[:,1]+4126, yerr=5.63*1000*vmT0[:,2], fmt='m>',markersize=12, label=r'$B_s^* T_0$')
plt.legend(fontsize=15)
plt.xlabel(r'T MeV', fontsize=20)
plt.ylabel(r'M(T) MeV', fontsize=20)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.axhline(y=5366, color='b', linestyle='--', label=r'$B^0_S(exp)$')
plt.axhline(y=5415, color='r', linestyle='dashed', label=r'$B^*_S(exp)$')
#plt.ylim(4600, 6200)
#plt.xlim(25,210)
plt.title('effmass - half')

####corr
with open(file_dir+"/corr_fits_ps_16_s0.csv") as file:
    psm16 = pd.read_csv(file)
    psm16 = psm16.to_numpy()
    file.close()

with open(file_dir+"/corr_fits_v_16_s0.csv") as file:
    vm16 = pd.read_csv(file)
    vm16 = vm16.to_numpy()
    file.close()

with open(file_dir+"/corr_fits_ps_24_s0.csv") as file:
    psm24 = pd.read_csv(file)
    psm24 = psm24.to_numpy()
    file.close()

with open(file_dir+"/corr_fits_v_24_s0.csv") as file:
    vm24 = pd.read_csv(file)
    vm24 = vm24.to_numpy()
    file.close()

with open(file_dir+"/corr_fits_ps_32_s0.csv") as file:
    psm32 = pd.read_csv(file)
    psm32 = psm32.to_numpy()
    file.close()

with open(file_dir+"/corr_fits_v_32_s0.csv") as file:
    vm32 = pd.read_csv(file)
    vm32 = vm32.to_numpy()
    file.close()

with open(file_dir+"/corr_fits_ps_T0_s0.csv") as file:
    psmT0 = pd.read_csv(file)
    psmT0 = psmT0.to_numpy()
    file.close()

with open(file_dir+"/corr_fits_v_T0_s0.csv") as file:
    vmT0 = pd.read_csv(file)
    vmT0 = vmT0.to_numpy()
    file.close()

#plot corr with Temp
plt.figure(figsize=(16,9))
plt.errorbar(5.63*1000/psm16[:,0], 5.63*1000*psm16[:,1]+4126, yerr=5.63*1000*psm16[:,2], fmt='bX',markersize=12)
plt.errorbar(5.63*1000/vm16[:,0], 5.63*1000*vm16[:,1]+4126, yerr=5.63*1000*vm16[:,2], fmt='r*',markersize=12)
plt.errorbar(5.63*1000/psm24[:,0], 5.63*1000*psm24[:,1]+4126, yerr=5.63*1000*psm24[:,2], fmt='bX',markersize=12, label=r'$B_s^0 Ns=24$')
plt.errorbar(5.63*1000/vm24[:,0], 5.63*1000*vm24[:,1]+4126, yerr=5.63*1000*vm24[:,2], fmt='r*',markersize=12, label=r'$B_s^* Ns=24$')
plt.errorbar(5.63*1000/psm32[:,0], 5.63*1000*psm32[:,1]+4126, yerr=5.63*1000*psm32[:,2], fmt='b+',markersize=12, label=r'$B_s^0 Ns=32$')
plt.errorbar(5.63*1000/vm32[:,0], 5.63*1000*vm32[:,1]+4126, yerr=5.63*1000*vm32[:,2], fmt='r.',markersize=12, label=r'$B_s^* Ns=32$')
plt.errorbar(5.63*1000/psmT0[:,0], 5.63*1000*psmT0[:,1]+4126, yerr=5.63*1000*psmT0[:,2], fmt='g^', markersize=12,label=r'$B_s^0 T_0$')
plt.errorbar(5.63*1000/vmT0[:,0], 5.63*1000*vmT0[:,1]+4126, yerr=5.63*1000*vmT0[:,2], fmt='m>',markersize=12, label=r'$B_s^* T_0$')
plt.legend(fontsize=15)
plt.xlabel(r'T MeV', fontsize=20)
plt.ylabel(r'M(T) MeV', fontsize=20)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.axhline(y=5366, color='b', linestyle='--', label=r'$B^0_S(exp)$')
plt.axhline(y=5415, color='r', linestyle='dashed', label=r'$B^*_S(exp)$')
#plt.ylim(4600, 6200)
plt.title('corr')
#plt.xlim(25,210)



with open(file_dir+"/corr_fits_ps_16_s0_half.csv") as file:
    psm16 = pd.read_csv(file)
    psm16 = psm16.to_numpy()
    file.close()

with open(file_dir+"/corr_fits_v_16_s0_half.csv") as file:
    vm16 = pd.read_csv(file)
    vm16 = vm16.to_numpy()
    file.close()

with open(file_dir+"/corr_fits_ps_24_s0_half.csv") as file:
    psm24 = pd.read_csv(file)
    psm24 = psm24.to_numpy()
    file.close()

with open(file_dir+"/corr_fits_v_24_s0_half.csv") as file:
    vm24 = pd.read_csv(file)
    vm24 = vm24.to_numpy()
    file.close()

with open(file_dir+"/corr_fits_ps_32_s0_half.csv") as file:
    psm32 = pd.read_csv(file)
    psm32 = psm32.to_numpy()
    file.close()

with open(file_dir+"/corr_fits_v_32_s0_half.csv") as file:
    vm32 = pd.read_csv(file)
    vm32 = vm32.to_numpy()
    file.close()
with open(file_dir+"/corr_fits_ps_T0_s0_half.csv") as file:
    psmT0 = pd.read_csv(file)
    psmT0 = psmT0.to_numpy()
    file.close()

with open(file_dir+"/corr_fits_v_T0_s0_half.csv") as file:
    vmT0 = pd.read_csv(file)
    vmT0 = vmT0.to_numpy()
    file.close()

#plot mass with Temp
plt.figure(figsize=(16,9))
plt.errorbar(5.63*1000/psm16[:,0], 5.63*1000*psm16[:,1]+4126, yerr=5.63*1000*psm16[:,2], fmt='bX',markersize=12)
plt.errorbar(5.63*1000/vm16[:,0], 5.63*1000*vm16[:,1]+4126, yerr=5.63*1000*vm16[:,2], fmt='r*',markersize=12)
plt.errorbar(5.63*1000/psm24[:,0], 5.63*1000*psm24[:,1]+4126, yerr=5.63*1000*psm24[:,2], fmt='bX',markersize=12, label=r'$B_s^0 Ns=24$')
plt.errorbar(5.63*1000/vm24[:,0], 5.63*1000*vm24[:,1]+4126, yerr=5.63*1000*vm24[:,2], fmt='r*',markersize=12, label=r'$B_s^* Ns=24$')
plt.errorbar(5.63*1000/psm32[:,0], 5.63*1000*psm32[:,1]+4126, yerr=5.63*1000*psm32[:,2], fmt='b+',markersize=12, label=r'$B_s^0 Ns=32$')
plt.errorbar(5.63*1000/vm32[:,0], 5.63*1000*vm32[:,1]+4126, yerr=5.63*1000*vm32[:,2], fmt='r.',markersize=12, label=r'$B_s^* Ns=32$')
plt.errorbar(5.63*1000/psmT0[:,0], 5.63*1000*psmT0[:,1]+4126, yerr=5.63*1000*psmT0[:,2], fmt='g^', markersize=12,label=r'$B_s^0 T_0$')
plt.errorbar(5.63*1000/vmT0[:,0], 5.63*1000*vmT0[:,1]+4126, yerr=5.63*1000*vmT0[:,2], fmt='m>',markersize=12, label=r'$B_s^* T_0$')
plt.legend(fontsize=15)
plt.xlabel(r'T MeV', fontsize=20)
plt.ylabel(r'M(T) MeV', fontsize=20)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.axhline(y=5366, color='b', linestyle='--', label=r'$B^0_S(exp)$')
plt.axhline(y=5415, color='r', linestyle='dashed', label=r'$B^*_S(exp)$')
#plt.ylim(4600, 6200)
#plt.xlim(25,210)
plt.title('corr - half')


plt.show()


