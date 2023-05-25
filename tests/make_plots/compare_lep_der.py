import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

me = 0.510998928 #MeV
mmu = 105.6583745 #MeV
MeV = 1.602176634e-6 
kB = 8.617333262145e-11

species = 1
if (species == 1):
    sp = 'electrons'
    n1 = 700
    mL = me
elif (species == 2):
    sp = 'muons'
    n1 = 750
    mL = mmu

abs_path = "/home/leonardo/Desktop/PhD_work/BNS_muons/EOS_module/"
table_folder = abs_path + "eos_table/"
plot_folder = abs_path + "tests/output/"

num_file = table_folder + sp + '/eos_'+sp+'_primitive_new_cs2_num.txt'
n1  = np.loadtxt(num_file,skiprows=0      ,max_rows=1,dtype=int)
n2  = np.loadtxt(num_file,skiprows=1      ,max_rows=1,dtype=int)
ne_array = 10.**(np.loadtxt(num_file,skiprows=2,max_rows=1,dtype=float))
t_array  = 10.**(np.loadtxt(num_file,skiprows=3,max_rows=1,dtype=float))
dPdn_num = np.loadtxt(num_file,skiprows=4+n1*9 ,max_rows=n1,unpack=True,dtype=float)
dsdn_num = np.loadtxt(num_file,skiprows=4+n1*10,max_rows=n1,unpack=True,dtype=float)
dPdt_num = np.loadtxt(num_file,skiprows=4+n1*11,max_rows=n1,unpack=True,dtype=float)
dsdt_num = np.loadtxt(num_file,skiprows=4+n1*12,max_rows=n1,unpack=True,dtype=float)

data_num = [dPdn_num, dsdn_num, dPdt_num, dsdt_num]

rec_file = table_folder + sp + '/eos_'+sp+'_primitive_new_cs2.txt'
dPdn_rec = np.loadtxt(rec_file,skiprows=4+n1*9 ,max_rows=n1,unpack=True,dtype=float)
dsdn_rec = np.loadtxt(rec_file,skiprows=4+n1*10,max_rows=n1,unpack=True,dtype=float)
dPdt_rec = np.loadtxt(rec_file,skiprows=4+n1*11,max_rows=n1,unpack=True,dtype=float)
dsdt_rec = np.loadtxt(rec_file,skiprows=4+n1*12,max_rows=n1,unpack=True,dtype=float)

data_rec = [dPdn_rec, dsdn_rec, dPdt_rec, dsdt_rec]

diff = [np.where(data_num[idx] == data_rec[idx], 1.e-10, abs(data_num[idx] - data_rec[idx]) / abs(data_num[idx])) for idx in range(len(data_num))]

title_list = [r'dPdn', r'dsdn', r'dPdt', r'dsdt']

nncols = len(diff)

# Plot figure
fig, axs = plt.subplots(nrows=1, ncols=4, figsize=(20,4))
plt.subplots_adjust(top=0.85,bottom=0.05,right=0.95,left=0.05) #,hspace=0.25,wspace=0.17)
for id_ax in range(nncols):
    id_col = int(id_ax%nncols)
    ax = axs[id_col]
    data_plot = diff[id_ax]
    cs = ax.pcolormesh(ne_array, t_array, data_plot, norm=colors.LogNorm(vmin=1.e-5,vmax=1.e0), shading="nearest")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_title(title_list[id_col]+" (Max = %.2e)" %np.amax(data_plot))
    if (id_col == nncols-1):    fig.colorbar(cs,ax=axs[:],pad=0.01,label=r"Relative difference")
    plt.suptitle(r"Finite difference vs recursive formula")
    ax.set_xlabel("Number density [fm-3]")
    if (id_col == 0):   ax.set_ylabel("Temperature [MeV]")
plt.savefig(plot_folder+sp+"/compare_lep_der.png", dpi=200, bbox_inches="tight")
plt.close()
