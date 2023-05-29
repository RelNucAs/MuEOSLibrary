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

num_file = table_folder + sp + '/eos_'+sp+'_primitive_new_cs2_num_fine_log.txt'
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

#for i in range(len(rec_file)):
#    print(np.where((data_num[i] == 0.) & (data_rec[i] != 0.)))

title_list = [r'dPdn', r'dsdn', r'dPdt', r'dsdt']

nnrows = 3
nncols = len(diff)

# Plot figure
fig, axs = plt.subplots(nrows=nnrows, ncols=nncols, figsize=(20,12))
plt.subplots_adjust(top=0.95,bottom=0.05,right=0.95,left=0.05) #,hspace=0.25,wspace=0.17)
for id_ax in range(nncols):
    id_col = int(id_ax%nncols)
    
    ax = axs[0,id_col]
    data_plot = abs(data_num[id_ax])
    data_plot = np.where(data_plot == 0., 1.e-10, data_plot)
    cs1 = ax.pcolormesh(ne_array, t_array, data_plot, norm=colors.LogNorm(vmin=np.amin(data_plot),vmax=np.amax(data_plot)), shading="nearest")
    ax.set_title(title_list[id_col]+" (Max = %.2e)" %np.amax(data_plot))
    
    ax = axs[1,id_col]
    data_plot = abs(data_rec[id_ax])
    data_plot = np.where(data_plot == 0., 1.e-10, data_plot)
    cs2 = ax.pcolormesh(ne_array, t_array, data_plot, norm=colors.LogNorm(vmin=np.amin(data_plot),vmax=np.amax(data_plot)), shading="nearest")
    ax.set_title(title_list[id_col]+" (Max = %.2e)" %np.amax(data_plot))
    
    ax = axs[2,id_col]
    data_plot = diff[id_ax]
    check_inf = np.transpose(np.argwhere(data_plot == np.inf))
    ax.scatter(ne_array[check_inf[1]], t_array[check_inf[0]], marker=".", facecolor="k", edgecolor=None, alpha=0.2)
    data_plot = np.where(data_plot == np.inf, 1.0, data_plot)
    cs3 = ax.pcolormesh(ne_array, t_array, data_plot, norm=colors.LogNorm(vmin=1.e-5,vmax=1.e0), shading="nearest")
    ax.set_title(title_list[id_col]+" (Max = %.2e)" %np.amax(data_plot))
    for i in range(nnrows):
        axs[i,id_col].set_xscale("log")
        axs[i,id_col].set_yscale("log")
        axs[i,id_col].set_xlabel("Number density [fm-3]")
    fig.colorbar(cs1,ax=axs[0,id_col],pad=0.01,label=r"Finite difference")
    fig.colorbar(cs2,ax=axs[1,id_col],pad=0.01,label=r"Recursion relations")
    fig.colorbar(cs3,ax=axs[2,id_col],pad=0.01,label=r"Relative difference")
    #plt.suptitle(r"Finite difference vs recursive formula")
    if (id_col == 0):
        for i in range(nnrows):
            axs[i,id_col].set_ylabel("Temperature [MeV]")
plt.savefig(plot_folder+sp+"/compare_lep_der_fine_log.png", dpi=200, bbox_inches="tight")
plt.close()
