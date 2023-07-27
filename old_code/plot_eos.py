import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

def col_levels(a,b):
    low = math.floor(np.log10(a))
    upp = math.floor(np.log10(b)) 
    if ((upp-low)%2 == 0):
        n = int((upp-low)/2)
    else:
        n = int((upp-low+1)/2)
        low -= 1
    return [low,upp,n]

def pad_to_divide(a,b):
    low = math.floor(np.log10(a))
    upp = math.floor(np.log10(b))
    i = 0
    while(i == 0):
        if ((upp-low)%4 == 0):
            n = int((upp-low)/4)
            i = 1
        else:
            upp += 1
    return [low,upp,n]

me = 0.510998928 #MeV
mmu = 105.6583745 #MeV
MeV = 1.602176634E-06

species = 2
if (species == 1):
    sp = 'electrons'
    n1 = 700
    mL = me
elif (species == 2):
    sp = 'muons'
    n1 = 750
    mL = mmu

HR = False
if HR: 
    res = "_HR"
else:
    res = ""

abs_path = "/home/leonardo/Desktop/PhD_work/BNS_muons/EOS_module/"
table_folder = abs_path + "eos_table/"
plot_folder = abs_path + "tests/output/"

complete_file = table_folder + sp + '/eos_'+sp+'_complete_leo.txt' #with_eta_ele.txt'
ne  = np.loadtxt(complete_file,skiprows=0,max_rows=1,dtype=int)
nt  = np.loadtxt(complete_file,skiprows=1,max_rows=1,dtype=int)
ne_array = 10.**(np.loadtxt(complete_file,skiprows=2,max_rows=1,dtype=float))
t_array  = 10.**(np.loadtxt(complete_file,skiprows=3,max_rows=1,dtype=float))
mu  = np.loadtxt(complete_file,skiprows=4+n1*0,max_rows=n1,unpack=True,dtype=float)
n   = np.loadtxt(complete_file,skiprows=4+n1*1,max_rows=n1,unpack=True,dtype=float)
a_n = np.loadtxt(complete_file,skiprows=4+n1*2,max_rows=n1,unpack=True,dtype=float)
P   = np.loadtxt(complete_file,skiprows=4+n1*3,max_rows=n1,unpack=True,dtype=float)
a_P = np.loadtxt(complete_file,skiprows=4+n1*4,max_rows=n1,unpack=True,dtype=float)
e   = np.loadtxt(complete_file,skiprows=4+n1*5,max_rows=n1,unpack=True,dtype=float)
a_e = np.loadtxt(complete_file,skiprows=4+n1*6,max_rows=n1,unpack=True,dtype=float)
s   = np.loadtxt(complete_file,skiprows=4+n1*7,max_rows=n1,unpack=True,dtype=float)
a_s = np.loadtxt(complete_file,skiprows=4+n1*8,max_rows=n1,unpack=True,dtype=float)

[X, Y] = np.meshgrid(ne_array,t_array)
therm_leo = np.array([P+a_P, e+a_e, s+a_s])
#therm_leo = np.where(therm_leo>0.,therm_leo,1.e-50)
#print(np.amin(therm_leo[3]))

titles = [r'Pressure $P_{L}$', r'Int. energy $e_{L}$', r'Entropy $s_{L}$']

#plot EOS comparison of sum particles + antiparticles
fig, axs = plt.subplots(1, 3, sharex='all', sharey='all', figsize=(15,3.5))
plt.suptitle('Plot summed EOS: '+sp)
for i in range(therm_leo.shape[0]):
    ax = axs[i]
    tmp = np.transpose(therm_leo[i])
    tmp = therm_leo[i]
    #cx = ax.pcolormesh(ne_rand,t_rand,tmp,norm=colors.LogNorm(vmin=1.e-5,vmax=1.e-1),shading='gouraud')
    cx = ax.contourf(ne_array,t_array,tmp,norm=colors.LogNorm(vmin=np.amin(tmp),vmax=np.amax(tmp)),extend='both') #,levels=np.logspace(-5,0,num=6),extend='both')
    plt.colorbar(cx,ax=ax)
    ax.set_title(titles[i])
    ax.set_xlabel(r'$n_L$ [fm$^{-3}$]')
axs[0].set_xscale('log')
axs[0].set_yscale('log')
axs[0].set_ylabel(r'$T$ [MeV]')
plt.subplots_adjust(top=0.85,bottom=0.1,right=0.8,left=0.1,wspace=0.1)
plt.savefig(plot_folder + sp + "/plot_eos_sum_"+sp+res+".png",dpi=200,bbox_inches='tight')
plt.close()

exit()

#plot for presentation
fig = plt.figure(figsize=(4.6,3.7))
tmp = therm_leo[0]*1.E-39/MeV
cx = plt.pcolormesh(ne_array,t_array,tmp,norm=colors.LogNorm(vmin=np.amin(tmp),vmax=np.amax(tmp)),shading='gouraud')
#cx = ax.contourf(ne_array,t_array,tmp,norm=colors.LogNorm(vmin=np.amin(tmp),vmax=np.amax(tmp)),extend='both') #,levels=np.logspace(-5,0,num=6),extend='both')
cbar = plt.colorbar(cx)
cbar.set_label(r'[${\rm MeV}\,{\rm fm}^{-3}$]')
#plt.title(titles[0])
plt.title(r"Pressure:$\quad\mu^+$ + $\mu^-$")
plt.xlabel(r'$n_L$ [fm$^{-3}$]')
plt.xscale('log')
plt.yscale('log')
plt.ylabel(r'$T$ [MeV]')
plt.savefig(plot_folder + sp + "/plot_pres_"+sp+res+".png",dpi=200,bbox_inches='tight')
plt.close()

