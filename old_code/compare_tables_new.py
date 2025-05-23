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

complete_file = table_folder + sp + '/eos_'+sp+'_primitive_new_cs2.txt' #with_eta_ele.txt'
n1  = np.loadtxt(complete_file,skiprows=0      ,max_rows=1,dtype=int)
n2  = np.loadtxt(complete_file,skiprows=1      ,max_rows=1,dtype=int)
ne_array = 10.**(np.loadtxt(complete_file,skiprows=2,max_rows=1,dtype=float))
t_array  = 10.**(np.loadtxt(complete_file,skiprows=3,max_rows=1,dtype=float))
n   = np.loadtxt(complete_file,skiprows=4+n1*0,max_rows=n1,unpack=True,dtype=float)
a_n = np.loadtxt(complete_file,skiprows=4+n1*1,max_rows=n1,unpack=True,dtype=float)
P   = np.loadtxt(complete_file,skiprows=4+n1*2,max_rows=n1,unpack=True,dtype=float)
a_P = np.loadtxt(complete_file,skiprows=4+n1*3,max_rows=n1,unpack=True,dtype=float)
e   = np.loadtxt(complete_file,skiprows=4+n1*4,max_rows=n1,unpack=True,dtype=float)
a_e = np.loadtxt(complete_file,skiprows=4+n1*5,max_rows=n1,unpack=True,dtype=float)
s   = np.loadtxt(complete_file,skiprows=4+n1*6,max_rows=n1,unpack=True,dtype=float)
a_s = np.loadtxt(complete_file,skiprows=4+n1*7,max_rows=n1,unpack=True,dtype=float)
mu  = np.loadtxt(complete_file,skiprows=4+n1*8,max_rows=n1,unpack=True,dtype=float)

[X, Y] = np.meshgrid(ne_array,t_array)

therm_leo = np.array([P, e, s*kB, a_P, a_e, a_s*kB]) * 1.e39 * MeV
#therm_leo = np.where(therm_leo>0.,therm_leo,1.e-50)
#print(np.amin(therm_leo[3]))

complete_file = table_folder + sp + '/eos_'+sp+'_complete_v4.txt' #_bis.txt
mu  = np.loadtxt(complete_file,skiprows=4,max_rows=700,unpack=True,dtype=float)
P   = np.loadtxt(complete_file,skiprows=4+n1*1,max_rows=n1,unpack=True,dtype=float)
a_P = np.loadtxt(complete_file,skiprows=4+n1*2,max_rows=n1,unpack=True,dtype=float)
e   = np.loadtxt(complete_file,skiprows=4+n1*3,max_rows=n1,unpack=True,dtype=float)
a_e = np.loadtxt(complete_file,skiprows=4+n1*4,max_rows=n1,unpack=True,dtype=float)
s   = np.loadtxt(complete_file,skiprows=4+n1*5,max_rows=n1,unpack=True,dtype=float)
a_s = np.loadtxt(complete_file,skiprows=4+n1*6,max_rows=n1,unpack=True,dtype=float)

therm_ele = np.array([P, e, s, a_P, a_e, a_s])
#therm_ele = np.where(therm_ele>0.,therm_ele,1.e-50)
#print(np.amin(therm_ele[3]))

diff = abs(therm_leo-therm_ele)/therm_ele
#diff = np.where(therm_leo > 0., diff, 1.e-10)
diff = np.where(diff>0.,diff,1.e-10)

print([np.amax(tmp) for tmp in diff])
titles = [r'Pressure $P_{L}$', r'Int. energy $e_{L}$', r'Entropy $s_{L}$', r'Anti pressure $P_{\bar{L}}$', r'Anti int. energy $e_{\bar{L}}$', r'Positron entropy $s_{\bar{L}}$']

#plot figure
for i in range(int(therm_leo.shape[0]/3)):
    fig, axs = plt.subplots(2, 3, sharex='all', sharey='all', figsize=(18,8))
    #plt.suptitle(titles[i])
    for j in range(int(therm_leo.shape[0]/2)):
        #c1 = axs[0][0].pcolormesh(X,Y,therm_leo[i]  ,norm=colors.LogNorm(vmin=np.amin(therm_leo[i]), vmax=np.amax(therm_leo[i])),shading='nearest')
        if (i==0):
            [low,upp,n] = col_levels(np.amin(therm_leo[i*3+j]),np.amax(therm_leo[i*3+j]))
            c1 = axs[0][j].contourf(X,Y,therm_leo[i*3+j],norm=colors.LogNorm(vmin=10**low,vmax=10**upp),levels=np.logspace(low,upp,num=n),extend='both')
            #c2 = axs[1][j].pcolormesh(X,Y,diff[i*3+j]     ,norm=colors.LogNorm(vmin=np.amin(diff[i*3+j]), vmax=np.amax(diff[i*3+j])),shading='nearest')
        else:
            [low,upp,n] = pad_to_divide(1.e-10,np.amax(therm_leo[i*3+j]))
            c1 = axs[0][j].contourf(X,Y,therm_leo[i*3+j],norm=colors.LogNorm(vmin=1.e-10,vmax=np.amax(therm_leo[i*3+j])),levels=10**np.arange(-10,math.floor(np.log10(np.amax(therm_leo[i*3+j]))),6,dtype=np.float),extend='both')

        c2 = axs[1][j].contourf(X,Y,diff[i*3+j]     ,norm=colors.LogNorm(vmin=1.e-8, vmax=1.e-4),levels=np.logspace(-8,-4,num=5),extend='both')
        print(np.amin(diff[i*3+j]))

        #c3 = axs[1][0].pcolormesh(X,Y,therm_leo[i+3],norm=colors.LogNorm(vmin=1.e-10, vmax=np.amax(therm_leo[i+3])),shading='nearest')
        #c4 = axs[1][1].pcolormesh(X,Y,diff[i+3]     ,norm=colors.LogNorm(vmin=np.amin(diff[i+3]), vmax=np.amax(diff[i+3])),shading='nearest')

        axs[1][j].set_xlabel(r'$n_L$ [fm$^{-3}$]')

        axs[0][j].set_title(titles[i*3+j])
        axs[1][j].set_title('Rel diff (max = %.2e)' %np.amax(diff[i*3+j]))
       
        plt.colorbar(c1,ax=axs[0][j])
        plt.colorbar(c2,ax=axs[1][j])
        
    axs[0][0].set_xscale('log')
    axs[0][0].set_yscale('log')
    axs[0][0].set_ylabel(r'$T$ [MeV]')
    axs[1][0].set_ylabel(r'$T$ [MeV]')
    
    plt.subplots_adjust(top=0.95,bottom=0.1,right=0.8,left=0.1,wspace=0.12)
    plt.savefig(plot_folder + sp + "/compare_eos_"+sp+"_v4_%d.png" %i,dpi=300,bbox_inches='tight')
    plt.close()

