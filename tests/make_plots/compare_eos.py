import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import argparse

parser = argparse.ArgumentParser(description="Choosing methods to compare")
parser.add_argument("-a", dest="method_a", required=True, type=int)
parser.add_argument("-b", dest="method_b", required=True, type=int)
args = parser.parse_args() #pass arguments

if ((args.method_a != 1) and (args.method_a != 2) and (args.method_a != 3)):
    print("-a must be 1, 2 or 3")
    exit()

if ((args.method_b != 1) and (args.method_b != 2) and (args.method_b != 3)):
    print("-b must be 1, 2 or 3")
    exit()

if (args.method_a == args.method_b):
    print("-a and -b must be different")
    exit()

HR = False

if HR:
    res = '_HR'
else:
    res = ''

species = 2
if (species==1):
    sp = 'electrons'
elif (species==2):
    sp = 'muons'

abs_path = "/home/leonardo/Desktop/PhD_work/BNS_muons/EOS_module/"
table_folder = abs_path + "eos_table/"
input_folder = abs_path + "tests/" + sp + "/"
plot_folder = abs_path + "tests/output/"

input_table = table_folder + sp + "/eos_"+sp+"_leo"+res+".txt"

#grid parameters
nne = np.loadtxt(input_table,skiprows=0,max_rows=1,dtype=int)
nt  = np.loadtxt(input_table,skiprows=1,max_rows=1,dtype=int)

#number density array
ne_edges = 10.**(np.loadtxt(input_table,skiprows=2,max_rows=1,dtype=float))
ne_array = np.array([0.5*(ne_edges[i]+ne_edges[i+1]) for i in range(ne_edges.size-1)])

#temperature array
t_edges  = 10.**(np.loadtxt(input_table,skiprows=3,max_rows=1,dtype=float))
t_array = np.array([0.5*(t_edges[i]+t_edges[i+1]) for i in range(t_edges.size-1)])

#eta_file = '../eos_table/eos_electrons_v2.txt'
#eta_ele = np.loadtxt(eta_file,skiprows=4,max_rows=700,unpack=True,dtype=float)

#import data
data = np.loadtxt(input_folder+"/test_accuracy_"+sp+res+".txt",unpack=True)
data_size = data.shape[0]

#reshape arrays
data = np.reshape(data,(data_size,nne-1,nt-1))

ne_rand = data[0]
t_rand  = data[1]

first  = data[2:11]
second = data[11:20]
third  = data[20:29]

eos = [first[2:8], second[2:8], third[2:8]]
#[p1, a_p1, e1, a_e1, s1, a_s1] = first[3:]

eos_sum = [np.array([tmp[2*i]+tmp[2*i+1] for i in range(int(len(tmp)/2))]) for tmp in eos]

titles = [r'Pressure $P_{L}$', r'Anti Pressure $P_{\bar{L}}$', r'Int. energy $e_{L}$', r'Anti Int. energy $e_{\bar{L}}$', r'Entropy $s_{L}$', r'Anti entropy $s_{\bar{L}}$']

id_a = args.method_a
id_b = args.method_b

ref = eos[0]
Y1  = eos[id_a-1]
Y2  = eos[id_b-1]

ref_sum = eos_sum[0]
Y1_sum  = eos_sum[id_a-1]
Y2_sum  = eos_sum[id_b-1]

thres = 1.e-10

Y = abs(Y1-Y2)/Y1
Y = np.where(Y>0., Y, 1.e-6) 
ref = np.where(ref>0.,ref,1.e-50)

Y_sum = abs(Y1_sum-Y2_sum)/Y1_sum
Y_sum = np.where(Y_sum>0., Y_sum, 1.e-6) 

#plot EOS (with method 1) for particles and antiparticles
fig, axs = plt.subplots(2, 3, sharex='all', sharey='all', figsize=(16,8))
plt.suptitle('Method 1')
for i in range(Y.shape[0]):
    ax = axs[i%2][int(i/2)]
    tmp = np.transpose(ref[i])
    #cx = ax.pcolormesh(ne_rand,t_rand,tmp,norm=colors.LogNorm(vmin=1.e-10,vmax=np.amax(tmp)),shading='gouraud')
    if (i%2==0):
        #cx = ax.pcolormesh(ne_rand,t_rand,tmp,norm=colors.LogNorm(vmin=np.amin(tmp),vmax=np.amax(tmp)),shading='gouraud')
        cx = ax.contourf(ne_array,t_array,tmp,norm=colors.LogNorm(vmin=np.amin(tmp),vmax=np.amax(tmp)),levels=10**np.arange(math.floor(np.log10(np.amin(tmp))),math.floor(np.log10(np.amax(tmp))),2,dtype=np.float),extend='both')
    else:
        #cx = ax.pcolormesh(ne_rand,t_rand,tmp,norm=colors.LogNorm(vmin=1.e-10,vmax=np.amax(tmp)),shading='gouraud')
        cx = ax.contourf(ne_array,t_array,tmp,norm=colors.LogNorm(vmin=1.e-10,vmax=np.amax(tmp)),levels=10**np.arange(-10,math.floor(np.log10(np.amax(tmp))),2,dtype=np.float),extend='both')
    plt.colorbar(cx,ax=ax)
    ax.set_title(titles[i])

axs[0][0].set_xscale('log')
axs[0][0].set_yscale('log')
for i in range(3):
    axs[1][i].set_xlabel(r'$n_L$ [fm$^{-3}$]')
for i in range(2):
    axs[i][0].set_ylabel(r'$T$ [MeV]')

plt.subplots_adjust(top=0.85,bottom=0.1,right=0.8,left=0.1,wspace=0.12)
plt.savefig(plot_folder + sp + "/plot_eos_onthefly_"+sp+res+".png",dpi=200,bbox_inches='tight')
plt.close()

#plot EOS comparison of particles and antiparticles separately 
fig, axs = plt.subplots(2, 3, sharex='all', sharey='all', figsize=(16,8))
plt.suptitle('Method %d vs method %d' %(id_a,id_b))
for i in range(Y.shape[0]):
    ax = axs[i%2][int(i/2)]
    tmp = np.transpose(Y[i])
    #cx = ax.pcolormesh(ne_rand,t_rand,tmp,norm=colors.LogNorm(vmin=1.e-5,vmax=1.e-1),shading='gouraud')
    cx = ax.contourf(ne_array,t_array,tmp,norm=colors.LogNorm(vmin=1.e-5,vmax=1.e0),levels=np.logspace(-5,0,num=6),extend='both')
    plt.colorbar(cx,ax=ax)
    ax.set_title(titles[i]+' (Max diff: %.2e)' %np.amax(tmp))
axs[0][0].set_xscale('log')
axs[0][0].set_yscale('log')
for i in range(3):
    axs[1][i].set_xlabel(r'$n_L$ [fm$^{-3}$]')
for i in range(2):
    axs[i][0].set_ylabel(r'$T$ [MeV]')

plt.subplots_adjust(top=0.85,bottom=0.1,right=0.8,left=0.1,wspace=0.12)
plt.savefig(plot_folder + sp + "/test_eos_"+sp+"_%d_vs_%d" %(id_a,id_b)+res+".png",dpi=200,bbox_inches='tight')
plt.close()




#plot EOS comparison of sum particles + antiparticles
fig, axs = plt.subplots(1, 3, sharex='all', sharey='all', figsize=(15,3.5))
plt.suptitle('Method %d vs method %d' %(id_a,id_b))
for i in range(int(Y.shape[0]/2)):
    ax = axs[i]
    tmp = np.transpose(Y_sum[i])
    #cx = ax.pcolormesh(ne_rand,t_rand,tmp,norm=colors.LogNorm(vmin=1.e-5,vmax=1.e-1),shading='gouraud')
    cx = ax.contourf(ne_array,t_array,tmp,norm=colors.LogNorm(vmin=1.e-5,vmax=1.e0),levels=np.logspace(-5,0,num=6),extend='both')
    plt.colorbar(cx,ax=ax)
    ax.set_title(titles[2*i]+' (Max diff: %.2e)' %np.amax(tmp))
    ax.set_xlabel(r'$n_L$ [fm$^{-3}$]')
axs[0].set_xscale('log')
axs[0].set_yscale('log')
axs[0].set_ylabel(r'$T$ [MeV]')
plt.subplots_adjust(top=0.85,bottom=0.1,right=0.8,left=0.1,wspace=0.1)
plt.savefig(plot_folder + sp + "/test_eos_sum_"+sp+"_%d_vs_%d" %(id_a,id_b)+res+".png",dpi=200,bbox_inches='tight')
plt.close()

