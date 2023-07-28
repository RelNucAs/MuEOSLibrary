import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler

plt.rcParams['axes.prop_cycle'] = cycler(linestyle=['-', ':'])*cycler(color=['#023eff','#ff7c00', '#1ac938', '#e8000b', '#8b2be2', '#9f4800', '#f14cc1', '#ffc400', '#00d7ff', 'k'])

input_file_cold = "diff_hist_DD2_ylmu_last.txt"
input_file_zero = "diff_hist_DD2.txt"

##############################
# Select correct folder name #
##############################

#folder_name = "data_ylmu_0"     # table with initial Ylmu = 0
#folder_name = "data_ylmu_cold"  # table with cold initial Ylmu
folder_name = "data_ylmu_tot"    # concatenation of the two tables

idx_1, nb_1, temp_1, ye_1, ym_1 = np.loadtxt(input_file_cold, comments='#', skiprows=1, usecols=(0,1,2,3,4), unpack=True)
idx_2, nb_2, temp_2, ye_2, ym_2 = np.loadtxt(input_file_zero, comments='#', skiprows=1, usecols=(0,1,2,3,4), unpack=True)

diff_list_cold = np.loadtxt(input_file_cold, comments='#', skiprows=1, usecols=(5,6,7,8,9,10,11,12,13,14,15), unpack=True)
diff_list_zero = np.loadtxt(input_file_zero, comments='#', skiprows=1, usecols=(5,6,7,8,9,10,11,12,13,14,15), unpack=True)

diff_list_tot = np.concatenate((diff_list_cold, diff_list_zero), axis=1)

if (folder_name == "data_ylmu_cold"):
    diff_list = diff_list_cold
elif (folder_name == "data_ylmu_0"):
    diff_list = diff_list_zero
elif (folder_name == "data_ylmu_tot"):
    diff_list = diff_list_tot
else:
    print("Wrong folder name!")
    exit()

leg_list = ["Total pressure", "Total entropy", "Y_nue", "Y_anue", "Y_num", "Y_anum", "Y_nux", "Neutron chemical potential", "Proton chemical potential", "Electron chemical potential", "Muon chemical potential"]

num_per_dec = 8
hist_edges = np.logspace(np.log10(1.0E-08), np.log10(1.0E+00), num=8*num_per_dec + 1)

for id_d, diff_data in enumerate(diff_list):
    fig = plt.figure()
    plt.hist(diff_data, bins=hist_edges, histtype='bar') #, label='$q=%.2lf$' %q)
    plt.xscale("log")
    plt.xticks(np.logspace(np.log10(1.0E-08), np.log10(1.0E+00), num=9))
    plt.gca().tick_params(which='both')
    plt.xlabel("Relative difference")
    plt.ylabel("Counts")
    plt.title(leg_list[id_d])
    fig.savefig("plots/"+folder_name+"/hist_%d.png" %id_d, dpi=300, bbox_inches='tight')
