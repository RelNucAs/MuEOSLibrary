import os
import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler

## Make histograms of relative differences in the comparison with Eleonora's reference table (with muons)

# Cycler for plot colors
plt.rcParams['axes.prop_cycle'] = cycler(linestyle=['-', ':'])*cycler(color=['#023eff','#ff7c00', '#1ac938', '#e8000b', '#8b2be2', '#9f4800', '#f14cc1', '#ffc400', '#00d7ff', 'k'])

# Define table names
input_file_cold = "../output/diff_hist_DD2_ylmu_last.txt"
input_file_zero = "../output/diff_hist_DD2.txt"

##############################
# Select correct folder name #
##############################

#folder_name = "data_ylmu_zero"     # table with initial Ylmu = 0
#folder_name = "data_ylmu_cold"  # table with cold initial Ylmu
folder_name = "data_ylmu_tot"    # concatenation of the two tables

# Read input arrays and EOS data from tables
try:
    idx_1, nb_1, temp_1, ye_1, ym_1 = np.loadtxt(input_file_cold, comments='#', skiprows=1, usecols=(0,1,2,3,4), unpack=True)
    diff_list_cold = np.loadtxt(input_file_cold, comments='#', skiprows=1, usecols=(5,6,7,8,9,10,11,12,13,14,15), unpack=True)
except:
    print(input_file_cold + " not found")
        
try:
    idx_2, nb_2, temp_2, ye_2, ym_2 = np.loadtxt(input_file_zero, comments='#', skiprows=1, usecols=(0,1,2,3,4), unpack=True)
    diff_list_zero = np.loadtxt(input_file_zero, comments='#', skiprows=1, usecols=(5,6,7,8,9,10,11,12,13,14,15), unpack=True)
except:
    print(input_file_zero + " not found")

# Concatenate the data of the two tables
try:
    diff_list_tot = np.concatenate((diff_list_cold, diff_list_zero), axis=1)
except:
    if ("diff_list_cold" in vars()):
        diff_list_tot = diff_list_cold
    elif ("diff_list_zero" in vars()):
        diff_list_tot = diff_list_zero


# Select data to be plotted
if (folder_name == "data_ylmu_cold"):
    diff_list = diff_list_cold
elif (folder_name == "data_ylmu_zero"):
    diff_list = diff_list_zero
elif (folder_name == "data_ylmu_tot"):
    diff_list = diff_list_tot
else:
    print("Wrong folder name!")
    exit()

# Define list of plot titles
leg_list = ["Total pressure", "Total entropy", "Y_nue", "Y_anue", "Y_num", "Y_anum", "Y_nux", "Neutron chemical potential", "Proton chemical potential", "Electron chemical potential", "Muon chemical potential"]

# Define plot folder and create it if not already existing
plot_dir = "../plots/comparison/with_mu/" + folder_name + "/"
cmd = "mkdir -p " + plot_dir
os.system(cmd)

# Define histogram edges
x_min = 1.0E-08  # minimum difference
x_max = 1.0E+00  # maximum difference
num_per_dec = 8  # number of edges per decade 
hist_edges = np.logspace(np.log10(x_min), np.log10(x_max), num=8*num_per_dec + 1)

# Make histogram plots
for id_d, diff_data in enumerate(diff_list):
    fig = plt.figure()
    plt.hist(diff_data, bins=hist_edges, histtype='bar') #, label='$q=%.2lf$' %q)
    plt.xscale("log")
    plt.xticks(np.logspace(np.log10(1.0E-08), np.log10(1.0E+00), num=9))
    plt.gca().tick_params(which='both')
    plt.xlabel("Relative difference")
    plt.ylabel("Counts")
    plt.title(leg_list[id_d])
    fig.savefig(plot_dir + "hist_comp_idx_%d.png" %id_d, dpi=300, bbox_inches='tight')
    plt.close()