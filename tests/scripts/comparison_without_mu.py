import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

## Comparison of EOS table without muons against Albino's Fortran reference table

# Table dimensions
n_nb = 60 # Number density
n_t  = 27 # Temperature
n_ye = 20 # Electron fraction

# Read and reshape output table
C_data = np.loadtxt("../output/table_comp_without_mu.txt")
C_data = C_data.reshape((n_nb,n_t,n_ye,21))

# Read and reshape reference table
F_data = np.loadtxt("../data/ic_file_table.txt", comments="#", skiprows=2)
F_data = F_data.reshape((n_nb,n_t,n_ye,21))

# Get input arrays
nb = F_data[:,0,0,3] * 1.0E-39 # Number density array [fm-3]
T  = F_data[0,:,0,4]           # Temperature [MeV]
ye = F_data[0,0,:,5]           # Electron fraction array

# Compute relative difference between tables
diff = np.where((C_data == 0.) & (F_data == 0.), 1.e-50, abs(C_data-F_data)/abs(F_data))

# Absoulute difference for nuclear abundances (yh, ya, yp, yn)
diff[...,8:12] = abs(C_data[...,8:12]-F_data[...,8:12])
diff[...,8:12] = np.where(diff[...,8:12] != 0., diff[...,8:12], 1.e-50)
diff = np.where(diff == 0., 1.e-50, diff)

# Define list of plot titles
tit_list = [r"$d~[{\rm g}~{\rm cm}^{-3}]$", r"$Y_{\rm l_{\rm e}}~[{\rm baryon}^{-1}]$", r"$u~[{\rm erg}~{\rm cm}^{-3}]$", r"$n_{\rm B}~[{\rm cm}^{-3}]$", r"$T~[{\rm MeV}]$", r"Y_{\rm e}~[{\rm baryon}^{-1}]$", r"$e~[{\rm erg}~{\rm cm}^{-3}]$", r"$P~[{\rm erg}~{\rm cm}^{-3}]$", r"$Y_{\rm H}~[{\rm baryon}^{-1}]$", r"$Y_\alpha~[{\rm baryon}^{-1}]$", r"$Y_{\rm p}~[{\rm baryon}^{-1}]$", r"$Y_{\rm n}~[{\rm baryon}^{-1}]$", r"$Y_{\nu_{\rm e}}~[{\rm baryon}^{-1}]$", r"$Y_{\bar{\nu}_{\rm e}}~[{\rm baryon}^{-1}]$", r"$Y_{\nu_{\rm x}}~[{\rm baryon}^{-1}]$", r"$Z_{\nu_{\rm e}}~[{\rm MeV}~{\rm baryon}^{-1}]$", r"$Z_{\bar{\nu}_{\rm e}}~[{\rm MeV}~{\rm baryon}^{-1}]$", r"$Z_{\nu_{\rm x}}~[{\rm MeV}~{\rm baryon}^{-1}]$", r"$\mu^0_{\rm p}~[{\rm MeV}]$", r"$\mu^0_{\rm n}~[{\rm MeV}]$", r"$\mu_{\rm e}~[{\rm MeV}]$"]

# Choose quantities to be plotted
diff_cut = diff[...,6:]
tit_list_cut = tit_list[6:]

data = diff_cut

# Define plot_folder
plot_dir = "../plots/comparison/without_mu/"

# Make plots about internal energy density comparison
nnrows = 2 # Number of rows
nncols = 5 # Number of columns

fig, axs = plt.subplots(nrows=nnrows, ncols=nncols, figsize=(20,6.8))
plt.subplots_adjust(top=0.92,bottom=0.05,right=0.95,left=0.05) #,hspace=0.25,wspace=0.17)
    
for id_ax in range(nnrows*nncols):
    id_row = int(id_ax/nncols)
    id_col = int(id_ax%nncols)
    ax = axs[id_row,id_col]
    id_ye = 2*id_ax
    print("Ye = %.2lf" %ye[id_ye])
    data_plot = np.transpose(data[:,:,id_ye,0]) # Internal energy density
    cs = ax.pcolormesh(nb, T, data_plot, norm=colors.LogNorm(vmin=1.e-5,vmax=1.e3), shading="nearest")
    (max_idt, max_idn)  = np.where(data_plot == np.amax(data_plot))
    #print(np.amin(data_plot))
    #print(np.amax(data_plot))
    #ax.scatter(nb[max_idn[0]], T[max_idt[0]], marker="x", color="k")
    diff_sign = np.where(C_data[:,:,id_ye,6]*F_data[:,:,id_ye,6]<0.)
    ax.scatter(nb[diff_sign[0]], T[diff_sign[1]], marker=".", edgecolor=None, facecolor="k", alpha=0.2)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_title(r"$Y_{\rm e} = %.2lf$ (%.1e)" %(ye[id_ye],np.amax(data_plot)))
    fig.colorbar(cs,ax=axs[id_row,id_col],pad=0.02) #,label="Relative difference")
    plt.suptitle("Internal energy density relative difference")
    if (id_row == nnrows-1):    ax.set_xlabel("Number density [fm-3]")
    if (id_col == 0):           ax.set_ylabel("Temperature [MeV]")
plt.savefig(plot_dir + "int_energy_comp.png", dpi=200, bbox_inches="tight")
plt.close()

# Make plots about internal energy density comparison
nnrows = 2
nncols = 2

for idx in range(0,data.shape[2],2):
    fig, axs = plt.subplots(nrows=nnrows, ncols=nncols, figsize=(8,6.8))
    plt.subplots_adjust(top=0.92,bottom=0.05,right=0.95,left=0.05) #,hspace=0.25,wspace=0.17)
    print("Ye = %.2lf" %ye[idx])
    cs1 = axs[0][0].pcolormesh(nb, T, np.transpose(C_data[:,:,idx,10]), norm=colors.Normalize(vmin=0.,vmax=1.), shading="nearest")
    cs2 = axs[0][1].pcolormesh(nb, T, np.transpose(F_data[:,:,idx,10]), norm=colors.Normalize(vmin=0.,vmax=1.), shading="nearest")
    cs3 = axs[1][0].pcolormesh(nb, T, np.transpose(C_data[:,:,idx,11]), norm=colors.Normalize(vmin=0.,vmax=1.), shading="nearest")
    cs4 = axs[1][1].pcolormesh(nb, T, np.transpose(F_data[:,:,idx,11]), norm=colors.Normalize(vmin=0.,vmax=1.), shading="nearest")
    for i in range(2):
        for j in range(2):
            axs[i][j].set_xscale("log")
            axs[i][j].set_yscale("log")
    fig.colorbar(cs2,ax=axs[0,:].ravel().tolist(),pad=0.02,label="Proton fraction")
    fig.colorbar(cs4,ax=axs[1,:].ravel().tolist(),pad=0.02,label="Neutron fraction")
    for i in range(2):
        axs[i][0].set_title("C++")
        axs[i][1].set_title("Fortran")
        axs[1][i].set_xlabel("Number density [fm-3]")
        axs[i][0].set_ylabel("Temperature [MeV]")
    plt.suptitle(r"$Y_{\rm e} = %.2lf$" %ye[idx])
    plt.savefig(plot_dir + "nucleon_frac_comp_idx_%d.png" %idx, dpi=200, bbox_inches="tight")
    plt.close()


# Make plots about general comparison
nnrows = 3
nncols = 5

for idx in range(data.shape[2]):
    print(idx)
    fig, axs = plt.subplots(nrows=nnrows, ncols=nncols, figsize=(20,10)) #sharex='all', sharey='all'
    plt.subplots_adjust(top=0.92,bottom=0.05,right=0.95,left=0.05) #,hspace=0.25,wspace=0.17)
    for id_ax in range(nnrows*nncols):
        id_row = int(id_ax/nncols)
        id_col = int(id_ax%nncols)
        ax = axs[id_row,id_col]
        data_plot = np.transpose(data[:,:,idx,id_ax])
        cs = ax.pcolormesh(nb, T, data_plot, norm=colors.LogNorm(vmin=1.e-5,vmax=1.e0), shading="nearest")
        (max_idt, max_idn)  = np.where(data_plot == np.amax(data_plot))
        if (id_ax == 0):      ax.scatter(nb[max_idn[0]], T[max_idt[0]], marker="x", color="k")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(tit_list_cut[id_ax]+" (%.1e)" %np.amax(data_plot))
        if (id_col == nncols-1):    fig.colorbar(cs,ax=axs[id_row,:].ravel().tolist(),pad=0.02,label="Relative difference")
        plt.suptitle(r"$Y_{\rm e} = %.2lf$" %ye[idx])
        if (id_row == nnrows-1):    ax.set_xlabel("Number density [fm-3]")
        if (id_col == 0):           ax.set_ylabel("Temperature [MeV]")
    plt.savefig(plot_dir + "general_comp_idx_%d.png" %idx, dpi=200, bbox_inches="tight")
    plt.close()

