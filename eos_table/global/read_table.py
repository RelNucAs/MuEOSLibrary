import numpy as np

#Dictionary of input and output files
files = {'postproc_table':'./data_DD2_ylmu_last.txt', #post-processing table
         'postproc_grid':'./grid_DD2.txt', #grid
        }

if __name__ == '__main__':

   f = files
   
   """Load post-processing grid getting tabulated density [g/cm^3], temperature [MeV] and net Ye. These are the initial input variables before postprocessing, spanning the entire thermodynamics space """
   
   print("Elaborating post-proc grid...")
   
   #load the postproc_grid table file with dens, initial temp, initial net Ye
   pp_grid = open(f['postproc_grid'],'r')

   #read the grid dimension
   line = pp_grid.readline()
   nd,nt,ny = line.split()
   
   nd = int(nd)
   nt = int(nt)
   ny = int(ny)

   #read dens, temp, Ye of grid in numpy arrays
   line = pp_grid.readline()
   d_array = np.fromstring(line, dtype=float, sep=' ')
   line = pp_grid.readline()
   t_array = np.fromstring(line, dtype=float, sep=' ')
   line = pp_grid.readline()
   ye_array = np.fromstring(line, dtype=float, sep=' ')
   
   pp_grid.close()
   
   #compute the log10 for density and temperature
   ld_array = np.log10(d_array)
   lt_array = np.log10(t_array)
   
   #compute minima and maxima
   ldmin = np.amin(ld_array)
   ldmax = np.amax(ld_array)
   ltmin = np.amin(lt_array)
   ltmax = np.amax(lt_array)
   ymin = np.amin(ye_array)
   ymax = np.amax(ye_array)
   
   print("Done!")
   
   """Load the post-processing data, define thermodynamic quantities and interpolating functions."""

   print("Elaborating post-proc data...")
   
   den,tin,yein,sin,pin,fyel,fypos,fymu,fyamu,pfin,sfin,tfin,yne,yane,ynmu,yanmu,yntau,mu_n,mu_pr,mu_el,mu_muons,flag = np.loadtxt(f['postproc_table'],dtype = 'float',unpack=True)
   
   #den: matter density [g/cm3]
   #tin: initial temperature [MeV] 
   #yein: initial net electron fraction []
   #sin: initial entropy [kb/baryon]
   #pin: initial pressure [erg/cm^3]
   #fyel: final fraction of electrons []
   #fypos: final fraction of positrons []
   #fymu: final fraction of muons []
   #fyamu: final fraction of anti-muons []
   #pfin: final pressure [erg/cm^3]
   #sfin: final entropy [kb/baryon]
   #tfin: final temperature [MeV]
   #yne: final electronic neutrino fraction []
   #yane: final electronic anti-neutrino fraction []
   #ynmu: final muonic neutrino fraction []
   #yanmu: final muonic anti-neutrino fraction []
   #yntau: final tauonic neutrino fraction []
   #mu_n: final non relat. neutrons chemic. potential [MeV]
   #mu_pr: final non relat. protons chemic. potential [MeV]
   #mu_el: final non relat. electrons chemic. potential [MeV]
   #mu_muons: final non relat. muons chemic. potential [MeV]
   #flag: ?

   #flip arrays to order them in growing dens, temp, Ye
   den = den[::-1]
   tin = tin[::-1]
   yein = yein[::-1]
   sin = sin[::-1]
   pin = pin[::-1]
   fyel = fyel[::-1]
   fypos = fypos[::-1]
   fymu = fymu[::-1]
   fyamu = fyamu[::-1]
   pfin = pfin[::-1]
   sfin = sfin[::-1]
   tfin = tfin[::-1]
   yne = yne[::-1]
   yane = yane[::-1]
   ynmu = ynmu[::-1]
   yanmu = yanmu[::-1]
   yntau = yntau[::-1]
   mu_n = mu_n[::-1]
   mu_pr = mu_pr[::-1]
   mu_el = mu_el[::-1]
   mu_muons = mu_muons[::-1]
   flag = flag[::-1]
