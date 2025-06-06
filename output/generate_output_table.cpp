#include <iostream>
#include <fstream>
#include <iomanip>
#include <array>

#include "../src/eos_species/eos_assembled.hpp"

/* Function printing the EOS output to file */
void print_EOS_output(FullEOSOutput *eos_out, std::ostream& os) {
  const double etot = eos_out->e + 1.0E+39 * eos_out->nb * MEOS_MeV2erg * eos_out->nuEOS.Z_tot;
  const double ptot = eos_out->P + 1.0E+39 * eos_out->nb * MEOS_MeV2erg * eos_out->nuEOS.Z_tot / 3.;
  const double stot = eos_out->s + eos_out->nuEOS.s_tot;

  os <<  eos_out->rho                          << "  ";  // Mass density [g/cm3]
  os <<  eos_out->Y_part.yle                   << "  ";  // Electron lepton fraction [#/baryon]
  os <<  eos_out->Y_part.ylm                   << "  ";  // Muon lepton fraction [#/baryon]
  os <<  etot                                  << "  ";  // Internal energy density (including neutrinos) [erg/cm3]
  os <<  eos_out->nb*1.0E+39                   << "  ";  // Number density [1/cm3]
  os <<  eos_out->T                            << "  ";  // Temperature [MeV]
  os <<  eos_out->Y_part.ye                    << "  ";  // Electron fraction [#/baryon]
  os <<  eos_out->Y_part.ym                    << "  ";  // Muon fraction [#/baryon]
  os <<  eos_out->e                            << "  ";  // Internal energy density (fluid only) [erg/cm3]
  os <<  ptot                                  << "  ";  // Pressure (including neutrinos) [erg/cm3]
  os <<  stot                                  << "  ";  // Entropy per baryon (including neutrinos) [#/baryon]
  os <<  eos_out->Y_part.yh                    << "  ";  // Heavy nuclei fraction [#/baryon]
  os <<  eos_out->Y_part.ya                    << "  ";  // Alpha particle fraction [#/baryon]
  os <<  eos_out->Y_part.yp                    << "  ";  // Proton fraction [#/baryon]
  os <<  eos_out->Y_part.yn                    << "  ";  // Neutron fraction [#/baryon]
  os <<  eos_out->Y_part.ynue                  << "  ";  // Electron neutrino fraction [#/baryon]
  os <<  eos_out->Y_part.yanue                 << "  ";  // Electron antineutrino fraction [#/baryon]
  os <<  eos_out->Y_part.ynum                  << "  ";  // Muon neutrino fraction [#/baryon]
  os <<  eos_out->Y_part.yanum                 << "  ";  // Muon antineutrino fraction [#/baryon]
  os <<  eos_out->Y_part.ynux                  << "  ";  // Tau (anti)neutrino fraction [#/baryon]
  os <<  eos_out->nuEOS.Z_nue                  << "  ";  // Electron neutrino energy per baryon [MeV/baryon]
  os <<  eos_out->nuEOS.Z_anue                 << "  ";  // Electron antineutrino energy per baryon [MeV/baryon]
  os <<  eos_out->nuEOS.Z_num                  << "  ";  // Muon neutrino energy per baryon [MeV/baryon]
  os <<  eos_out->nuEOS.Z_anum                 << "  ";  // Muon antineutrino energy per baryon [MeV/baryon]
  os <<  eos_out->nuEOS.Z_nux                  << "  ";  // Tau (anti)neutrino energy per baryon [MeV/baryon]
  os <<  eos_out->chem_pot.mu_p - eos_out->mb  << "  ";  // Proton chemical potential wrt baryon (neutron) mass [MeV]
  os <<  eos_out->chem_pot.mu_n - eos_out->mb  << "  ";  // Neutron chemical potential wrt baryon (neutron) mass [MeV]
  os <<  eos_out->chem_pot.mu_e                << "  ";  // Electron chemical potential (with rest mass) [MeV]
  os <<  eos_out->chem_pot.mu_m                << "  ";  // Muon chemical potential (with rest mass) [MeV]
  os << std::endl;

  return;
}

int main () {
  /* Boolean variable for inclusion of muons */
  const bool with_mu = true; // false
	
  std::cout << std::scientific; // scientific format for std output
    	
  /* Print to std output if EOS contains muons contribution or not */
  std::cout << "#####################" << std::endl;
  if (with_mu == true) {
    std::cout << "#   EOS with muons  #" << std::endl;
  } else {
    std::cout << "# EOS without muons #" << std::endl;
  }
  std::cout << "#####################" << std::endl;
  std::cout << std::endl;

  /* Name of baryon EOS table */
  std::string BarTableName = "DD2_bar.h5";  // baryon table

  /* Initialize global EOS class

  Constructor -> EOS_assembled(const int id_eos, const bool el_flag, const bool mu_flag, std::string BarTableName)

  Inputs:
   - id_EOS: method for EOS computation (1: interpolation, 2: on-the-fly)
   - el_flag: flag for activating electrons
   - mu_flag: flag for activating muons
   - BarTableName: path of baryon EOS table  */
  EOS_assembled eos(2, true, with_mu, BarTableName);
  
  /* Define name of output table */
  std::string table_name;
  if (with_mu == true) {
    table_name = "eos_global_w_mu.txt";
  } else {
    table_name = "eos_global_wo_mu.txt";
  }

  /* Output stream */
  std::ofstream Iout("output/data/"+table_name);
  Iout << std::scientific << std::setprecision(10); // scientific format for output table

  /* Define table input arrays 
	
    Input parameters:
      - nb: baryon number density [1/fm3]
      -  T: temperature [MeV]
	    - yq: charge fraction [#/baryon] (ye = yq if muons are not included)
	    - ym: muon fraction   [#/baryon] (only needed if muons are included, the electron fraction is then
                                        obtained from charge neutrality, ye = yq - ym, and discarded in the
                                        loop in the case it is negative)

    For (nb, T, yq) we take the same grid as in the input baryon EOS CompOSE table, while for ym we define
    a uniformly log-spaced grid
  */

  /* Table dimensions */
  const int n_nb = eos.m_nn; // Number density
  const int n_t  = eos.m_nt; // Temperature
  const int n_yq = eos.m_ny; // Charge fraction
  const int n_ym = 60;       // Muon fraction

  /* Arrays of input variables */
  // Number density
  double const * n_array  = eos.GetRawLogNumberDensity(); // log spaced
	
  // Temperature
  double const * t_array  = eos.GetRawLogTemperature(); // log spaced
    	
  // Charge fraction
  double const * yq_array = eos.GetRawYq(); // lin spaced
	
  // Muon fraction
  std::vector<double> ym_array; // log spaced
	
  const double ymmin = 2.0e-05;
  const double ymmax = 2.0e-01;

  const double log_ymmin = log10(ymmin);
  const double log_ymmax = log10(ymmax);
        
  for (int i=0; i<n_ym; i++) ym_array.push_back(log_ymmin + static_cast<double>(i) * (log_ymmax-log_ymmin) / static_cast<double>(n_ym));

  /* Output stream -> legend */
  Iout << "# d [g/cm3], yle [#/baryon], ylm [#/baryon], u [erg/cm3], nb [1/cm3], T [MeV], ye [#/baryon], ";
  Iout <<   "e [erg/cm3], P [erg/cm3], s [#/baryon], yh [#/baryon], ya [#/baryon], yp [#/baryon], yn [#/baryon], ";
  Iout <<   "ynue [#/baryon], yanue [#/baryon], ynum [#/baryon], yanum [#/baryon], ynux [#/baryon], ";
  Iout <<   "znue [MeV/baryon], zanue [MeV/baryon], znum [MeV/baryon], zanum [MeV/baryon], znux [MeV/baryon], ";
  Iout <<   "mup wrt mn [MeV], mun wrt mn [MeV], mue w me [MeV]" << std::endl;

  double nb, T;
  double yq, ye, ym;
  double Y[2] = {0.0};
	
  /* Choose step size for iteration over input variable arrays */
  const int di = 10, dj = 10, dk = 10, dl = 10; //di: nb step, dj: T step, dk: yq step, dl: ym step

  FullEOSOutput eos_out;

  if (with_mu == true) {
    /* Write table WITH MUONS to output file */
    std::cout << "# id_n, id_t, id_ye, id_ym,  nb [1/fm3], T [MeV], Yq [#/baryon], Ym [#/baryon]" << std::endl;
    for (int i=0; i<n_nb; i+=di) {
      nb = exp(n_array[i]);
      for (int j=0; j<n_t; j+=dj) {
        T = exp(t_array[j]);
        for (int k=0; k<n_yq; k+=dk) {
          yq = yq_array[k];
          for (int l=0; l<n_ym; l+=dl) {
            ym = pow(10.,ym_array[l]);
            ye = yq - ym;

            if (ye <= 0.) continue; // exclude EOS point if Ye is negative

            Y[0] = ye;
            Y[1] = ym;

            // Print input quantities to standard output
            std::cout << i  << "  " << j << "  " << k    << "  " << l    << "  ";
            std::cout << nb << "  " << T << "  " << Y[0] << "  " << Y[1] << std::endl;

            eos_out = eos.compute_full_EOS(nb, T, Y); // compute full EOS output
            print_EOS_output(&eos_out, Iout);  // print output to file
          }
        }
      }
    }
  } else {
    /* Write table WITHOUT MUONS to output file */
    std::cout << "# id_n, id_t, id_ye, nb [1/fm3], T [MeV], Yq [#/baryon]" << std::endl;
    for (int i=0; i<n_nb; i+=di) {
      nb = exp(n_array[i]);
      for (int j=0; j<n_t; j+=dj) {
        T = exp(t_array[j]);
        for (int k=0; k<n_yq; k+=dk) {
          Y[0] = yq_array[k];

          // Print input quantities to standard output
          std::cout << i  << "  " << j << "  " << k    << "  ";
          std::cout << nb << "  " << T << "  " << Y[0] << std::endl;

          eos_out = eos.compute_full_EOS(nb, T, Y); // compute full EOS output
          print_EOS_output(&eos_out, Iout);   // print output to file
        }  
      }
    }
  }

  Iout.close(); // close output stream

  return 0;
}
