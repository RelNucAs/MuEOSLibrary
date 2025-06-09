#include <iostream>
#include <fstream>
#include <iomanip>
#include <array>
#include <cmath>

#include "eos_species/eos_species.hpp"

/* Function printing to file the EOS output */
void print_EOS_output(const double nb, const double T, const double *Y, EOSstruct *out, std::ostream& os) {
  const double etot = out->e + 1.0E+39 * nb * MEOS_MeV2erg * out->nuEOS.Z_tot;
  const double ptot = out->P + 1.0E+39 * nb * MEOS_MeV2erg * out->nuEOS.Z_tot / 3.;
  const double stot = out->s + out->nuEOS.s_tot;

  //os << out->rho                     << " ";  // Mass density [g/cm3]
  //os << out->Y_part.yle              << " ";  // Net electron lepton fraction [#/baryon]
  os << etot                         << " ";  // Internal energy density (including neutrinos) [erg/cm3]
  os << 1.0E+39 * nb            << " ";  // Number density [1/cm3]
  os << T                       << " ";  // Temperature [MeV]
  os << out->comp[4]               << " ";  // Electron fraction [#/baryon]
  os << out->e                       << " ";  // Internal energy density (without neutrinos) [erg/cm3]
  os << ptot                         << " ";  // Pressure (including neutrinos) [erg/cm3]
  //os << stot                         << " " << "#/baryon"   << std::endl;
  os << out->comp[1]               << " ";  // Heavy nuclei fraction [#/baryon]
  os << out->comp[0]               << " ";  // Alpha particle fraction [#/baryon]
  os << out->comp[3]               << " ";  // Proton fraction [#/baryon]
  os << out->comp[2]               << " ";  // Neutron fraction [#/baryon]
  os << out->nuEOS.Y_nu[0]             << " ";  // Electron neutrino fraction [#/baryon]
  os << out->nuEOS.Y_nu[1]           << " ";  // Electron antineutrino fraction [#/baryon]
  os << out->nuEOS.Y_nu[4]             << " ";  // Heavy (anti)neutrino fraction [#/baryon]
  os << out->nuEOS.Z_nue             << " ";  // Electron neutrino energy per baryon [MeV/baryon]
  os << out->nuEOS.Z_anue            << " ";  // Electron antineutrino energy per baryon [MeV/baryon]
  os << out->nuEOS.Z_nux             << " ";  // Heavy (anti)neutrino energy per baryon [MeV/baryon]
  os << out->chem_pot[0] - out->mb << " ";  // Proton chemical potential wrt baryon (neutron) mass [MeV]
  os << out->chem_pot[1] - out->mb << " ";  // Neutron chemical potential wrt baryon (neutron) mass [MeV]
  os << out->chem_pot[2]           << " ";  // Electron chemical potential (with rest mass) [MeV]
  os << std::endl;
  
  return;
}

int main () {
  std::cout << std::scientific; // scientific format for std output
    	
  std::cout << "#####################" << std::endl;
  std::cout << "# EOS without muons #" << std::endl;
  std::cout << "#####################" << std::endl;
  std::cout << std::endl;

  /* Name of baryon EOS table */
  std::string BarTableName = "DD2_bar.h5";  // baryon table

  /* Initialize global EOS class

  Constructor -> MuEOSClass(const int id_eos, const bool el_flag, const bool mu_flag, std::string BarTableName)

  Inputs:
   - id_EOS: method for EOS computation (1: interpolation, 2: on-the-fly)
   - el_flag: flag for activating electrons
   - mu_flag: flag for activating muons
   - BarTableName: path of baryon EOS table  */
  MuEOSClass eos(id_test, false, BarTableName);      

  /* Define name of output table */
  std::string out_table = "table_comp_without_mu.txt";

  /* Output stream */
  std::ofstream Iout("tests/output/" + out_table);
  Iout << std::scientific << std::setprecision(10); // scientific format for output table

  /* Define table input paramters
	
	Input parameters:
	- nb: baryon number density [1/fm3]
	-  T: temperature [MeV]
	- ye: electron fraction [#/baryon]

  */

  // The following parameters are calibrated on the shape of the grid in "tests/data/ic_file_table.txt"

  /* Table dimensions */
  const int n_nb = 60; // Number density
  const int n_t  = 27; // Temperature
  const int n_ye = 57; // Electron fraction

  /* Arrays of input variables */
  // Number density
  const double nbmin = 7.58576778E-07;
  const double nbmax = 9.12011723E+00;

  const double d_nb = pow(nbmax / nbmin, 1./(n_nb-1.));

  std::vector<double> n_array; // log spaced
  for (int i=0; i<n_nb; i++) n_array.push_back(nbmin * pow(d_nb,static_cast<double>(i)));
	
  // Temperature
  const double tmin = 1.00000000E-01;
  const double tmax = 1.31825639E+02;

  const double d_t = pow(tmax / tmin, 1./(n_t-1.));
	
  std::vector<double> t_array; // log spaced
  for (int i=0; i<n_t; i++) t_array.push_back(tmin * pow(d_t,static_cast<double>(i)));
	
  // Electron fraction
  const double yemin = 9.99999978E-03;
  const double yemax = 5.79999983E-01;

  std::vector<double> ye_array; // lin spaced
  for (int i=0; i<n_ye+1; i++) ye_array.push_back(yemin + static_cast<double>(i) * (yemax-yemin) / static_cast<double>(n_ye));


  /* Output stream -> legend */
  Iout << "# d [g/cm3], yle [#/baryon], u [erg/cm3], nb [1/cm3], T [MeV], ye [#/baryon], ";
  Iout <<   "e [erg/cm3], P [erg/cm3], yh [#/baryon], ya [#/baryon], yp [#/baryon], yn [#/baryon], ";
  Iout <<   "ynue [#/baryon], yanue [#/baryon], ynux [#/baryon], ";
  Iout <<   "znue [MeV/baryon], zanue [MeV/baryon], znux [MeV/baryon], ";
  Iout <<   "mup wrt mn [MeV], mun wrt mn [MeV], mue w me [MeV]" << std::endl;

  double nb, T;
  double Y[2];
	
  Y[1] = 0.0; // set muon fraction equal to zero

  /* Choose step size for iteration over input variable arrays */
  const int di = 1, dj = 1, dk = 3; //di: nb step, dj: T step, dk: ye step

  EOSstruct eos_out;

  /* Write table WITHOUT MUONS to output file */
  std::cout << "# id_n, id_t, id_ye, nb [1/fm3], T [MeV], Ye [#/baryon]" << std::endl;
  for (int i=0; i<n_nb; i+=di) {
    nb = n_array[i];
    for (int j=0; j<n_t; j+=dj) {
      T = t_array[j];
      for (int k=0; k<n_ye+1; k+=dk) {
        Y[0] = ye_array[k];

        std::cout << i  << "  " << j << "  " << k    << "  ";
        std::cout << nb << "  " << T << "  " << Y[0] << std::endl;

        eos_out = eos.compute_full_EOS(nb, T, Y);
		    print_EOS_output(nb, T, Y, &eos_out, Iout);
      }
    }
  }

  Iout.close(); // close output stream

  return 0;
}
