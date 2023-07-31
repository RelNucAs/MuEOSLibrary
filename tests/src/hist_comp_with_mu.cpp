
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <array>
#include <cmath>

#include "../../src/eos_species/eos_assembled.hpp"

#define n_var 22

#define ym_min   2.0E-05
#define ym_max   2.0E-01
#define m_e      0.51099895         // electron mass in MeV
#define m_mu     105.65837          // muon mass in MeV
#define m_n      939.56542052       // neutron mass in MeV
#define m_p      938.27208816       // proton mass in MeV
#define m_b      1.674E-24          // baryonic mass in MeV
#define amu_g    1.66054E-24        // atomic mass unit in g
#define de_lim   1.0E11             // exponential fading of electron (anti)neutrinos in g/cm3
#define dx_lim   1.0E12             // exponential fading of heavy (anti)neutrinos in g/cm3
#define nb_thr   0.01               // minimum number density for neutrino trapping in fm-3
#define corner_check true           // boolean variable for checking corner cases

/* Function computing the EOS output */
FullEOSOutput compute_EOS_nu_thres(EOS_assembled* eos, double nb, double T, double *Y) {
  FullEOSOutput eos_out;

  eos_out = eos->compute_full_EOS(nb, T, Y);

  const double rho = 1.0E+39 * nb * amu_g;   // use amu_g for conversion instead of CompOSE baryon
                                             // mass to be consistent with Eleonora's choice
  eos_out.rho = rho;  // Mass density [g/cm3]
  
  eos_out.Y_part.ynue  = exp(-de_lim/rho) * eos_out.Y_part.ynue;  // Electron neutrino fraction [#/baryon]
  eos_out.Y_part.yanue = exp(-de_lim/rho) * eos_out.Y_part.yanue; // Electron antineutrino fraction [#/baryon]
  eos_out.Y_part.ynum  = exp(-dx_lim/rho) * eos_out.Y_part.ynum;  // Muon neutrino fraction [#/baryon]
  eos_out.Y_part.yanum = exp(-dx_lim/rho) * eos_out.Y_part.yanum; // Muon antineutrino fraction [#/baryon]
  eos_out.Y_part.ynux  = exp(-dx_lim/rho) * eos_out.Y_part.ynux;  // Tau (anti)neutrino fraction [#/baryon]

  eos_out.Y_part.yle = eos_out.Y_part.ye + eos_out.Y_part.ynue - eos_out.Y_part.yanue; // Net electron lepton fraction [#/baryon]
  eos_out.Y_part.ylm = eos_out.Y_part.ym + eos_out.Y_part.ynum - eos_out.Y_part.yanum; // Net muon lepton fraction [#/baryon]

  eos_out.nuEOS.Y_nu = eos_out.Y_part;

  eos_out.nuEOS.Z_nue  = exp(-de_lim/rho) * eos_out.nuEOS.Z_nue;   // Electron neutrino energy per baryon [MeV/baryon]
  eos_out.nuEOS.Z_anue = exp(-de_lim/rho) * eos_out.nuEOS.Z_anue;  // Electron antineutrino energy per baryon [MeV/baryon]
  eos_out.nuEOS.Z_num  = exp(-dx_lim/rho) * eos_out.nuEOS.Z_num;   // Muon neutrino energy per baryon [MeV/baryon]
  eos_out.nuEOS.Z_anum = exp(-dx_lim/rho) * eos_out.nuEOS.Z_anum;  // Muon antineutrino energy per baryon [MeV/baryon]
  eos_out.nuEOS.Z_nux  = exp(-dx_lim/rho) * eos_out.nuEOS.Z_nux;   // Tau (anti)neutrino energy per baryon [MeV/baryon]

  eos_out.nuEOS.Z_tot = eos_out.nuEOS.Z_nue + eos_out.nuEOS.Z_anue + eos_out.nuEOS.Z_num + eos_out.nuEOS.Z_anum + 2. * eos_out.nuEOS.Z_nux;

  eos_out.nuEOS.s_nue  = exp(-de_lim/rho) * eos_out.nuEOS.s_nue;
  eos_out.nuEOS.s_anue = exp(-de_lim/rho) * eos_out.nuEOS.s_anue;
  eos_out.nuEOS.s_num  = exp(-dx_lim/rho) * eos_out.nuEOS.s_num;
  eos_out.nuEOS.s_anum = exp(-dx_lim/rho) * eos_out.nuEOS.s_anum;
  eos_out.nuEOS.s_nux  = exp(-dx_lim/rho) * eos_out.nuEOS.s_nux;

  eos_out.nuEOS.s_tot = eos_out.nuEOS.s_nue + eos_out.nuEOS.s_anue + eos_out.nuEOS.s_num + eos_out.nuEOS.s_anum + 2. * eos_out.nuEOS.s_nux;       

  return eos_out;
}

double compare_single_point(std::vector<std::array<double,n_var>>* tab_in,
  const int idx_p, const int idx_ref, const double q_eos) {
  const double q_ref = (*tab_in)[idx_p][idx_ref];
  double q_diff = fabs(q_eos-q_ref) / fabs(q_ref);
  return q_diff;
}

bool IsValid(EOS_assembled* eos, std::vector<std::array<double,n_var>>* tab_in, const int idx_p) {
  double nb = (*tab_in)[idx_p][0] * 1.0E-39 / amu_g;
  double  T = (*tab_in)[idx_p][11];
  double ye = (*tab_in)[idx_p][5] - (*tab_in)[idx_p][6];
  double ym = (*tab_in)[idx_p][7] - (*tab_in)[idx_p][8];
  double yq = ye + ym;
    
  const double * Cnb = eos->GetRawLogNumberDensity();
  const double *  CT = eos->GetRawLogTemperature();
  const double * Cyq = eos->GetRawYq();

  // Check if point is valid
  if ((nb < nb_thr))                                      return false;
  if ((nb < exp(Cnb[0])) || (nb > exp(Cnb[eos->m_nn-1]))) return false;
  if (( T < exp( CT[0])) || ( T > exp( CT[eos->m_nt-1]))) return false;
  if ((yq < Cyq[0]) || (yq > Cyq[eos->m_ny-1]))           return false;
  if ( ye < 0. )                                          return false;
  if ((ym < ym_min) || (ym > ym_max))                     return false;

  return true;
}


void print_rel_diff(EOS_assembled* eos, std::vector<std::array<double,n_var>>* tab_in, const int idx_p, std::ostream& os) {
  double nb = (*tab_in)[idx_p][0] * 1.0E-39 / amu_g;
  double  T = (*tab_in)[idx_p][11];
  double ye = (*tab_in)[idx_p][5] - (*tab_in)[idx_p][6];
  double ym = (*tab_in)[idx_p][7] - (*tab_in)[idx_p][8];
  //double yq = ye + ym;
  double Y[2] = {0.0, 0.0};
  Y[0] = ye;
  Y[1] = ym;

  // Print input quantities
  os << idx_p << "   ";
  os << nb    << "   ";
  os << T     << "   ";
  os << ye    << "   ";
  os << ym    << "   ";
    
  // Compute EOS
  FullEOSOutput eos_out = compute_EOS_nu_thres(eos, nb, T, Y);

  const double ptot = eos_out.P + 1.0E+39 * nb * MeV * eos_out.nuEOS.Z_tot / 3.;
  const double stot = eos_out.s + eos_out.nuEOS.s_tot;

  // Compare EOS against reference table
  os << compare_single_point(tab_in, idx_p,  9, ptot)                        << "   ";  // Total Pressure
  os << compare_single_point(tab_in, idx_p, 10, stot)                        << "   ";  // Total Entropy
  os << compare_single_point(tab_in, idx_p, 12, eos_out.nuEOS.Y_nu.ynue)     << "   ";  // nu_e fraction  
  os << compare_single_point(tab_in, idx_p, 13, eos_out.nuEOS.Y_nu.yanue)    << "   ";  // anu_e fraction
  os << compare_single_point(tab_in, idx_p, 14, eos_out.nuEOS.Y_nu.ynum)     << "   ";  // nu_mu  fraction
  os << compare_single_point(tab_in, idx_p, 15, eos_out.nuEOS.Y_nu.yanum)    << "   ";  // anu_mu fraction
  os << compare_single_point(tab_in, idx_p, 16, eos_out.nuEOS.Y_nu.ynux)     << "   ";  // nu_tau fraction  
  os << compare_single_point(tab_in, idx_p, 17, eos_out.chem_pot.mu_n - m_n) << "   ";  // Neutrons  NR chem potential 
  os << compare_single_point(tab_in, idx_p, 18, eos_out.chem_pot.mu_p - m_p) << "   ";  // Protons   NR chem potential 
  os << compare_single_point(tab_in, idx_p, 19, eos_out.chem_pot.mu_e - m_e) << "   ";  // Electrons NR chem potential
  os << compare_single_point(tab_in, idx_p, 20, eos_out.chem_pot.mu_m -m_mu) << "   ";  // Muons     NR chem potential 
  os << std::endl;

  return;
}




int main () {
  /* Name of baryon EOS table */
  std::string BarTableName = "eos_table/baryons/DD2_bar.h5";  // baryon table
  
  /* Initialize global EOS class

  Constructor -> EOS_assembled(const int id_eos, const bool el_bool, const bool mu_bool, std::string BarTableName)

  Inputs:
   - id_EOS: method for EOS computation (1: interpolation, 2: on-the-fly)
   - el_bool: flag for activating electrons
   - mu_bool: flag for activating muons
   - BarTableName: path of baryon EOS table  */
  EOS_assembled eos(2, true, true, BarTableName);

  /* Input stream */
  //std::string ref_table = "data_DD2_ylmu_last.txt"; // initial Yl_mu = Ym_cold
  std::string ref_table = "data_DD2.txt";           // initial Yl_mu = 0.
  std::ifstream fin("tests/data/" + ref_table);
    
  /* Print header */
  std::cout << "#####################" << std::endl;
  std::cout << "#   EOS with muons  #" << std::endl;
  std::cout << "#####################" << std::endl;
  std::cout << std::endl;
   
  /* Read Baryon EOS table input arrays */
  const double * nb_bar = eos.GetRawLogNumberDensity();
  const double *  T_bar = eos.GetRawLogTemperature();
  const double * yq_bar = eos.GetRawYq();

  /* Print conditions for validity of comparison */
  std::cout << "Comparing global EOS with muons with Eleonora's table: " << ref_table << std::endl;
  std::cout << std::endl;
  std::cout << "Points must lie within the range of validity of the EOS:"     << std::endl;
  std::cout << " - Number density    : [" << exp(nb_bar[0]) << ", " << exp(nb_bar[eos.m_nn-1]) << "]" << " " << "fm-3" << std::endl;
  std::cout << " - Temperature       : [" << exp( T_bar[0]) << ", " << exp( T_bar[eos.m_nt-1]) << "]" << " " <<  "MeV" << std::endl;
  std::cout << " - Charge fraction   : [" <<      yq_bar[0] << ", " <<     yq_bar[eos.m_ny-1]  << "]" << std::endl;
  std::cout << " - Muon fraction     : [" <<         ym_min << ", " <<                  ym_max << "]" << std::endl;
  std::cout << " - Electron fraction : positive"                                                      << std::endl;
  std::cout << std::endl;
  std::cout << "Number density must be also > " << nb_thr << " fm-3 " << "to justify neutrino trapped" << std::endl;
  std::cout << std::endl;

  /* Input vector */
  std::array<double,n_var> data_in;
  std::vector<std::array<double,n_var>> table_data;

  std::string dummyLine;
    
  /* Read data from input table */
  int np = 0;
  std::cout << "Reading data from the input table!" << std::endl;
  while (!fin.fail()) {
    std::getline(fin, dummyLine);
    std::stringstream ss(dummyLine);
        
    for (int i=0; i<n_var; i++) ss >> data_in[i];
    table_data.push_back(data_in);
    np += 1;
  }    

  fin.close();
  std::cout << "Done!" << std::endl;
  std::cout << std::endl;

   
  /* Define input paramters for EOS calculation
  *
  * Pick n_rand random EOS points from the input table
  * and read the input parameters
  *
  * Input parameters:
  * - nb: baryon number density [1/fm3]
  * -  T: temperature [MeV]
  * - ye: electron fraction [#/baryon]
  * - ym: muon fraction [#/baryon]
  *
  */
    
  /* Output stream */
  std::string hist_file;
  if (ref_table == "data_DD2_ylmu_last.txt") {
    hist_file = "diff_hist_DD2_ylmu_last.txt";
  } else if (ref_table == "data_DD2.txt") {
    hist_file = "diff_hist_DD2.txt";
  } else {
    std::cout << "ERROR: wrong reference table name!" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  std::ofstream Fout("tests/output/" + hist_file);
  Fout << std::scientific << std::setprecision(8); // scientific format for std output

  /* Output stream -> legend */
  Fout << "# id, nb, temp, Ye, Ymu, Ptot, Stot, Ynue, Yanue, Ynum, Yanum, Ynux, mu_n, mu_p, mu_e, mu_m" << std::endl;

  int count = 0;

  /* Loop over input table */
  for (std::array<double,n_var> tmp : table_data) {
    if (IsValid(&eos, &table_data, count)) {
      print_rel_diff(&eos, &table_data, count, Fout);
    }
    std::cout << "Count = " << count << std::endl;
    count++;
  }

  return 0;
}