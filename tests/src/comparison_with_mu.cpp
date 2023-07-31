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

double print_single_diff(std::vector<std::array<double,n_var>>* tab_in,
    const int idx_p, const int idx_ref, const double q_eos,
    const std::string str,  std::ostream& os) {
    
    const double q_ref = (*tab_in)[idx_p][idx_ref];
    const double q_diff = fabs(q_eos-q_ref) / fabs(q_ref);

    os << str << ": " << q_diff << "\t" << q_ref << "\t" << q_eos << std::endl;
    return q_diff;
}

void print_pressure_spec(EOS_assembled* eos, double nb, double temp, double *Y, FullEOSOutput *eos_out, std::ostream& os) {
    std::array<double,10> tmp;
    
    tmp[0] = eos->BarPressure(nb, temp, Y);                     // Baryon pressure
    tmp[1] = eos->EOS_leptons<0>::LepPressure<2>(nb, temp, Y);  // Electron pressure
    tmp[2] = eos->EOS_leptons<1>::LepPressure<2>(nb, temp, Y);  // Muon pressure
    tmp[3] = eos->RadPressure(temp);                            // Photon pressure
    tmp[4] = nb * eos_out->nuEOS.Z_nue  / 3.;  // Electron neutrino pressure
    tmp[5] = nb * eos_out->nuEOS.Z_anue / 3.;  // Electron antineutrino pressure
    tmp[6] = nb * eos_out->nuEOS.Z_num  / 3.;  // Muon neutrino pressure
    tmp[7] = nb * eos_out->nuEOS.Z_anum / 3.;  // Muon antineutrino pressure
    tmp[8] = nb * eos_out->nuEOS.Z_nux  / 3.;  // Tau neutrino pressure
    tmp[9] = nb * eos_out->nuEOS.Z_nux  / 3.;  // Tau antineutrino pressure
    
    double Ptot = 0.;
    for (int i=0; i<10; i++) {
        tmp[i] = 1.0E+39 * MeV * tmp[i]; // convert from MeV/fm3 to erg/cm3
        Ptot += tmp[i];
    }

    os << std::endl;
    os << "Pressure contribution from the various species:" << std::endl;
    os << "Baryons        : " << tmp[0] << " erg/cm3" << std::endl;
    os << "Electrons      : " << tmp[1] << " erg/cm3" << std::endl;
    os << "Muons          : " << tmp[2] << " erg/cm3" << std::endl;
    os << "Photons        : " << tmp[3] << " erg/cm3" << std::endl;
    os << "Electronic nu  : " << tmp[4] << " erg/cm3" << std::endl;
    os << "Electronic anu : " << tmp[5] << " erg/cm3" << std::endl;
    os << "Electronic nu  : " << tmp[6] << " erg/cm3" << std::endl;
    os << "Muonic anu     : " << tmp[7] << " erg/cm3" << std::endl;
    os << "Tauonic (a)nu  : " << tmp[8] << " erg/cm3" << std::endl;
    os << std::endl;
    os << "Total          : " <<  Ptot  << " erg/cm3" << std::endl;
    os << std::endl;

    return;
}


void print_entropy_spec(EOS_assembled* eos, double nb, double temp, double *Y, FullEOSOutput *eos_out, std::ostream& os) {
    std::array<double,10> tmp;
    
    tmp[0] = eos->BarEntropy(nb, temp, Y);                          // Baryon entropy
    tmp[1] = eos->EOS_leptons<0>::LepEntropy<2>(nb, temp, Y) / nb;  // Electron entropy
    tmp[2] = eos->EOS_leptons<1>::LepEntropy<2>(nb, temp, Y) / nb;  // Muon entropy
    tmp[3] = eos->RadEntropy(temp) / nb;                            // Photon entropy
    tmp[4] = eos_out->nuEOS.s_nue;                                  // Electron neutrino entropy
    tmp[5] = eos_out->nuEOS.s_anue;                                 // Electron antineutrino entropy
    tmp[6] = eos_out->nuEOS.s_num;                                  // Muon neutrino entropy
    tmp[7] = eos_out->nuEOS.s_anum;                                 // Muon antineutrino entropy
    tmp[8] = eos_out->nuEOS.s_nux;                                  // Tau neutrino entropy
    tmp[9] = eos_out->nuEOS.s_nux;                                  // Tau antineutrino entropy
    
    double stot = 0.;
    for (int i=0; i<10; i++) stot += tmp[i];

    os << std::endl;
    os << "Specific entropy contribution from the various species:" << std::endl;
    os << "Baryons        : " << tmp[0] << " kB/baryon" << std::endl;
    os << "Electrons      : " << tmp[1] << " kB/baryon" << std::endl;
    os << "Muons          : " << tmp[2] << " kB/baryon" << std::endl;
    os << "Photons        : " << tmp[3] << " kB/baryon" << std::endl;
    os << "Electronic nu  : " << tmp[4] << " kB/baryon" << std::endl;
    os << "Electronic anu : " << tmp[5] << " kB/baryon" << std::endl;
    os << "Electronic nu  : " << tmp[6] << " kB/baryon" << std::endl;
    os << "Muonic anu     : " << tmp[7] << " kB/baryon" << std::endl;
    os << "Tauonic (a)nu  : " << tmp[8] << " kB/baryon" << std::endl;
    os << std::endl;
    os << "Total          : " <<  stot  << " kB/baryon" << std::endl;
    os << std::endl;

    return;
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

void compare_EOS_single(EOS_assembled* eos, std::vector<std::array<double,n_var>>* tab_in, const int idx_p, const int n, std::ostream& os) {
    double nb = (*tab_in)[idx_p][0] * 1.0E-39 / amu_g;
    double  T = (*tab_in)[idx_p][11];
    double ye = (*tab_in)[idx_p][5] - (*tab_in)[idx_p][6];
    double ym = (*tab_in)[idx_p][7] - (*tab_in)[idx_p][8];
    //double yq = ye + ym;
    double Y[2] = {0.0, 0.0};
    Y[0] = ye;
    Y[1] = ym;

    os << "===================================================" << std::endl;

    os << std::endl;
    os << "###############################" << std::endl;
    os << "#    EOS point number " << std::setw(2) << n << "      #" << std::endl;
    os << "###############################" << std::endl;
    os << std::endl;

    // Print input quantities
    os << std::endl;
    os << "Index of EOS point: " << idx_p << std::endl;
    os << "Input values: " << std::endl;
    os << "nb = " <<   nb << " " << "fm-3" << std::endl;
    os << "T  = " <<    T << " " <<  "MeV" << std::endl;
    os << "Ye = " <<   ye                  << std::endl;
    os << "Ym = " <<   ym                  << std::endl;
    os << std::endl;
    
    // Compute EOS
    FullEOSOutput eos_out = compute_EOS_nu_thres(eos, nb, T, Y);

    const double ptot = eos_out.P + 1.0E+39 * nb * MeV * eos_out.nuEOS.Z_tot / 3.;
    const double stot = eos_out.s + eos_out.nuEOS.s_tot;

    std::array<double,11> tmp; // array to store relative difference values

    os << "                            ";    
    os << "  Relative diff   ";    
    os << "    Ref value     ";
    os << "  EOS output      ";
    os << std::endl;    
    os << std::endl;
    
    // Compare EOS against reference table
    tmp[0]  = print_single_diff(tab_in, idx_p,  9, ptot                       , "Total Pressure              ", os);
    tmp[1]  = print_single_diff(tab_in, idx_p, 10, stot                       , "Total Entropy               ", os);
    tmp[2]  = print_single_diff(tab_in, idx_p, 12, eos_out.nuEOS.Y_nu.ynue    , "nu_e fraction               ", os);
    tmp[3]  = print_single_diff(tab_in, idx_p, 13, eos_out.nuEOS.Y_nu.yanue   , "anu_e fraction              ", os);
    tmp[4]  = print_single_diff(tab_in, idx_p, 14, eos_out.nuEOS.Y_nu.ynum    , "nu_mu  fraction             ", os);
    tmp[5]  = print_single_diff(tab_in, idx_p, 15, eos_out.nuEOS.Y_nu.yanum   , "anu_mu fraction             ", os);
    tmp[6]  = print_single_diff(tab_in, idx_p, 16, eos_out.nuEOS.Y_nu.ynux    , "nu_tau fraction             ", os);
    tmp[7]  = print_single_diff(tab_in, idx_p, 17, eos_out.chem_pot.mu_n - m_n, "Neutrons  NR chem potential ", os);
    tmp[8]  = print_single_diff(tab_in, idx_p, 18, eos_out.chem_pot.mu_p - m_p, "Protons   NR chem potential ", os);
    tmp[9]  = print_single_diff(tab_in, idx_p, 19, eos_out.chem_pot.mu_e - m_e, "Electrons NR chem potential ", os);
    tmp[10] = print_single_diff(tab_in, idx_p, 20, eos_out.chem_pot.mu_m -m_mu, "Muons     NR chem potential ", os);

    if (tmp[0] > 1.0E-02) print_pressure_spec(eos, nb, T, Y, &eos_out, os);
    if (tmp[1] > 1.0E-02) print_entropy_spec(eos , nb, T, Y, &eos_out,  os);

    double avg = 0., max = 0.;

    for (int j=0; j<11; j++) {
        avg += tmp[j] / 11.;
        max = std::max(max,tmp[j]);
    }
        
    os << std::endl;
    os << "Average difference: " << avg << std::endl;
    os << "Maximum difference: " << max << std::endl;
    os << std::endl;

    os << "===================================================" << std::endl;

    return;
}

void FindCornerCases(double var, double *lim, int *idx, int count) {
    if (var < lim[0]) {
        lim[0] = var;
        idx[0] = count;
    }
    if (var > lim[1]) {
        lim[1] = var;
        idx[1] = count; 
    }
    return;
}

void PrintHeader(EOS_assembled* eos, std::string table_name, std::ostream& os) {
    os << std::scientific << std::setprecision(8); // scientific format for std output
        
    os << "#####################" << std::endl;
    os << "#   EOS with muons  #" << std::endl;
    os << "#####################" << std::endl;
    os << std::endl;
   
    /* Read Baryon EOS table input arrays */
    const double * nb_bar = eos->GetRawLogNumberDensity();
    const double *  T_bar = eos->GetRawLogTemperature();
    const double * yq_bar = eos->GetRawYq();

    os << "Comparing global EOS with muons with Eleonora's table: " << table_name << std::endl;
    os << std::endl;
    os << "Points must lie within the range of validity of the EOS:"     << std::endl;
    os << " - Number density    : [" << exp(nb_bar[0]) << ", " << exp(nb_bar[eos->m_nn-1]) << "]" << " " << "fm-3" << std::endl;
    os << " - Temperature       : [" << exp( T_bar[0]) << ", " << exp( T_bar[eos->m_nt-1]) << "]" << " " <<  "MeV" << std::endl;
    os << " - Charge fraction   : [" <<      yq_bar[0] << ", " <<     yq_bar[eos->m_ny-1]  << "]" << std::endl;
    os << " - Muon fraction     : [" <<         ym_min << ", " <<                  ym_max  << "]" << std::endl;
    os << " - Electron fraction : positive"                                                       << std::endl;
    os << std::endl;
    os << "Number density must be also > " << nb_thr << " fm-3 " << "to justify neutrino trapped" << std::endl;
    os << std::endl;

    return;
}

int main () {
    /* Name of baryon EOS table */
    std::string BarTableName = "DD2_bar.h5";  // baryon table

    /* Initialize global EOS class

    Constructor -> EOS_assembled(const int id_eos, const bool el_bool, const bool mu_bool, std::string BarTableName)

    Inputs:
     - id_EOS: method for EOS computation (1: interpolation, 2: on-the-fly)
     - el_bool: flag for activating electrons
     - mu_bool: flag for activating muons
     - BarTableName: path of baryon EOS table  */
    EOS_assembled eos(2, true, true, BarTableName);

    /* Input stream */
    std::string ref_table = "data_DD2_ylmu_last.txt"; // initial Yl_mu = Ym_cold
    //std::string ref_table = "data_DD2.txt";           // initial Yl_mu = 0.
    std::ifstream fin("tests/data/" + ref_table);
    
    /* Output stream */
	  std::ofstream Rout("tests/output/random_comp_with_mu.txt");

    PrintHeader(&eos, ref_table, std::cout);
    PrintHeader(&eos, ref_table, Rout);

    /* Number of random EOS points picked for the comparison */
    const int n_rand = 10;


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

    Rout << std::endl;
    Rout << "Pick " << n_rand << " random points from Eleonora's table"    << std::endl;
    Rout << std::endl;

    /* Providing a seed value */
    srand((unsigned)time(NULL));

    /* Loop to get n_rand random numbers */
    for (int i=0; i<n_rand; i++) {
      int random = rand() % np; // Random number between 0 and np

      // std::cout << "Random number: " << random << std::endl;
      if (IsValid(&eos, &table_data, random)) {
        compare_EOS_single(&eos, &table_data, random, i + 1, Rout);
      } else {
        i -= 1;
      }
    }

    if (not corner_check) return 0;

    /* Check EOS in the corner cases */

    /* Output stream */
	  std::ofstream Cout("tests/output/corner_comp_with_mu.txt");
    PrintHeader(&eos, ref_table, Cout);

    Cout << std::endl;
    Cout << "Check corner cases:" << std::endl;
    Cout << std::endl;

    double nb_lim[2] = {+numeric_limits<double>::max(), -numeric_limits<double>::max()};
    double  T_lim[2] = {+numeric_limits<double>::max(), -numeric_limits<double>::max()};
    double ye_lim[2] = {+numeric_limits<double>::max(), -numeric_limits<double>::max()};
    double ym_lim[2] = {+numeric_limits<double>::max(), -numeric_limits<double>::max()};
    
    int id_nb[2] = {0};
    int id_T[2]  = {0};
    int id_ye[2] = {0};
    int id_ym[2] = {0};
    
    int count = 0;

    /* Find corner cases */
    for (std::array<double,n_var> tmp : table_data) {
      if (IsValid(&eos, &table_data, count)) {
        double nb   = tmp[0] * 1.0E-39 / amu_g;
        double T    = tmp[11];
        double ye   = tmp[5] - tmp[6];
        double ym   = tmp[7] - tmp[8];

        FindCornerCases(nb, nb_lim, id_nb, count);
        FindCornerCases(T ,  T_lim, id_T , count);
        FindCornerCases(ye, ye_lim, id_ye, count);
        FindCornerCases(ym, ym_lim, id_ym, count);
      }
      count++;
    }
   
    compare_EOS_single(&eos, &table_data, id_nb[0], 1, Cout);
    compare_EOS_single(&eos, &table_data, id_nb[1], 2, Cout);
    compare_EOS_single(&eos, &table_data,  id_T[0], 3, Cout);
    compare_EOS_single(&eos, &table_data,  id_T[1], 4, Cout);
    compare_EOS_single(&eos, &table_data, id_ye[0], 5, Cout);
    compare_EOS_single(&eos, &table_data, id_ye[1], 6, Cout);
    compare_EOS_single(&eos, &table_data, id_ym[0], 7, Cout);
    compare_EOS_single(&eos, &table_data, id_ym[1], 8, Cout);
    return 0;
}
