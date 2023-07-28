#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <array>

#include "parameters.hpp"
#include "eos_photons.hpp"
#include "eos_neutrinos.hpp"
#include "eos_assembled.hpp"

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
#define corner_check true          // boolean variable for checking corner cases

using namespace parameters;

/* Boolean variable for including of muons */
const bool with_mu = true;

/* Function computing the EOS output */
std::array<double,29> compute_EOS(EOS_assembled* eos, double nb, double T, double *Y) {

    std::array<double,29> EOS_output;

    const double mb   = eos->GetBaryonMass();
    const double d    = 1.0E39*nb*amu_g;
    const double mu_n = eos->NeutronChemicalPotential(nb, T, Y);
    const double mu_p = eos->ProtonChemicalPotential(nb, T, Y);
    const double mu_e = eos->EOS_leptons<0>::LepChemicalPotential<2>(nb, T, Y);
    const double mu_m = eos->EOS_leptons<1>::LepChemicalPotential<2>(nb, T, Y);

    const double mu_nue = mu_p - mu_n + mu_e;
    double mu_num = 0.;
    if (with_mu == true) mu_num = mu_p - mu_n + mu_m;

    const double e = MeV*1.e39*(eos->Energy(nb, T, Y) - mb*nb); // + 1.e39*8.265e18*mb/(c*c)*MeV*nb; //erg/cm3
    const double p = MeV*1.e39*eos->Pressure(nb, T, Y);
    const double s = eos->Entropy(nb, T, Y); // kB=1

    const double ye = Y[0];
    const double ym = Y[1];

    const double ya = eos->AlphaFraction(nb, T, Y);
    const double yh = eos->HeavyFraction(nb, T, Y);
    const double yn = eos->NeutronFraction(nb, T, Y);
    const double yp = eos->ProtonFraction(nb, T, Y);

    const double ynue  = exp(-de_lim/d)*eos->nu_fraction(nb, T, mu_nue/T);
    const double yanue = exp(-de_lim/d)*eos->anu_fraction(nb, T, mu_nue/T);
    const double ynum  = exp(-dx_lim/d)*eos->nu_fraction(nb, T, mu_num/T);
    const double yanum = exp(-dx_lim/d)*eos->anu_fraction(nb, T, mu_num/T);
    const double ynux  = exp(-dx_lim/d)*eos->nu_fraction(nb, T, 0.);

    const double yle = ye + ynue - yanue;
    const double ylm = ym + ynum - yanum;
    
    const double znue  = exp(-de_lim/d)*eos->nu_energy(nb, T, mu_nue/T);
    const double zanue = exp(-de_lim/d)*eos->anu_energy(nb, T, mu_nue/T);
    const double znum  = exp(-dx_lim/d)*eos->nu_energy(nb, T, mu_num/T);
    const double zanum = exp(-dx_lim/d)*eos->anu_energy(nb, T, mu_num/T);
    const double znux  = exp(-dx_lim/d)*eos->nu_energy(nb, T, 0.);

    const double snue  = exp(-de_lim/d)*eos->nu_entropy(nb, T, mu_nue/T);
    const double sanue = exp(-de_lim/d)*eos->anu_entropy(nb, T, mu_nue/T);
    const double snum  = exp(-dx_lim/d)*eos->nu_entropy(nb, T, mu_num/T);
    const double sanum = exp(-dx_lim/d)*eos->anu_entropy(nb, T, mu_num/T);
    const double snux  = exp(-dx_lim/d)*eos->nu_entropy(nb, T, 0.);
        
    const double u = e + 1.e39*nb * MeV * (znue + zanue + znum + zanum + 2. * znux);
    const double ptot = p + 1.e39*nb * MeV * (znue + zanue + znum + zanum + 2. * znux) / 3.;
    const double stot = s + (snue + sanue + snum + sanum + 2. * snux);

    EOS_output[0] = d;          // Mass density [g/cm3]
    EOS_output[1] = yle;        // Electron lepton fraction [#/baryon]
    EOS_output[2] = ylm;        // Muon lepton fraction [#/baryon]
    EOS_output[3] = u;          // Internal energy density (including neutrinos) [erg/cm3]
    EOS_output[4] = nb*1.e39;   // Number density [1/cm3]
    EOS_output[5] = T;          // Temperature [MeV]
    EOS_output[6] = ye;         // Electron fraction [#/baryon]
    EOS_output[7] = ym;         // Muon fraction [#/baryon]
    EOS_output[8] = e;          // Internal energy density (including neutrinos) [erg/cm3]
    EOS_output[9] = ptot;       // Pressure (including neutrinos) [erg/cm3]
    EOS_output[10] = stot;      // Entropy per baryon (including neutrinos) [#/baryon]
    EOS_output[11] = yh;        // Heavy nuclei fraction [#/baryon]
    EOS_output[12] = ya;        // Alpha particle fraction [#/baryon]
    EOS_output[13] = yp;        // Proton fraction [#/baryon]
    EOS_output[14] = yn;        // Neutron fraction [#/baryon]
    EOS_output[15] = ynue;      // Electron neutrino fraction [#/baryon]
    EOS_output[16] = yanue;     // Electron antineutrino fraction [#/baryon]
    EOS_output[17] = ynum;      // Muon neutrino fraction [#/baryon]
    EOS_output[18] = yanum;     // Muon antineutrino fraction [#/baryon]
    EOS_output[19] = ynux;      // Tau (anti)neutrino fraction [#/baryon]
    EOS_output[20] = znue;      // Electron neutrino energy per baryon [MeV/baryon]
    EOS_output[21] = zanue;     // Electron antineutrino energy per baryon [MeV/baryon]
    EOS_output[22] = znum;      // Muon neutrino energy per baryon [MeV/baryon]
    EOS_output[23] = zanum;     // Muon antineutrino energy per baryon [MeV/baryon]
    EOS_output[24] = znux;      // Tau (anti)neutrino energy per baryon [MeV/baryon]
    EOS_output[25] = mu_p - mb; // Proton chemical potential wrt baryon (neutron) mass [MeV]
    EOS_output[26] = mu_n - mb; // Neutron chemical potential wrt baryon (neutron) mass [MeV]
    EOS_output[27] = mu_e;      // Electron chemical potential (with rest mass) [MeV]
    EOS_output[28] = mu_m;      // Muon chemical potential (with rest mass) [MeV]

    return EOS_output;
}

double compare_single_point(std::vector<std::array<double,n_var>>* tab_in,
    std::array<double,29>* eos_arr,
    const int idx_p, const int idx_ref, const int idx_eos,
    const double conv_fact, const double shift,
    const std::string str,  std::ostream& os) {
    const double q_ref = (*tab_in)[idx_p][idx_ref];
    const double q_eos = (*eos_arr)[idx_eos] * conv_fact + shift;

    double q_diff = fabs(q_eos-q_ref) / fabs(q_ref);
    os << str << ": " << q_diff << "\t" << q_ref << "\t" << q_eos << std::endl;
    return q_diff;
}

void print_pressure_spec(EOS_assembled* eos, double nb, double temp, double *Y, std::array<double,29>* eos_arr, std::ostream& os) {
    std::array<double,10> tmp;
    
    tmp[0] = 1.e39 * MeV * eos->BarPressure(nb, temp, Y);                     // Baryon pressure
    tmp[1] = 1.e39 * MeV * eos->EOS_leptons<0>::LepPressure<2>(nb, temp, Y);  // Electron pressure
    tmp[2] = 1.e39 * MeV * eos->EOS_leptons<1>::LepPressure<2>(nb, temp, Y);  // Muon pressure
    tmp[3] = 1.e39 * MeV * eos->RadPressure(temp);                            // Photon pressure
    tmp[4] = 1.e39 * MeV * nb * (*eos_arr)[20] / 3.;  // Electron neutrino pressure
    tmp[5] = 1.e39 * MeV * nb * (*eos_arr)[21] / 3.;  // Electron antineutrino pressure
    tmp[6] = 1.e39 * MeV * nb * (*eos_arr)[22] / 3.;  // Muon neutrino pressure
    tmp[7] = 1.e39 * MeV * nb * (*eos_arr)[23] / 3.;  // Muon antineutrino pressure
    tmp[8] = 1.e39 * MeV * nb * (*eos_arr)[24] / 3.;  // Tau neutrino pressure
    tmp[9] = 1.e39 * MeV * nb * (*eos_arr)[24] / 3.;  // Tau antineutrino pressure
    double Ptot = 0.;
    for (int i=0; i<10; i++) Ptot += tmp[i];

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


void print_entropy_spec(EOS_assembled* eos, double nb, double temp, double *Y, std::array<double,29>* eos_arr, std::ostream& os) {
    const double d    = 1.0E39*nb*amu_g;
    const double mu_n = eos->NeutronChemicalPotential(nb, temp, Y);
    const double mu_p = eos->ProtonChemicalPotential(nb, temp, Y);
    const double mu_e = eos->EOS_leptons<0>::LepChemicalPotential<2>(nb, temp, Y);
    const double mu_m = eos->EOS_leptons<1>::LepChemicalPotential<2>(nb, temp, Y);
    const double mu_nue = mu_p - mu_n + mu_e;
    double mu_num = 0.;
    if (with_mu == true) mu_num = mu_p - mu_n + mu_m;

    std::array<double,10> tmp;
    
    tmp[0] = eos->BarEntropy(nb, temp, Y);                          // Baryon entropy
    tmp[1] = eos->EOS_leptons<0>::LepEntropy<2>(nb, temp, Y) / nb;  // Electron entropy
    tmp[2] = eos->EOS_leptons<1>::LepEntropy<2>(nb, temp, Y) / nb;  // Muon entropy
    tmp[3] = eos->RadEntropy(temp) / nb;                            // Photon entropy
    tmp[4] = exp(-de_lim/d)*eos->nu_entropy(nb, temp, mu_nue/temp);      // Electron neutrino entropy
    tmp[5] = exp(-de_lim/d)*eos->anu_entropy(nb, temp, mu_nue/temp);     // Electron antineutrino entropy
    tmp[6] = exp(-dx_lim/d)*eos->nu_entropy(nb, temp, mu_num/temp);      // Muon neutrino entropy
    tmp[7] = exp(-dx_lim/d)*eos->anu_entropy(nb, temp, mu_num/temp);     // Muon antineutrino entropy
    tmp[8] = exp(-dx_lim/d)*eos->nu_energy(nb, temp, 0.);                // Tau neutrino entropy
    tmp[9] = exp(-dx_lim/d)*eos->nu_energy(nb, temp, 0.);                // Tau antineutrino entropy
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

void compare_EOS(EOS_assembled* eos, std::vector<std::array<double,n_var>>* tab_in, const int idx_p, const int n, std::ostream& os) {
    double nb = (*tab_in)[idx_p][0] * 1.0E-39 / amu_g;
    double  T = (*tab_in)[idx_p][11];
    double ye = (*tab_in)[idx_p][5] - (*tab_in)[idx_p][6];
    double ym = (*tab_in)[idx_p][7] - (*tab_in)[idx_p][8];
    //double yq = ye + ym;
    double Y[2] = {0.0, 0.0};
    Y[0] = ye;
    Y[1] = ym;

    const double  MBar = eos->GetBaryonMass();

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
    std::array<double,29> eos_array = compute_EOS(eos, nb, T, Y);

    std::array<double,11> tmp; // array to store relative difference values

    os << "                            ";    
    os << "  Relative diff   ";    
    os << "    Ref value     ";
    os << "  EOS output      ";
    os << std::endl;    
    os << std::endl;
    
    // Compare EOS against reference table
    tmp[0]  = compare_single_point(tab_in, &eos_array, idx_p,  9,  9, 1.,       0., "Total Pressure              ", os);
    tmp[1]  = compare_single_point(tab_in, &eos_array, idx_p, 10, 10, 1.,       0., "Total Entropy               ", os);
    tmp[2]  = compare_single_point(tab_in, &eos_array, idx_p, 12, 15, 1.,       0., "nu_e fraction               ", os);
    tmp[3]  = compare_single_point(tab_in, &eos_array, idx_p, 13, 16, 1.,       0., "anu_e fraction              ", os);
    tmp[4]  = compare_single_point(tab_in, &eos_array, idx_p, 14, 17, 1.,       0., "nu_mu  fraction             ", os);
    tmp[5]  = compare_single_point(tab_in, &eos_array, idx_p, 15, 18, 1.,       0., "anu_mu fraction             ", os);
    tmp[6]  = compare_single_point(tab_in, &eos_array, idx_p, 16, 19, 1.,       0., "nu_tau fraction             ", os);
    tmp[7]  = compare_single_point(tab_in, &eos_array, idx_p, 17, 26, 1., MBar-m_n, "Neutrons  NR chem potential ", os);
    tmp[8]  = compare_single_point(tab_in, &eos_array, idx_p, 18, 25, 1., MBar-m_p, "Protons   NR chem potential ", os);
    tmp[9]  = compare_single_point(tab_in, &eos_array, idx_p, 19, 27, 1.,     -m_e, "Electrons NR chem potential ", os);
    tmp[10] = compare_single_point(tab_in, &eos_array, idx_p, 20, 28, 1.,    -m_mu, "Muons     NR chem potential ", os);

    if (tmp[0] > 1.0E-02) print_pressure_spec(eos, nb, T, Y, &eos_array, os);
    if (tmp[1] > 1.0E-02) print_entropy_spec(eos, nb, T, Y, &eos_array,  os);

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
    /* Global EOS class */
    EOS_assembled eos;

    /* Read EOS tables */
    eos.ReadBarTableFromFile("eos_table/baryons/DD2_bar.h5");
    eos.EOS_leptons<0>::m_lep_active = true;
    if (with_mu == true) eos.EOS_leptons<1>::m_lep_active = true;

    /* Input stream */
    std::string table_name = "data_DD2_ylmu_last.txt"; // initial Yl_mu = Ym_cold
    //std::string table_name = "data_DD2.txt";           // initial Yl_mu = 0.
    std::ifstream fin("eos_table/global/"+table_name);
    
    /* Output stream */
	std::ofstream Rout("eos_table/random_comp_ele.txt");

    PrintHeader(&eos, table_name, std::cout);
    PrintHeader(&eos, table_name, Rout);

    //* Number of random EOS points picked for the comparison */
    const int n_rand = 10;


    /* Input vector */
    std::array<double,n_var> tmp;
    std::vector<std::array<double,n_var>> ref_table;

    std::string dummyLine;
    
    /* Read data from input table */
    int np = 0;
    std::cout << "Reading data from the input table!" << std::endl;
    while (!fin.fail()) {
        std::getline(fin, dummyLine);
        std::stringstream ss(dummyLine);
        
        for (int i=0; i<n_var; i++) ss >> tmp[i];
        ref_table.push_back(tmp);
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
    for (int i = 0; i < n_rand; i++) {
      int random = rand() % np; // Random number between 0 and np

      // std::cout << "Random number: " << random << std::endl;
      if (IsValid(&eos, &ref_table, random)) {
        compare_EOS(&eos, &ref_table, random, i + 1, Rout);
      } else {
        i -= 1;
      }
    }

    if (not corner_check) return 0;

    /* Check EOS in the corner cases */

    /* Output stream */
	std::ofstream Cout("eos_table/corner_comp_ele.txt");

    PrintHeader(&eos, table_name, Cout);

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
    for (std::array<double,n_var> tmp : ref_table) {
      if (IsValid(&eos, &ref_table, count)) {
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
   
    compare_EOS(&eos, &ref_table, id_nb[0], 1, Cout);
    compare_EOS(&eos, &ref_table, id_nb[1], 2, Cout);
    compare_EOS(&eos, &ref_table,  id_T[0], 3, Cout);
    compare_EOS(&eos, &ref_table,  id_T[1], 4, Cout);
    compare_EOS(&eos, &ref_table, id_ye[0], 5, Cout);
    compare_EOS(&eos, &ref_table, id_ye[1], 6, Cout);
    compare_EOS(&eos, &ref_table, id_ym[0], 7, Cout);
    compare_EOS(&eos, &ref_table, id_ym[1], 8, Cout);
    return 0;
}
