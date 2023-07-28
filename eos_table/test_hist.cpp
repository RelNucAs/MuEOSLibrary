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
    const double conv_fact, const double shift) {
    const double q_ref = (*tab_in)[idx_p][idx_ref];
    const double q_eos = (*eos_arr)[idx_eos] * conv_fact + shift;

    return fabs(q_eos-q_ref) / fabs(q_ref);
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

void compare_EOS(EOS_assembled* eos, std::vector<std::array<double,n_var>>* tab_in, const int idx_p, std::ostream& os) {
    double nb = (*tab_in)[idx_p][0] * 1.0E-39 / amu_g;
    double  T = (*tab_in)[idx_p][11];
    double ye = (*tab_in)[idx_p][5] - (*tab_in)[idx_p][6];
    double ym = (*tab_in)[idx_p][7] - (*tab_in)[idx_p][8];
    //double yq = ye + ym;
    double Y[2] = {0.0, 0.0};
    Y[0] = ye;
    Y[1] = ym;

    const double  MBar = eos->GetBaryonMass();

    // Print input quantities
    os << idx_p << "   ";
    os << nb    << "   ";
    os << T     << "   ";
    os << ye    << "   ";
    os << ym    << "   ";
    
    // Compute EOS
    std::array<double,29> eos_array = compute_EOS(eos, nb, T, Y);
   
    // Compare EOS against reference table
    os << compare_single_point(tab_in, &eos_array, idx_p,  9,  9, 1.,       0.) << "   "; // Total Pressure              
    os << compare_single_point(tab_in, &eos_array, idx_p, 10, 10, 1.,       0.) << "   "; // Total Entropy               
    os << compare_single_point(tab_in, &eos_array, idx_p, 12, 15, 1.,       0.) << "   "; // nu_e fraction               
    os << compare_single_point(tab_in, &eos_array, idx_p, 13, 16, 1.,       0.) << "   "; // anu_e fraction              
    os << compare_single_point(tab_in, &eos_array, idx_p, 14, 17, 1.,       0.) << "   "; // nu_mu  fraction             
    os << compare_single_point(tab_in, &eos_array, idx_p, 15, 18, 1.,       0.) << "   "; // anu_mu fraction             
    os << compare_single_point(tab_in, &eos_array, idx_p, 16, 19, 1.,       0.) << "   "; // nu_tau fraction             
    os << compare_single_point(tab_in, &eos_array, idx_p, 17, 26, 1., MBar-m_n) << "   "; // Neutrons  NR chem potential 
    os << compare_single_point(tab_in, &eos_array, idx_p, 18, 25, 1., MBar-m_p) << "   "; // Protons   NR chem potential 
    os << compare_single_point(tab_in, &eos_array, idx_p, 19, 27, 1.,     -m_e) << "   "; // Electrons NR chem potential
    os << compare_single_point(tab_in, &eos_array, idx_p, 20, 28, 1.,    -m_mu) << "   "; // Muons     NR chem potential 
    os << std::endl;

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
    //std::string table_name = "data_DD2_ylmu_last.txt"; // initial Yl_mu = Ym_cold
    std::string table_name = "data_DD2.txt";           // initial Yl_mu = 0.
    std::ifstream fin("eos_table/global/"+table_name);
    
    PrintHeader(&eos, table_name, std::cout);
    
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


    /* Output stream */
	//std::ofstream Fout("eos_table/diff_hist_DD2_ylmu_last.txt");	
	std::ofstream Fout("eos_table/diff_hist_DD2.txt");

    Fout << std::scientific << std::setprecision(8); // scientific format for std output

    /* Output stream -> legend */
	Fout << "# id, nb, temp, Ye, Ymu, Ptot, Stot, Ynue, Yanue, Ynum, Yanum, Ynux, mu_n, mu_p, mu_e, mu_m" << std::endl;

    int count = 0;

    /* Loop over input table */
    for (std::array<double,n_var> tmp : ref_table) {
      if (IsValid(&eos, &ref_table, count)) {
        compare_EOS(&eos, &ref_table, count, Fout);
      }
      std::cout << "Count = " << count << std::endl;
      count++;
    }
   
    return 0;
}
