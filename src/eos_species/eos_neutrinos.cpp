/* NEUTRINO EOS AT EQUILIBRIUM WITH RADIATION (neutrinos degeneracy parameter is equal to minus anti-neutrinos degeneracy parameter) */

#include "eos_species/eos_species.hpp"
#include "constants.hpp"
#include "fermi_integrals/fermi_integrals.hpp"

/* Input parameters:
 *    - number density                nb [fm-3]
 *    - temperature                   t  [MeV]
 *    - neutrino degeneracy parameter eta
*/

constexpr double fdi_p2_0 = 1.803085354739391; // Fermi_integral_p2(0) = 1.5 * zeta(3.)
constexpr double fdi_p3_0 = 5.682196976983475; // Fermi_integral_p3(0) = 7. * MEOS_pi**4. / 60.

constexpr double twice_fdi_p3_0 = 2. * fdi_p3_0; // 7. * MEOS_pi**4. / 60.

void EOS_neutrinos::NeutrinoEOS(const double nb, const double T, const double* chem_pot) {
 
  const double mu_nue = chem_pot[0] - chem_pot[1] + chem_pot[2]; // Electron neutrino chemical potential [MeV]
  double mu_num = 0.;
  if (m_mu_active) mu_num = chem_pot[0] - chem_pot[1] + chem_pot[3];// Muon neutrino chemical potential [MeV]

  const double eta_nue = mu_nue / T; //mu_nue/T
  const double eta_num = mu_num / T; //mu_num/T

  eta_nu[0] =  eta_nue;
  eta_nu[1] = -eta_nue;
  eta_nu[2] =  eta_num;
  eta_nu[3] = -eta_num;

  const double eta_nue_sq = POW2(eta_nue);
  const double eta_num_sq = POW2(eta_num);


  const double ynu_factor = MEOS_fourpi_hc3 * POW3(T) / nb;
  const double znu_factor = ynu_factor * T;

  const double fdi_p2_nue = Fermi_integral_p2(eta_nue);
  const double fdi_p2_num = Fermi_integral_p2(eta_num);
  const double fdi_p3_nue = Fermi_integral_p3(eta_nue);
  const double fdi_p3_num = Fermi_integral_p3(eta_num);

  const double fdi_p2_anue = fdi_p2_nue - eta_nue * (MEOS_pi_sq + eta_nue_sq) / 3.;
  const double fdi_p2_anum = fdi_p2_num - eta_num * (MEOS_pi_sq + POW2(eta_num)) / 3.;
  const double fdi_p3_anue = twice_fdi_p3_0 + 0.5 * eta_nue_sq * (MEOS_pi_sq + 0.5 * eta_nue_sq) - fdi_p3_nue;
  const double fdi_p3_anum = twice_fdi_p3_0 + 0.5 * eta_num_sq * (MEOS_pi_sq + 0.5 * eta_num_sq) - fdi_p3_num;

  // Number fractions
  Y_nu[0] = ynu_factor * fdi_p2_nue;  //Y_nue
  Y_nu[1] = ynu_factor * fdi_p2_anue;   //Y_anue
  Y_nu[2] = ynu_factor * fdi_p2_num;   //Y_num
  Y_nu[3] = ynu_factor * fdi_p2_anum;   //Y_anum
  Y_nu[4] = ynu_factor * fdi_p2_0;        //Y_nux

// Energy per baryon [MeV/baryon]
  Z_nu[0] = znu_factor * fdi_p3_nue;   // Z_nue
  Z_nu[1] = znu_factor * fdi_p3_anue;   // Z_anue
  Z_nu[2] = znu_factor * fdi_p3_num;    // Z_num
  Z_nu[3] = znu_factor * fdi_p3_anum;   // Z_anum
  Z_nu[5] = znu_factor * fdi_p3_0;         // Z_nux

  // Pressure per baryon [MeV/baryon]
  P_nu[0] = Z_nu[0] / 3.; // P_nue
  P_nu[1] = Z_nu[1] / 3.; // P_nue
  P_nu[2] = Z_nu[2] / 3.; // P_nue
  P_nu[3] = Z_nu[3] / 3.; // P_nue
  P_nu[4] = Z_nu[4] / 3.; // P_nue

  // Entropy per baryon [1/baryon (kB=1)]
  s_nu[0] = ynu_factor * (1.3333333333333333 * fdi_p3_nue - eta_nue * fdi_p2_nue);   // s_nue
  s_nu[1] = ynu_factor * (1.3333333333333333 * fdi_p3_anue + eta_nue * fdi_p2_anue);   // s_anue
  s_nu[2] = ynu_factor * (1.3333333333333333 * fdi_p3_num - eta_num * fdi_p2_num);    // s_num
  s_nu[3] = ynu_factor * (1.3333333333333333 * fdi_p3_anum + eta_num * fdi_p2_anum);   // s_anum
  s_nu[4] = ynu_factor * (1.3333333333333333 * fdi_p3_0);         // s_nux

  return;
}

// Number fraction of neutrinos
double EOS_neutrinos::GetNuNumberFraction(const int idx) {
  return Y_nu[idx];
}

// Energy per baryon of neutrinos
double EOS_neutrinos::GetNuEnergyFraction(const int idx) {
  return Z_nu[idx]; // MeV/baryon
}

// Pressure per baryon of neutrinos
double EOS_neutrinos::GetNuPressure(const int idx) {
  return P_nu[idx]; // MeV/baryon
}

// Entropy per baryon of electron neutrinos
double EOS_neutrinos::GetNuEntropy(const int idx) {
  return s_nu[idx];
}

// Degeneracy parameter of neutrinos
double EOS_neutrinos::GetNuDegeneracyParameter(const int idx) {
	return eta_nu[idx]; // 1/baryon (kB=1)
}
