#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "eos_species/eos_species.hpp"

std::string EFullTableName = "eos_table/electrons/eos_electrons_table.txt"; // electron full table
std::string MFullTableName = "eos_table/muons/eos_muons_table.txt"; // muon full table

std::string EEtaTableName = "eos_table/electrons/eos_electrons_eta.txt"; // electron full table
std::string MEtaTableName = "eos_table/muons/eos_muons_eta.txt"; // muon full table

/*
  Inputs:
   - id_EOS: method for EOS computation (0: full interpolation, 1: eta interpolation, 2: fully on-the-fly)
   - mu_flag: flag for activating muons
   - BarTableName: path of baryon EOS table
*/
MuEOSClass::MuEOSClass(const int id_eos, const bool mu_flag, std::string BarTableName) {
  ReadBarTableFromFile("eos_table/baryons/" + BarTableName);
  EOS_leptons<0>::m_lep_active = true;
  EOS_leptons<1>::m_lep_active = mu_flag;
  if (id_eos == 0) {
    EOS_leptons<0>::ReadFullLepTableFromFile(EFullTableName);
    if (mu_flag == true) EOS_leptons<1>::ReadFullLepTableFromFile(MFullTableName);
  } else if (id_eos == 1) {
    EOS_leptons<0>::ReadEtaLepTableFromFile(EEtaTableName);
    if (mu_flag == true) EOS_leptons<1>::ReadEtaLepTableFromFile(MEtaTableName);
  } else if (id_eos != 2) {
    std::cout << "ERROR: wrong index for EOS computation in InitalizeEOS" << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

MuEOSClass::~MuEOSClass() {};


/* double MuEOSClass::Energy(double n, double T, double *Y) {
  return BarEnergy(n, T, Y) + EOS_leptons<0>::LepEnergy<id_test>(n, T, Y) + EOS_leptons<1>::LepEnergy<id_test>(n, T, Y) + RadEnergy(T);
}

double MuEOSClass::Pressure(double n, double T, double *Y) {
  return BarPressure(n, T, Y) + EOS_leptons<0>::LepPressure<id_test>(n, T, Y) + EOS_leptons<1>::LepPressure<id_test>(n, T, Y) + RadPressure(T);
}

double MuEOSClass::Entropy(double n, double T, double *Y) {
  return BarEntropy(n, T, Y) + (EOS_leptons<0>::LepEntropy<id_test>(n, T, Y) + EOS_leptons<1>::LepEntropy<id_test>(n, T, Y) + RadEntropy(T)) / n;
}

double MuEOSClass::Enthalpy(double n, double T, double *Y) {
  double const P = Pressure(n, T, Y);
  double const e = Energy(n, T, Y);
  return (P + e)/n;
} */

/* double MuEOSClass::SoundSpeed(double n, double T, double *Y) {

  double const dPdn = BardPdn(n, T, Y) + EOS_leptons<0>::LdPdn<id_test>(n, T, Y) + EOS_leptons<1>::LdPdn<id_test>(n, T, Y);
  double const dsdn = Bardsdn(n, T, Y) + EOS_leptons<0>::Ldsdn<id_test>(n, T, Y) + EOS_leptons<1>::Ldsdn<id_test>(n, T, Y) - (RadEntropy(T) / (n*n));
  double const dPdt = BardPdT(n, T, Y) + EOS_leptons<0>::LdPdt<id_test>(n, T, Y) + EOS_leptons<1>::LdPdt<id_test>(n, T, Y) +  RaddPdT(T);
  double const dsdt = BardsdT(n, T, Y) + EOS_leptons<0>::Ldsdt<id_test>(n, T, Y) + EOS_leptons<1>::Ldsdt<id_test>(n, T, Y) + (RaddsdT(T)/n);

  double const cs2 = sqrt(std::max(1.e-6,(dPdn - dsdn/dsdt*dPdt) / Enthalpy(n, T, Y)));
  //if (cs2 >= 1.) cout << "cs2 > 1!" << endl;
  //cout << "dPdn = " << dPdn << "\t" << "dsdn = " << dsdn << "\t" << "dPdt = " << dPdt << "\t" << "dsdt = " << dsdt << endl << endl; 
  return cs2;
} */

/* ChemPotentials MuEOSClass::GetChemicalPotentials(double n, double T, double *Y) {
  ChemPotentials out;
  out.mu_n = NeutronChemicalPotential(n, T, Y);                                  // Neutron chemical potential [MeV]
  out.mu_p = ProtonChemicalPotential(n, T, Y);                                   // Proton chemical potential [MeV]
  out.mu_e = EOS_leptons<0>::LepChemicalPotential<id_test>(n, T, Y);             // Electron chemical potential [MeV]
  out.mu_m = EOS_leptons<1>::LepChemicalPotential<id_test>(n, T, Y);             // Muon chemical potential [MeV]
  if (EOS_leptons<0>::m_lep_active == true) {
    out.mu_nue = out.mu_p - out.mu_n + out.mu_e; // Electron neutrino chemical potential [MeV]
  } else {
    out.mu_nue = 0.;
  }
  if (EOS_leptons<1>::m_lep_active == true) {
    out.mu_num = out.mu_p - out.mu_n + out.mu_m; // Muon neutrino chemical potential [MeV]
  } else {
    out.mu_num = 0.;
  }
  return out;
}

ParticleFractions MuEOSClass::GetParticleFractions(double n, double T, double *Y) {
  ParticleFractions out;

  double *YY = Y;
  if (EOS_leptons<0>::m_lep_active == true) {
    out.ye = Y[0]; // Electron fraction [#/baryon]
  } else {
    YY[0] = 0.;
    out.ye = 0.;
  }
  if (EOS_leptons<1>::m_lep_active == true) {
    out.ym = Y[1]; // Muon fraction [#/baryon]
  } else {
    YY[1] = 0.;
    out.ym = 0.;
  }

  // Baryon fractions
  ParticleFractions Y_bar = BaryonFractions(n, T, YY);
  out.yh = Y_bar.yh;    // Heavy nuclei fraction [#/baryon]
  out.ya = Y_bar.ya;    // Alpha particle fraction [#/baryon]
  out.yn = Y_bar.yn;    // Neutron fraction [#/baryon]
  out.yp = Y_bar.yp;    // Proton fraction [#/baryon]
  
  ChemPotentials chem_pot = GetChemicalPotentials(n, T, Y);

  // Neutrino fractions
  ParticleFractions Y_nu = NeutrinoFractions(n, T, Y, &chem_pot);
  out.ynue  = Y_nu.ynue;  // Electron neutrino fraction [#/baryon]
  out.yanue = Y_nu.yanue; // Electron antineutrino fraction [#/baryon]
  out.ynum  = Y_nu.ynum;  // Muon neutrino fraction [#/baryon]
  out.yanum = Y_nu.yanum; // Muon antineutrino fraction [#/baryon]
  out.ynux  = Y_nu.ynux;  // Tau (anti)neutrino fraction [#/baryon]
  
  out.yle = out.ye + out.ynue - out.yanue;
  out.ylm = out.ym + out.ynum - out.yanum;

  return out;
}

ParticleFractions MuEOSClass::BaryonFractions(double n, double T, double *Y) {
  ParticleFractions out;

  out.yh = HeavyFraction(n, T, Y);   // Heavy nuclei fraction
  out.ya = AlphaFraction(n, T, Y);   // Alpha fraction
  out.yn = NeutronFraction(n, T, Y); // Neutron fracton
  out.yp = ProtonFraction(n, T, Y);  // Proton fraction

  return out;
} */

ParticleFractions MuEOSClass::NeutrinoFractions(double n, double T, double *Y, ChemPotentials *chem_pot) {
  const double eta_nue = chem_pot->mu_nue / T; //mu_nue/T
  const double eta_num = chem_pot->mu_num / T; //mu_num/T

  ParticleFractions out;

  out.ynue  = nu_fraction(n, T,  eta_nue);   //Y_nue
  out.yanue = anu_fraction(n, T, eta_nue);   //Y_anue
  out.ynum  = nu_fraction(n, T, eta_num);    //Y_num
  out.yanum = anu_fraction(n, T, eta_num);   //Y_anum
  out.ynux  = nu_fraction(n, T, 0.);         //Y_nux

  return out;
}

// @TODO: add possibility to set a density threshold for neutrinos
NeutrinoEOSOutput MuEOSClass::compute_neutrino_EOS(double n, double T, double *Y, double chem_pot[4]) {
  NeutrinoEOSOutput out;
  
  out.chem_pot[0] = chem_pot[0] - chem_pot[1] + chem_pot[2]; // Electron neutrino chemical potential [MeV]
  if (EOS_leptons<1>::m_lep_active == true) {
    out.chem_pot[1] = chem_pot[0] - chem_pot[1] + chem_pot[3];// Muon neutrino chemical potential [MeV]
  } else {
    out.chem_pot[1] = 0.;
  }

  const double eta_nue = out.chem_pot[0] / T; //mu_nue/T
  const double eta_num = out.chem_pot[1] / T; //mu_num/T

  out.Y_nu[0] =  nu_fraction(n, T,  eta_nue);  //Y_nue
  out.Y_nu[1] = anu_fraction(n, T, eta_nue);   //Y_anue
  out.Y_nu[2] =  nu_fraction(n, T, eta_num);   //Y_num
  out.Y_nu[3] = anu_fraction(n, T, eta_num);   //Y_anum
  out.Y_nu[4] =  nu_fraction(n, T, 0.);        //Y_nux


  out.Z_nue  = nu_energy(n, T,  eta_nue);   // Z_nue
  out.Z_anue = anu_energy(n, T, eta_nue);   // Z_anue
  out.Z_num  = nu_energy(n, T, eta_num);    // Z_num
  out.Z_anum = anu_energy(n, T, eta_num);   // Z_anum
  out.Z_nux  = nu_energy(n, T, 0.);         // Z_nux

  out.Z_tot = out.Z_nue + out.Z_anue + out.Z_num + out.Z_anum + 2. * out.Z_nux;
  
  out.s_nue  = nu_entropy(n, T,  eta_nue);   // s_nue
  out.s_anue = anu_entropy(n, T, eta_nue);   // s_anue
  out.s_num  = nu_entropy(n, T, eta_num);    // s_num
  out.s_anum = anu_entropy(n, T, eta_num);   // s_anum
  out.s_nux  = nu_entropy(n, T, 0.);         // s_nux

  out.s_tot = out.s_nue + out.s_anue + out.s_num + out.s_anum + 2. * out.s_nux;

  return out;
}

EOSstruct MuEOSClass::compute_full_EOS(double n, double T, double *Y) {
  EOSstruct eos_tot = {0};
  
  //out.nb = n;  // Baryon number density [1/fm3]
  //out.T  = T;  // Temperature [MeV]
 
  eos_tot.mb  = GetBaryonMass();
  //out.rho = 1.0E+39 * n * out.mb * MEOS_MeV2erg / (MEOS_c*MEOS_c); // Mass density [g/cm3]

  EOSstruct eos_bar = BaryonEOS(n, T, Y); // baryons
  EOSstruct eos_lep_0 = EOS_leptons<0>::LeptonEOS<id_test>(n, T, Y); // electrons and positrons
  EOSstruct eos_lep_1 = EOS_leptons<1>::LeptonEOS<id_test>(n, T, Y); // muons and antimuons

  eos_tot.comp[4] = Y[0];
  if (EOS_leptons<1>::m_lep_active == true) {
    eos_tot.comp[5] = Y[1];
  } else {
    eos_tot.comp[5] = 0;
  }

  eos_tot.e = eos_bar.e + eos_lep_0.e + eos_lep_1.e + RadEnergy(T);
  eos_tot.P = eos_bar.P + eos_lep_0.P + eos_lep_1.P + RadPressure(T);
  eos_tot.s = eos_bar.s + (eos_lep_0.s + eos_lep_1.s + RadEntropy(T)) / n;

  double enthalpy = (eos_tot.P + eos_tot.e) / n;

  eos_tot.e = MEOS_MeV2erg * 1.0E+39 * (eos_tot.e - eos_tot.mb * n); // + 1.e39*8.265e18*mb/(c*c)*MeV*nb; //erg/cm3
  eos_tot.P = MEOS_MeV2erg * 1.0E+39 * eos_tot.P;
  
  double const dPdn = eos_bar.eos_der[0] + eos_lep_0.eos_der[0] + eos_lep_1.eos_der[0];
  double const dsdn = eos_bar.eos_der[1] + eos_lep_0.eos_der[1] + eos_lep_1.eos_der[1] - (RadEntropy(T) / (n*n));
  double const dPdt = eos_bar.eos_der[2] + eos_lep_0.eos_der[2] + eos_lep_1.eos_der[2] +  RaddPdT(T);
  double const dsdt = eos_bar.eos_der[3] + eos_lep_0.eos_der[3] + eos_lep_1.eos_der[3] + (RaddsdT(T)/n);

  eos_tot.cs2 = sqrt(std::max(1.e-6,(dPdn - dsdn/dsdt*dPdt) / enthalpy));
  
  eos_tot.chem_pot[0] = eos_bar.chem_pot[0];
  eos_tot.chem_pot[1] = eos_bar.chem_pot[1];
  eos_tot.chem_pot[2] = eos_lep_0.chem_pot[2];
  eos_tot.chem_pot[3] = eos_lep_0.chem_pot[3];

  eos_tot.comp[0] = eos_bar.comp[0];
  eos_tot.comp[1] = eos_bar.comp[1];
  eos_tot.comp[2] = eos_bar.comp[2];
  eos_tot.comp[3] = eos_bar.comp[3];
  eos_tot.comp[4] = Y[0];
  eos_tot.comp[5] = Y[1];
  eos_tot.comp[6] = eos_lep_0.comp[6];
  eos_tot.comp[7] = eos_lep_1.comp[7];

  //if (cs2 >= 1.) cout << "cs2 > 1!" << endl;
  //cout << "dPdn = " << dPdn << "\t" << "dsdn = " << dsdn << "\t" << "dPdt = " << dPdt << "\t" << "dsdt = " << dsdt << endl << endl; 
  
  eos_tot.nuEOS = compute_neutrino_EOS(n, T, Y, eos_tot.chem_pot);

  return eos_tot;
}
