//! \file eos_assembled.cpp
//  \brief Implementation of EOS_assembled

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

#include "eos_assembled.hpp"

std::string ETableName = "eos_table/electrons/eos_electrons_table.txt"; // electron table
std::string MTableName = "eos_table/muons/eos_muons_table.txt"; // muon table

/*
  Inputs:
   - id_EOS: method for EOS computation (1: interpolation, 2: on-the-fly)
   - el_flag: flag for activating electrons
   - mu_flag: flag for activating muons
   - BarTableName: path of baryon EOS table
*/
EOS_assembled::EOS_assembled(const int id_eos, const bool el_flag, const bool mu_flag, std::string BarTableName) {
  ReadBarTableFromFile("eos_table/baryons/" + BarTableName);
  EOS_leptons<0>::m_lep_active = el_flag;
  EOS_leptons<1>::m_lep_active = mu_flag;
  if (id_eos == 1) {
    if (el_flag == true) EOS_leptons<0>::ReadLepTableFromFile(ETableName);
    if (mu_flag == true) EOS_leptons<1>::ReadLepTableFromFile(MTableName);
  } else if (id_eos != 2) {
    std::cout << "ERROR: wrong index for EOS computation in InitalizeEOS" << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

double EOS_assembled::Energy(double n, double T, double *Y) {
  return BarEnergy(n, T, Y) + EOS_leptons<0>::LepEnergy<id_test>(n, T, Y) + EOS_leptons<1>::LepEnergy<id_test>(n, T, Y) + RadEnergy(T);
}

double EOS_assembled::Pressure(double n, double T, double *Y) {
  return BarPressure(n, T, Y) + EOS_leptons<0>::LepPressure<id_test>(n, T, Y) + EOS_leptons<1>::LepPressure<id_test>(n, T, Y) + RadPressure(T);
}

double EOS_assembled::Entropy(double n, double T, double *Y) {
  return BarEntropy(n, T, Y) + (EOS_leptons<0>::LepEntropy<id_test>(n, T, Y) + EOS_leptons<1>::LepEntropy<id_test>(n, T, Y) + RadEntropy(T)) / n;
}

double EOS_assembled::Enthalpy(double n, double T, double *Y) {
  double const P = Pressure(n, T, Y);
  double const e = Energy(n, T, Y);
  return (P + e)/n;
}

double EOS_assembled::SoundSpeed(double n, double T, double *Y) {

  double const dPdn = BardPdn(n, T, Y) + EOS_leptons<0>::LdPdn<id_test>(n, T, Y) + EOS_leptons<1>::LdPdn<id_test>(n, T, Y);
  double const dsdn = Bardsdn(n, T, Y) + EOS_leptons<0>::Ldsdn<id_test>(n, T, Y) + EOS_leptons<1>::Ldsdn<id_test>(n, T, Y) - (RadEntropy(T) / (n*n));
  double const dPdt = BardPdT(n, T, Y) + EOS_leptons<0>::LdPdt<id_test>(n, T, Y) + EOS_leptons<1>::LdPdt<id_test>(n, T, Y) +  RaddPdT(T);
  double const dsdt = BardsdT(n, T, Y) + EOS_leptons<0>::Ldsdt<id_test>(n, T, Y) + EOS_leptons<1>::Ldsdt<id_test>(n, T, Y) + (RaddsdT(T)/n);

  double const cs2 = sqrt(std::max(1.e-6,(dPdn - dsdn/dsdt*dPdt) / Enthalpy(n, T, Y)));
  //if (cs2 >= 1.) cout << "cs2 > 1!" << endl;
  //cout << "dPdn = " << dPdn << "\t" << "dsdn = " << dsdn << "\t" << "dPdt = " << dPdt << "\t" << "dsdt = " << dsdt << endl << endl; 
  return cs2;
}

ChemPotentials EOS_assembled::GetChemicalPotentials(double n, double T, double *Y) {
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

ParticleFractions EOS_assembled::GetParticleFractions(double n, double T, double *Y) {
  ParticleFractions out;

  if (EOS_leptons<0>::m_lep_active == true) {
    out.ye = Y[0]; // Electron fraction [#/baryon]
  } else {
    out.ye = 0.;
  }
  if (EOS_leptons<1>::m_lep_active == true) {
    out.ym = Y[1]; // Muon fraction [#/baryon]
  } else {
    out.ym = 0.;
  }

  // Baryon fractions
  ParticleFractions Y_bar = BaryonFractions(n, T, Y);
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

ParticleFractions EOS_assembled::BaryonFractions(double n, double T, double *Y) {
  ParticleFractions out;

  out.yh = HeavyFraction(n, T, Y);   // Heavy nuclei fraction
  out.ya = AlphaFraction(n, T, Y);   // Alpha fraction
  out.yn = NeutronFraction(n, T, Y); // Neutron fracton
  out.yp = ProtonFraction(n, T, Y);  // Proton fraction

  return out;
}

ParticleFractions EOS_assembled::NeutrinoFractions(double n, double T, double *Y, ChemPotentials *chem_pot) {
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
NeutrinoEOSOutput EOS_assembled::compute_neutrino_EOS(double n, double T, double *Y) {
  NeutrinoEOSOutput out;
  
  ChemPotentials chem_pot = GetChemicalPotentials(n, T, Y);

  out.Y_nu = NeutrinoFractions(n, T, Y, &chem_pot);

  const double eta_nue = chem_pot.mu_nue / T; //mu_nue/T
  const double eta_num = chem_pot.mu_num / T; //mu_num/T

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

FullEOSOutput EOS_assembled::compute_full_EOS(double n, double T, double *Y) {
  FullEOSOutput out;
  
  out.nb = n;  // Baryon number density [1/fm3]
  out.T  = T;  // Temperature [MeV]
 
  out.mb  = GetBaryonMass();
  out.rho = 1.0E+39 * n * out.mb * MeV / (c*c); // Mass density [g/cm3]

  out.chem_pot = GetChemicalPotentials(n, T, Y);
  out.Y_part   = GetParticleFractions(n, T, Y);

  out.e = MeV * 1.0E+39 * (Energy(n, T, Y) - out.mb * n); // + 1.e39*8.265e18*mb/(c*c)*MeV*nb; //erg/cm3
  out.P = MeV * 1.0E+39 * Pressure(n, T, Y);
  out.s = Entropy(n, T, Y);

  out.nuEOS = compute_neutrino_EOS(n, T, Y);

  return out;	
}