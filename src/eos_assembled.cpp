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

#include <eos_assembled.hpp>

//EOS_assembled::EOS_assembled():{};
  //EOS_leptons<0>::m_lep_active(true),
  //EOS_leptons<1>::m_lep_active(false) {};
  //m_eos_initialized(false) {};

//EOS_assembled::~EOS_assembled() {
//}


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

void EOS_assembled::BaryonFractions(double n, double T, double *Y, double *out) {
  out[0] = HeavyFraction(n, T, Y);   // Heavy nuclei fraction
  out[1] = AlphaFraction(n, T, Y);   // Alpha fraction
  out[2] = NeutronFraction(n, T, Y); // Neutron fracton
  out[3] = ProtonFraction(n, T, Y);  // Proton fraction   
  return;
}

void EOS_assembled::NuParticleFractions(double n, double T, double *Y, double *out) {
  //double chem_pot[6] = {0.0};
  ChemPotentials chem_pot = GetChemicalPotentials(n, T, Y);
  
  const double eta_nue = chem_pot.mu_nue / T; //mu_nue/T
  const double eta_num = chem_pot.mu_num / T; //mu_num/T

  out[0] = nu_fraction(n, T,  eta_nue);   //Y_nue
  out[1] = anu_fraction(n, T, eta_nue);   //Y_anue
  out[2] = nu_fraction(n, T, eta_num);    //Y_num
  out[3] = anu_fraction(n, T, eta_num);   //Y_anum
  out[4] = nu_fraction(n, T, 0.);         //Y_nux
  return;
}

void EOS_assembled::NuEnergyFractions(double n, double T, double *Y, double *out) {
  //double chem_pot[6] = {0.0};
  ChemPotentials chem_pot = GetChemicalPotentials(n, T, Y);
  
  const double eta_nue = chem_pot.mu_nue / T; //mu_nue/T
  const double eta_num = chem_pot.mu_num / T; //mu_num/T

  out[0] = nu_energy(n, T,  eta_nue);   //Z_nue
  out[1] = anu_energy(n, T, eta_nue);   //Z_anue
  out[2] = nu_energy(n, T, eta_num);    //Z_num
  out[3] = anu_energy(n, T, eta_num);   //Z_anum
  out[4] = nu_energy(n, T, 0.);         //Z_nux
  return;
}

void EOS_assembled::NuEntropies(double n, double T, double *Y, double *out) {
  //double chem_pot[6] = {0.0};
  ChemPotentials chem_pot = GetChemicalPotentials(n, T, Y);
  
  const double eta_nue = chem_pot.mu_nue / T; //mu_nue/T
  const double eta_num = chem_pot.mu_num / T; //mu_num/T

  out[0] = nu_entropy(n, T,  eta_nue);   //Z_nue
  out[1] = anu_entropy(n, T, eta_nue);   //Z_anue
  out[2] = nu_entropy(n, T, eta_num);    //Z_num
  out[3] = anu_entropy(n, T, eta_num);   //Z_anum
  out[4] = nu_entropy(n, T, 0.);         //Z_nux
  return;
}
void EOS_assembled::ReadTables(std::string BarTableName,
                               std::string ETableName,
                               std::string MTableName) {
  ReadBarTableFromFile(BarTableName);
  EOS_leptons<0>::ReadLepTableFromFile(ETableName);
  EOS_leptons<1>::ReadLepTableFromFile(MTableName);
  return;
}
