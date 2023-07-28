#ifndef EOS_ASSEMBLED_H
#define EOS_ASSEMBLED_H

//! \file eos_assembled.hpp

#include <cstddef>
#include <string>

#include <eos_photons.hpp>
#include <eos_leptons.hpp>
#include <eos_baryons.hpp>
#include <eos_neutrinos.hpp>

struct ChemPotentials {
  double mu_n;    // Neutron chemical potential [MeV]
  double mu_p;    // Proton chemical potential [MeV]
  double mu_e;    // Electron chemical potential [MeV]
  double mu_m;    // Muon chemical potential [MeV]
  double mu_nue;  // Electron neutrino chemical potential [MeV] 
  double mu_num;  // Muon neutrino chemical potential [MeV]
  double mu_nux;  // Tau neutrino chemical potential [MeV]
};
typedef struct ChemPotentials ChemPotentials;

struct ParticleFractions {
  double yh;    // Heavy nuclei fraction [#/baryon]
  double ya;    // Alpha particle fraction [#/baryon]
  double yn;    // Neutron fraction [#/baryon]
  double yp;    // Muon chemical potential [MeV]
  double ye;    // Electron fraction [#/baryon]
  double ym;    // Muon fraction [#/baryon]
  double ynue;  // Electron neutrino fraction [#/baryon]
  double yanue; // Electron antineutrino fraction [#/baryon]
  double ynum;  // Muon neutrino fraction [#/baryon]
  double yanum; // Muon antineutrino fraction [#/baryon]
  double ynux;  // Tau (anti)neutrino fraction [#/baryon]};
};
typedef struct ParticleFrac ParticleFrac;

class EOS_assembled : public EOS_baryons, public EOS_leptons<0>, public EOS_leptons<1>, public EOS_photons, public EOS_neutrinos {
    
    public:
    /// Constructor
    //EOS_assembled();

    /// Destructor
    //~EOS_assembled();

    /// Calculate the energy density using.
    double Energy(double n, double T, double *Y);

    /// Calculate the pressure using.
    double Pressure(double n, double T, double *Y);

    /// Calculate the entropy per baryon using.
    double Entropy(double n, double T, double *Y);

    /// Calculate the enthalpy per baryon using.
    double Enthalpy(double n, double T, double *Y);

    /// Calculate the sound speed using.
    double SoundSpeed(double n, double T, double *Y);

    /// Calculate the chemical potentials using.
    ChemPotentials GetChemicalPotentials(double n, double T, double *Y);
    
    /// Calculate the particle fractions using.
    void GetParticleFractions(double n, double T, double *Y, double *out);
    
    /// Calculate the fractions of baryons using.
    void BaryonFractions(double n, double T, double *Y, double *out);

    /// Calculate the neutrino particle fractions using.
    void NuParticleFractions(double n, double T, double *Y, double *out);

    /// Calculate the neutrino energy fractions using.
    void NuEnergyFractions(double n, double T, double *Y, double *out);

    /// Calculate the neutrino entropies per baryon using.
    void NuEntropies(double n, double T, double *Y, double *out);

    /// Read in tables
    void ReadTables(std::string BarTableName,
                    std::string ETableName,
                    std::string MTableName);

};

#endif
