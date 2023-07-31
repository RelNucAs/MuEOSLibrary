#ifndef EOS_ASSEMBLED_HPP
#define EOS_ASSEMBLED_HPP

#include "eos_species.hpp"
#include "eos_leptons.hpp"

/*============================================================================*/

// file: eos_assembled.cpp

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
  double yp;    // Proton fraction [#/baryon]
  double ye;    // Electron fraction [#/baryon]
  double ym;    // Muon fraction [#/baryon]
  double ynue;  // Electron neutrino fraction [#/baryon]
  double yanue; // Electron antineutrino fraction [#/baryon]
  double ynum;  // Muon neutrino fraction [#/baryon]
  double yanum; // Muon antineutrino fraction [#/baryon]
  double ynux;  // Tau (anti)neutrino fraction [#/baryon]
  double yle;   // Net electronic lepton fraction [#/baryons]
  double ylm;   // Net muonic lepton fraction [#/baryon]
};
typedef struct ParticleFractions ParticleFractions;

struct NeutrinoEOSOutput {
  ParticleFractions Y_nu;
  double Z_nue;
  double Z_anue;
  double Z_num;
  double Z_anum;
  double Z_nux;
  double Z_tot;
  double s_nue;
  double s_anue;
  double s_num;
  double s_anum;
  double s_nux;
  double s_tot;
};
typedef struct NeutrinoEOSOutput NeutrinoEOSOutput;

// @TODO: store detailed EOS information for each species
struct FullEOSOutput {
  double rho;                // Mass density [g/cm3]
  double nb;                 // Baryon number density [1/cm3]
  double T;                  // Temperature [MeV]
  double mb;                 // Baryon mass [MeV]
  double e;                  // Internal energy density [erg/cm3]
  double P;                  // Pressure [erg/cm3]
  double s;                  // Entropy per baryon [#/baryon]
  ChemPotentials chem_pot;   // Relativistic chemical potentials [MeV]
  ParticleFractions Y_part;  // Particle fractions [#/baryons]
  NeutrinoEOSOutput nuEOS;   // Neutrino EOS quantities
};
typedef struct FullEOSOutput FullEOSOutput;

class EOS_assembled : public EOS_baryons, public EOS_leptons<0>, public EOS_leptons<1>, public EOS_photons, public EOS_neutrinos {

    public:
    /// Constructor
    EOS_assembled(const int id_eos, const bool el_bool, const bool mu_bool, std::string BarTableName);

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
    ParticleFractions GetParticleFractions(double n, double T, double *Y);

    /// Calculate the fractions of baryons using.
    ParticleFractions BaryonFractions(double n, double T, double *Y);

    /// Calculate the neutrino particle fractions using.
    ParticleFractions NeutrinoFractions(double n, double T, double *Y, ChemPotentials *chem_pot);

    /// Calculate the neutrino EOS quantities using.
    NeutrinoEOSOutput compute_neutrino_EOS(double n, double T, double *Y);

    /// Calculate the full EOS output using.
    FullEOSOutput compute_full_EOS(double n, double T, double *Y);
};

#endif //EOS_ASSEMBLED_HPP
