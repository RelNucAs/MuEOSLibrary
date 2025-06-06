#ifndef MUEOSLIBRARY_SRC_MUEOSCLASS_HPP_
#define MUEOSLIBRARY_SRC_MUEOSCLASS_HPP_

//#ifndef REAL_TYPE
//#define REAL_TYPE double
//#define REAL_TYPE_IS_DOUBLE
//#warning "REAL_TYPE for MuEOSLibrary has been set to double (default)."
//#endif

//typedef REAL_TYPE MEOS_REAL;

#define POW0(X) ((1))
#define POW1(X) ((X))
#define POW2(X) ((X) * (X))
#define POW3(X) ((X) * (X) * (X))
#define POW4(X) ((X) * (X) * (X) * (X))
#define POW5(X) ((X) * (X) * (X) * (X) * (X))
#define POW6(X) ((X) * (X) * (X) * (X) * (X) * (X))
#define POW7(X) ((X) * (X) * (X) * (X) * (X) * (X) * (X))


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
  double Y_nu[5];
  double chem_pot[2];
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

struct EOSstruct {
//  double rho;                // Mass density [g/cm3]
//  double T;                  // Temperature [MeV]
  double mb;                 // Baryon mass [MeV]
  double e;                  // Internal energy density [erg/cm3]
  double P;                  // Pressure [erg/cm3]
  double s;                  // Entropy per baryon [#/baryon]
  double chem_pot[4];        // Relativistic chemical potentials [MeV]
  double comp[8];            // Particle fractions [#/baryons]
  double eos_der[4];         // EOS derivatives
  double cs2;                // Sound speed squared
  NeutrinoEOSOutput nuEOS;   // Neutrino EOS quantities
};
typedef struct EOSstruct EOSstruct;

#endif // MUEOSLIBRARY_SRC_MUEOSCLASS_HPP_
