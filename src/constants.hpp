#ifndef MUEOSLIBRARY_SRC_CONSTANTS_HPP_
#define MUEOSLIBRARY_SRC_CONSTANTS_HPP_

// Unit conversions
constexpr double MEOS_MeV2erg = 1.602176634e-6; // MeV to erg (CGS) conversion
constexpr double MEOS_cm2fm   = 1.0e13;         // cm to fm conversion 
 	
// Particle masses [MeV]
constexpr double MEOS_me  = 0.51099895;        // electron
constexpr double MEOS_mmu = 105.65837;         // muon
constexpr double MEOS_mn  = 939.5654133;       // neutron
constexpr double MEOS_mp  = 938.2720813;       // proton
constexpr double MEOS_Q   = MEOS_mn - MEOS_mp; // neutron-proton mass difference
	
// Physical constants
constexpr double MEOS_h    = 4.135667696e-21;    // Planck constant [MeV s]
constexpr double MEOS_hbar = 6.582119569e-22;    // reduced Planck constant [MeV s]
constexpr double MEOS_c    = 2.99792458e+10;     // speed of light in vacuum [cm s-1]
constexpr double MEOS_hc   = 1.2398419840550368e+3; // h*c [MeV fm]
constexpr double MEOS_hc3  = 1.905895196929953e+9;  // (h*c)**3 [MeV3 fm3]
constexpr double MEOS_fourpi_hc3 = 6.593421629164755e-09; // (4*pi)/(h*c)**3 [MeV-3 fm-3]
constexpr double MEOS_GF   = 8.957e-44;          // MeV*cm^3
constexpr double MEOS_pi   = 3.141592653589793;  // pi constant
constexpr double MEOS_kB   = 8.617333262145e-11; // Boltzmann constant [MeV K-1]
constexpr double MEOS_mb   = 1.674e-24;          // reference baryonic mass [g]
constexpr double MEOS_mu   = 1.66054e-24;        // atomic mass unit [g]

// Coupling constants
constexpr double MEOS_gA = 1.23;
constexpr double MEOS_gV = 1.;
constexpr double MEOS_gS = 0.;
constexpr double MEOS_sinsqthetaw = 0.2325;

// Neutral current nucleon form factors (Q^2=0)
constexpr double MEOS_hnv = -0.5;
constexpr double MEOS_hna = -0.5*MEOS_gA;
constexpr double MEOS_hpv =  0.5-2.*MEOS_sinsqthetaw;
constexpr double MEOS_hpa =  0.5*MEOS_gA;

#endif // MUEOSLIBRARY_SRC_CONSTANTS_HPP_
