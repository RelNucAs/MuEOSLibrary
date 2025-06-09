#ifndef MUEOSLIBRARY_SRC_HELMHOLTZ_EOS_HPP_
#define MUEOSLIBRARY_SRC_HELMHOLTZ_EOS_HPP_

#include <iostream>
#include <array>
#include <cmath>

#include "constants.hpp"
#include "fermi_integrals/fermi_integrals.hpp"

enum EOSQuantities {
	IL_N  = 0, //! Number density of leptons
	IA_N  = 1, //! Number density of anti-leptons
	IL_P  = 2, //! Pressure of leptons
	IA_P  = 3, //! Pressure of anti-leptons
	IL_E  = 4, //! Internal energy density of leptons
	IA_E  = 5, //! Internal energy density of anti-leptons
	IL_S  = 6, //! Entropy density of leptons
	IA_S  = 7, //! Entropy density of anti-leptons
	IL_MU = 8, //! Chemical potential of leptons
	NVARS = 9
};


// Definition of constants
const double K0    = 8. * sqrt(2.) * MEOS_pi / pow(MEOS_h*MEOS_c*1.e13,3.);
const double mL[2] = {MEOS_me, MEOS_mmu};
const double K[2]  = {K0*pow(mL[0],3.), K0*pow(mL[1],3.)}; // 31217845.162531383*mLep**3 in cm-3
const double K3[2] = {K[0]/3., K[1]/3.};

struct HelmEOSOutput {
  double nl;   // number density of leptons
  double a_nl; // number density of anti-leptons
  double pl;   // pressure of leptons
  double a_pl; // pressure of anti-leptons
  double el;   // internal energy density of leptons
  double a_el; // internal energy density of anti-leptons
  double sl;   // entropy density of leptons
  double a_sl; // entropy density of anti-leptons
  double mul;  // chemical potential of leptons
};
typedef struct HelmEOSOutput HelmEOSOutput;

struct HelmEOSDer {
  double dPdn;
  double dsdn;
  double dPdt;
  double dsdt;
};
typedef struct HelmEOSDer HelmEOSDer;

HelmEOSOutput eos_helm_from_eta(const double eta, const double T, const int id_L, GFDs *FD);

HelmEOSOutput eos_helm_full(double nLep, double temp, const int id_L);

HelmEOSDer der_cs2(double nLep, double temp, const int id_L);
HelmEOSDer der_cs2_from_eta(double eta, double temp, const int id_L);
HelmEOSDer der_cs2_num(double nLep, double temp, double id_L);

#endif // MUEOSLIBRARY_SRC_HELMHOLTZ_EOS_HPP_
