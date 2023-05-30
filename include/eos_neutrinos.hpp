#ifndef EOS_NU_H
#define EOS_NU_H

/* NEUTRINO EOS AT EQUILIBRIUM WITH RADIATION (neutrinos degeneracy parameter is equal to minus anti-neutrinos degeneracy parameter) */

#include "constants.hpp"
#include "fermi_integrals.h"

using namespace constants;

/* Compute constant in front of neutrino EOS once for all */
const double nu_const = 4. * pi / pow(h*c*1.e13,3.); // MeV-3 fm-3

/* Input parameters:
 *    - number density                nb [fm-3]
 *    - temperature                   t  [MeV]
 *    - neutrino degeneracy parameter eta
*/

// Number fraction of neutrinos (yn)
double nu_fraction(const double nb, const double T, const double eta) {
	return nu_const * pow(T,3.) * Fermi_integral_p2(eta) / nb; // 1/baryon
}

// Number fraction of antineutrinos (yan)
double anu_fraction(const double nb, const double T, const double eta) {
	return nu_fraction(nb, T, -eta); // 1/baryon
	//return nu_const * pow(T,3.) * Fermi_integral_p2(-eta) / nb; // 1/baryon
}

// Energy per baryon of neutrinos (zn)
double nu_energy(const double nb, const double T, const double eta) {
        return nu_const * pow(T,4.) * Fermi_integral_p3(eta) / nb; // MeV/baryon
}

// Energy per baryon of antineutrinos (zan)
double anu_energy(const double nb, const double T, const double eta) {
        return nu_energy(nb, T, -eta); // MeV/baryon
        //return nu_const * pow(T,4.) * Fermi_integral_p3(-eta) / nb; // MeV/baryon
}

// Pressure per baryon of neutrinos (pn)
double nu_pressure(const double nb, const double T, const double eta) {
        return nu_energy(nb, T, eta) / 3.; // MeV/baryon
	//return nu_const * pow(T,4.) * Fermi_integral_p3(eta) / (3. * nb); // MeV/baryon
}

// Pressure per baryon of antineutrinos (pan)
double anu_pressure(const double nb, const double T, const double eta) {
        return nu_pressure(nb, T, -eta); // MeV/baryon
        //return anu_energy(nb, T, eta) / 3.; // MeV/baryon
        //return nu_const * pow(T,4.) * Fermi_integral_p3(-eta) / (3. * nb); // MeV/baryon
}

// Entropy per baryon of neutrinos (sn)
double nu_entropy(const double nb, const double T, const double eta) {
	return nu_const * pow(T,3.) * (4. * Fermi_integral_p3(eta) / 3. - eta * Fermi_integral_p2(eta)) / nb; // 1/baryon (kB=1)
	//return (4.*nu_energy(nb, T, eta) /(3.*T) - eta*nu_fraction(nb, T, eta)); // equivalent expression (yet not tested)
}

// Entropy per baryon of antineutrinos (san)
double anu_entropy(const double nb, const double T, const double eta) {
	return nu_entropy(nb, T, -eta); // 1/baryon (kB=1)
	//return (4.*anu_energy(nb, T, eta)/(3.*T) + eta*anu_fraction(nb, T, eta)); // equivalent expression (yet not tested)
}

#endif
