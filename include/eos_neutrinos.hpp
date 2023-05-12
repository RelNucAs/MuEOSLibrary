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
	return nu_const * pow(T,3.) * Fermi_integral_p2(eta) / nb;
}

// Number fraction of antineutrinos (yan)
double anu_fraction(const double nb, const double T, const double eta) {
	return nu_fraction(nb, T, -eta);
	//return nu_const * pow(T,3.) * Fermi_integral_p2(-eta) / nb;
}

// Energy per baryon of neutrinos (en)
double nu_energy(const double nb, const double T, const double eta) {
        return nu_const * pow(T,4.) * Fermi_integral_p3(eta) / nb; // MeV/baryon
}

// Energy per baryon of antineutrinos (ean)
double anu_energy(const double nb, const double T, const double eta) {
        return nu_energy(nb, T, -eta); // MeV/baryon
        //return nu_const * pow(T,4.) * Fermi_integral_p3(-eta) / nb; // MeV/baryon
}

#endif
