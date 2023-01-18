/* NEUTRINO EOS AT EQUILIBRIUM WITH RADIATION (neutrinos degeneracy parameter is equal to minus anti-neutrinos degeneracy parameter) */
#pragma once

//fraction of neutrinos (yn)
double nu_fraction(double rho, double T, double eta);

//fraction of antineutrinos (yan)
double anu_fraction(double rho, double T, double eta);

//energy per baryon of neutrinos (en)
double nu_energy(double rho, double T, double eta);

//energy per baryon of antineutrinos (ean)
double anu_energy(double rho, double T, double eta);
