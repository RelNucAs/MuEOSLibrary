/* RADIATION EOS */
#pragma once

#include <cmath>

double zeta3 = std::riemann_zeta(3.);

//number density
double n_rad(double T); //T in MeV

//internal energy (per unit volume)
double e_rad(double T); //T in MeV

//pressure 
double P_rad(double T); //T in MeV
	
//entropy
double s_rad(double T); //T in MeV

