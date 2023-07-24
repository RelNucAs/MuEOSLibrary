#ifndef ROOT_FINDING_HPP
#define ROOT_FINDING_HPP

#include "../../fermi_integrals/fermi_integrals.hpp"

/*============================================================================*/

// file: newt_raphson.cpp

double n_net_f(const double eta, const double T, const int id_L, GFDs *FD);

double n_net_df(const double eta, const double T, const int id_L);

double find_guess_eta(double nLep, double T, const int id_L);
 
double rtsafe(const double nLep, const double T, const int id_L, GFDs *FD);

#endif //ROOT_FINDING_HPP