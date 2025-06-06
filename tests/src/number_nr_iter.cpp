#include <iostream>
#include <fstream>
#include <iomanip>
#include <array>
#include <cmath>

#include "../../src/eos_species/eos_assembled.hpp"
#include "../../src/num_tools/root_finding/root_finding.hpp"

int main () {
  const double nb_min[2] = {1.0000000e-12, 1.00e-12}; // {5.0e+03/amu_g*1.e-39, ...}
  const double nb_max[2] = {1.9054607e+00, 1.90e+00}; // {2.0e+16/amu_g*1.e-39, ...}
  const double  t_min[2] = {1.0000000e-01, 1.00e-01}; // {5.00e-02, ...}
  const double  t_max[2] = {1.4454398e+02, 1.58e+02}; // {1.05e+02, ...}  
  const double  y_min[2] = {1.00e-02, 5.0e-08}; // {5.0e-03, ...}
  const double  y_max[2] = {6.00e-01, 5.0e-01}; // {5.5e-01, ...}
  
  const int n_nl[2] = {700, 750}; // {1400, 1500}
  const int  n_t[2] = {150, 150}; // {300, 300}

  const int id_L = 1; // 0: electrons, 1: muons

  const double nl_min = y_min[id_L] * nb_min[id_L];
  const double nl_max = y_max[id_L] * nb_max[id_L];

  const double log_nmin = log10(nl_min);
  const double log_nmax = log10(nl_max);

  const double log_tmin = log10(t_min[id_L]);
  const double log_tmax = log10(t_max[id_L]);

  const int n1 = n_nl[id_L];
  const int n2 = n_t[id_L];

  std::vector<double>  n_array, t_array;

  for (int i=0; i<n1; i++) n_array.push_back(log_nmin + static_cast<double>(i) * (log_nmax-log_nmin) / static_cast<double>(n1-1));
  for (int i=0; i<n2 ;i++) t_array.push_back(log_tmin + static_cast<double>(i) * (log_tmax-log_tmin) / static_cast<double>(n2-1));

  double nLep, temp, eta;
  double r1, r2, dn, dt;

  GFDs FD_integ;
	
  //Initialization of random number genLrator
  srand(time(NULL));

  for (int i=0; i<n1-1; i++) {
    //std::cout << "i = " << i << std::endl;
     for (int j=0; j<n2-1; j++) {

	for (int k=0; k<100; k++) {
  	  r1 = rand() / (double) RAND_MAX; //Returns a pseudo-random number between 0 and 1
	  r2 = rand() / (double) RAND_MAX;
          
	  dn = n_array[i+1] - n_array[i];
          dt = t_array[j+1] - t_array[j];

	  nLep = pow(10.,n_array[i] + r1*dn);
          temp = pow(10.,t_array[j] + r2*dt);

	  // std::cout << nLep << " " << temp << std::endl;

	  eta = rtsafe(1.e39*nLep, temp, id_L, &FD_integ);
       }
     }
  }

  return 0;
}
