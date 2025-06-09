#include <iostream>
#include <fstream>
#include <iomanip>
#include <array>
#include <cmath>
#include <chrono>

#include "eos_species/eos_assembled.hpp"

using namespace std::chrono;

typedef std::chrono::time_point<std::chrono::high_resolution_clock> chrono_t; //Chrono type.
typedef std::chrono::duration<long int, std::ratio<1, 1000> > chrono_int;

int main () {
  /* Name of baryon EOS table */
  std::string BarTableName = "DD2_bar.h5";  // baryon table

  const int niter = 10;

  double r1, r2, r3, r4;



  chrono_t begin, stop;
  chrono_int duration;

  /* Initialize global EOS class

  Constructor -> EOS_assembled(const int id_eos, const bool el_flag, const bool mu_flag, std::string BarTableName)

  Inputs:
   - id_EOS: method for EOS computation (1: interpolation, 2: on-the-fly)
   - el_flag: flag for activating electrons
   - mu_flag: flag for activating muons
   - BarTableName: path of baryon EOS table  */
  EOS_assembled eos(1, true, false, BarTableName);      

  /* Arrays of input variables */
  // Number density
  const double nbmin = 7.58576778E-07;
  const double nbmax = 9.12011723E+00;
  const double dn = log10(nbmax) - log10(nbmin);

  // Temperature
  const double tmin = 1.00000000E-01;
  const double tmax = 1.31825639E+02;
  const double dt = log10(tmax) - log10(tmin);

  // Electron fraction
  const double yemin = 9.99999978E-03;
  const double yemax = 5.79999983E-01;
  const double dye = yemax - yemin;

  // Muon fraction
  const double ymmin = 2.0E-05;
  const double ymmax = 2.0E-01;
  const double dym = ymmax - ymmin;

  //Initialization of random number generator
  srand(time(NULL));

  std::vector<double> n_array; // log spaced
  std::vector<double> t_array; // log spaced
  std::vector<double> ye_array; // log spaced
  // std::vector<double> ym_array; // log spaced

  //Loop over (ne-t) grid
  for (int i=0; i<niter; i++) {
      r1 = rand() / (double) RAND_MAX; //Returns a pseudo-random number between 0 and 1
      r2 = rand() / (double) RAND_MAX;
      r3 = rand() / (double) RAND_MAX;
      r4 = rand() / (double) RAND_MAX;

      n_array.push_back(pow(10., log10(nbmin) + r1*dn));
      t_array.push_back(pow(10., log10(tmin)  + r2*dt));
      ye_array.push_back(yemin  + r3*dye);
      // ym_array.push_back(ymmin  + r4*dym);
  }


  double nb, T;
  double Y[2];
	
  Y[1] = 0.0; // set muon fraction equal to zero

  FullEOSOutput eos_out;

  begin = high_resolution_clock::now();
  for (int i=0; i<niter; i++) {
    nb = n_array[i];
    T  = t_array[i];
    Y[0] = ye_array[i];

    eos_out = eos.compute_full_EOS(nb, T, Y);  
  }

  stop = high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - begin);
  printf("%.5e  s\n", duration.count()*1.e-3);

  return 0;
}
