#include <math.h>

#include "root_finding.hpp"
#include "../../constants.hpp"
#include "../../helmholtz_eos/helmholtz_eos.hpp"
#include "../../fermi_integrals/fermi_integrals.hpp"

/* 0: electrons, 1: muons */
const double eta_min[2] = {-2.3e+04, -2.3e+04}; // minimum eta for root-finding bracketing
const double eta_max[2] = {+2.3e+04, +5.0e+04}; //1.3e+04 // maximum eta for root-finding bracketing


// Given the density rho, the non relativistic degeneracy eta, the temperature T
// and the particle species, the subroutine computes the net fraction ynet and the first eta derivative.
//.......UNITS: rho in g/cm**3; T in MeV
double n_net_f(const double eta, const double T, const int id_L, GFDs *FD) {
  const double theta = T / mL[id_L];
  const double a_eta = - eta - 2. / theta;
  
  double f1, f2, f3, f4;
    
  f1 = compute_res(  eta, theta, 0.5);
  f2 = compute_res(a_eta, theta, 0.5);
  f3 = compute_res(  eta, theta, 1.5);
  f4 = compute_res(a_eta, theta, 1.5);
  
  FD->f12   = f1;
  FD->a_f12 = f2;
  FD->f32   = f3;
  FD->a_f32 = f4;

  // @TODO: fix units
  return 1.0e+39 * K[id_L] * pow(theta,1.5) * (f1 - f2 + theta * (f3 - f4));
}


double n_net_df(const double eta, const double T, const int id_L) {
  const double theta = T / mL[id_L];
  const double a_eta = - eta - 2. / theta;

  double f1, f2, f3, f4;

  f1 = compute_res_ed(  eta, theta, 0.5);
  f2 = compute_res_ed(a_eta, theta, 0.5);
  f3 = compute_res_ed(  eta, theta, 1.5);
  f4 = compute_res_ed(a_eta, theta, 1.5);

  // @TODO: fix units
  return 1.0e+39 * K[id_L] * pow(theta,1.5) * (f1 - f2 + theta * (f3 - f4));
}

// @TODO: use constexpr from C++17 on to avoid if statement at execution time
//template<int species>
double find_guess_eta(double nLep, double T, const int id_L) {
  double nq;

  // @TODO: fix units and computation of constants
  /* Electrons */
  if (id_L == 0) { //if constexpr(species == 0) {
    const double t = T/MEOS_kB; //temp in K
    const double nfm = nLep*1.e-39; //number density in fm^{-3}
    const double pF = MEOS_h*MEOS_c*pow(3.*nLep/(8.*MEOS_pi),1./3.); //[MeV]

    // Find the guess
    if ((t<=1.e10) && (nfm <= 1.e-17)) { //classical NR
      nq = pow(2.*MEOS_pi*MEOS_me*T/(MEOS_h*MEOS_h*MEOS_c*MEOS_c),1.5); //cm^{-3}
      return -log(fabs(2.*nq/nLep));
    } else if ((t<=1.e10) && (nfm>1.e-17) && (nfm <= 1.e-10)) { //degenerate NR
      return (pF*pF/(2.*MEOS_me*T)) - (MEOS_me/T);
    } else if ((t<=1.e10) && (nfm>1.e-10)) { //degenerate UR
      return (pF/T) - (MEOS_me/T);
    } else if (t>1.e10) {
      nq = 8.*MEOS_pi*pow(T/(MEOS_h*MEOS_c),3.); //classical UR
      return -log(fabs(2.*nq/nLep)) - (MEOS_me/T);
    }

  /* Muons */
  } else if (id_L == 1) { //else if constexpr(species == 1) {
    // Find the guess
    nq = pow(2.*MEOS_pi*MEOS_mmu*T/(MEOS_h*MEOS_h*MEOS_c*MEOS_c),1.5); //cm^{-3}
    return -log(fabs(2.*nq/nLep)); 
  }
}


double rtsafe(const double nLep, const double T, const int id_L, GFDs *FD) {
//Using a combination of Newton-Raphson and bisection, return the root of a function bracketed
//between x1 and x2. The root will be refined until its accuracy is known within xacc.
  const int MAXIT=150; // Maximum allowed number of iterations.
  const double xacc = 1.e-7; // set the accuracy for Newton Raphson
  double x1, x2;
  double xh, xl;

  x1 = eta_min[id_L];
  x2 = eta_max[id_L];


  double fl = n_net_f(x1, T, id_L, FD) - nLep;
  double fh = n_net_f(x2, T, id_L, FD) - nLep;

  if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)) {
    printf("nLep = %.3e cm-3, T = %.3e MeV, fl = %.3e cm-3, fh = %.3e cm-3\n", nLep, T, fl, fh);
    throw("Root must be bracketed in rtsafe");
  }
  if (fl == 0.0) return x1;
  if (fh == 0.0) return x2;
  if (fl < 0.0) { // Orient the search so that f(xl) < 0.
    xl = x1;
    xh = x2;
  } else {
    xh = x1;
    xl = x2;
  }


  double rts = find_guess_eta(nLep, T, id_L); //0.5*(x1+x2);  // Initialize the guess for root,

  double dxold = fabs(x2-x1);         // the “stepsize before last,”
  double dx = dxold;                  // and the last step.
  double f  = n_net_f(rts, T, id_L, FD) - nLep;
  double df = n_net_df(rts, T, id_L);
  for (int j=0;j<MAXIT;j++) { // Loop over allowed iterations.
    if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0) || (fabs(2.0*f) > fabs(dxold*df))) { // Bisect if Newton out of range, or not decreasing fast enough.
      dxold = dx;
      dx = 0.5 * (xh - xl);
      rts = xl + dx;
      if (xl == rts) {
        // std::cout << j << std::endl;
	return rts;
      }
    } else { // Change in root is negligible. Newton step acceptable. Take it.
      dxold = dx;
      dx = f/df;
      double temp = rts;
      rts -= dx;
      if (temp == rts) {
        // std::cout << j << std::endl;
	return rts;
      }
    }
    if (fabs(dx) < xacc) {
      //std::cout << std::endl << "Num steps: " << j << std::endl;
      // std::cout << j << std::endl;
      return rts; // Convergence criterion.
    }

    f  = n_net_f(rts, T, id_L, FD) - nLep;
    df = n_net_df(rts, T, id_L); // The one new function evaluation per iteration.
    if (f < 0.0) { //Maintain the bracket on the root.
      xl = rts;
    } else {
      xh = rts;
    }
  }

  throw("Maximum number of iterations exceeded in rtsafe");
}
