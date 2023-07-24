#include <vector>
#include <array>
#include <sstream>

#include "helmholtz_eos.hpp"
#include "../constants.hpp"
#include "../fermi_integrals/fermi_integrals.hpp"
#include "../num_tools/root_finding/root_finding.hpp"

const double fivethirds = 5./3.;
const double fourthirds = 4./3.;

HelmEOSOutput eos_helm_from_eta(const double eta, const double T, const int id_L, GFDs *FD) {
  // Define theta (relativity parameter)
  const double theta = T / mL[id_L];

  // Degeneracy parameter of anti-leptons
  const double a_eta = - (eta + 2. / theta);
 
  // Compute missing Generalized Fermi-Dirac integrals
  FD->f52   = compute_res(eta,   theta, 2.5); // k = 5/2
  FD->a_f52 = compute_res(a_eta, theta, 2.5); // k = 5/2
  
  HelmEOSOutput out;

  const double K_theta32 = K[id_L]  * pow(theta, 1.5);
  const double K_mL_theta32 = K_theta32 * mL[id_L];
  const double K3_mL_theta52 = K3[id_L] * mL[id_L] * pow(theta, 2.5);
  const double thetasq = theta * theta;

  // Number density of leptons
  out.nl = K_theta32 * (FD->f12 + theta * FD->f32);

  // Number density of anti-leptons
  out.a_nl = K_theta32 * (FD->a_f12 + theta * FD->a_f32);

  // Pressure of leptons
  out.pl = K3_mL_theta52 * (2. * FD->f32 + theta * FD->f52); //MeV/fm^3

  // Pressure of anti-leptons
  out.a_pl = K3_mL_theta52 * (2. * FD->a_f32 + theta * FD->a_f52); //MeV/fm^3

  // @TODO: compute 2*theta only once
  // Internal energy density of leptons
  out.el = K_mL_theta32 * (FD->f12 + 2. * theta * FD->f32 + thetasq * FD->f52); //MeV/fm^3

  // Internal energy density of anti-leptons
  out.a_el = K_mL_theta32 * (FD->a_f12 + 2. * theta * FD->a_f32 + thetasq * FD->a_f52); //MeV/fm^3
  
  // @TODO: compute 4./3.*theta only once
  // Entropy density of leptons
  out.sl = K_theta32 * (- eta * FD->f12 + (fivethirds - eta * theta) * FD->f32 + fourthirds * theta * FD->f52); //1/fm^3

  // Entropy density of anti-leptons
  out.a_sl = K_theta32 * (- a_eta * FD->a_f12 + (fivethirds - a_eta * theta) * FD->a_f32 + fourthirds * theta * FD->a_f52); //1/fm^3

  // Chemical potential of leptons
  out.mul = mL[id_L] + eta * T;

  return out;
}

HelmEOSOutput eos_helm_full(double nLep, double temp, const int id_L) {
  GFDs FD_integ;
  const double eta = rtsafe(1.e39*nLep, temp, id_L, &FD_integ);

  return eos_helm_from_eta(eta, temp, id_L, &FD_integ);
}

HelmEOSDer der_cs2(double nLep, double temp, const int id_L) {
  // Define theta (relativity parameter)
  const double theta = temp / mL[id_L];

  GFDs FD_integ;
  const double eta = rtsafe(1.e39*nLep, temp, id_L, &FD_integ);

  // @TODO: optimize computation of Fermi-Dirac integrals here

  // Compute missing Generalized Fermi-Dirac integrals
  const double f12 = FD_integ.f12;
  const double f32 = FD_integ.f32;
  const double f52 = compute_res(eta,   theta, 2.5); // k = 5/2
 
  const double f12_dn = compute_res_ed(eta, theta, 0.5);
  const double f32_dn = compute_res_ed(eta, theta, 1.5);

  const double f12_dT = (f32_dn - 1.5*f12) / theta;
  const double f32_dT = (f32 - 4.*f12_dT) / (2.*theta);
  const double f52_dT = (f52 - 4.*f32_dT) / (2.*theta);
  const double f52_dn = theta*f32_dT + 2.5*f32;

  // Degeneracy parameter of anti-leptons
  const double a_eta = - (eta + 2. / theta);

  const double a_f12 = FD_integ.a_f12;
  const double a_f32 = FD_integ.a_f32;
  const double a_f52 = compute_res(a_eta, theta, 2.5); // k = 5/2

  const double a_f12_dn = compute_res_ed(a_eta, theta, 0.5);
  const double a_f32_dn = compute_res_ed(a_eta, theta, 1.5);

  const double a_f12_dT = (a_f32_dn - 1.5*a_f12) / theta;
  const double a_f32_dT = (a_f32 - 4.*a_f12_dT) / (2.*theta);
  const double a_f52_dT = (a_f52 - 4.*a_f32_dT) / (2.*theta);
  const double a_f52_dn = theta*a_f32_dT + 2.5*a_f32;

  const double dn = f12_dn+a_f12_dn + theta*(f32_dn+a_f32_dn);

  HelmEOSDer out;

  // @TODO: optimize calculation of constants
  out.dPdn = mL[id_L]/3. * theta * (2.*(f32_dn-a_f32_dn) + theta*(f52_dn-a_f52_dn)) / dn;

  out.dsdn = (-f12+a_f12-eta*f12_dn+a_eta*a_f12_dn+5./3.*(f32_dn-a_f32_dn) - theta*(f32-a_f32+eta*f32_dn-a_eta*a_f32_dn-4./3.*(f52_dn-a_f52_dn))) / dn;

  out.dPdt = K3[id_L] * pow(theta,1.5) * (5.*(f32+a_f32) + theta*(3.5*(f52+a_f52)+2.*(f32_dT+a_f32_dT)) + theta*theta*(f52_dT+a_f52_dT));
        
  out.dsdt = K[id_L]/mL[id_L] * pow(theta,0.5) * (-1.5*(eta*f12+a_eta*a_f12) + 2.5*(f32+a_f32) + theta*(-2.5*(eta*f32+a_eta*a_f32) + 10./3.*(f52+a_f52) - (eta*f12_dT+a_eta*a_f12_dT) + 5./3.*(f32_dT+a_f32_dT)) + theta*theta * (-(eta*f32_dT+a_eta*a_f32_dT) + 4./3.*(f52_dT+a_f52_dT)));

  return out;
}


HelmEOSDer der_cs2_num(double nLep, double temp, double id_L) {
  GFDs FD_integ;
  const double eta = rtsafe(1.e39*nLep, temp, id_L, &FD_integ);
 
  HelmEOSOutput tmp = eos_helm_from_eta(eta, temp, id_L, &FD_integ);
  
  const double P = tmp.pl + tmp.a_pl;
  const double s = tmp.sl + tmp.a_sl;

  const double eps_t = temp*0.005;
  const double eps_n = nLep*0.005;

  const double temp_1 = temp - eps_t;
  const double temp_2 = temp + eps_t;

  const double nLep_1 = nLep - eps_n;
  const double nLep_2 = nLep + eps_n;

  // N.B.: order here is important because of FD_integ
  double eta_1 = rtsafe(1.e39*nLep_1, temp, id_L, &FD_integ);
  HelmEOSOutput tmp_1 = eos_helm_from_eta(eta_1, temp, id_L, &FD_integ);

  double eta_2 = rtsafe(1.e39*nLep_2, temp, id_L, &FD_integ);
  HelmEOSOutput tmp_2 = eos_helm_from_eta(eta_2, temp, id_L, &FD_integ);

  HelmEOSDer out;

  // @TODO: decide between different kind of derivatives (lin-lin, lin-log, log-log)
  out.dPdn = P * (log(tmp_2.pl + tmp_2.a_pl) - log(tmp_1.pl + tmp_1.a_pl)) / (nLep_2-nLep_1);
  out.dsdn = s * (log(tmp_2.sl + tmp_2.a_sl) - log(tmp_1.sl + tmp_1.a_sl)) / (nLep_2-nLep_1);

  eta_1 = rtsafe(1.e39*nLep, temp_1, id_L, &FD_integ);
  tmp_1 = eos_helm_from_eta(eta, temp_1, id_L, &FD_integ);

  eta_2 = rtsafe(1.e39*nLep, temp_2, id_L, &FD_integ);
  tmp_2 = eos_helm_from_eta(eta, temp_2, id_L, &FD_integ);

  out.dPdt = P * (log(tmp_2.pl + tmp_2.a_pl) - log(tmp_1.pl + tmp_1.a_pl)) / (temp_2-temp_1);
  out.dsdt = s * (log(tmp_2.sl + tmp_2.a_sl) - log(tmp_1.sl + tmp_1.a_sl)) / (temp_2-temp_1);

  return out;
}