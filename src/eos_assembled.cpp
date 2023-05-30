//! \file eos_assembled.cpp
//  \brief Implementation of EOS_assembled

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <iomanip>

#include <eos_assembled.hpp>
#include <eos_photons.hpp>

using namespace std;
using namespace radiation;

//EOS_assembled::EOS_assembled():
//  m_id_log_nb(numeric_limits<double>::quiet_NaN()),
//  m_bar_id_log_t(numeric_limits<double>::quiet_NaN()),
//  m_id_yq(numeric_limits<double>::quiet_NaN()),
//  m_id_log_nl(numeric_limits<double>::quiet_NaN()),
//  m_lep_id_log_t(numeric_limits<double>::quiet_NaN()),
//  m_nn(0), m_bar_nt(0), m_ny(0), m_nl(0), m_lep_nt(0),
//  m_min_h(numeric_limits<double>::max()),
//  m_log_nb(nullptr),
//  m_bar_log_t(nullptr),
//  m_yq(nullptr),
//  m_log_nl(nullptr),
//  m_lep_log_t(nullptr),
//  m_bar_table(nullptr),
//  m_lep_table(nullptr),
//  m_bar_initialized(false),
//  m_lep_initialized(false) {
//  n_species = 1;
//  eos_units = &Nuclear;
//{ };

//EOS_assembled::~EOS_assembled() {
//  if (m_bar_initialized) {
//    delete[] m_log_nb;
//    delete[] m_bar_log_t;
//    delete[] m_yq;
//    delete[] m_bar_table;
//  }
//  if (m_lep_initialized) {
//    delete[] m_log_nl;
//    delete[] m_lep_log_t;
//    delete[] m_lep_table;
//  }
//}

double EOS_assembled::TemperatureFromE(double n, double e, double *Y) {
  return temperature_from_e(log(e), n, Y);
}

double EOS_assembled::TemperatureFromP(double n, double p, double *Y) {
  return temperature_from_p(log(p), n, Y);
}

double EOS_assembled::Energy(double n, double T, double *Y) {
  return BarEnergy(n, T, Y) + EOS_leptons<0>::LepEnergy<id_test>(n, T, Y) + EOS_leptons<1>::LepEnergy<id_test>(n, T, Y) + RadEnergy(T);
}

double EOS_assembled::Pressure(double n, double T, double *Y) {
  return BarPressure(n, T, Y) + EOS_leptons<0>::LepPressure<id_test>(n, T, Y) + EOS_leptons<1>::LepPressure<id_test>(n, T, Y) + RadPressure(T);
}

double EOS_assembled::Entropy(double n, double T, double *Y) {
  return BarEntropy(n, T, Y) + (EOS_leptons<0>::LepEntropy<id_test>(n, T, Y) + EOS_leptons<1>::LepEntropy<id_test>(n, T, Y) + RadEntropy(T)) / n;
}

double EOS_assembled::Enthalpy(double n, double T, double *Y) {
  double const P = Pressure(n, T, Y);
  double const e = Energy(n, T, Y);
  return (P + e)/n;
}

double EOS_assembled::SoundSpeed(double n, double T, double *Y) {

  double const dPdn = BardPdn(n, T, Y) + EOS_leptons<0>::LdPdn<id_test>(n, T, Y) + EOS_leptons<1>::LdPdn<id_test>(n, T, Y);
  double const dsdn = Bardsdn(n, T, Y) + EOS_leptons<0>::Ldsdn<id_test>(n, T, Y) + EOS_leptons<1>::Ldsdn<id_test>(n, T, Y) - (RadEntropy(T) / (n*n));
  double const dPdt = BardPdT(n, T, Y) + EOS_leptons<0>::LdPdt<id_test>(n, T, Y) + EOS_leptons<1>::LdPdt<id_test>(n, T, Y) +  RaddPdT(T);
  double const dsdt = BardsdT(n, T, Y) + EOS_leptons<0>::Ldsdt<id_test>(n, T, Y) + EOS_leptons<1>::Ldsdt<id_test>(n, T, Y) + (RaddsdT(T)/n);

  double const cs2 = sqrt(std::max(1.e-6,(dPdn - dsdn/dsdt*dPdt) / Enthalpy(n, T, Y)));
  //if (cs2 >= 1.) cout << "cs2 > 1!" << endl;
  //cout << "dPdn = " << dPdn << "\t" << "dsdn = " << dsdn << "\t" << "dPdt = " << dPdt << "\t" << "dsdt = " << dsdt << endl << endl; 
  return cs2;
}


double EOS_assembled::temperature_from_e(double var, double n, double *Y) { //const {
  const int iv = EOS_baryons::BLOGE;
  int in, iy;
  double wn0, wn1, wy0, wy1;

  EOS_baryons::weight_idx_ln(&wn0, &wn1, &in, log(n));
  EOS_baryons::weight_idx_yq(&wy0, &wy1, &iy, Y[0]+Y[1]);

  auto f = [=](int it){
    double var_pt =
      log(exp(wn0 * (wy0 * m_table[EOS_baryons::index(iv, in+0, iy+0, it)]  +
                     wy1 * m_table[EOS_baryons::index(iv, in+0, iy+1, it)]) +
              wn1 * (wy0 * m_table[EOS_baryons::index(iv, in+1, iy+0, it)]  +
                     wy1 * m_table[EOS_baryons::index(iv, in+1, iy+1, it)])) +
          EOS_leptons<0>::LepEnergy<id_test>(n, exp(EOS_baryons::m_log_t[it]), Y) + EOS_leptons<1>::LepEnergy<id_test>(n, exp(EOS_baryons::m_log_t[it]), Y) +  RadEnergy(exp(EOS_baryons::m_log_t[it])));

    return var - var_pt;
  };

  int ilo = 0;
  int ihi = m_nt-1;
  double flo = f(ilo);
  double fhi = f(ihi);
  assert(flo*fhi <= 0);
  while (ihi - ilo > 1) {
    int ip = ilo + (ihi - ilo)/2;
    double fp = f(ip);
    if (fp*flo <= 0) {
      ihi = ip;
      fhi = fp;
    }
    else {
      ilo = ip;
      flo = fp;
    }
  }
  assert(ihi - ilo == 1);
  double lthi = m_log_t[ihi];
  double ltlo = m_log_t[ilo];

  if (flo == 0) {
    return exp(ltlo);
  }
  if (fhi == 0) {
    return exp(lthi);
  }

  double lt = m_log_t[ilo] - flo*(lthi - ltlo)/(fhi - flo);
  return exp(lt);
}


double EOS_assembled::temperature_from_p(double var, double n, double *Y) { //const {
  const int iv = EOS_baryons::BLOGP;
  int in, iy;
  double wn0, wn1, wy0, wy1;

  EOS_baryons::weight_idx_ln(&wn0, &wn1, &in, log(n));
  EOS_baryons::weight_idx_yq(&wy0, &wy1, &iy, Y[0]+Y[1]);

  auto f = [=](int it){
    double temp = exp(EOS_baryons::m_log_t[it]);
    double var_pt = exp(wn0 * (wy0 * EOS_baryons::m_table[EOS_baryons::index(iv, in+0, iy+0, it)]   +
                               wy1 * EOS_baryons::m_table[EOS_baryons::index(iv, in+0, iy+1, it)])  +
                        wn1 * (wy0 * EOS_baryons::m_table[EOS_baryons::index(iv, in+1, iy+0, it)]   +
                               wy1 * EOS_baryons::m_table[EOS_baryons::index(iv, in+1, iy+1, it)])) -
	            Pmin + EOS_leptons<0>::LepPressure<id_test>(n, temp, Y) + EOS_leptons<1>::LepPressure<id_test>(n, temp, Y) + RadPressure(temp);
    assert(var_pt > 0.);  
    return var - log(var_pt);
  };

  int ilo = 0;
  int ihi = m_nt-1;
  double flo = f(ilo);
  double fhi = f(ihi);
  //if (!(flo*fhi <= 0)) {
	//cout << iv << "\t" << ilo << "\t" << ihi << endl;
  	//cout << in << "\t" << iy << endl;
  	//cout << wn0 << "\t" << wn1 << "\t" << wy0 << "\t" << wy1 << endl;
  	//cout << "Number density = " << n << ", Electron fraction = " << Y[0] << ": flo = " << var-flo << ", fhi = " << var-fhi << endl;
  	//cout << var << endl;
	//cout << "Baryon flo = " << f1 << endl;;
  	//assert(flo*fhi <= 0);
	//return 0;
  //}
  assert(flo*fhi <= 0);
  while (ihi - ilo > 1) {
    int ip = ilo + (ihi - ilo)/2;
    double fp = f(ip);
    if (fp*flo <= 0) {
      ihi = ip;
      fhi = fp;
    }
    else {
      ilo = ip;
      flo = fp;
    }
  }
  assert(ihi - ilo == 1);
  double lthi = m_log_t[ihi];
  double ltlo = m_log_t[ilo];

  if (flo == 0) {
    return exp(ltlo);
  }
  if (fhi == 0) {
    return exp(lthi);
  }

  double lt = m_log_t[ilo] - flo*(lthi - ltlo)/(fhi - flo);
  return exp(lt);
}

