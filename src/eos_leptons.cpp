//! \file eos_leptons.cpp
//  \brief Implementation of EOS_leptons

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <fstream>

#include <eos_leptons.hpp>

using namespace std;

#define muon_check(check, value) (check ? value : 0)

EOS_leptons::EOS_leptons():
  m_id_log_ne(numeric_limits<double>::quiet_NaN()),
  m_id_log_te(numeric_limits<double>::quiet_NaN()),
  m_id_log_nm(numeric_limits<double>::quiet_NaN()),
  m_id_log_tm(numeric_limits<double>::quiet_NaN()),
  m_ne(0), m_nte(0), m_nm(0), m_ntm(0),
  m_log_ne(nullptr),
  m_log_te(nullptr),
  m_log_nm(nullptr),
  m_log_tm(nullptr),
  m_el_table(nullptr),
  m_mu_table(nullptr),
  m_el_initialized(false),
  m_mu_initialized(false) {};

EOS_leptons::~EOS_leptons() {
  if (m_el_initialized) {
    delete[] m_log_ne;
    delete[] m_log_te;
    delete[] m_el_table;
  }
  if (m_mu_initialized) {
    delete[] m_log_nm;
    delete[] m_log_tm;
    delete[] m_mu_table;
  }
}

double EOS_leptons::mu_check(const bool check, const double value) {
	if (check) {
		return value;
	} else {
		return 0.;
	}
}

double EOS_leptons::ElectronNumberDensity(double n, double T, double *Y) {
  assert(m_el_initialized);
  return exp(eval_el_at_nty(IL_N, n, T, Y[id_e]));
}

double EOS_leptons::PositronNumberDensity(double n, double T, double *Y) {
  assert(m_el_initialized);
  //return exp(eval_el_at_nty(IA_N, n, T, Y));
  return ElectronNumberDensity(n, T, Y) - n*Y[id_e];
}

double EOS_leptons::MuonNumberDensity(double n, double T, double *Y) {
  return mu_check(m_mu_initialized, exp(eval_mu_at_nty(IL_N, n, T, Y[id_mu])));
  //assert(m_mu_initialized);
  //return exp(eval_mu_at_nty(IL_N, n, T, Y[id_mu]));
}

double EOS_leptons::AntiMuonNumberDensity(double n, double T, double *Y) {
  return mu_check(m_mu_initialized, MuonNumberDensity(n, T, Y) - n*Y[id_mu]);
  //assert(m_mu_initialized);
  //return MuonNumberDensity(n, T, Y) - n*Y[id_mu];
}

double EOS_leptons::ElectronEnergy(double n, double T, double Y) {
  assert(m_el_initialized);
  return exp(eval_el_at_nty(IL_E, n, T, Y)) + exp(eval_el_at_nty(IA_E, n, T, Y));
}

double EOS_leptons::MuonEnergy(double n, double T, double Y) {
  //muon_check(m_mu_initialized,exp(eval_mu_at_nty(IL_E, n, T, Y)) + exp(eval_mu_at_nty(IA_E, n, T, Y)));
  //assert(m_mu_initialized);
  return mu_check(m_mu_initialized,exp(eval_mu_at_nty(IL_E, n, T, Y)) + exp(eval_mu_at_nty(IA_E, n, T, Y)));
}

double EOS_leptons::LepEnergy(double n, double T, double *Y) {
  return ElectronEnergy(n, T, Y[id_e]) + MuonEnergy(n, T, Y[id_mu]);
}

double EOS_leptons::ElectronPressure(double n, double T, double *Y) {
  assert(m_el_initialized);
  return exp(eval_el_at_nty(IL_P, n, T, Y[id_e])) + exp(eval_el_at_nty(IA_P, n, T, Y[id_e]));
}

double EOS_leptons::MuonPressure(double n, double T, double *Y) {
  return mu_check(m_mu_initialized,exp(eval_mu_at_nty(IL_P, n, T, Y[id_mu])) + exp(eval_mu_at_nty(IA_P, n, T, Y[id_mu])));
  //assert(m_mu_initialized);
  //return exp(eval_mu_at_nty(IL_P, n, T, Y[id_mu])) + exp(eval_mu_at_nty(IA_P, n, T, Y[id_mu]));
}

double EOS_leptons::LepPressure(double n, double T, double *Y) {
  return ElectronPressure(n, T, Y) + MuonPressure(n, T, Y);
}

double EOS_leptons::ElectronEntropy(double n, double T, double *Y) {
  assert(m_el_initialized);
  return exp(eval_el_at_nty(IL_S, n, T, Y[id_e])) + exp(eval_el_at_nty(IA_S, n, T, Y[id_e]));
}

double EOS_leptons::MuonEntropy(double n, double T, double *Y) {
  return mu_check(m_mu_initialized,exp(eval_mu_at_nty(IL_S, n, T, Y[id_mu])) + exp(eval_mu_at_nty(IA_S, n, T, Y[id_mu])));
  //assert(m_mu_initialized);
  //return exp(eval_mu_at_nty(IL_S, n, T, Y[id_mu])) + exp(eval_mu_at_nty(IA_S, n, T, Y[id_mu]));
}

double EOS_leptons::LepEntropy(double n, double T, double *Y) {
  return ElectronEntropy(n, T, Y) + MuonEntropy(n, T, Y);
}

double EOS_leptons::ElectronChemicalPotential(double n, double T, double *Y) {
  assert(m_el_initialized);
  return eval_el_at_nty(IL_MU, n, T, Y[id_e]);
}

double EOS_leptons::MuonChemicalPotential(double n, double T, double *Y) {
  return mu_check(m_mu_initialized,eval_mu_at_nty(IL_MU, n, T, Y[id_mu]));
  //assert(m_mu_initialized);
  //return eval_mu_at_nty(IL_MU, n, T, Y[id_mu]);
}

double EOS_leptons::EdPdn(double n, double T, double *Y) {
  assert(m_el_initialized);
  return eval_el_at_nty(ID_DPDN, n, T, Y[id_e])*Y[id_e];
}

double EOS_leptons::Edsdn(double n, double T, double *Y) {
  assert(m_el_initialized);
  return eval_el_at_nty(ID_DSDN, n, T, Y[id_e]) * Y[id_e] / n - ElectronEntropy(n, T, Y) / (n*n);
}

double EOS_leptons::EdPdt(double n, double T, double *Y) {
  assert(m_el_initialized);
  return eval_el_at_nty(ID_DPDT, n, T, Y[id_e]);
}

double EOS_leptons::Edsdt(double n, double T, double *Y) {
  assert(m_el_initialized);
  return eval_el_at_nty(ID_DSDT, n, T, Y[id_e]);
}

double EOS_leptons::MdPdn(double n, double T, double *Y) {
  return mu_check(m_mu_initialized,eval_mu_at_nty(ID_DPDN, n, T, Y[id_mu])*Y[id_mu]);
  //assert(m_mu_initialized);
  //return eval_mu_at_nty(ID_DPDN, n, T, Y[id_mu]);
}

double EOS_leptons::Mdsdn(double n, double T, double *Y) {
  return mu_check(m_mu_initialized,eval_mu_at_nty(ID_DSDN, n, T, Y[id_mu]) * Y[id_mu] / n - MuonEntropy(n, T, Y) / (n*n));
  //assert(m_mu_initialized);
  //return eval_mu_at_nty(ID_DSDN, n, T, Y[id_mu]);
}

double EOS_leptons::MdPdt(double n, double T, double *Y) {
  return mu_check(m_mu_initialized,eval_mu_at_nty(ID_DPDT, n, T, Y[id_mu]));
  //assert(m_mu_initialized);
  //return eval_mu_at_nty(ID_DPDT, n, T, Y[id_mu]);
}

double EOS_leptons::Mdsdt(double n, double T, double *Y) {
  return mu_check(m_mu_initialized, eval_mu_at_nty(ID_DSDT, n, T, Y[id_mu]));
  //assert(m_mu_initialized);
  //return eval_mu_at_nty(ID_DSDT, n, T, Y[id_mu]);
}

std::array<double,EOS_leptons::NVARS> EOS_leptons::ComputeFullElectronEOS(double n, double T, double *Y) {
  return eval_all_el_at_lnt(log(n*Y[id_e]),log(T));
}

std::array<double,EOS_leptons::NVARS> EOS_leptons::ComputeFullMuonEOS(double n, double T, double *Y) {
  assert(m_mu_initialized);
  return eval_all_mu_at_lnt(log(n*Y[id_mu]),log(T));
}

void EOS_leptons::ReadETableFromFile(std::string fname) {
  // Open input file
  // -------------------------------------------------------------------------
  ifstream EOSinput;
  EOSinput.open(fname);
  
  if (!EOSinput) {
    cout << "Electron EOS table not found!" << endl;
    exit(EXIT_FAILURE);
  }

  // Get dataset sizes
  // -------------------------------------------------------------------------
  string EOSline;
  getline(EOSinput, EOSline);
  stringstream s1(EOSline);
  s1 >> m_ne;

  getline(EOSinput, EOSline);
  stringstream s2(EOSline);
  s2 >> m_nte;

  // Allocate memory
  // -------------------------------------------------------------------------
  m_log_ne   = new double[m_ne];
  m_log_te   = new double[m_nte];
  m_el_table = new double[NTOT*m_ne*m_nte];

  double number;

  // Read nb, t, yq
  // -------------------------------------------------------------------------
  getline(EOSinput, EOSline);
  stringstream s3(EOSline);
  for (int in = 0; in < m_ne; ++in) {
    s3 >> number;
    m_log_ne[in] = log(pow(10.,number));
  }
  m_id_log_ne = 1.0/(m_log_ne[1] - m_log_ne[0]);

  getline(EOSinput, EOSline);
  stringstream s4(EOSline);
  for (int it = 0; it < m_nte; ++it) {
    s4 >> number;
    m_log_te[it] = log(pow(10.,number));
  }
  m_id_log_te = 1.0/(m_log_te[1] - m_log_te[0]);

  // Read other thermodynamics quantities
  // -------------------------------------------------------------------------
  for (int in = 0; in < m_ne; ++in) {
    getline(EOSinput, EOSline);
    stringstream s5(EOSline);
      for (int it = 0; it < m_nte; ++it) {
        s5 >> number;
        m_el_table[el_index(IL_MU, in, it)] = log(number);
  }}

  for (int in = 0; in < m_ne; ++in) {
    getline(EOSinput, EOSline);
    stringstream s5(EOSline);
      for (int it = 0; it < m_nte; ++it) {
        s5 >> number;
        m_el_table[el_index(IL_N, in, it)] = log(number);
  }}

  for (int in = 0; in < m_ne; ++in) {
    getline(EOSinput, EOSline);
    stringstream s5(EOSline);
      for (int it = 0; it < m_nte; ++it) {
        s5 >> number;
        m_el_table[el_index(IA_N, in, it)] = log(number);
  }}

  for (int in = 0; in < m_ne; ++in) {
    getline(EOSinput, EOSline);
    stringstream s5(EOSline);
      for (int it = 0; it < m_nte; ++it) {
        s5 >> number;
        m_el_table[el_index(IL_P, in, it)] = log(number);
  }}

  for (int in = 0; in < m_ne; ++in) {
    getline(EOSinput, EOSline);
    stringstream s5(EOSline);
      for (int it = 0; it < m_nte; ++it) {
        s5 >> number;
        m_el_table[el_index(IA_P, in, it)] = log(number);
  }}

  for (int in = 0; in < m_ne; ++in) {
    getline(EOSinput, EOSline);
    stringstream s5(EOSline);
      for (int it = 0; it < m_nte; ++it) {
        s5 >> number;
        m_el_table[el_index(IL_E, in, it)] = log(number);
  }}

  for (int in = 0; in < m_ne; ++in) {
    getline(EOSinput, EOSline);
    stringstream s5(EOSline);
      for (int it = 0; it < m_nte; ++it) {
        s5 >> number;
        m_el_table[el_index(IA_E, in, it)] = log(number);
  }}

  for (int in = 0; in < m_ne; ++in) {
    getline(EOSinput, EOSline);
    stringstream s5(EOSline);
      for (int it = 0; it < m_nte; ++it) {
        s5 >> number;
        m_el_table[el_index(IL_S, in, it)] = log(number);
  }}
  
  for (int in = 0; in < m_ne; ++in) {
    getline(EOSinput, EOSline);
    stringstream s5(EOSline);
      for (int it = 0; it < m_nte; ++it) {
        s5 >> number;
        m_el_table[el_index(IA_S, in, it)] = number;
  }}

  for (int in = 0; in < m_ne; ++in) {
    getline(EOSinput, EOSline);
    stringstream s5(EOSline);
      for (int it = 0; it < m_nte; ++it) {
        s5 >> number;
        m_el_table[el_index(ID_DPDN, in, it)] = number;
  }}

  for (int in = 0; in < m_ne; ++in) {
    getline(EOSinput, EOSline);
    stringstream s5(EOSline);
      for (int it = 0; it < m_nte; ++it) {
        s5 >> number;
        m_el_table[el_index(ID_DSDN, in, it)] = number;
  }}
  
  for (int in = 0; in < m_ne; ++in) {
    getline(EOSinput, EOSline);
    stringstream s5(EOSline);
      for (int it = 0; it < m_nte; ++it) {
        s5 >> number;
        m_el_table[el_index(ID_DPDT, in, it)] = number;
  }}

  for (int in = 0; in < m_ne; ++in) {
    getline(EOSinput, EOSline);
    stringstream s5(EOSline);
      for (int it = 0; it < m_nte; ++it) {
        s5 >> number;
        m_el_table[el_index(ID_DSDT, in, it)] = number;
  }}

  // Cleanup
  // -------------------------------------------------------------------------
  //delete[] s5;
  EOSinput.close();

  m_el_initialized = true;
}


void EOS_leptons::ReadMTableFromFile(std::string fname) {
  // Open input file
  // -------------------------------------------------------------------------
  ifstream EOSinput;
  EOSinput.open(fname);
  
  if (!EOSinput) {
    cout << "Muon EOS table not found!" << endl;
    exit(EXIT_FAILURE);
  }

  // Get dataset sizes
  // -------------------------------------------------------------------------
  string EOSline;
  getline(EOSinput, EOSline);
  stringstream s1(EOSline);
  s1 >> m_nm;

  getline(EOSinput, EOSline);
  stringstream s2(EOSline);
  s2 >> m_ntm;

  // Allocate memory
  // -------------------------------------------------------------------------
  m_log_nm   = new double[m_nm];
  m_log_tm   = new double[m_ntm];
  m_mu_table = new double[NTOT*m_nm*m_ntm];


  double number;

  // Read nb, t, yq
  // -------------------------------------------------------------------------
  getline(EOSinput, EOSline);
  stringstream s3(EOSline);
  for (int in = 0; in < m_nm; ++in) {
    s3 >> number;
    m_log_nm[in] = log(pow(10.,number));
  }
  m_id_log_nm = 1.0/(m_log_nm[1] - m_log_nm[0]);

  getline(EOSinput, EOSline);
  stringstream s4(EOSline);
  for (int it = 0; it < m_ntm; ++it) {
    s4 >> number;
    m_log_tm[it] = log(pow(10.,number));
  }
  m_id_log_tm = 1.0/(m_log_tm[1] - m_log_tm[0]);

  // Read other thermodynamics quantities
  // -------------------------------------------------------------------------
  for (int in = 0; in < m_nm; ++in) {
    getline(EOSinput, EOSline);
    stringstream s5(EOSline);
      for (int it = 0; it < m_ntm; ++it) {
        s5 >> number;
        m_mu_table[mu_index(IL_MU, in, it)] = log(number);
  }}

  for (int in = 0; in < m_nm; ++in) {
    getline(EOSinput, EOSline);
    stringstream s5(EOSline);
      for (int it = 0; it < m_ntm; ++it) {
        s5 >> number;
        m_mu_table[mu_index(IL_N, in, it)] = log(number);
  }}

  for (int in = 0; in < m_nm; ++in) {
    getline(EOSinput, EOSline);
    stringstream s5(EOSline);
      for (int it = 0; it < m_ntm; ++it) {
        s5 >> number;
        m_mu_table[mu_index(IA_N, in, it)] = log(number);
  }}

  for (int in = 0; in < m_nm; ++in) {
    getline(EOSinput, EOSline);
    stringstream s5(EOSline);
      for (int it = 0; it < m_ntm; ++it) {
        s5 >> number;
        m_mu_table[mu_index(IL_P, in, it)] = log(number);
  }}

  for (int in = 0; in < m_nm; ++in) {
    getline(EOSinput, EOSline);
    stringstream s5(EOSline);
      for (int it = 0; it < m_ntm; ++it) {
        s5 >> number;
        m_mu_table[mu_index(IA_P, in, it)] = log(number);
  }}

  for (int in = 0; in < m_nm; ++in) {
    getline(EOSinput, EOSline);
    stringstream s5(EOSline);
      for (int it = 0; it < m_ntm; ++it) {
        s5 >> number;
        m_mu_table[mu_index(IL_E, in, it)] = log(number);
  }}

  for (int in = 0; in < m_nm; ++in) {
    getline(EOSinput, EOSline);
    stringstream s5(EOSline);
      for (int it = 0; it < m_ntm; ++it) {
        s5 >> number;
        m_mu_table[mu_index(IA_E, in, it)] = log(number);
  }}

  for (int in = 0; in < m_nm; ++in) {
    getline(EOSinput, EOSline);
    stringstream s5(EOSline);
      for (int it = 0; it < m_ntm; ++it) {
        s5 >> number;
        m_mu_table[mu_index(IL_S, in, it)] = log(number);
  }}
  
  for (int in = 0; in < m_nm; ++in) {
    getline(EOSinput, EOSline);
    stringstream s5(EOSline);
      for (int it = 0; it < m_ntm; ++it) {
        s5 >> number;
        m_mu_table[mu_index(IA_S, in, it)] = number;
  }}

  for (int in = 0; in < m_nm; ++in) {
    getline(EOSinput, EOSline);
    stringstream s5(EOSline);
      for (int it = 0; it < m_ntm; ++it) {
        s5 >> number;
        m_mu_table[mu_index(ID_DPDN, in, it)] = number;
  }}

  for (int in = 0; in < m_nm; ++in) {
    getline(EOSinput, EOSline);
    stringstream s5(EOSline);
      for (int it = 0; it < m_ntm; ++it) {
        s5 >> number;
        m_mu_table[mu_index(ID_DSDN, in, it)] = number;
  }}
  
  for (int in = 0; in < m_nm; ++in) {
    getline(EOSinput, EOSline);
    stringstream s5(EOSline);
      for (int it = 0; it < m_ntm; ++it) {
        s5 >> number;
        m_mu_table[mu_index(ID_DPDT, in, it)] = number;
  }}

  for (int in = 0; in < m_nm; ++in) {
    getline(EOSinput, EOSline);
    stringstream s5(EOSline);
      for (int it = 0; it < m_ntm; ++it) {
        s5 >> number;
        m_mu_table[mu_index(ID_DSDT, in, it)] = number;
  }}

  // Cleanup
  // -------------------------------------------------------------------------
  //delete[] s5;
  EOSinput.close();

  m_mu_initialized = true;
}

double EOS_leptons::eval_el_at_nty(int vi, double n, double T, double Yl) const {
  return eval_el_at_lnty(vi, log(n*Yl), log(T));
}

double EOS_leptons::eval_mu_at_nty(int vi, double n, double T, double Yl) const {
  return eval_mu_at_lnty(vi, log(n*Yl), log(T));
}


void EOS_leptons::weight_idx_lne(double *w0, double *w1, int *in, double log_n) const {
  *in = (log_n - m_log_ne[0])*m_id_log_ne;
  *w1 = (log_n - m_log_ne[*in])*m_id_log_ne;
  *w0 = 1.0 - (*w1);
}

void EOS_leptons::weight_idx_lnm(double *w0, double *w1, int *in, double log_n) const {
  *in = (log_n - m_log_nm[0])*m_id_log_nm;
  *w1 = (log_n - m_log_nm[*in])*m_id_log_nm;
  *w0 = 1.0 - (*w1);
}

void EOS_leptons::weight_idx_lte(double *w0, double *w1, int *it, double log_t) const {
  *it = (log_t - m_log_te[0])*m_id_log_te;
  *w1 = (log_t - m_log_te[*it])*m_id_log_te;
  *w0 = 1.0 - (*w1);
}

void EOS_leptons::weight_idx_ltm(double *w0, double *w1, int *it, double log_t) const {
  *it = (log_t - m_log_tm[0])*m_id_log_tm;
  *w1 = (log_t - m_log_tm[*it])*m_id_log_tm;
  *w0 = 1.0 - (*w1);
}

double EOS_leptons::eval_el_at_lnty(int iv, double log_n, double log_t) const {
  int in, it;
  double wn0, wn1, wt0, wt1;

  weight_idx_lne(&wn0, &wn1, &in, log_n);
  weight_idx_lte(&wt0, &wt1, &it, log_t);

  return
    wn0 * (wt0 * m_el_table[el_index(iv, in+0, it+0)]   +
           wt1 * m_el_table[el_index(iv, in+0, it+1)])  +
    wn1 * (wt0 * m_el_table[el_index(iv, in+1, it+0)]   +
           wt1 * m_el_table[el_index(iv, in+1, it+1)]);
}

double EOS_leptons::eval_mu_at_lnty(int iv, double log_n, double log_t) const {
  int in, it;
  double wn0, wn1, wt0, wt1;

  weight_idx_lnm(&wn0, &wn1, &in, log_n);
  weight_idx_ltm(&wt0, &wt1, &it, log_t);

  return
    wn0 * (wt0 * m_mu_table[mu_index(iv, in+0, it+0)]   +
           wt1 * m_mu_table[mu_index(iv, in+0, it+1)])  +
    wn1 * (wt0 * m_mu_table[mu_index(iv, in+1, it+0)]   +
           wt1 * m_mu_table[mu_index(iv, in+1, it+1)]);
}

std::array<double,EOS_leptons::NVARS> EOS_leptons::eval_all_el_at_lnt(double log_n, double log_t) const {
  int in, it;
  double wn0, wn1, wt0, wt1;
  std::array<double,EOS_leptons::NVARS> eos_out;

  weight_idx_lne(&wn0, &wn1, &in, log_n);
  weight_idx_lte(&wt0, &wt1, &it, log_t);
  
  for (int iv=0; iv<NVARS; iv++) {
     eos_out[iv] = wn0 * (wt0 * m_el_table[el_index(iv, in+0, it+0)]   +
                          wt1 * m_el_table[el_index(iv, in+0, it+1)])  +
                   wn1 * (wt0 * m_el_table[el_index(iv, in+1, it+0)]   +
                          wt1 * m_el_table[el_index(iv, in+1, it+1)]);
  }
  return eos_out;
}

std::array<double,EOS_leptons::NVARS> EOS_leptons::eval_all_mu_at_lnt(double log_n, double log_t) const {
  int in, it;
  double wn0, wn1, wt0, wt1;
  std::array<double,EOS_leptons::NVARS> eos_out;

  weight_idx_lnm(&wn0, &wn1, &in, log_n);
  weight_idx_ltm(&wt0, &wt1, &it, log_t);

  for (int iv=0; iv<NVARS; iv++) {
     eos_out[iv] = wn0 * (wt0 * m_mu_table[mu_index(iv, in+0, it+0)]   +
                          wt1 * m_mu_table[mu_index(iv, in+0, it+1)])  +
                   wn1 * (wt0 * m_mu_table[mu_index(iv, in+1, it+0)]   +
                          wt1 * m_mu_table[mu_index(iv, in+1, it+1)]);
  }
  return eos_out;
}



