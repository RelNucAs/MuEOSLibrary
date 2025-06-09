#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <iostream>

#include <hdf5.h>
#include <hdf5_hl.h>

#include "eos_species/eos_species.hpp"

using namespace std;

#define MYH5CHECK(ierr) \
  if(ierr < 0) { \
    stringstream ss; \
    ss << __FILE__ << ":" << __LINE__ << " error reading EOS table!"; \
    throw runtime_error(ss.str().c_str()); \
  }

EOS_baryons::EOS_baryons():
  m_id_log_nb(numeric_limits<double>::quiet_NaN()),
  m_id_log_t(numeric_limits<double>::quiet_NaN()),
  m_id_yq(numeric_limits<double>::quiet_NaN()),
  m_nn(0), m_nt(0), m_ny(0),
  Pmin(0),
  m_log_nb(nullptr),
  m_log_t(nullptr),
  m_yq(nullptr),
  m_table(nullptr),
  m_initialized(false) {};

EOS_baryons::~EOS_baryons() {
  if (m_initialized) {
    delete[] m_log_nb;
    delete[] m_log_t;
    delete[] m_yq;
    delete[] m_table;
  }
}

EOSstruct EOS_baryons::BaryonEOS(double n, double T, double *Y) {
  assert (m_initialized);

  EOSstruct eos_out = {0};

  weight_idx_ln(log(n));
  weight_idx_yq(Y[0] + Y[1]);
  weight_idx_lt(log(T));

  eos_out.mb = m_bar;
  
  eos_out.e = exp(interp_3d(BLOGE));
  eos_out.P = exp(interp_3d(BLOGP)) - Pmin;
  eos_out.s = interp_3d(BENT);
  eos_out.chem_pot[0] = interp_3d(BMUB) + interp_3d(BMUQ);
  eos_out.chem_pot[1] = interp_3d(BMUB);
  eos_out.comp[0] = interp_3d(BYALP);
  eos_out.comp[1] = interp_3d(BYNUC);
  eos_out.comp[2] = interp_3d(BYNTR);
  eos_out.comp[3] = interp_3d(BYPTN);
  eos_out.eos_der[0] = interp_3d(BDPDN);
  eos_out.eos_der[1] = interp_3d(BDSDN);
  eos_out.eos_der[2] = interp_3d(BDPDT);
  eos_out.eos_der[3] = interp_3d(BDSDT);

  return eos_out;
}

/* double EOS_baryons::BarEnergy(double n, double T, double *Y) {
  assert (m_initialized);
  return exp(eval_at_nty(BLOGE, n, T, Y[0]+Y[1]));
}

double EOS_baryons::BarPressure(double n, double T, double *Y) {
  assert (m_initialized);
  return exp(eval_at_nty(BLOGP, n, T, Y[0]+Y[1])) - Pmin;
}

double EOS_baryons::BarEntropy(double n, double T, double *Y) {
  assert (m_initialized);
  return eval_at_nty(BENT, n, T, Y[0]+Y[1]);
}

double EOS_baryons::ProtonChemicalPotential(double n, double T, double *Y) {
  assert (m_initialized);
  return eval_at_nty(BMUB, n, T, Y[0]+Y[1]) + eval_at_nty(BMUQ, n, T, Y[0]+Y[1]);
}

double EOS_baryons::NeutronChemicalPotential(double n, double T, double *Y) {
  assert (m_initialized);
  return eval_at_nty(BMUB, n, T, Y[0]+Y[1]);
}

double EOS_baryons::AlphaFraction(double n, double T, double *Y) {
  assert (m_initialized);
  return eval_at_nty(BYALP, n, T, Y[0]+Y[1]);
}

double EOS_baryons::HeavyFraction(double n, double T, double *Y) {
  assert (m_initialized);
  return eval_at_nty(BYNUC, n, T, Y[0]+Y[1]);
}

double EOS_baryons::NeutronFraction(double n, double T, double *Y) {
  assert (m_initialized);
  return eval_at_nty(BYNTR, n, T, Y[0]+Y[1]);
}

double EOS_baryons::ProtonFraction(double n, double T, double *Y) {
  assert (m_initialized);
  return eval_at_nty(BYPTN, n, T, Y[0]+Y[1]);
}

double EOS_baryons::BardPdn(double n, double T, double *Y) {
  assert (m_initialized);
  return eval_at_nty(BDPDN, n, T, Y[0]+Y[1]);
}

double EOS_baryons::Bardsdn(double n, double T, double *Y) {
  assert (m_initialized);
  return eval_at_nty(BDSDN, n, T, Y[0]+Y[1]);
}

double EOS_baryons::BardPdT(double n, double T, double *Y) {
  assert (m_initialized);
  return eval_at_nty(BDPDT, n, T, Y[0]+Y[1]);
}

double EOS_baryons::BardsdT(double n, double T, double *Y) {
  assert (m_initialized);
  return eval_at_nty(BDSDT, n, T, Y[0]+Y[1]);
}
 */

void EOS_baryons::ReadBarTableFromFile(std::string fname) {
  herr_t ierr;
  hid_t file_id;
  hsize_t snb, st, syq;

  // Open input file
  // -------------------------------------------------------------------------
  file_id = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    MYH5CHECK(file_id);

  // Get dataset sizes
  // -------------------------------------------------------------------------
  ierr = H5LTget_dataset_info(file_id, "nb", &snb, NULL, NULL);
    MYH5CHECK(ierr);
  ierr = H5LTget_dataset_info(file_id, "t", &st, NULL, NULL);
    MYH5CHECK(ierr);
  ierr = H5LTget_dataset_info(file_id, "yq", &syq, NULL, NULL);
    MYH5CHECK(ierr);
  m_nn = snb;
  m_nt = st;
  m_ny = syq;

  // Allocate memory
  // -------------------------------------------------------------------------
  m_log_nb = new double[m_nn];
  m_log_t = new double[m_nt];
  m_yq = new double[m_ny];
  m_table = new double[BNVARS*m_nn*m_ny*m_nt];
  double * scratch = new double[m_nn*m_ny*m_nt];

  // Read nb, t, yq
  // -------------------------------------------------------------------------
  ierr = H5LTread_dataset_double(file_id, "nb", scratch);
    MYH5CHECK(ierr);
  //min_n = scratch[0];
  //max_n = scratch[m_nn-1];
  for (int in = 0; in < m_nn; ++in) {
    m_log_nb[in] = log(scratch[in]);
  }
  m_id_log_nb = 1.0/(m_log_nb[1] - m_log_nb[0]);

  ierr = H5LTread_dataset_double(file_id, "t", scratch);
    MYH5CHECK(ierr);
  //min_T = scratch[0];
  //max_T = scratch[m_nt-1];
  for (int it = 0; it < m_nt; ++it) {
    m_log_t[it] = log(scratch[it]);
  }
  m_id_log_t = 1.0/(m_log_t[1] - m_log_t[0]);

  ierr = H5LTread_dataset_double(file_id, "yq", scratch);
    MYH5CHECK(ierr);
  //min_Y[0] = scratch[0];
  //max_Y[0] = scratch[m_ny-1];
  for (int iy = 0; iy < m_ny; ++iy) {
    m_yq[iy] = scratch[iy];
  }
  m_id_yq = 1.0/(m_yq[1] - m_yq[0]);

  // the neutron mass is used as the baryon mass in CompOSE
  ierr = H5LTread_dataset_double(file_id, "mn", scratch);
    MYH5CHECK(ierr);
  m_bar = scratch[0];

  // Read other thermodynamics quantities
  // -------------------------------------------------------------------------
  ierr = H5LTread_dataset_double(file_id, "Q1", scratch);
    MYH5CHECK(ierr);
  for (int inb = 0; inb < m_nn; ++inb) {
  for (int iyq = 0; iyq < m_ny; ++iyq) {
  for (int it = 0; it < m_nt; ++it) {
    Pmin = min(Pmin,exp(m_log_nb[inb])*scratch[index(0, inb, iyq, it)]);
  }}}
  Pmin = fabs(Pmin) + 1.e-10;

  for (int inb = 0; inb < m_nn; ++inb) {
  for (int iyq = 0; iyq < m_ny; ++iyq) {
  for (int it = 0; it < m_nt; ++it) {
    m_table[index(BLOGP, inb, iyq, it)] =
        //log(scratch[index(0, inb, iyq, it)]) + m_log_nb[inb];
        log(exp(m_log_nb[inb])*scratch[index(0, inb, iyq, it)] + Pmin);
    //cout << m_table[index(BLOGP, inb, iyq, it)] << endl;
  }}}

  ierr = H5LTread_dataset_double(file_id, "Q2", scratch);
    MYH5CHECK(ierr);
  copy(&scratch[0], &scratch[m_nn*m_ny*m_nt], &m_table[index(BENT, 0, 0, 0)]);

  ierr = H5LTread_dataset_double(file_id, "Q3", scratch);
    MYH5CHECK(ierr);
  for (int in = 0; in < m_nn; ++in) {
  for (int iy = 0; iy < m_ny; ++iy) {
  for (int it = 0; it < m_nt; ++it) {
    m_table[index(BMUB, in, iy, it)] =
      m_bar*(scratch[index(0, in, iy, it)] + 1);
  }}}

  ierr = H5LTread_dataset_double(file_id, "Q4", scratch);
    MYH5CHECK(ierr);
  for (int in = 0; in < m_nn; ++in) {
  for (int iy = 0; iy < m_ny; ++iy) {
  for (int it = 0; it < m_nt; ++it) {
    m_table[index(BMUQ, in, iy, it)] = m_bar*scratch[index(0, in, iy, it)];
  }}}

  //ierr = H5LTread_dataset_double(file_id, "Q5", scratch);
  //  MYH5CHECK(ierr);
  //for (int in = 0; in < m_nn; ++in) {
  //for (int iy = 0; iy < m_ny; ++iy) {
  //for (int it = 0; it < m_nt; ++it) {
  //  m_table[index(ECMUL, in, iy, it)] = m_bar*scratch[index(0, in, iy, it)];
  //}}}

  ierr = H5LTread_dataset_double(file_id, "Q7", scratch);
    MYH5CHECK(ierr);
  for (int in = 0; in < m_nn; ++in) {
  for (int iy = 0; iy < m_ny; ++iy) {
  for (int it = 0; it < m_nt; ++it) {
    m_table[index(BLOGE, in, iy, it)] =
      log(m_bar*(scratch[index(0, in, iy, it)] + 1)) + m_log_nb[in];
  }}}

  ierr = H5LTread_dataset_double(file_id, "Y[He4]", scratch);
    MYH5CHECK(ierr);
  for (int in = 0; in < m_nn; ++in) {
  for (int iy = 0; iy < m_ny; ++iy) {
  for (int it = 0; it < m_nt; ++it) {
    m_table[index(BYALP, in, iy, it)] = scratch[index(0, in, iy, it)];
  }}}

  ierr = H5LTread_dataset_double(file_id, "Y[N]", scratch);
    MYH5CHECK(ierr);
  for (int in = 0; in < m_nn; ++in) {
  for (int iy = 0; iy < m_ny; ++iy) {
  for (int it = 0; it < m_nt; ++it) {
    m_table[index(BYNUC, in, iy, it)] = scratch[index(0, in, iy, it)];
  }}}

  ierr = H5LTread_dataset_double(file_id, "Y[n]", scratch);
    MYH5CHECK(ierr);
  for (int in = 0; in < m_nn; ++in) {
  for (int iy = 0; iy < m_ny; ++iy) {
  for (int it = 0; it < m_nt; ++it) {
    m_table[index(BYNTR, in, iy, it)] = scratch[index(0, in, iy, it)];
  }}}

  ierr = H5LTread_dataset_double(file_id, "Y[p]", scratch);
    MYH5CHECK(ierr);
  for (int in = 0; in < m_nn; ++in) {
  for (int iy = 0; iy < m_ny; ++iy) {
  for (int it = 0; it < m_nt; ++it) {
    m_table[index(BYPTN, in, iy, it)] = scratch[index(0, in, iy, it)];
  }}}

  ierr = H5LTread_dataset_double(file_id, "dPdn", scratch);
    MYH5CHECK(ierr);
  for (int in = 0; in < m_nn; ++in) {
  for (int iy = 0; iy < m_ny; ++iy) {
  for (int it = 0; it < m_nt; ++it) {
    m_table[index(BDPDN, in, iy, it)] = scratch[index(0, in, iy, it)];
  }}}

  ierr = H5LTread_dataset_double(file_id, "dSdn", scratch);
    MYH5CHECK(ierr);
  for (int in = 0; in < m_nn; ++in) {
  for (int iy = 0; iy < m_ny; ++iy) {
  for (int it = 0; it < m_nt; ++it) {
    m_table[index(BDSDN, in, iy, it)] = scratch[index(0, in, iy, it)];
  }}}

  ierr = H5LTread_dataset_double(file_id, "dPdt", scratch);
    MYH5CHECK(ierr);
  for (int in = 0; in < m_nn; ++in) {
  for (int iy = 0; iy < m_ny; ++iy) {
  for (int it = 0; it < m_nt; ++it) {
    m_table[index(BDPDT, in, iy, it)] = scratch[index(0, in, iy, it)];
  }}}

  ierr = H5LTread_dataset_double(file_id, "dSdt", scratch);
    MYH5CHECK(ierr);
  for (int in = 0; in < m_nn; ++in) {
  for (int iy = 0; iy < m_ny; ++iy) {
  for (int it = 0; it < m_nt; ++it) {
    m_table[index(BDSDT, in, iy, it)] = scratch[index(0, in, iy, it)];
  }}}

  // Cleanup
  // -------------------------------------------------------------------------
  delete[] scratch;
  H5Fclose(file_id);

  m_initialized = true;
}

double EOS_baryons::eval_at_nty(int vi, double n, double T, double Yq) {
  return eval_at_lnty(vi, log(n), log(T), Yq);
}

void EOS_baryons::weight_idx_ln(double log_n) {
  in = (log_n - m_log_nb[0])*m_id_log_nb;
  wn1 = (log_n - m_log_nb[in])*m_id_log_nb;
  wn0 = 1.0 - wn1;
}

void EOS_baryons::weight_idx_yq(double yq) {
  iy = (yq - m_yq[0])*m_id_yq;
  wy1 = (yq - m_yq[iy])*m_id_yq;
  wy0 = 1.0 - wy1;
}

void EOS_baryons::weight_idx_lt(double log_t) {
  it = (log_t - m_log_t[0])*m_id_log_t;
  wt1 = (log_t - m_log_t[it])*m_id_log_t;
  wt0 = 1.0 - wt1;
}

double EOS_baryons::eval_at_lnty(int iv, double log_n, double log_t, double yq) {

  weight_idx_ln(log_n);
  weight_idx_yq(yq);
  weight_idx_lt(log_t);

  return
    wn0 * (wy0 * (wt0 * m_table[index(iv, in+0, iy+0, it+0)]   +
                  wt1 * m_table[index(iv, in+0, iy+0, it+1)])  +
           wy1 * (wt0 * m_table[index(iv, in+0, iy+1, it+0)]   +
                  wt1 * m_table[index(iv, in+0, iy+1, it+1)])) +
    wn1 * (wy0 * (wt0 * m_table[index(iv, in+1, iy+0, it+0)]   +
                  wt1 * m_table[index(iv, in+1, iy+0, it+1)])  +
           wy1 * (wt0 * m_table[index(iv, in+1, iy+1, it+0)]   +
                  wt1 * m_table[index(iv, in+1, iy+1, it+1)]));
}

double EOS_baryons::interp_3d(int iv) const {
   return
    wn0 * (wy0 * (wt0 * m_table[index(iv, in+0, iy+0, it+0)]   +
                  wt1 * m_table[index(iv, in+0, iy+0, it+1)])  +
           wy1 * (wt0 * m_table[index(iv, in+0, iy+1, it+0)]   +
                  wt1 * m_table[index(iv, in+0, iy+1, it+1)])) +
    wn1 * (wy0 * (wt0 * m_table[index(iv, in+1, iy+0, it+0)]   +
                  wt1 * m_table[index(iv, in+1, iy+0, it+1)])  +
           wy1 * (wt0 * m_table[index(iv, in+1, iy+1, it+0)]   +
                  wt1 * m_table[index(iv, in+1, iy+1, it+1)]));
}


/* double EOS_baryons::eval_at_nty_new(int iv, double n, double t, double yq) const {
  int in, iy, it;
  double wn0, wn1, wy0, wy1, wt0, wt1;

  weight_idx_ln(&wn0, &wn1, &in, log(n));
  weight_idx_yq(&wy0, &wy1, &iy, yq);
  weight_idx_lt(&wt0, &wt1, &it, log(t));

  return
    pow(m_table[index(iv, in+0, iy+0, it+0)], wn0*wy0*wt0) *
    pow(m_table[index(iv, in+0, iy+0, it+1)], wn0*wy0*wt1) *
    pow(m_table[index(iv, in+0, iy+1, it+0)], wn0*wy1*wt0) *
    pow(m_table[index(iv, in+0, iy+1, it+1)], wn0*wy1*wt1) *
    pow(m_table[index(iv, in+1, iy+0, it+0)], wn1*wy0*wt0) *
    pow(m_table[index(iv, in+1, iy+0, it+1)], wn1*wy0*wt1) *
    pow(m_table[index(iv, in+1, iy+1, it+0)], wn1*wy1*wt0) *
    pow(m_table[index(iv, in+1, iy+1, it+1)], wn1*wy1*wt1);
} */
