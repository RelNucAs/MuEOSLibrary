//! \file eos_compose.cpp
//  \brief Implementation of EOSCompose

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <limits>
#include <sstream>
#include <stdexcept>

#include <hdf5.h>
#include <hdf5_hl.h>

#include "eos_baryons.hpp"

using namespace std;

#define MYH5CHECK(ierr) \
  if(ierr < 0) { \
    stringstream ss; \
    ss << __FILE__ << ":" << __LINE__ << " error reading EOS table!"; \
    throw runtime_error(ss.str().c_str()); \
  }

void ReadTableFromFile(std::string fname, int *m_nn, int *m_nt, int *m_ny, double *m_log_nb, double *m_log_t, double *m_yq, double *m_table) {
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
  int mm_nn = snb;
  int mm_nt = st;
  int mm_ny = syq;

  // Allocate memory
  // -------------------------------------------------------------------------
  m_log_nb = new double[mm_nn];
  m_log_t = new double[mm_nt];
  m_yq = new double[mm_ny];
  m_table = new double[ECNVARS*mm_nn*mm_ny*mm_nt];
  double * scratch = new double[mm_nn*mm_ny*mm_nt];

  double mb;  
  double min_n, max_n, m_id_log_nb;
  double min_T, max_T, m_id_log_t;
  double min_Y[MAX_SPECIES], max_Y[MAX_SPECIES];
  double m_id_yq;

  // Read nb, t, yq
  // -------------------------------------------------------------------------
  ierr = H5LTread_dataset_double(file_id, "nb", scratch);
    MYH5CHECK(ierr);
  min_n = scratch[0];
  max_n = scratch[mm_nn-1];
  for (int in = 0; in < mm_nn; ++in) {
    m_log_nb[in] = log(scratch[in]);
  }
  m_id_log_nb = 1.0/(m_log_nb[1] - m_log_nb[0]);

  ierr = H5LTread_dataset_double(file_id, "t", scratch);
    MYH5CHECK(ierr);
  min_T = scratch[0];
  max_T = scratch[mm_nt-1];
  for (int it = 0; it < mm_nt; ++it) {
    m_log_t[it] = log(scratch[it]);
  }
  m_id_log_t = 1.0/(m_log_t[1] - m_log_t[0]);

  ierr = H5LTread_dataset_double(file_id, "yq", scratch);
    MYH5CHECK(ierr);
  min_Y[0] = scratch[0];
  max_Y[0] = scratch[mm_ny-1];
  for (int iy = 0; iy < mm_ny; ++iy) {
    m_yq[iy] = scratch[iy];
  }
  m_id_yq = 1.0/(m_yq[1] - m_yq[0]);

  // the neutron mass is used as the baryon mass in CompOSE
  ierr = H5LTread_dataset_double(file_id, "mn", scratch);
    MYH5CHECK(ierr);
  mb = scratch[0];

  // Read other thermodynamics quantities
  // -------------------------------------------------------------------------
  ierr = H5LTread_dataset_double(file_id, "Q1", scratch);
    MYH5CHECK(ierr);
  for (int inb = 0; inb < mm_nn; ++inb) {
  for (int iyq = 0; iyq < mm_ny; ++iyq) {
  for (int it = 0; it < mm_nt; ++it) {
    m_table[index(mm_nn, mm_nt, mm_ny, ECLOGP, inb, iyq, it)] =
        log(scratch[index(mm_nn, mm_nt, mm_ny, 0, inb, iyq, it)]) + m_log_nb[inb];
  }}}

  ierr = H5LTread_dataset_double(file_id, "Q2", scratch);
    MYH5CHECK(ierr);
  copy(&scratch[0], &scratch[mm_nn*mm_ny*mm_nt], &m_table[index(mm_nn, mm_nt, mm_ny, ECENT, 0, 0, 0)]);

  ierr = H5LTread_dataset_double(file_id, "Q3", scratch);
    MYH5CHECK(ierr);
  for (int in = 0; in < mm_nn; ++in) {
  for (int iy = 0; iy < mm_ny; ++iy) {
  for (int it = 0; it < mm_nt; ++it) {
    m_table[index(mm_nn, mm_nt, mm_ny, ECMUB, in, iy, it)] =
      mb*(scratch[index(mm_nn, mm_nt, mm_ny, 0, in, iy, it)] + 1);
  }}}

  ierr = H5LTread_dataset_double(file_id, "Q4", scratch);
    MYH5CHECK(ierr);
  for (int in = 0; in < mm_nn; ++in) {
  for (int iy = 0; iy < mm_ny; ++iy) {
  for (int it = 0; it < mm_nt; ++it) {
    m_table[index(mm_nn, mm_nt, mm_ny, ECMUQ, in, iy, it)] = mb*scratch[index(mm_nn, mm_nt, mm_ny, 0, in, iy, it)];
  }}}

  ierr = H5LTread_dataset_double(file_id, "Q5", scratch);
    MYH5CHECK(ierr);
  for (int in = 0; in < mm_nn; ++in) {
  for (int iy = 0; iy < mm_ny; ++iy) {
  for (int it = 0; it < mm_nt; ++it) {
    m_table[index(mm_nn, mm_nt, mm_ny, ECMUL, in, iy, it)] = mb*scratch[index(mm_nn, mm_nt, mm_ny, 0, in, iy, it)];
  }}}

  ierr = H5LTread_dataset_double(file_id, "Q7", scratch);
    MYH5CHECK(ierr);
  for (int in = 0; in < mm_nn; ++in) {
  for (int iy = 0; iy < mm_ny; ++iy) {
  for (int it = 0; it < mm_nt; ++it) {
    m_table[index(mm_nn, mm_nt, mm_ny, ECLOGE, in, iy, it)] =
      log(mb*(scratch[index(mm_nn, mm_nt, mm_ny, 0, in, iy, it)] + 1)) + m_log_nb[in];
  }}}

  ierr = H5LTread_dataset_double(file_id, "cs2", scratch);
    MYH5CHECK(ierr);
  for (int in = 0; in < mm_nn; ++in) {
  for (int iy = 0; iy < mm_ny; ++iy) {
  for (int it = 0; it < mm_nt; ++it) {
    m_table[index(mm_nn, mm_nt, mm_ny, ECCS, in, iy, it)] = sqrt(scratch[index(mm_nn, mm_nt, mm_ny, 0, in, iy, it)]);
  }}}

  // Cleanup
  // -------------------------------------------------------------------------
  delete[] scratch;
  H5Fclose(file_id);

  m_nn = &mm_nn;
  m_nt = &mm_nt;
  m_ny = &mm_ny;
  //m_initialized = true;

  // Compute minimum enthalpy
  // -------------------------------------------------------------------------
  //for (int in = 0; in < mm_nn; ++in) {
  //  double const nb = exp(m_log_nb[in]);
  //  for (int it = 0; it < mm_nt; ++it) {
  //    double const t = exp(m_log_t[it]);
  //    for (int iy = 0; iy < mm_ny; ++iy) {
  //      m_min_h = min(m_min_h, Enthalpy(nb, t, &m_yq[iy]));
  //    }
  //  }
  //}
}
