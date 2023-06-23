#ifndef EOS_LEPTONS_H
#define EOS_LEPTONS_H

//! \file eos_leptons.hpp
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
#include <array>

#include "eos_fermions.hpp"

using namespace std; 

template <int idx_lep>
class EOS_leptons {
  public:
    enum EOSQuantities {
        IL_N  = 0, //! Number density of leptons
        IA_N  = 1, //! Number density of anti-leptons
        IL_P  = 2, //! Pressure of leptons
        IA_P  = 3, //! Pressure of anti-leptons
        IL_E  = 4, //! Internal energy density of leptons
        IA_E  = 5, //! Internal energy density of anti-leptons
        IL_S  = 6, //! Entropy density of leptons
        IA_S  = 7, //! Entropy density of anti-leptons
        IL_MU = 8, //! Chemical potential of leptons
        NVARS = 9,
	ID_DPDN = 9,
	ID_DSDN = 10,
	ID_DPDT = 11,
	ID_DSDT = 12,
	NTOT = 13
    };
    
    static const int nEOS = NVARS;

  private:
    // Inverse of table spacing
    double m_id_log_nl, m_id_log_tl;

  public:
    // Table storage, care should be made to store these data on the GPU later
    double * m_log_nl;
    double * m_log_tl;
    double * m_lep_table;

    //const int idx_lep  = 0;
    //const int idx_lep = 1;

    // Table size
    int m_nl, m_ntl;

  private:    
    bool m_lep_initialized;

  public:
    bool m_lep_active;

  public:
    /// Constructor
    EOS_leptons():
      m_id_log_nl(std::numeric_limits<double>::quiet_NaN()),
      m_id_log_tl(std::numeric_limits<double>::quiet_NaN()),
      m_log_nl(nullptr),
      m_log_tl(nullptr),
      m_lep_table(nullptr),
      m_nl(0), m_ntl(0),
      m_lep_initialized(false),
      m_lep_active(false) {};

  
    /// Destructor
    ~EOS_leptons() {
      if (m_lep_initialized) {
        delete[] m_log_nl;
        delete[] m_log_tl;
        delete[] m_lep_table;
      }
    }

    /// Calculate the number density of leptons using.
    template <int id_EOS>
    double LepNumberDensity(double n, double T, double *Y) {
      if (m_lep_active) {
	if constexpr(id_EOS == 1) {
          assert(m_lep_initialized);
          return exp(eval_lep_at_nty(IL_N, n, T, Y[idx_lep]));
        } else if constexpr(id_EOS == 2) {
          return eos_ferm_single<id_EOS,idx_lep>(n*Y[idx_lep], T)[IL_N];
        }
      } else {
	return 0.;
      }
    }

    /// Calculate the number density of antileptons using.
    template <int id_EOS>
    double AntiLepNumberDensity(double n, double T, double *Y) {
      if (m_lep_active) {
        //assert(m_lep_initialized);
        //return exp(eval_lep_at_nty(IA_N, n, T, Y[idx_lep]));
        return LepNumberDensity<id_EOS>(n, T, Y) - n*Y[idx_lep];
      } else {
	return 0.;
      }
    }

    /// Calculate the internal energy density of leptons using.
    template <int id_EOS>
    double LepEnergy(double n, double T, double *Y) {
      if (m_lep_active) {
	if constexpr(id_EOS == 1) {
            assert(m_lep_initialized);
            return exp(eval_lep_at_nty(IL_E, n, T, Y[idx_lep])) + exp(eval_lep_at_nty(IA_E, n, T, Y[idx_lep]));
        } else if constexpr(id_EOS == 2) {
	    return eos_ferm_single<id_EOS,idx_lep>(n*Y[idx_lep], T)[IL_E] + eos_ferm_single<id_EOS,idx_lep>(n*Y[idx_lep], T)[IA_E]; 
        }
      } else {
	return 0.;
      }
    }

    /// Calculate the pressure of leptons using.
    template <int id_EOS>
    double LepPressure(double n, double T, double *Y) {
      if (m_lep_active) {
        if constexpr(id_EOS == 1) {
            assert(m_lep_initialized);
            return exp(eval_lep_at_nty(IL_P, n, T, Y[idx_lep])) + exp(eval_lep_at_nty(IA_P, n, T, Y[idx_lep]));
        } else if constexpr(id_EOS == 2) {
            return eos_ferm_single<id_EOS,idx_lep>(n*Y[idx_lep], T)[IL_P] + eos_ferm_single<id_EOS,idx_lep>(n*Y[idx_lep], T)[IA_P];
        }
      } else {
	return 0.;
      }
    }

    /// Calculate the entropy density of leptons using.
    template <int id_EOS>
    double LepEntropy(double n, double T, double *Y) {
      if (m_lep_active) {
        if constexpr(id_EOS == 1) {
            assert(m_lep_initialized);
            return exp(eval_lep_at_nty(IL_S, n, T, Y[idx_lep])) + exp(eval_lep_at_nty(IA_S, n, T, Y[idx_lep]));
        } else if constexpr(id_EOS == 2) {
            return eos_ferm_single<id_EOS,idx_lep>(n*Y[idx_lep], T)[IL_S] + eos_ferm_single<id_EOS,idx_lep>(n*Y[idx_lep], T)[IA_S];
        }
      } else {
	return 0.;
      }
    }

    /// Calculate the chemical potential of leptons using.
    template <int id_EOS>
    double LepChemicalPotential(double n, double T, double *Y) {
      if (m_lep_active) {
        if constexpr(id_EOS == 1) {
            assert(m_lep_initialized);
            return eval_lep_at_nty(IL_MU, n, T, Y[idx_lep]);
        } else if constexpr(id_EOS == 2) {
            return eos_ferm_single<id_EOS,idx_lep>(n*Y[idx_lep], T)[IL_MU];
        }
      } else {
	return 0.;
      }
    }

    /// Calculate the sound speed derivatives using.
    template <int id_EOS>
    double LdPdn(double n, double T, double *Y) {
      if (m_lep_active) {
        if constexpr(id_EOS == 1) {
            assert(m_lep_initialized);
            return eval_lep_at_nty(ID_DPDN, n, T, Y[idx_lep])*Y[idx_lep];
        } else if constexpr(id_EOS == 2) {
            return der_cs2_num<idx_lep>(n*Y[idx_lep], T)[0]*Y[idx_lep];
	}
      } else {
	return 0.;
      }
    }

    template <int id_EOS>
    double Ldsdn(double n, double T, double *Y) {
      if (m_lep_active) {
        if constexpr(id_EOS == 1) {
            assert(m_lep_initialized);
            return eval_lep_at_nty(ID_DSDN, n, T, Y[idx_lep]) * Y[idx_lep] / n - LepEntropy<1>(n, T, Y) / (n*n);
        } else if constexpr(id_EOS == 2) {
            return der_cs2_num<idx_lep>(n*Y[idx_lep], T)[1] * Y[idx_lep] / n - LepEntropy<2>(n, T, Y) / (n*n);
	}
      } else {
	return 0.;
      }
    }
    
    template <int id_EOS>
    double LdPdt(double n, double T, double *Y) {
      if (m_lep_active) {
        if constexpr(id_EOS == 1) {
            assert(m_lep_initialized);
            return eval_lep_at_nty(ID_DPDT, n, T, Y[idx_lep]);
        } else if constexpr(id_EOS == 2) {
            return der_cs2_num<idx_lep>(n*Y[idx_lep], T)[2];
	}
      } else {
	return 0.;
      }
    }
    
    template <int id_EOS>
    double Ldsdt(double n, double T, double *Y) {
      if (m_lep_active) {
        if constexpr(id_EOS == 1) {
            assert(m_lep_initialized);
            return eval_lep_at_nty(ID_DSDT, n, T, Y[idx_lep]) / n;
        } else if constexpr(id_EOS == 2) {
            return der_cs2_num<idx_lep>(n*Y[idx_lep], T)[3] / n;
	}
      } else {
	return 0.;
      }
    }

    /// Calculate the full leptonic EOS using.
    std::array<double,NVARS> ComputeFullLepEOS(double n, double T, double *Y) {
      assert(m_lep_initialized);
      return eval_all_lep_at_lnt(log(n*Y[idx_lep]),log(T));
    }

  public:
    /// Reads the table file.
    void ReadLepTableFromFile(std::string fname) {
      // Open input file
      // -------------------------------------------------------------------------
      ifstream EOSinput;
      EOSinput.open(fname);

      if (!EOSinput) {
        cout << "Lepton EOS table not found!: id_lep = " << idx_lep << endl;
        exit(EXIT_FAILURE);
      }

      // Get dataset sizes
      // -------------------------------------------------------------------------
      string EOSline;
      getline(EOSinput, EOSline);
      stringstream s1(EOSline);
      s1 >> m_nl;

      getline(EOSinput, EOSline);
      stringstream s2(EOSline);
      s2 >> m_ntl;

      // Allocate memory
      // -------------------------------------------------------------------------
      m_log_nl   = new double[m_nl];
      m_log_tl   = new double[m_ntl];
      m_lep_table = new double[NTOT*m_nl*m_ntl];

      double number;

      // Read nb, t, yq
      // -------------------------------------------------------------------------
      getline(EOSinput, EOSline);
      stringstream s3(EOSline);
      for (int in = 0; in < m_nl; ++in) {
        s3 >> number;
        m_log_nl[in] = log(pow(10.,number));
      }
      m_id_log_nl = 1.0/(m_log_nl[1] - m_log_nl[0]);

      getline(EOSinput, EOSline);
      stringstream s4(EOSline);
      for (int it = 0; it < m_ntl; ++it) {
        s4 >> number;
        m_log_tl[it] = log(pow(10.,number));
      }
      m_id_log_tl = 1.0/(m_log_tl[1] - m_log_tl[0]);


      // Read other thermodynamics quantities
      // -------------------------------------------------------------------------
      for (int in = 0; in < m_nl; ++in) {
        getline(EOSinput, EOSline);
        stringstream s5(EOSline);
          for (int it = 0; it < m_ntl; ++it) {
            s5 >> number;
            m_lep_table[lep_index(IL_N, in, it)] = log(number);
      }}
    
      for (int in = 0; in < m_nl; ++in) {
        getline(EOSinput, EOSline);
        stringstream s5(EOSline);
          for (int it = 0; it < m_ntl; ++it) {
            s5 >> number;
            m_lep_table[lep_index(IA_N, in, it)] = log(number);
      }}
    
      for (int in = 0; in < m_nl; ++in) {
        getline(EOSinput, EOSline);
        stringstream s5(EOSline);
          for (int it = 0; it < m_ntl; ++it) {
            s5 >> number;
            m_lep_table[lep_index(IL_P, in, it)] = log(number);
      }}
    
      for (int in = 0; in < m_nl; ++in) {
        getline(EOSinput, EOSline);
        stringstream s5(EOSline);
          for (int it = 0; it < m_ntl; ++it) {
            s5 >> number;
            m_lep_table[lep_index(IA_P, in, it)] = log(number);
      }}


      for (int in = 0; in < m_nl; ++in) {
        getline(EOSinput, EOSline);
        stringstream s5(EOSline);
          for (int it = 0; it < m_ntl; ++it) {
            s5 >> number;
            m_lep_table[lep_index(IL_E, in, it)] = log(number);
      }}
    
      for (int in = 0; in < m_nl; ++in) {
        getline(EOSinput, EOSline);
        stringstream s5(EOSline);
          for (int it = 0; it < m_ntl; ++it) {
            s5 >> number;
            m_lep_table[lep_index(IA_E, in, it)] = log(number);
      }}
    
      for (int in = 0; in < m_nl; ++in) {
        getline(EOSinput, EOSline);
        stringstream s5(EOSline);
          for (int it = 0; it < m_ntl; ++it) {
            s5 >> number;
            m_lep_table[lep_index(IL_S, in, it)] = log(number);
      }}
    
      for (int in = 0; in < m_nl; ++in) {
        getline(EOSinput, EOSline);
        stringstream s5(EOSline);
          for (int it = 0; it < m_ntl; ++it) {
            s5 >> number;
            m_lep_table[lep_index(IA_S, in, it)] = log(number);
      }}
    
      for (int in = 0; in < m_nl; ++in) {
        getline(EOSinput, EOSline);
        stringstream s5(EOSline);
          for (int it = 0; it < m_ntl; ++it) {
            s5 >> number;
            m_lep_table[lep_index(IL_MU, in, it)] = number;
      }}
    
      for (int in = 0; in < m_nl; ++in) {
        getline(EOSinput, EOSline);
        stringstream s5(EOSline);
          for (int it = 0; it < m_ntl; ++it) {
            s5 >> number;
            m_lep_table[lep_index(ID_DPDN, in, it)] = number;
      }}
    
      for (int in = 0; in < m_nl; ++in) {
        getline(EOSinput, EOSline);
        stringstream s5(EOSline);
          for (int it = 0; it < m_ntl; ++it) {
            s5 >> number;
            m_lep_table[lep_index(ID_DSDN, in, it)] = number;
      }}
    
      for (int in = 0; in < m_nl; ++in) {
        getline(EOSinput, EOSline);
        stringstream s5(EOSline);
          for (int it = 0; it < m_ntl; ++it) {
            s5 >> number;
            m_lep_table[lep_index(ID_DPDT, in, it)] = number;
      }}
    
      for (int in = 0; in < m_nl; ++in) {
        getline(EOSinput, EOSline);
        stringstream s5(EOSline);
          for (int it = 0; it < m_ntl; ++it) {
            s5 >> number;
            m_lep_table[lep_index(ID_DSDT, in, it)] = number;
      }}
    
      // Cleanup
      // -------------------------------------------------------------------------
      //delete[] s5;
      EOSinput.close();
    
      m_lep_initialized = true;
    }
    
    /// Get the raw number density
    double const * GetRawLepLogNumberDensity() const {
      return m_log_nl;
    }

    /// Get the raw table temperature
    double const * GetRawLepLogTemperature() const {
      return m_log_tl;
    }

    /// Get the raw table data
    double const * GetLepTable() const {
      return m_lep_table;
    }

    // Indexing used to access the lepton data
    inline ptrdiff_t lep_index(int iv, int in, int it) const {
      return it + m_ntl*(in + m_nl*iv);
    }

    /// Check if the EOS has been initialized properly.
    inline bool IsLepInitialized() const {
      return m_lep_initialized;
    }

    /// Check if the lepton species is active.
    inline bool IsLepActive() const {
      return m_lep_active;
    }

  private:
    /// Low level evaluation function, not intended for outside use
    double eval_lep_at_nty(int vi, double n, double T, double Yl) const {
      // @FIXME
      if (Yl == 0.) return 0.;
      return eval_lep_at_lnty(vi, log(n*Yl), log(T));
    }
    
    /// Evaluate interpolation weight for lepton density
    void weight_idx_lnl(double *w0, double *w1, int *in, double log_n) const {
      *in = (log_n - m_log_nl[0])*m_id_log_nl;
      *w1 = (log_n - m_log_nl[*in])*m_id_log_nl;
      *w0 = 1.0 - (*w1);
    }
    
    /// Evaluate interpolation weight for lepton temperature
    void weight_idx_ltl(double *w0, double *w1, int *it, double log_t) const {
      *it = (log_t - m_log_tl[0])*m_id_log_tl;
      *w1 = (log_t - m_log_tl[*it])*m_id_log_tl;
      *w0 = 1.0 - (*w1);
    }
    
    /// Low level evaluation function, not intended for outside use
    double eval_lep_at_lnty(int iv, double log_n, double log_t) const {
      int in, it;
      double wn0, wn1, wt0, wt1;
    
      weight_idx_lnl(&wn0, &wn1, &in, log_n);
      weight_idx_ltl(&wt0, &wt1, &it, log_t);
    
      return
        wn0 * (wt0 * m_lep_table[lep_index(iv, in+0, it+0)]   +
               wt1 * m_lep_table[lep_index(iv, in+0, it+1)])  +
        wn1 * (wt0 * m_lep_table[lep_index(iv, in+1, it+0)]   +
               wt1 * m_lep_table[lep_index(iv, in+1, it+1)]);
    }


    /// Low level evaluation function for entire EOS, not intended for outside use
    std::array<double,NVARS> eval_all_lep_at_lnt(double log_n, double log_t) const {
      int in, it;
      double wn0, wn1, wt0, wt1;
      std::array<double,EOS_leptons::NVARS> eos_out;

      weight_idx_lnl(&wn0, &wn1, &in, log_n);
      weight_idx_ltl(&wt0, &wt1, &it, log_t);

      for (int iv=0; iv<NVARS; iv++) {
         eos_out[iv] = wn0 * (wt0 * m_lep_table[lep_index(iv, in+0, it+0)]   +
                              wt1 * m_lep_table[lep_index(iv, in+0, it+1)])  +
                       wn1 * (wt0 * m_lep_table[lep_index(iv, in+1, it+0)]   +
                              wt1 * m_lep_table[lep_index(iv, in+1, it+1)]);
      }
      return eos_out;
    }

 
    double mu_check(const bool check, const double value) {
      if (check) {
        return value;
      } else {
        return 0.;
      }
    }
};


#endif
