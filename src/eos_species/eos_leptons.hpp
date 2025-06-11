#ifndef MUEOSLIBRARY_SRC_EOS_LEPTONS_HPP_
#define MUEOSLIBRARY_SRC_EOS_LEPTONS_HPP_

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

#include "mueosclass.hpp"
#include "helmholtz_eos/helmholtz_eos.hpp"

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


    /// Calculate the full leptonic EOS using.
    template <int id_EOS>
    void LeptonEOS(double n, double T, double *Y) {
      if (!m_lep_active) {
        n_lep = 0.;
        e_lep = 0.;
        P_lep = 0.;
        s_lep = 0.;
        mu_lep = 0.;
        der_lep[0] = 0.;
        der_lep[1] = 0.;
        der_lep[2] = 0.;
        der_lep[3] = 0.;
        return;
       }

      if constexpr(id_EOS == 0) {
        assert(m_lep_initialized);
        weight_idx_lnl(log(n * Y[idx_lep]));
        weight_idx_ltl(log(T));

        n_lep = exp(interp_2d(IL_N));
        e_lep = exp(interp_2d(IL_E)) + exp(interp_2d(IA_E));
        P_lep = exp(interp_2d(IL_P)) + exp(interp_2d(IA_P));
        s_lep = exp(interp_2d(IL_S)) + exp(interp_2d(IA_S));
        mu_lep = interp_2d(IL_MU);
        der_lep[0] = interp_2d(ID_DPDN) * Y[idx_lep];
        der_lep[1] = interp_2d(ID_DSDN) * Y[idx_lep] / n - s_lep / (n*n);
        der_lep[2] = interp_2d(ID_DPDT);
        der_lep[3] = interp_2d(ID_DSDT) / n;


      } else if constexpr(id_EOS == 1) {
        assert(m_lep_initialized);
        
        const double eta = eval_lep_at_nty(0, n, T, Y[idx_lep]);

        const double theta = T / mL[idx_lep];
        const double a_eta = - (eta + 2. / theta);

        GFDs FD;

        FD.f12   = compute_res(eta,   theta, 0.5); // k = 1/2
        FD.f32   = compute_res(eta,   theta, 1.5); // k = 3/2      
        FD.a_f12 = compute_res(a_eta, theta, 0.5); // k = 1/2
        FD.a_f32 = compute_res(a_eta, theta, 1.5); // k = 3/2

        HelmEOSOutput tmp = eos_helm_from_eta(eta, T, idx_lep, &FD);
        HelmEOSDer tmp_der = der_cs2_from_eta(eta, T, idx_lep);

        n_lep = tmp.nl;
        e_lep = tmp.el + tmp.a_el;
        P_lep = tmp.pl + tmp.a_pl;
        s_lep = tmp.sl + tmp.a_sl;
        mu_lep = tmp.mul;
        der_lep[0] = tmp_der.dPdn * Y[idx_lep];
        der_lep[1] = tmp_der.dsdn * Y[idx_lep] / n - s_lep / (n*n);
        der_lep[2] = tmp_der.dPdt;
        der_lep[3] = tmp_der.dsdt / n;
       
      } else if constexpr(id_EOS == 2) {

        HelmEOSOutput tmp = eos_helm_full(n*Y[idx_lep], T, idx_lep);

        double eta = (tmp.mul - mL[idx_lep]) / T;

        HelmEOSDer tmp_der = der_cs2_from_eta(eta, T, idx_lep);

        n_lep = tmp.nl;
        e_lep = tmp.el + tmp.a_el;
        P_lep = tmp.pl + tmp.a_pl;
        s_lep = tmp.sl + tmp.a_sl;
        mu_lep = tmp.mul;
        der_lep[0] = tmp_der.dPdn * Y[idx_lep];
        der_lep[1] = tmp_der.dsdn * Y[idx_lep] / n - s_lep / (n*n);
        der_lep[2] = tmp_der.dPdt;
        der_lep[3] = tmp_der.dsdt / n;

      }
      return;
    }

// Calculate the number density of leptons using.
    double GetLeptonNumberDensity() {
      return n_lep;
    }

    // Calculate the number density of antileptons using.
    double GetAntiLeptonNumberDensity(double n, double* Y) {
      return n_lep - n * Y[idx_lep];
    } 

    /// Calculate the internal energy density of leptons using.
    double GetLeptonEnergy() {
      return e_lep;
    }

    /// Calculate the pressure of leptons using.
    double GetLeptonPressure() {
      return P_lep;
    }

    /// Calculate the entropy density of leptons using.
    double GetLeptonEntropy() {
      return s_lep;
    }

    /// Calculate the chemical potential of leptons using.
    double GetLeptonChemicalPotential() {
      return mu_lep;
    }

    /// Calculate the sound speed derivatives using.
    double GetLeptondPdn() {
      return der_lep[0];
    }

    double GetLeptondsdn() {
  	      return der_lep[1];
    }
    
    double GetLeptondPdt() {
	      return der_lep[2];
      }
    
    double GetLeptondsdt() {
	      return der_lep[3];
    }

/*     /// Calculate the full leptonic EOS using.
    template <int id_EOS>
    HelmEOSOutput ComputeFullLepEOS(double n, double T, double *Y) {
      HelmEOSOutput out = {0.0};
      if (m_lep_active) {
        if constexpr(id_EOS == 1) {
          assert(m_lep_initialized);
          // @TODO: fix problem when Yl = 0
          return eval_all_lep_at_lnty(log(n*Y[idx_lep]), log(T));
        } else if constexpr(id_EOS == 2) {
          return eos_helm_full(n*Y[idx_lep], T, idx_lep);
        }
      }
    } 
  */

 
  public:
    /// Reads the table file.
    void ReadFullLepTableFromFile(std::string fname) {
      // Open input file
      // -------------------------------------------------------------------------
      std::ifstream EOSinput;
      EOSinput.open(fname);

      if (!EOSinput) {
        std::cout << "Lepton EOS table not found!: id_lep = " << idx_lep << std::endl;
        std::cout << "Lepton EOS tables can be generated by running 'make lep_table' in the parent directory" << std::endl;
        exit(EXIT_FAILURE);
      }

      // Get dataset sizes
      // -------------------------------------------------------------------------
      std::string EOSline;
      getline(EOSinput, EOSline);
      std::stringstream s1(EOSline);
      s1 >> m_nl;

      getline(EOSinput, EOSline);
      std::stringstream s2(EOSline);
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
      std::stringstream s3(EOSline);
      for (int in = 0; in < m_nl; ++in) {
        s3 >> number;
        m_log_nl[in] = log(pow(10.,number));
      }
      m_id_log_nl = 1.0/(m_log_nl[1] - m_log_nl[0]);

      getline(EOSinput, EOSline);
      std::stringstream s4(EOSline);
      for (int it = 0; it < m_ntl; ++it) {
        s4 >> number;
        m_log_tl[it] = log(pow(10.,number));
      }
      m_id_log_tl = 1.0/(m_log_tl[1] - m_log_tl[0]);


      // Read other thermodynamics quantities
      // -------------------------------------------------------------------------
      for (int in = 0; in < m_nl; ++in) {
        getline(EOSinput, EOSline);
        std::stringstream s5(EOSline);
          for (int it = 0; it < m_ntl; ++it) {
            s5 >> number;
            m_lep_table[lep_index(IL_N, in, it)] = log(number);
      }}
    
      for (int in = 0; in < m_nl; ++in) {
        getline(EOSinput, EOSline);
        std::stringstream s5(EOSline);
          for (int it = 0; it < m_ntl; ++it) {
            s5 >> number;
            m_lep_table[lep_index(IA_N, in, it)] = log(number);
      }}
    
      for (int in = 0; in < m_nl; ++in) {
        getline(EOSinput, EOSline);
        std::stringstream s5(EOSline);
          for (int it = 0; it < m_ntl; ++it) {
            s5 >> number;
            m_lep_table[lep_index(IL_P, in, it)] = log(number);
      }}
    
      for (int in = 0; in < m_nl; ++in) {
        getline(EOSinput, EOSline);
        std::stringstream s5(EOSline);
          for (int it = 0; it < m_ntl; ++it) {
            s5 >> number;
            m_lep_table[lep_index(IA_P, in, it)] = log(number);
      }}


      for (int in = 0; in < m_nl; ++in) {
        getline(EOSinput, EOSline);
        std::stringstream s5(EOSline);
          for (int it = 0; it < m_ntl; ++it) {
            s5 >> number;
            m_lep_table[lep_index(IL_E, in, it)] = log(number);
      }}
    
      for (int in = 0; in < m_nl; ++in) {
        getline(EOSinput, EOSline);
        std::stringstream s5(EOSline);
          for (int it = 0; it < m_ntl; ++it) {
            s5 >> number;
            m_lep_table[lep_index(IA_E, in, it)] = log(number);
      }}
    
      for (int in = 0; in < m_nl; ++in) {
        getline(EOSinput, EOSline);
        std::stringstream s5(EOSline);
          for (int it = 0; it < m_ntl; ++it) {
            s5 >> number;
            m_lep_table[lep_index(IL_S, in, it)] = log(number);
      }}
    
      for (int in = 0; in < m_nl; ++in) {
        getline(EOSinput, EOSline);
        std::stringstream s5(EOSline);
          for (int it = 0; it < m_ntl; ++it) {
            s5 >> number;
            m_lep_table[lep_index(IA_S, in, it)] = log(number);
      }}
    
      for (int in = 0; in < m_nl; ++in) {
        getline(EOSinput, EOSline);
        std::stringstream s5(EOSline);
          for (int it = 0; it < m_ntl; ++it) {
            s5 >> number;
            m_lep_table[lep_index(IL_MU, in, it)] = number;
      }}
    
      for (int in = 0; in < m_nl; ++in) {
        getline(EOSinput, EOSline);
        std::stringstream s5(EOSline);
          for (int it = 0; it < m_ntl; ++it) {
            s5 >> number;
            m_lep_table[lep_index(ID_DPDN, in, it)] = number;
      }}
    
      for (int in = 0; in < m_nl; ++in) {
        getline(EOSinput, EOSline);
        std::stringstream s5(EOSline);
          for (int it = 0; it < m_ntl; ++it) {
            s5 >> number;
            m_lep_table[lep_index(ID_DSDN, in, it)] = number;
      }}
    
      for (int in = 0; in < m_nl; ++in) {
        getline(EOSinput, EOSline);
        std::stringstream s5(EOSline);
          for (int it = 0; it < m_ntl; ++it) {
            s5 >> number;
            m_lep_table[lep_index(ID_DPDT, in, it)] = number;
      }}
    
      for (int in = 0; in < m_nl; ++in) {
        getline(EOSinput, EOSline);
        std::stringstream s5(EOSline);
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
    

      /// Reads the table file.
    void ReadEtaLepTableFromFile(std::string fname) {
      // Open input file
      // -------------------------------------------------------------------------
      std::ifstream EOSinput;
      EOSinput.open(fname);

      if (!EOSinput) {
        std::cout << "Lepton EOS table not found!: id_lep = " << idx_lep << std::endl;
        std::cout << "Lepton EOS tables can be generated by running 'make lep_table' in the parent directory" << std::endl;
        exit(EXIT_FAILURE);
      }

      // Get dataset sizes
      // -------------------------------------------------------------------------
      std::string EOSline;
      getline(EOSinput, EOSline);
      std::stringstream s1(EOSline);
      s1 >> m_nl;

      getline(EOSinput, EOSline);
      std::stringstream s2(EOSline);
      s2 >> m_ntl;

      // Allocate memory
      // -------------------------------------------------------------------------
      m_log_nl   = new double[m_nl];
      m_log_tl   = new double[m_ntl];
      m_lep_table = new double[m_nl*m_ntl];

      double number;

      // Read nb, t, yq
      // -------------------------------------------------------------------------
      getline(EOSinput, EOSline);
      std::stringstream s3(EOSline);
      for (int in = 0; in < m_nl; ++in) {
        s3 >> number;
        m_log_nl[in] = log(pow(10.,number));
      }
      m_id_log_nl = 1.0/(m_log_nl[1] - m_log_nl[0]);

      getline(EOSinput, EOSline);
      std::stringstream s4(EOSline);
      for (int it = 0; it < m_ntl; ++it) {
        s4 >> number;
        m_log_tl[it] = log(pow(10.,number));
      }
      m_id_log_tl = 1.0/(m_log_tl[1] - m_log_tl[0]);


      // Read other thermodynamics quantities
      // -------------------------------------------------------------------------   
      for (int in = 0; in < m_nl; ++in) {
        getline(EOSinput, EOSline);
        std::stringstream s5(EOSline);
          for (int it = 0; it < m_ntl; ++it) {
            s5 >> number;
            m_lep_table[lep_index(0, in, it)] = number;
        }
      }
    
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
    double eval_lep_at_nty(int vi, double n, double T, double Yl) {
      // @FIXME
      if (Yl == 0.) return 0.;
      return eval_lep_at_lnty(vi, log(n*Yl), log(T));
    }
    
    /// Evaluate interpolation weight for lepton density
    void weight_idx_lnl(double log_n) {
      in = (log_n - m_log_nl[0])*m_id_log_nl;
      wn1 = (log_n - m_log_nl[in])*m_id_log_nl;
      wn0 = 1.0 - wn1;
    }
    
    /// Evaluate interpolation weight for lepton temperature
    void weight_idx_ltl(double log_t) {
      it = (log_t - m_log_tl[0])*m_id_log_tl;
      wt1 = (log_t - m_log_tl[it])*m_id_log_tl;
      wt0 = 1.0 - wt1;
    }

    // @TODO: add action for interpolation outside the table range
    /// Low level evaluation function, not intended for outside use
    double eval_lep_at_lnty(int iv, double log_n, double log_t) {   
      weight_idx_lnl(log_n);
      weight_idx_ltl(log_t);
    
      return
        wn0 * (wt0 * m_lep_table[lep_index(iv, in+0, it+0)]   +
               wt1 * m_lep_table[lep_index(iv, in+0, it+1)])  +
        wn1 * (wt0 * m_lep_table[lep_index(iv, in+1, it+0)]   +
               wt1 * m_lep_table[lep_index(iv, in+1, it+1)]);
    }

   // @TODO: add action for interpolation outside the table range
    /// Low level evaluation function, not intended for outside use
    double interp_2d(int iv) const {
      return
        wn0 * (wt0 * m_lep_table[lep_index(iv, in+0, it+0)]   +
               wt1 * m_lep_table[lep_index(iv, in+0, it+1)])  +
        wn1 * (wt0 * m_lep_table[lep_index(iv, in+1, it+0)]   +
               wt1 * m_lep_table[lep_index(iv, in+1, it+1)]);
    }

    int in, it;
    double wn0, wn1, wt0, wt1;

    /// Low level evaluation function for entire EOS, not intended for outside use
    HelmEOSOutput eval_all_lep_at_lnty(double log_n, double log_t) const {


      std::array<double,EOS_leptons::NVARS> tmp;

      HelmEOSOutput out;

      weight_idx_lnl(log_n);
      weight_idx_ltl(log_t);

      for (int iv=0; iv<NVARS; iv++) {
        tmp[iv] = wn0 * (wt0 * m_lep_table[lep_index(iv, in+0, it+0)]   +
                         wt1 * m_lep_table[lep_index(iv, in+0, it+1)])  +
                  wn1 * (wt0 * m_lep_table[lep_index(iv, in+1, it+0)]   +
                         wt1 * m_lep_table[lep_index(iv, in+1, it+1)]);
      }

      out.nl   = tmp[0];
      out.a_nl = tmp[1];
      out.pl   = tmp[2];
      out.a_pl = tmp[3];
      out.el   = tmp[4];
      out.a_el = tmp[5];
      out.sl   = tmp[6];
      out.a_sl = tmp[7];
      out.mul  = tmp[8];

      return out;
    }

private:
  double e_lep, s_lep, P_lep;
  double n_lep;
  double mu_lep;
  double der_lep[4];

};

#endif // MUEOSLIBRARY_SRC_EOS_LEPTONS_HPP_
