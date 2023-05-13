#ifndef EOS_LEPTONS_H
#define EOS_LEPTONS_H

//! \file eos_leptons.hpp
//  \brief Defines EOSTable, which stores information from a tabulated
//         equation of state

#include <array>
#include <cassert>
#include "eos_fermions.hpp"

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
  
  public:
    /// Constructor
    EOS_leptons();

    /// Destructor
    ~EOS_leptons();

    /// Calculate the number density of electrons using.
    template <int id_EOS>
    double ElectronNumberDensity(double n, double T, double *Y) {
        if constexpr(id_EOS == 1) {
            assert(m_el_initialized);
            return exp(eval_el_at_nty(IL_N, n, T, Y[id_e]));
        } else if constexpr(id_EOS == 2) {
            return eos_ferm_single<id_EOS,0>(n*Y[id_e], T)[IL_N];
        }
    }

    /// Calculate the number density of positrons using.
    template <int id_EOS>
    double PositronNumberDensity(double n, double T, double *Y) {
      //assert(m_el_initialized);
      //return exp(eval_el_at_nty(IA_N, n, T, Y[id_e]));
      return ElectronNumberDensity<id_EOS>(n, T, Y) - n*Y[id_e];
    }

    /// Calculate the number density of muons using.
    template <int id_EOS>
    double MuonNumberDensity(double n, double T, double *Y) {
      if (m_mu_initialized) {
        if constexpr(id_EOS == 1) {
            assert(m_mu_initialized);
            return exp(eval_mu_at_nty(IL_N, n, T, Y[id_mu]));
        } else if constexpr(id_EOS == 2) {
	    return eos_ferm_single<id_EOS,1>(n*Y[id_mu], T)[IL_N];
        }
      } else {
        return 0.;
      }
    }

    /// Calculate the number density of anti-muons using.
    template <int id_EOS>
    double AntiMuonNumberDensity(double n, double T, double *Y) {
      //assert(m_mu_initialized);
      //return exp(eval_mu_at_nty(IA_N, n, T, Y[id_mu]));
      if (m_mu_initialized) {
      return MuonNumberDensity<id_EOS>(n, T, Y) - n*Y[id_mu];
      } else {
        return 0.;
      }
    }
    
    /// Calculate the internal energy density of electrons + positrons using.
    template <int id_EOS>
    double ElectronEnergy(double n, double T, double *Y) {
	if constexpr(id_EOS == 1) {
            assert(m_el_initialized);
            return exp(eval_el_at_nty(IL_E, n, T, Y[id_e])) + exp(eval_el_at_nty(IA_E, n, T, Y[id_e]));
        } else if constexpr(id_EOS == 2) {
	    return eos_ferm_single<id_EOS,0>(n*Y[id_e], T)[IL_E] + eos_ferm_single<id_EOS,0>(n*Y[id_e], T)[IA_E]; 
        }
    }

    /// Calculate the internal energy density of muons + anti-muons using.
    template <int id_EOS>
    double MuonEnergy(double n, double T, double *Y) {
      if (m_mu_initialized) {
        if constexpr(id_EOS == 1) {
            assert(m_mu_initialized);
            return exp(eval_mu_at_nty(IL_E, n, T, Y[id_mu])) + exp(eval_mu_at_nty(IA_E, n, T, Y[id_mu]));
        } else if constexpr(id_EOS == 2) {
            return eos_ferm_single<id_EOS,1>(n*Y[id_mu], T)[IL_E] + eos_ferm_single<id_EOS,1>(n*Y[id_mu], T)[IA_E];
        }
      } else {
        return 0.;
      }
    }
    
    /// Calculate the internal energy density of leptons using.
    template <int id_EOS>
    double LepEnergy(double n, double T, double *Y) {
        return ElectronEnergy<id_EOS>(n, T, Y) + MuonEnergy<id_EOS>(n, T, Y);	    
    }

    /// Calculate the pressure of electrons + positrons using.
    template <int id_EOS>
    double ElectronPressure(double n, double T, double *Y) {
        if constexpr(id_EOS == 1) {
            assert(m_el_initialized);
            return exp(eval_el_at_nty(IL_P, n, T, Y[id_e])) + exp(eval_el_at_nty(IA_P, n, T, Y[id_e]));
        } else if constexpr(id_EOS == 2) {
            return eos_ferm_single<id_EOS,0>(n*Y[id_e], T)[IL_P] + eos_ferm_single<id_EOS,0>(n*Y[id_e], T)[IA_P];
        }
    }

    /// Calculate the pressure of muons + anti-muons using.
    template <int id_EOS>
    double MuonPressure(double n, double T, double *Y) {
      if (m_mu_initialized) {
        if constexpr(id_EOS == 1) {
            assert(m_mu_initialized);
            return exp(eval_mu_at_nty(IL_P, n, T, Y[id_mu])) + exp(eval_mu_at_nty(IA_P, n, T, Y[id_mu]));
        } else if constexpr(id_EOS == 2) {
            return eos_ferm_single<id_EOS,1>(n*Y[id_mu], T)[IL_P] + eos_ferm_single<id_EOS,1>(n*Y[id_mu], T)[IA_P];
        }
      } else {
        return 0.;
      }
    }
    
    /// Calculate the pressure of leptons using.
    template <int id_EOS>
    double LepPressure(double n, double T, double *Y) {
        return ElectronPressure<id_EOS>(n, T, Y) + MuonPressure<id_EOS>(n, T, Y);
    }

    /// Calculate the entropy density of electrons + positrons using.
    template <int id_EOS>
    double ElectronEntropy(double n, double T, double *Y) {
        if constexpr(id_EOS == 1) {
            assert(m_el_initialized);
            return exp(eval_el_at_nty(IL_S, n, T, Y[id_e])) + exp(eval_el_at_nty(IA_S, n, T, Y[id_e]));
        } else if constexpr(id_EOS == 2) {
            return eos_ferm_single<id_EOS,0>(n*Y[id_e], T)[IL_S] + eos_ferm_single<id_EOS,0>(n*Y[id_e], T)[IA_S];
        }
    }

    /// Calculate the entropy density of muons + anti-muons using.
    template <int id_EOS>
    double MuonEntropy(double n, double T, double *Y) {
      if (m_mu_initialized) {
        if constexpr(id_EOS == 1) {
            assert(m_mu_initialized);
            return exp(eval_mu_at_nty(IL_S, n, T, Y[id_mu])) + exp(eval_mu_at_nty(IA_S, n, T, Y[id_mu]));
        } else if constexpr(id_EOS == 2) {
            return eos_ferm_single<id_EOS,1>(n*Y[id_mu], T)[IL_S] + eos_ferm_single<id_EOS,1>(n*Y[id_mu], T)[IA_S];
        }
      } else {
        return 0.;
      }
    }

    /// Calculate the entropy density of leptons using.
    template <int id_EOS>
    double LepEntropy(double n, double T, double *Y) {
        return ElectronEntropy<id_EOS>(n, T, Y) + MuonEntropy<id_EOS>(n, T, Y);
    }

    /// Calculate the chemical potential of electrons using.
    template <int id_EOS>
    double ElectronChemicalPotential(double n, double T, double *Y) {
        if constexpr(id_EOS == 1) {
            assert(m_el_initialized);
            return eval_el_at_nty(IL_MU, n, T, Y[id_e]);
        } else if constexpr(id_EOS == 2) {
            return eos_ferm_single<id_EOS,0>(n*Y[id_e], T)[IL_MU];
        }
    }

    /// Calculate the chemical potential of muons using.
    template <int id_EOS>
    double MuonChemicalPotential(double n, double T, double *Y) {
      if (m_mu_initialized) {
        if constexpr(id_EOS == 1) {
            assert(m_mu_initialized);
            return eval_mu_at_nty(IL_MU, n, T, Y[id_mu]);
        } else if constexpr(id_EOS == 2) {
            return eos_ferm_single<id_EOS,1>(n*Y[id_mu], T)[IL_MU];
        }
      } else {
        return 0.;
      }
    }
    /// Calculate the sound speed derivatives using.
    double EdPdn(double n, double T, double *Y);
    double Edsdn(double n, double T, double *Y);
    double EdPdt(double n, double T, double *Y);
    double Edsdt(double n, double T, double *Y);

    double MdPdn(double n, double T, double *Y);
    double Mdsdn(double n, double T, double *Y);
    double MdPdt(double n, double T, double *Y);
    double Mdsdt(double n, double T, double *Y);
    
    /// Calculate the full leptonic EOS using.
    std::array<double,NVARS> ComputeFullElectronEOS(double n, double T, double *Y);
    std::array<double,NVARS> ComputeFullMuonEOS(double n, double T, double *Y);

  public:
    /// Reads the table file.
    void ReadETableFromFile(std::string fname);
    void ReadMTableFromFile(std::string fname);

    /// Get the raw number density
    double const * GetRawELogNumberDensity() const {
      return m_log_ne;
    }
    double const * GetRawMLogNumberDensity() const {
      return m_log_nm;
    }

    /// Get the raw table temperature
    double const * GetRawELogTemperature() const {
      return m_log_te;
    }
    double const * GetRawMLogTemperature() const {
      return m_log_tm;
    }

    /// Get the raw table data
    double const * GetElectronTable() const {
      return m_el_table;
    }

    double const * GetMuonTable() const {
      return m_mu_table;
    }

    // Indexing used to access the lepton data
    inline ptrdiff_t el_index(int iv, int in, int it) const {
      return it + m_nte*(in + m_ne*iv);
    }
    
    inline ptrdiff_t mu_index(int iv, int in, int it) const {
      return it + m_ntm*(in + m_nm*iv);
    }


    /// Check if the EOS has been initialized properly.
    inline bool IsElectronInitialized() const {
      return m_el_initialized;
    }

    inline bool IsMuonInitialized() const {
      return m_mu_initialized;
    }

  private:
    /// Low level evaluation function, not intended for outside use
    double eval_el_at_nty(int vi, double n, double T, double Yl) const;
    double eval_mu_at_nty(int vi, double n, double T, double Yl) const;
    /// Low level evaluation function, not intended for outside use
    double eval_el_at_lnty(int vi, double ln, double lT) const;
    double eval_mu_at_lnty(int vi, double ln, double lT) const;

    /// Evaluate interpolation weight for lepton density
    void weight_idx_lne(double *w0, double *w1, int *in, double log_n) const;
    void weight_idx_lnm(double *w0, double *w1, int *in, double log_n) const;
    /// Evaluate interpolation weight for lepton temperature
    void weight_idx_lte(double *w0, double *w1, int *it, double log_t) const;
    void weight_idx_ltm(double *w0, double *w1, int *it, double log_t) const;

    /// Low level evaluation function for entire EOS, not intended for outside use
    std::array<double,NVARS> eval_all_el_at_lnt(double ln, double lT) const;
    std::array<double,NVARS> eval_all_mu_at_lnt(double ln, double lT) const;
 
    double mu_check(const bool check, const double value);

  private:
    // Inverse of table spacing
    double m_id_log_ne, m_id_log_te;
    double m_id_log_nm, m_id_log_tm;

    // Table storage, care should be made to store these data on the GPU later
    double * m_log_ne;
    double * m_log_te;
    double * m_el_table;
    
    double * m_log_nm;
    double * m_log_tm;
    double * m_mu_table;
    
    const int id_e  = 0;
    const int id_mu = 1;

    bool m_el_initialized;
    bool m_mu_initialized;

    bool m_mu_active;
  public:
    // Table size
    int m_ne, m_nte;
    int m_nm, m_ntm;

};


#endif
