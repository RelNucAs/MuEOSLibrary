#ifndef EOS_COMPOSE_H
#define EOS_COMPOSE_H

//! \file eos_compose.hpp
//  \brief Defines EOSTable, which stores information from a tabulated
//         equation of state in CompOSE format.
//
//  Tables should be generated using
//  <a href="https://bitbucket.org/dradice/pycompose">PyCompOSE</a>

///  \warning This code assumes the table to be uniformly spaced in
///           log nb, log t, and yq

#include <cstddef>
#include <string>

class EOS_baryons {
  public:
    enum TableVariables {
      BLOGP  = 0,  //! log (pressure / 1 MeV fm^-3)
      BENT   = 1,  //! entropy per baryon [kb]
      BMUB   = 2,  //! baryon chemical potential [MeV]
      BMUQ   = 3,  //! charge chemical potential [MeV]
      //ECMUL   = 4,  //! lepton chemical potential [MeV]
      BLOGE  = 4,  //! log (total energy density / 1 MeV fm^-3)
      BDPDN  = 5,  //! dP/dn [MeV]
      BDSDN  = 6,  //! ds/dn [fm^3]
      BDPDT  = 7,  //! dP/dT [fm^-3]
      BDSDT  = 8,  //! ds/dT [MeV^-1]
      BNVARS = 9
    };

  protected:
    /// Constructor
    EOS_baryons();

    /// Destructor
    ~EOS_baryons();

    /// Calculate the energy density using.
    double BarEnergy(double n, double T, double *Y);

    /// Calculate the pressure using.
    double BarPressure(double n, double T, double *Y);

    /// Calculate the entropy per baryon using.
    double BarEntropy(double n, double T, double *Y);

    /// Calculate the sound speed derivatives using.
    double BardPdn(double n, double T, double *Y);
    double Bardsdn(double n, double T, double *Y);
    double BardPdT(double n, double T, double *Y);
    double BardsdT(double n, double T, double *Y);

  public:
    /// Reads the table file.
    void ReadBarTableFromFile(std::string fname);

    /// Get the raw number density
    double const * GetRawLogNumberDensity() const {
      return m_log_nb;
    }
    double const * GetRawYq() const {
      return m_yq;
    }
    /// Get the raw number density
    double const * GetRawLogTemperature() const {
      return m_log_t;
    }
    /// Get the raw table data
    double const * GetRawTable() const {
      return m_table;
    }

    /// Get the baryon mass
    double const GetBaryonMass() const {
      return m_bar;
    }

    // Indexing used to access the data
    inline ptrdiff_t index(int iv, int in, int iy, int it) const {
      return it + m_nt*(iy + m_ny*(in + m_nn*iv));
    }

    /// Check if the EOS has been initialized properly.
    inline bool IsBarInitialized() const {
      return m_initialized;
    }

  private:
    /// Low level evaluation function, not intended for outside use
    double eval_at_nty(int vi, double n, double T, double Yq) const;
    /// Low level evaluation function, not intended for outside use
    double eval_at_lnty(int vi, double ln, double lT, double Yq) const;

  public:
    /// Evaluate interpolation weight for density
    void weight_idx_ln(double *w0, double *w1, int *in, double log_n) const;
    /// Evaluate interpolation weight for composition
    void weight_idx_yq(double *w0, double *w1, int *iy, double yq) const;
    /// Evaluate interpolation weight for temperature
    void weight_idx_lt(double *w0, double *w1, int *it, double log_t) const;

  public:
    // Inverse of table spacing
    double m_id_log_nb, m_id_log_t, m_id_yq;
    // Table size
    int m_nn, m_nt, m_ny;

    // Table storage, care should be made to store these data on the GPU later
    double * m_log_nb;
    double * m_log_t;
    double * m_yq;
    double * m_table;

  private:
    double m_bar;

    bool m_initialized;
};


#endif
