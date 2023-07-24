#ifndef EOS_SPECIES_HPP
#define EOS_SPECIES_HPP

#include <cstddef>
#include <string>

/*============================================================================*/

// file: eos_neutrinos.cpp

class EOS_neutrinos {
  public:
    
    // Number fraction of neutrinos (yn)
    double nu_fraction(const double nb, const double T, const double eta);

    // Number fraction of antineutrinos (yan)
    double anu_fraction(const double nb, const double T, const double eta);

    // Energy per baryon of neutrinos (zn)
    double nu_energy(const double nb, const double T, const double eta);

    // Energy per baryon of antineutrinos (zan)
    double anu_energy(const double nb, const double T, const double eta);

    // Pressure per baryon of neutrinos (pn)
    double nu_pressure(const double nb, const double T, const double eta);

    // Pressure per baryon of antineutrinos (pan)
    double anu_pressure(const double nb, const double T, const double eta);

    // Entropy per baryon of neutrinos (sn)
    double nu_entropy(const double nb, const double T, const double eta);

    // Entropy per baryon of antineutrinos (san)
    double anu_entropy(const double nb, const double T, const double eta);
};

/*============================================================================*/

// file: eos_photons.cpp
class EOS_photons {
  protected:
    //number density
    double RadNumberDensity(double T);

    //internal energy (per unit volume)
    double RadEnergy(double T);

    //pressure 
    double RadPressure(double T);

    //entropy
    double RadEntropy(double T);

    //pressure derivative wrt temperature
    double RaddPdT(double T);

    //entropy derivative wrt temperature
    double RaddsdT(double T);
};

/*============================================================================*/

// file: eos_baryons.cpp

class EOS_baryons {
  public:
    enum TableVariables {
      BLOGP  = 0,  //! log (pressure / 1 MeV fm^-3)
      BENT   = 1,  //! entropy per baryon [kb]
      BMUB   = 2,  //! baryon chemical potential [MeV]
      BMUQ   = 3,  //! charge chemical potential [MeV]
      //ECMUL   = 4,  //! lepton chemical potential [MeV]
      BLOGE  = 4,  //! log (total energy density / 1 MeV fm^-3)
      BYALP  = 5,  //! alpha fraction 
      BYNUC  = 6,  //! heavy nuclei fraction 
      BYNTR  = 7,  //! neutron fraction 
      BYPTN  = 8,  //! proton fraction 
      BDPDN  = 9,  //! dP/dn [MeV]
      BDSDN  = 10, //! ds/dn [fm^3]
      BDPDT  = 11, //! dP/dT [fm^-3]
      BDSDT  = 12, //! ds/dT [MeV^-1]
      BNVARS = 13
    };

  public:
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

    /// Calculate the relativistic proton chemical potential using.
    double ProtonChemicalPotential(double n, double T, double *Y);

    /// Calculate the relativistic neutron chemical potential using.
    double NeutronChemicalPotential(double n, double T, double *Y);

    /// Calculate the fraction of alpha particles using.
    double AlphaFraction(double n, double T, double *Y);

    /// Calculate the fraction of heavy nuclei using.
    double HeavyFraction(double n, double T, double *Y);

    /// Calculate the fraction of free neutrons using.
    double NeutronFraction(double n, double T, double *Y);

    /// Calculate the fraction of free protons using.
    double ProtonFraction(double n, double T, double *Y);

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
    /// Low level evaluation function, not intended for outside use
    double eval_at_nty_new(int vi, double n, double T, double Yq) const;

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

    // Minimum Pressure
    double Pmin;

    // Table storage, care should be made to store these data on the GPU later
    double * m_log_nb;
    double * m_log_t;
    double * m_yq;
    double * m_table;

  private:
    double m_bar;

    bool m_initialized;
};

/*============================================================================*/

#endif //EOS_SPECIES_HPP