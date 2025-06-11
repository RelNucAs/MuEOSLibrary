#ifndef MUEOSLIBRARY_SRC_EOS_SPECIES_HPP_
#define MUEOSLIBRARY_SRC_EOS_SPECIES_HPP_

#include <cstddef>
#include <string>

#include "mueosclass.hpp"
#include "eos_species/eos_leptons.hpp"

/*============================================================================*/

// file: eos_neutrinos.cpp

class EOS_neutrinos {
  private:
    double Y_nu[5], Z_nu[5];
    double P_nu[5], s_nu[5];
    double eta_nu[5];

  public:
    EOS_neutrinos() {
      eta_nu[4] = 0;
    };

    void NeutrinoEOS(const double nb, const double T, const double* chem_pot);

    // Number fraction of neutrinos
    double GetNuNumberFraction(const int idx);

    // Energy per baryon of neutrinos
    double GetNuEnergyFraction(const int idx);

    // Pressure per baryon of neutrinos
    double GetNuPressure(const int idx);

    // Entropy per baryon of neutrinos
    double GetNuEntropy(const int idx);

    // Degeneracy parameter of neutrinos
    double GetNuDegeneracyParameter(const int idx);

  public:
    bool m_mu_active;
};

/*============================================================================*/

// file: eos_photons.cpp
class EOS_photons {
  public:
    void PhotonEOS(const double T);

    //number density
    double GetPhotonNumberDensity();

    //internal energy (per unit volume)
    double GetPhotonEnergy();

    //pressure 
    double GetPhotonPressure();

    //entropy
    double GetPhotonEntropy();

    //pressure derivative wrt temperature
    double GetPhotondPdT();

    //entropy derivative wrt temperature
    double GetPhotondsdT();

   private:
    double n_ph, e_ph, P_ph, s_ph;
    double dPdt_ph, dsdt_ph;
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

    /// Calculate full baryonic EOS using.
    void BaryonEOS(double n, double T, double *Y);

     /// Calculate the energy density using.
    double GetBaryonEnergy();

    /// Calculate the pressure using.
    double GetBaryonPressure();

    /// Calculate the entropy per baryon using.
    double GetBaryonEntropy();

    /// Calculate the relativistic proton chemical potential using.
    double GetProtonChemicalPotential();

    /// Calculate the relativistic neutron chemical potential using.
    double GetNeutronChemicalPotential();

    /// Calculate the fraction of alpha particles using.
    double GetAlphaFraction();

    /// Calculate the fraction of heavy nuclei using.
    double GetHeavyFraction();

    /// Calculate the fraction of free neutrons using.
    double GetNeutronFraction();

    /// Calculate the fraction of free protons using.
    double GetProtonFraction();

    /// Calculate the sound speed derivatives using.
    double GetBaryondPdn();
    double GetBaryondsdn();
    double GetBaryondPdT();
    double GetBaryondsdT(); 

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
    // @TODO: add action for interpolation outside the table range
    /// Low level evaluation function, not intended for outside use
    double eval_at_nty(int vi, double n, double T, double Yq);
    /// Low level evaluation function, not intended for outside use
    double eval_at_lnty(int vi, double ln, double lT, double Yq);
    /// Low level evaluation function, not intended for outside use
    double eval_at_nty_new(int vi, double n, double T, double Yq) const;
    double interp_3d(int iv) const;

  private:
    /// Evaluate interpolation weight for density
    void weight_idx_ln(double log_n);
    /// Evaluate interpolation weight for composition
    void weight_idx_yq(double yq);
    /// Evaluate interpolation weight for temperature
    void weight_idx_lt(double log_t);

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
    int in, iy, it;
    double wn0, wn1, wy0, wy1, wt0, wt1;

    double e_bar, P_bar, s_bar;
    double mu_bar[2];
    double Y_comp[4];
    double der_bar[4];

    bool m_initialized;
};

/*============================================================================*/

class MuEOSClass : public EOS_baryons, public EOS_leptons<0>, public EOS_leptons<1>, public EOS_photons, public EOS_neutrinos {

    public:
    /// Constructor
    MuEOSClass(const int id_eos, const bool mu_flag, std::string BarTableName);

    /// Destructor
    ~MuEOSClass();

/*     /// Calculate the energy density using.
    double Energy(double n, double T, double *Y);

    /// Calculate the pressure using.
    double Pressure(double n, double T, double *Y);

    /// Calculate the entropy per baryon using.
    double Entropy(double n, double T, double *Y);

    /// Calculate the enthalpy per baryon using.
    double Enthalpy(double n, double T, double *Y); 
    
    /// Calculate the sound speed using.
    double SoundSpeed(double n, double T, double *Y);

     /// Calculate the chemical potentials using.
    ChemPotentials GetChemicalPotentials(double n, double T, double *Y);

    /// Calculate the particle fractions using.
    ParticleFractions GetParticleFractions(double n, double T, double *Y);

    /// Calculate the fractions of baryons using.
    ParticleFractions BaryonFractions(double n, double T, double *Y); */

    /// Calculate the neutrino particle fractions using.
    ParticleFractions NeutrinoFractions(double n, double T, double *Y, ChemPotentials *chem_pot);

    /// Calculate the neutrino EOS quantities using.
    NeutrinoEOSOutput compute_neutrino_EOS(double n, double T, double *Y, double chem_pot[4]);

    /// Calculate the full EOS output using.
    EOSstruct compute_full_EOS(double n, double T, double *Y);
};

#endif // MUEOSLIBRARY_SRC_EOS_SPECIES_HPP_
