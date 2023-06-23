#ifndef EOS_ASSEMBLED_H
#define EOS_ASSEMBLED_H

//! \file eos_assembled.hpp

#include <cstddef>
#include <string>

#include <eos_photons.hpp>
#include <eos_leptons.hpp>
#include <eos_baryons.hpp>

class EOS_assembled : public EOS_baryons, public EOS_leptons<0>, public EOS_leptons<1>, public EOS_photons {
    
    public:
    /// Constructor
    //EOS_assembled();

    /// Destructor
    //~EOS_assembled();

    /// Calculate the energy density using.
    double Energy(double n, double T, double *Y);

    /// Calculate the pressure using.
    double Pressure(double n, double T, double *Y);

    /// Calculate the entropy per baryon using.
    double Entropy(double n, double T, double *Y);

    /// Calculate the enthalpy per baryon using.
    double Enthalpy(double n, double T, double *Y);

    /// Calculate the sound speed.
    double SoundSpeed(double n, double T, double *Y);

    /// Read in tables
    void ReadTables(std::string BarTableName,
                    std::string ETableName,
                    std::string MTableName);

};

#endif
