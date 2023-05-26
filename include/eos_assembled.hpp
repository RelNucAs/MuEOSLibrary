#ifndef EOS_ASSEMBLED_H
#define EOS_ASSEMBLED_H

//! \file eos_assembled.hpp

#include <cstddef>
#include <string>

#include <eos_leptons.hpp>
#include <eos_baryons.hpp>

constexpr int id_test = 2;

class EOS_assembled : public EOS_baryons, public EOS_leptons {
    
    public:
    /// Constructor
    //EOS_assembled();

    /// Destructor
    //~EOS_assembled();

    /// Temperature from energy density
    double TemperatureFromE(double n, double e, double *Y);

    /// Temperature from pressure
    double TemperatureFromP(double n, double p, double *Y);

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

    public:
    /// Low level function, not intended for outside use
    double temperature_from_e(double var, double n, double *Y); //const;
    double temperature_from_p(double var, double n, double *Y); //const;
};

#endif
