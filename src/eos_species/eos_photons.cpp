#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

#include <tr1/cmath>

#include "eos_species.hpp"
#include "../constants.hpp"

using namespace constants;

const double zeta3 = std::tr1::riemann_zeta(3.);
const double c_fms = c * 1.0E+13;

//number density
double EOS_photons::RadNumberDensity(double T){ //T in MeV
  return 16.*pi*zeta3*pow(T/(h*c_fms),3.); //fm^-3
}		

//internal energy (per unit volume)
double EOS_photons::RadEnergy(double T){ //T in MeV
	return (8.*pow(pi,5.)*pow(T,4.)) / (15*pow(h*c_fms,3.)); //MeV fm^-3	
}

//pressure 
double EOS_photons::RadPressure(double T){ //T in MeV
	return RadEnergy(T)/3.; //MeV cm^-3
}
	
//entropy
double EOS_photons::RadEntropy(double T){ //T in MeV
	return (32.*pow(pi,5.)/45.) * pow(T/(h*c_fms),3.); //fm^-3
}

//pressure derivative wrt temperature
double EOS_photons::RaddPdT(double T){ //T in MeV
	return 4. * RadPressure(T) / T;
}

//entropy derivative wrt temperature
double EOS_photons::RaddsdT(double T){ //T in MeV
	return 3. * RadEntropy(T) / T;
}
