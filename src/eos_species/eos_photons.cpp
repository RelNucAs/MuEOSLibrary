#include "eos_species/eos_species.hpp"
#include "constants.hpp"

constexpr double zeta3 = 1.202056903159594; //std::tr1::riemann_zeta(3.);

//number density
double EOS_photons::RadNumberDensity(double T){ //T in MeV
  constexpr double nph_const = 4. * MEOS_fourpi_hc3 * zeta3; // 16 * pi * zeta(3) / (hc)**3 [MeV-3 fm-3]
  return nph_const * POW3(T); //fm^-3
}		

//internal energy (per unit volume)
double EOS_photons::RadEnergy(double T){ //T in MeV
	constexpr double eph_const = 12.987878804533656 * MEOS_fourpi_hc3; // 8 * pi**5 / (15 * (hc)**3) [MeV-3 fm-3]
	return eph_const * POW4(T); //MeV fm^-3	
}

//pressure 
double EOS_photons::RadPressure(double T){ //T in MeV
	return RadEnergy(T) / 3.; //MeV cm^-3
}
	
//entropy
double EOS_photons::RadEntropy(double T){ //T in MeV
	constexpr double sph_const = 17.317171739378207 * MEOS_fourpi_hc3; // 32 * pi**5 / (45 * (hc)**3) [MeV-3 fm-3]
	return sph_const * POW3(T); //fm^-3
}

//pressure derivative wrt temperature
double EOS_photons::RaddPdT(double T){ //T in MeV
	return 4. * RadPressure(T) / T;
}

//entropy derivative wrt temperature
double EOS_photons::RaddsdT(double T){ //T in MeV
	return 3. * RadEntropy(T) / T;
}
