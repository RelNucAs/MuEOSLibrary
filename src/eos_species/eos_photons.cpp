#include "eos_species/eos_species.hpp"
#include "constants.hpp"

constexpr double zeta3 = 1.202056903159594; //std::tr1::riemann_zeta(3.);

void EOS_photons::PhotonEOS(const double T) { //T in MeV
  constexpr double nph_const = 4. * MEOS_fourpi_hc3 * zeta3; // 16 * pi * zeta(3) / (hc)**3 [MeV-3 fm-3]
  constexpr double eph_const = 12.987878804533656 * MEOS_fourpi_hc3; // 8 * pi**5 / (15 * (hc)**3) [MeV-3 fm-3]
  constexpr double sph_const = 17.317171739378207 * MEOS_fourpi_hc3; // 32 * pi**5 / (45 * (hc)**3) [MeV-3 fm-3]

  
  const double T_third =  POW3(T);
  const double T_fourth = T * T_third;

  n_ph = nph_const * T_third; // number density [fm^-3]
  e_ph = eph_const * T_fourth; // internal energy (per unit volume) [MeV fm^-3]	
  P_ph = e_ph / 3.; // pressure [MeV cm^-3]
  s_ph = sph_const * T_third; // entropy [fm^-3]
  
  dPdt_ph = 4. * P_ph / T; //pressure derivative wrt temperature
  dsdt_ph = 3. * s_ph / T; //entropy derivative wrt temperature

  return;
  
}

//number density
double EOS_photons::GetPhotonNumberDensity() {
	return n_ph;
}		

//internal energy (per unit volume)
double EOS_photons::GetPhotonEnergy(){ //T in MeV
	return e_ph;
}

//pressure 
double EOS_photons::GetPhotonPressure(){ //T in MeV
	return P_ph; //MeV cm^-3
}
	
//entropy
double EOS_photons::GetPhotonEntropy(){ //T in MeV
	return s_ph; //fm^-3
}

//pressure derivative wrt temperature
double EOS_photons::GetPhotondPdT(){ //T in MeV
	return dPdt_ph;
}

//entropy derivative wrt temperature
double EOS_photons::GetPhotondsdT(){ //T in MeV
	return dsdt_ph;
}
