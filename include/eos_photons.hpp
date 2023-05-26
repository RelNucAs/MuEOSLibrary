#ifndef EOS_RADIATION_H
#define EOS_RADIATION_H

#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

#include <tr1/cmath>

namespace radiation {

const double zeta3 = std::tr1::riemann_zeta(3.);
const double pi = std::acos(-1.);
const double h = 4.1356943e-21; // MeV*s
const double c = 2.997924562e+23; // fm/s

//number density
double RadNumberDensity(double T){ //T in MeV
      return 16.*pi*zeta3*pow(T/(h*c),3.); //fm^-3
}

//internal energy (per unit volume)
double RadEnergy(double T){ //T in MeV
	return (8.*pow(pi,5.)*pow(T,4.)) / (15*pow(h*c,3.)); //MeV fm^-3	
}

//pressure 
double RadPressure(double T){ //T in MeV
	return RadEnergy(T)/3.; //MeV cm^-3
}
	
//entropy
double RadEntropy(double T){ //T in MeV
	return (32.*pow(pi,5.)/45.) * pow(T/(h*c),3.); //fm^-3
}

//pressure derivative wrt temperature
double RaddPdT(double T){ //T in MeV
	return 4. * RadPressure(T) / T;
}

//entropy derivative wrt temperature
double RaddsdT(double T){ //T in MeV
	return 3. * RadEntropy(T) / T;
}

}
#endif
