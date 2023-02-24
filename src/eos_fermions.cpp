#include <vector>
#include <array>
#include <sstream>

#include "complete_FG.hpp"
#include "find_eta.hpp"
#include "interp.hpp"
#include "eos_fermions.hpp"

std::array<double,9> eos_ferm_onthefly(const double eta, const double T, const int id_L) {
        std::array<double,9> eos_array;

//.......Define theta (relativity parameter)
        const double theta = T/mL[id_L];

//.......Degeneracy parameter of anti-leptons
	const double a_eta = -(eta+2./theta);

//.......Compute Generalized Fermi-Dirac integrals for leptons
	const double f12 = compute_res(eta, theta, 0.5); //k = 1/2
	const double f32 = compute_res(eta, theta, 1.5); //k = 3/2
	const double f52 = compute_res(eta, theta, 2.5); //k = 5/2

//.......Compute Generalized Fermi-Dirac integrals for anti-leptons
	const double a_f12 = compute_res(a_eta, theta, 0.5); //k = 1/2
	const double a_f32 = compute_res(a_eta, theta, 1.5); //k = 3/2
	const double a_f52 = compute_res(a_eta, theta, 2.5); //k = 5/2
	
//.......Number density of leptons
	eos_array[IL_N] = K[id_L]*pow(theta,1.5)*(f12+theta*f32);

//.......Number density of anti-leptons
	eos_array[IA_N] = K[id_L]*pow(theta,1.5)*(a_f12+theta*a_f32);

//.......Pressure of leptons
	eos_array[IL_P] = K3[id_L]*mL[id_L]*pow(theta,2.5)*(2.*f32+theta*f52); //MeV/fm^3

//.......Pressure of anti-leptons
	eos_array[IA_P] = K3[id_L]*mL[id_L]*pow(theta,2.5)*(2.*a_f32+theta*a_f52); //MeV/cm^3

//.......Internal energy density of leptons
	eos_array[IL_E] = K[id_L]*mL[id_L]*pow(theta,1.5)*(f12+2.*theta*f32+theta*theta*f52); //MeV/fm^3

//.......Internal energy density of anti-leptons
	eos_array[IA_E] = K[id_L]*mL[id_L]*pow(theta,1.5)*(a_f12+2.*theta*a_f32+theta*theta*a_f52); //MeV/fm^3

//.......Entropy density of leptons
	eos_array[IL_S] = K[id_L]*pow(theta,1.5)*(-eta*f12+(5./3.-eta*theta)*f32+4./3.*theta*f52); //1/fm^3

//.......Entropy density of anti-leptons
	eos_array[IA_S] = K[id_L]*pow(theta,1.5)*(-a_eta*a_f12+(5./3.-a_eta*theta)*a_f32+4./3.*theta*a_f52); //1/fm^3

//.......Chemical potential of leptons rnal energy density of anti-leptons
	eos_array[IL_MU] = mL[id_L] + eta*T;

	return eos_array;
}
