#include "constants.hpp"
#include "parameters.hpp"
#include "complete_FG.hpp"
#include "find_eta.hpp"

using namespace constants;
using namespace parameters;

//.......Given the density rho, the non relativistic degeneracy eta, the temperature T
//.......and the particle species, the subroutine computes the net fraction ynet and the first eta derivative.
//.......UNITS: rho in g/cm**3; T in MeV
double n_net_f(const double eta, const double T, const double mLep) {
        const double K = 8.*sqrt(2.)*pi*pow(mLep/(h*c),3.); //31217845.162531383*mLep**3
        const double theta = T/mLep;
        double f1, f2, f3, f4;

        f1 = compute_res( eta,          theta, 0.5);
        f2 = compute_res(-eta-2./theta, theta, 0.5);
        f3 = compute_res( eta,          theta, 1.5);
        f4 = compute_res(-eta-2./theta, theta, 1.5);

        return  K*pow(theta,1.5)*(f1-f2+theta*(f3-f4));
}


double n_net_df(const double eta, const double T, const double mLep) {
        const double K = 8.*sqrt(2.)*pi*pow(mLep/(h*c),3.); //31217845.162531383*mLep**3
        const double theta = T/mLep;
        double f1, f2, f3, f4;

        f1 = compute_res_ed( eta,          theta, 0.5);
        f2 = compute_res_ed(-eta-2./theta, theta, 0.5);
        f3 = compute_res_ed( eta,          theta, 1.5);
        f4 = compute_res_ed(-eta-2./theta, theta, 1.5);

        return  K*pow(theta,1.5)*(f1-f2+theta*(f3-f4));
}

void n_net_fdf(const double eta, const double T, const double mLep, struct n_net_FD &FD) {
	const double K = 8.*sqrt(2.)*pi*pow(mLep/(h*c),3.);
	const double theta = T/mLep;
	double f1, f2, f3, f4;

	FD.f12   = compute_res( eta,          theta, 0.5);
	FD.f32   = compute_res( eta,          theta, 1.5);
	FD.a_f12 = compute_res(-eta-2./theta, theta, 0.5);
	FD.a_f32 = compute_res(-eta-2./theta, theta, 1.5);

        f1 = compute_res_ed( eta,          theta, 0.5);
        f2 = compute_res_ed(-eta-2./theta, theta, 0.5);
        f3 = compute_res_ed( eta,          theta, 1.5);
        f4 = compute_res_ed(-eta-2./theta, theta, 1.5);

	FD.n_net  = K*pow(theta,1.5)*(FD.f12-FD.a_f12+theta*(FD.f32-FD.a_f32));
	FD.dn_net = K*pow(theta,1.5)*(f1-f2+theta*(f3-f4));
	return;
}
