#ifndef EOS_FERM_H
#define EOS_FERM_H

#include <iostream>

#include "find_eta.hpp"
#include "interp.hpp"
#include "constants.hpp"
#include "complete_FG.hpp"
#include "eos_leptons.hpp"

using namespace constants;

enum EOSQuantities {
	IL_N  = 0, //! Number density of leptons
	IA_N  = 1, //! Number density of anti-leptons
	IL_P  = 2, //! Pressure of leptons
	IA_P  = 3, //! Pressure of anti-leptons
	IL_E  = 4, //! Internal energy density of leptons
	IA_E  = 5, //! Internal energy density of anti-leptons
	IL_S  = 6, //! Entropy density of leptons
	IA_S  = 7, //! Entropy density of anti-leptons
	IL_MU = 8, //! Chemical potential of leptons
	NVARS = 9
};

const double K0    = 8.*sqrt(2.)*pi/pow(h*c*1.e13,3.);
const double mL[2] = {me, mmu};
const double K[2]  = {K0*pow(mL[0],3.), K0*pow(mL[1],3.)};
const double K3[2] = {K[0]/3., K[1]/3.};


std::array<double,9> eos_ferm_onthefly(const double eta, const double T, const int species);

template <int species>
std::array<double,9> eos_ferm_fromNR(const double nLep, const double temp) {
        struct n_net_FD FDout;
        const double guess = find_guess_eta<species>(1.e39*nLep, temp);
        const double eta = rtsafe_mod<species>(1.e39*nLep, temp, guess, FDout);

        double f12, f32, f52, a_f12, a_f32, a_f52;

//.......Define theta (relativity parameter)
        const double theta = temp/mL[species];

        f12 = FDout.f12; //compute_res(eta, theta, 0.5);
        f32 = FDout.f32; //compute_res(eta, theta, 1.5);
        f52 = compute_res(eta, theta, 2.5);

        const double a_eta = -(eta+2./theta);

        a_f12 = FDout.a_f12; //compute_res(a_eta, theta, 0.5);
        a_f32 = FDout.a_f32; //compute_res(a_eta, theta, 1.5);
        a_f52 = compute_res(a_eta, theta, 2.5);

        std::array<double,9> eos_array;

        //.......Compute eos_array[0]
        //eos_array[0] = n_part(T, eta, mLep);
        eos_array[0] = K[species]*pow(theta,1.5)*(f12+theta*f32);

        //.......Compute eos_array[1]
        //eos_array[1] = n_antipart(T, eta, mLep);
        eos_array[1] = K[species]*pow(theta,1.5)*(a_f12+theta*a_f32);

        //.......Compute eos_array[2]
        //eos_array[2] = P_part(T, eta, mLep);
        eos_array[2] = K3[species]*mL[species]*MeV*pow(theta,2.5)*(2.*f32+theta*f52); //erg/cm^3

        //.......Compute eos_array[3]
        //eos_array[3] = P_antipart(T, eta, mLep)
        eos_array[3] = K3[species]*mL[species]*MeV*pow(theta,2.5)*(2.*a_f32+theta*a_f52); //erg/cm^3

        //.......Compute eos_array[4]
        //eos_array[4] = e_part(T, eta, mLep);
        eos_array[4] = K[species]*mL[species]*MeV*pow(theta,1.5)*(f12+2.*theta*f32+theta*theta*f52); //erg/cm^3

        //.......Compute eos_array[5]
        //eos_array[5] = e_antipart(T, eta, mLep);
        eos_array[5] = K[species]*mL[species]*MeV*pow(theta,1.5)*(a_f12+2.*theta*a_f32+theta*theta*a_f52); //erg/cm^3

        //.......Compute eos_array[6]
        //eos_array[6] = s_part(T, eta, mLep);
        eos_array[6] = K[species]*kB*MeV*pow(theta,1.5)*(-eta*f12+(5./3.-eta*theta)*f32+4./3.*theta*f52); //erg/K/cm^3

        //.......Compute eos_array[7]
        //eos_array[7] = s_antipart(T, eta, mLep);
        eos_array[7] = K[species]*kB*MeV*pow(theta,1.5)*(-a_eta*a_f12+(5./3.-a_eta*theta)*a_f32+4./3.*theta*a_f52); //erg/K/cm^3

        //.......Compute eos_array[8]
        eos_array[8] = eta;

        return eos_array;
}

template<int eos_method, int species>
std::array<double,9> eos_ferm_array(const double nLep, const double temp, struct EOSeta &eta_table, struct EOScomplete &EOS_table) {
	//!........................INPUT........................................................
//!...... matter density rho [g/cm^3]; net particle fraction y_tilde; temperature T [MeV];
//!.......species 1 for electrons, 2 for muons; interval for non relat. degen. parameter [eta1,eta2]
//!.......el_interp = .true. to interpolate electrons, el_interp = .false. to use newrap 1D
//!.......mu_interp = .true. to interpolate muons, mu_interp = .false. to use newrap 1D
//!.....................................................................................
//!.......................OUTPUT........................................................
//!.......eos_array(1) = number density of particles [1/cm^3]...........................
//!.......eos_array(2) = number density of antiparticles [1/cm^3].......................
//!.......eos_array(3) = pressure of particles [erg/cm^3]...............................
//!.......eos_array(4) = pressure of antiparticles [erg/cm^3]...........................
//!.......eos_array(5) = energy of particles [erg/g]....................................
//!.......eos_array(6) = energy of antiparticles [erg/g]................................
//!.......eos_array(7) = entropy of particles [erg/K/g].................................
//!.......eos_array(8) = entropy of antiparticles [erg/K/g].............................
//!.......eos_array(9) = non relat. degeneracy parameter of particles...................
//!.......For a complete description of output variables, see subroutines above.........

        //double eta1_Lep, eta2_Lep;
        //bool guess_from_interp = false;


        if constexpr(eos_method == 3) {
                std::array<double,13> tmp = eos_interp(nLep, temp, EOS_table);
                tmp[8] = (tmp[8] - mL[species]) / temp;
		std::array<double,9> eos_array;
		for (int i=0;i<9;i++) eos_array[i] = tmp[i];
                return eos_array;
        } else if constexpr(eos_method == 1) {
                double guess = find_guess_eta<species>(1.e39*nLep, temp);
		std::cout << 1.e39*nLep << std::endl;
		double eta = rtsafe<species>(1.e39*nLep, temp, guess);
                return eos_ferm_onthefly(eta, temp, species);
                //return eos_ferm_fromNR<species>(nLep, t);
        } else if constexpr(eos_method == 2) {
                double eta = eos_tintep(nLep, temp, eta_table.nL, eta_table.t, eta_table.eta);
                //if (guess_from_interp == true) {
                //      guess = eta;
                //      eta = rtsafe(1.e39*nLep, T, mLep, guess, eta1_Lep, eta2_Lep);i
                //}
                return eos_ferm_onthefly(eta, temp, species);
        }
}


template<int species>
std::array<double,4> der_cs2(const double nLep, const double temp) {
        std::array<double,4> der_array;

//.......Define theta (relativity parameter)
       	const double theta = temp/mL[species];

	const double guess = find_guess_eta<species>(1.e39*nLep, temp);
	const double eta   = rtsafe<species>(1.e39*nLep, temp, guess);

        const double f12 = compute_res(eta, theta, 0.5);
        const double f32 = compute_res(eta, theta, 1.5);
        const double f52 = compute_res(eta, theta, 2.5);
	const double f12_dn = compute_res_ed(eta, theta, 0.5);
	const double f32_dn = compute_res_ed(eta, theta, 1.5);

	const double f12_dT = (f32_dn - 1.5*f12) / theta;
	const double f32_dT = (f32 - 4.*f12_dT) / (2.*theta);
	const double f52_dT = (f52 - 4.*f32_dT) / (2.*theta);
	const double f52_dn = theta*f32_dT + 2.5*f32;

        const double a_eta = -(eta+2./theta);

        const double a_f12 = compute_res(a_eta, theta, 0.5);
        const double a_f32 = compute_res(a_eta, theta, 1.5);
        const double a_f52 = compute_res(a_eta, theta, 2.5);
	const double a_f12_dn = compute_res_ed(a_eta, theta, 0.5);
	const double a_f32_dn = compute_res_ed(a_eta, theta, 1.5);

	const double a_f12_dT = (a_f32_dn - 1.5*a_f12) / theta;
	const double a_f32_dT = (a_f32 - 4.*a_f12_dT) / (2.*theta);
	const double a_f52_dT = (a_f52 - 4.*a_f32_dT) / (2.*theta);
	const double a_f52_dn = theta*a_f32_dT + 2.5*a_f32;

	const double s = -eta*f12-a_eta*a_f12+(5./3.-eta*theta)*f32+(5./3.-a_eta*theta)*a_f32+4./3.*theta*(f52+a_f52); //erg/K/cm^3
	const double n = f12+a_f12+theta*(f32+a_f32);

	const double dn = f12_dn+a_f12_dn + theta*(f32_dn+a_f32_dn);
	der_array[0] = mL[species]/3. * theta * (2.*(f32_dn-a_f32_dn) + theta*(f52_dn-a_f52_dn)) / dn;
	der_array[1] = (-f12+a_f12-eta*f12_dn+a_eta*a_f12_dn+5./3.*(f32_dn-a_f32_dn) - theta*(f32-a_f32+eta*f32_dn-a_eta*a_f32_dn-4./3.*(f52_dn-a_f52_dn))) / dn;
	//der_array[1] = der_array[1] - s/n;
	der_array[2] = K3[species] * pow(theta,1.5) * (5.*(f32+a_f32) + theta*(3.5*(f52+a_f52)+2.*(f32_dT+a_f32_dT)) + theta*theta*(f52_dT+a_f52_dT));
	der_array[3] = K[species]/mL[species] * pow(theta,0.5) * (-1.5*(eta*f12+a_eta*a_f12) + 2.5*(f32+a_f32)
						+ theta*(-2.5*(eta*f32+a_eta*a_f32) + 10./3.*(f52*a_f52) - (eta*f12_dT+a_eta*a_f12_dT) 
							+ 5./3.*(f32_dT+a_f32_dT))
						+ theta*theta * (-(eta*f32_dT+a_eta*a_f32_dT) + 4./3.*(f52_dT+a_f52_dT)));
	return der_array;
}


template<int species>
std::array<double,4> der_cs2_num(const double nLep, const double temp) {
        std::array<double,4> der_array;

//.......Define theta (relativity parameter)
       	const double theta = temp/mL[species];

	const double g   = find_guess_eta<species>(1.e39*nLep, temp);
	const double eta = rtsafe<species>(1.e39*nLep, temp, g);
	
	const double eps_t = temp*0.02;
	const double eps_n = nLep*0.02;

	const double temp_1 = temp - eps_t;
	const double temp_2 = temp + eps_t;

	const double nLep_1 = nLep - eps_n;
	const double nLep_2 = nLep + eps_n;

	double g1    = find_guess_eta<species>(1.e39*(nLep_1), temp);
	double eta_1 = rtsafe<species>(1.e39*nLep_1, temp, g1);
	
	double g2    = find_guess_eta<species>(1.e39*(nLep_2), temp);
	double eta_2 = rtsafe<species>(1.e39*nLep_2, temp, g2);

	std::array<double,9> tmp_1 = eos_ferm_onthefly(eta_1, temp, species);
	std::array<double,9> tmp_2 = eos_ferm_onthefly(eta_2, temp, species);

        der_array[0] = (tmp_2[2]+tmp_2[3] - (tmp_1[2]+tmp_1[3])) / (nLep_2-nLep_1);
        der_array[1] = (tmp_2[6]+tmp_2[7] - (tmp_1[6]+tmp_1[7])) / (nLep_2-nLep_1);

	g1    = find_guess_eta<species>(1.e39*nLep, temp_1);
	eta_1 = rtsafe<species>(1.e39*nLep, temp_1, g1);
	
	g2    = find_guess_eta<species>(1.e39*nLep, temp_2);
	eta_2 = rtsafe<species>(1.e39*nLep, temp_2, g2);

	tmp_1 = eos_ferm_onthefly(eta, temp_1, species);
	tmp_2 = eos_ferm_onthefly(eta, temp_2, species);

        der_array[2] = (tmp_2[2]+tmp_2[3] - (tmp_1[2]+tmp_1[3])) / (temp_2-temp_1);
        der_array[3] = (tmp_2[6]+tmp_2[7] - (tmp_1[6]+tmp_1[7])) / (temp_2-temp_1);

        return der_array;
}

template <int id_EOS, int species>
std::array<double,9> eos_ferm_single(double nLep, double temp) {
        std::array<double,9> eos_array;
	const double guess = find_guess_eta<species+1>(1.e39*nLep, temp);
	const double eta = rtsafe<species+1>(1.e39*nLep, temp, guess);
        const double theta = temp/mL[species];

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
        eos_array[IL_N] = K[species]*pow(theta,1.5)*(f12+theta*f32);

//.......Number density of anti-leptons
        eos_array[IA_N] = K[species]*pow(theta,1.5)*(a_f12+theta*a_f32);

//.......Pressure of leptons
        eos_array[IL_P] = K3[species]*mL[species]*pow(theta,2.5)*(2.*f32+theta*f52); //MeV/fm^3

//.......Pressure of anti-leptons
        eos_array[IA_P] = K3[species]*mL[species]*pow(theta,2.5)*(2.*a_f32+theta*a_f52); //MeV/fm^3

//.......Internal energy density of leptons
        eos_array[IL_E] = K[species]*mL[species]*pow(theta,1.5)*(f12+2.*theta*f32+theta*theta*f52); //MeV/fm^3

//.......Internal energy density of anti-leptons
        eos_array[IA_E] = K[species]*mL[species]*pow(theta,1.5)*(a_f12+2.*theta*a_f32+theta*theta*a_f52); //MeV/fm^3

//.......Entropy density of leptons
        eos_array[IL_S] = K[species]*pow(theta,1.5)*(-eta*f12+(5./3.-eta*theta)*f32+4./3.*theta*f52); //1/fm^3

//.......Entropy density of anti-leptons
        eos_array[IA_S] = K[species]*pow(theta,1.5)*(-a_eta*a_f12+(5./3.-a_eta*theta)*a_f32+4./3.*theta*a_f52); //1/fm^3

//.......Chemical potential of leptons rnal energy density of anti-leptons
        eos_array[IL_MU] = mL[species] + eta*temp;

        return eos_array;
}


#endif
