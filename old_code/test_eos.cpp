#include <stdio.h>
#include <math.h>
#include <iostream> //Input&Output stream model based library
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <array>

#include "constants.hpp"
#include "parameters.hpp"
#include "interp.hpp"
#include "eos_fermions.hpp"
#include "find_eta.hpp"
#include "eos_assembled.hpp"

using namespace constants;
using namespace std;

void eos_cout(const double theta, const std::array<double,9> &eos_array) {
	cout << "Particle: " << endl;
	cout << "eta = " << eos_array[8] << endl;
	cout << "n = "   << eos_array[0] << " cm-3"         << endl;
	cout << "P = "   << eos_array[2] << " erg cm-3"     << endl;
	cout << "e = "   << eos_array[4] << " erg cm-3"     << endl;
	cout << "s = "   << eos_array[6] << " erg cm-3 K-1" << endl;
	cout << "Antiparticle: " << endl;
	cout << "eta = " << -(eos_array[8]+2./theta) << endl;
	cout << "n = "   << eos_array[1] << " cm-3"         << endl;
	cout << "P = "   << eos_array[3] << " erg cm-3"     << endl;
	cout << "e = "   << eos_array[5] << " erg cm-3"     << endl;
	cout << "s = "   << eos_array[7] << " erg cm-3 K-1" << endl << endl;
	return;
}

void compare_cout(const std::array<double,9> &eos_1, const std::array<double,9> &eos_2) {
	std::array<double,9> eos_diff;
	for (int i=0;i<9;i++) eos_diff[i] = fabs(eos_1[i]-eos_2[i]) / fabs(eos_1[i]);

	cout << "d_eta = " <<  eos_diff[8] << endl;
	cout << "Particle: ";
	cout << "d_n = " <<  eos_diff[0] << ", ";
	cout << "d_P = " <<  eos_diff[2] << ", ";
	cout << "d_e = " <<  eos_diff[4] << ", ";
	cout << "d_s = " <<  eos_diff[6] << endl;
	cout << "Antiparticle: ";
	cout << "d_n = " <<  eos_diff[1] << ", ";
	cout << "d_P = " <<  eos_diff[3] << ", ";
	cout << "d_e = " <<  eos_diff[5] << ", ";
	cout << "d_s = " <<  eos_diff[7] << endl << endl;

	return;
}

template<int species>
void test_lep_eos(const double nLep, const double temp) {
	/*Read EOS tables*/
        struct EOSeta eta_table = read_eta_table(species);
        struct EOScomplete EOS_table = read_total_table(species);
       
        double mLep;

        //EOS arrays
	std::array<double,9> eos_1, eos_2, eos_3;

	if constexpr(species == 1) {
		mLep = MEOS_me;
	} else if constexpr(species == 2) {
		mLep = mmu;
	}

	const double theta = temp/mLep;

	//First method: completely on-the-fly (both eta and e.g. energy)
	cout << "First method: EOS calculation completely on-the-fly" << endl;
	eos_1 = eos_ferm_array<1,species>(nLep, temp, eta_table, EOS_table);
	eos_cout(theta, eos_1);

	//Second method: eta interpolation + e.g. energy on-the-fly
	cout << "Second method: eta interpolation, then calculation on-the-fly" << endl;
	eos_2 = eos_ferm_array<2,species>(nLep, temp, eta_table, EOS_table);
	eos_cout(theta, eos_2);

	//Third method: direct interpolation of e.g. energy
	cout << "Third method: direct interpolation" << endl;
	eos_3 = eos_ferm_array<3,species>(nLep, temp, eta_table, EOS_table);
	eos_cout(theta, eos_3);

	cout << "Comparison between first and second method...." << endl;
	compare_cout(eos_1, eos_2);

	cout << "Comparison between first and third  method...." << endl;
	compare_cout(eos_1, eos_3);
	return;
}


int main (){
	double n_e, n_mu, temp;
	//EOS_leptons eos;	
	//eos.ReadETableFromFile("/home/leonardo/Desktop/PhD_work/BNS_muons/EOS_module/eos_table/electrons/eos_electrons_complete_leo.txt");
	//const double * test = eos.GetRawELogNumberDensity();
	//for (int i=0; i<eos.m_ne; i++) {
	//	cout << pow(10.,test[i]) << endl;
	
	n_e  = 1.54223401E-14; //fm-3
	temp = 5.12913428E-02; //MeV
	//double Ye[1] = {1.};	

        //double guess = find_guess_eta<1>(1.e39*n_e, temp);
        //double eta = rtsafe<1>(1.e39*n_e, temp, guess);
	//std::array<double,9> tmp = eos_ferm_onthefly(eta, temp, 0);
	//for (int i=0;i<9;i++) cout << tmp[i]*1.e39*MeV << endl;
        	
                //return eos_ferm_fromNR<species>(nLep, t);



	//cout << eos.ElectronPressure(n_e, temp, Ye) * 1.e39 / MeV << endl;

	//exit(EXIT_SUCCESS);
	//cout << "Testing accuracy of electron/positron EOS in a single (n_e,T) point:" << endl;
	//cout << "n_e = " << n_e << " fm-3, T = " << temp << " MeV (theta = " << temp/me << ")" << endl << endl;
	//test_lep_eos<1>(n_e, temp); //1: electrons


	//n_mu = 1.e-05; //fm-3
	//temp = 1.e+00; //MeV

	//cout << "Testing accuracy of (anti)muon EOS in a single (n_mu,T) point:" << endl;
	//cout << "n_mu = " << n_mu << " fm-3, T = " << temp << " MeV (theta = " << temp/mmu << ")" << endl << endl;
	//test_lep_eos<2>(n_mu, temp); //2: muons
	
	std::array<double,4> tmp1 = der_cs2<1>(n_e, temp);
	std::array<double,4> tmp2 = der_cs2_num<1>(n_e, temp);
	for (int i=0;i<4;i++) cout << "der_" << i << " = " << tmp1[i] << endl;
	for (int i=0;i<4;i++) cout << "der_" << i << " = " << tmp2[i] << endl;
	return 0;
}


//double h1 = compute_res(eta, theta, 0.5);
//double h2 = compute_res(eta, theta, 1.5);
//double h3 = compute_res(eta, theta, 2.5);
//double h2 = compute_res_ed(eta, theta, k);

//double test_1 = exp(500.);
//long double test_2 = exp((long double) 500.);
//long double test_3 = expl((long double) 500.);

//long double tmp = -25.;
//long double tt = 1. - exp(tmp)/(1.+exp(tmp));
//long double tt1 = pow(1.+exp(tmp),-1.);
//long double tt2 = 1. - exp(tmp);


