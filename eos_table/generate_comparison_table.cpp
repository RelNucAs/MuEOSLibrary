#include <iostream>
#include <fstream>
#include <iomanip>
#include <array>

#include "parameters.hpp"
#include "eos_neutrinos.hpp"
#include "eos_assembled.hpp"

using namespace parameters;

/* Boolean variable for including of muons */
const bool with_mu = false;

/* Function computing the EOS output */
std::array<double,21> compute_EOS(EOS_assembled* eos, double nb, double T, double *Y) {

	std::array<double,21> EOS_output;

	const double mb   = eos->GetBaryonMass();
        const double d    = 1.e39*nb*mb*MeV/(c*c);
        const double mu_n = eos->NeutronChemicalPotential(nb, T, Y);
        const double mu_p = eos->ProtonChemicalPotential(nb, T, Y);
        const double mu_e = eos->EOS_leptons<0>::LepChemicalPotential<id_test>(nb, T, Y);
        const double mu_m = 0.; //eos->EOS_leptons<1>LepChemicalPotential<id_test>(nb, T, Y);

        const double mu_nue = mu_p - mu_n + mu_e;
        const double mu_num = 0.;

        const double e = MeV*1.e39*(eos->Energy(nb, T, Y) - mb*nb); // + 1.e39*8.265e18*mb/(c*c)*MeV*nb; //erg/cm3
        const double p = MeV*1.e39*eos->Pressure(nb, T, Y);

        const double ye = Y[0];
        const double ym = Y[1];

        const double ya = eos->AlphaFraction(nb, T, Y);
        const double yh = eos->HeavyFraction(nb, T, Y);
        const double yn = eos->NeutronFraction(nb, T, Y);
        const double yp = eos->ProtonFraction(nb, T, Y);

        const double ynue  = nu_fraction(nb, T, mu_nue/T);
        const double yanue = anu_fraction(nb, T, mu_nue/T);
        const double ynum  = nu_fraction(nb, T, mu_num/T);
        const double yanum = anu_fraction(nb, T, mu_num/T);
        const double ynux  = nu_fraction(nb, T, 0.);

        const double yle = ye + ynue - yanue;
        const double ylm = ym + ynum - yanum;

        const double znue  = nu_energy(nb, T, mu_nue/T);
        const double zanue = anu_energy(nb, T, mu_nue/T);
        const double znum  = nu_energy(nb, T, mu_num/T);
        const double zanum = anu_energy(nb, T, mu_num/T);
        const double znux  = nu_energy(nb, T, 0.);

        const double u = e + 1.e39*nb * MeV * (znue + zanue + znum + zanum + 2. * znux);
        const double ptot = p + 1.e39*nb * MeV * (znue + zanue + znum + zanum + 2. * znux) / 3.;

        EOS_output[0] = d;          // Mass density [g/cm3]
	EOS_output[1] = yle;        // Electron lepton fraction [#/baryon]
	EOS_output[2] = u;          // Internal energy density (including neutrinos) [erg/cm3]
	EOS_output[3] = nb*1.e39;   // Number density [1/cm3]
	EOS_output[4] = T;          // Temperature [MeV]
	EOS_output[5] = ye;         // Electron fraction [#/baryon]
	EOS_output[6] = e;          // Internal energy density (including neutrinos) [erg/cm3]
	EOS_output[7] = ptot;       // Pressure (including neutrinos) [erg/cm3]
        //std::cout << "Entropy      = " << eos.Entropy(nb, T, Y)        << std::endl;
	EOS_output[8] = yh;         // Heavy nuclei fraction [#/baryon]
	EOS_output[9] = ya;         // Alpha particle fraction [#/baryon]
	EOS_output[10] = yp;        // Proton fraction [#/baryon]
	EOS_output[11] = yn;        // Neutron fraction [#/baryon]
	EOS_output[12] = ynue;      // Electron neutrino fraction [#/baryon]
	EOS_output[13] = yanue;     // Electron antineutrino fraction [#/baryon]
	EOS_output[14] = ynux;      // Heavy (anti)neutrino fraction [#/baryon]
	EOS_output[15] = znue;      // Electron neutrino energy per baryon [MeV/baryon]
	EOS_output[16] = zanue;     // Electron antineutrino energy per baryon [MeV/baryon]
	EOS_output[17] = znux;      // Heavy (anti)neutrino energy per baryon [MeV/baryon]
	EOS_output[18] = mu_p - mb; // Proton chemical potential wrt baryon (neutron) mass [MeV]
	EOS_output[19] = mu_n - mb; // Neutron chemical potential wrt baryon (neutron) mass [MeV]
	EOS_output[20] = mu_e;      // Electron chemical potential (with rest mass) [MeV]

        return EOS_output;
}

int main () {
	std::cout << std::scientific; // scientific format for std output
    	
	std::cout << "#####################" << std::endl;
	std::cout << "# EOS without muons #" << std::endl;
        std::cout << "#####################" << std::endl;
        std::cout << std::endl;
        
	/* Global EOS class */
	EOS_assembled eos;

	/* Read EOS tables (electron table is read just to initialize the class but is not used) */
        eos.ReadBarTableFromFile("/home/leonardo/Desktop/PhD_work/BNS_muons/EOS_module/eos_table/baryons/DD2_bar.h5");
        eos.EOS_leptons<0>::ReadLepTableFromFile("/home/leonardo/Desktop/PhD_work/BNS_muons/EOS_module/eos_table/electrons/eos_electrons_primitive_new.txt");
	eos.EOS_leptons<0>::m_lep_active = true;

	/* Define name of output table */
	std::string table_name = "eos_comparison_wo_mu.txt";

	/* Output stream */
	std::ofstream Iout("/home/leonardo/Desktop/PhD_work/BNS_muons/EOS_module/eos_table/global/"+table_name);
	
	/* Define table input paramters
	
	Input parameters:
	- nb: baryon number density [1/fm3]
	-  T: temperature [MeV]
	- ye: electron fraction [#/baryon]

	*/

	/* Table dimensions */
	const int n_nb = 60; // Number density
	const int n_t  = 27; // Temperature
	const int n_ye = 57; // Electron fraction

	/* Arrays of input variables */
	// Number density
	std::vector<double> n_array; // log spaced
	
	const double nbmin = 7.58576778E-07;
	const double nbmax = 9.12011723E+00;

	const double d_nb = pow(nbmax / nbmin, 1./(n_nb-1.));
	//const double log_nbmin = log10(nbmin);
	//const double log_nbmax = log10(nbmax);

	//for (int i=0; i<n_nb; i++) n_array.push_back(log_nbmin + static_cast<double>(i) * (log_nbmax-log_nbmin) / static_cast<double>(n_nb));
	for (int i=0; i<n_nb; i++) n_array.push_back(nbmin * pow(d_nb,static_cast<double>(i)));
	
	// Temperature
	std::vector<double> t_array; // log spaced
	
	const double tmin = 1.00000000E-01;
	const double tmax = 1.31825639E+02;

	const double d_t = pow(tmax / tmin, 1./(n_t-1.));
	//const double log_tmin = log10(tmin);
	//const double log_tmax = log10(tmax);

	//for (int i=0; i<n_t; i++) t_array.push_back(log_tmin + static_cast<double>(i) * (log_tmax-log_tmin) / static_cast<double>(n_t));
	for (int i=0; i<n_t; i++) t_array.push_back(tmin * pow(d_t,static_cast<double>(i)));
	
	// Electron fraction
	std::vector<double> ye_array; // lin spaced
	
	const double yemin = 9.99999978E-03;
	const double yemax = 5.79999983E-01;

	for (int i=0; i<n_ye+1; i++) ye_array.push_back(yemin + static_cast<double>(i) * (yemax-yemin) / static_cast<double>(n_ye));

	Iout << std::scientific << std::setprecision(10); // scientific format for output table


	/* Output stream -> legend */
	Iout << "# d [g/cm3], yle [#/baryon], u [erg/cm3], nb [1/cm3], T [MeV], ye [#/baryon], ";
	Iout <<   "e [erg/cm3], P [erg/cm3], yh [#/baryon], ya [#/baryon], yp [#/baryon], yn [#/baryon], ";
	Iout <<   "ynue [#/baryon], yanue [#/baryon], ynux [#/baryon], ";
	Iout <<   "znue [MeV/baryon], zanue [MeV/baryon], znux [MeV/baryon], ";
	Iout <<   "mup wrt mn [MeV], mun wrt mn [MeV], mue w me [MeV]" << std::endl;

	double nb, T;
	double yq, ye;
	double Y[2];
	
	Y[1] = 0.0; // set muon fraction equal to zero

	/* Choose step size for iteration over input variable arrays */
	const int di = 1, dj = 1, dk = 3; //di: nb step, dj: T step, dk: ye step

	std::array<double,21> eos_array;

	/* Write table WITHOUT MUONS to output file */
	std::cout << "# id_n, id_t, id_ye, nb [1/fm3], T [MeV], Ye [#/baryon]" << std::endl;
	for (int i=0; i<n_nb; i+=di) {
		nb = n_array[i];
		for (int j=0; j<n_t; j+=dj) {
			T = t_array[j];
			for (int k=0; k<n_ye+1; k+=dk) {
				Y[0] = ye_array[k];

				std::cout << i  << "  " << j << "  " << k    << "  ";
				std::cout << nb << "  " << T << "  " << Y[0] << std::endl;

				eos_array = compute_EOS(&eos, nb, T, Y);
				for (int i=0; i<21; i++) Iout << eos_array[i] << "  ";
				Iout << std::endl;
			}
		}
	}

	Iout.close(); // close output stream

	return 0;
}
