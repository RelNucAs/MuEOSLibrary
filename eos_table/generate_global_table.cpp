#include <iostream>
#include <fstream>
#include <iomanip>
#include <array>

#include "parameters.hpp"
#include "eos_neutrinos.hpp"
#include "eos_assembled.hpp"

using namespace parameters;

/* Boolean variable for including of muons */
const bool with_mu = false; // true

/* Function computing the EOS output */
std::array<double,29> compute_EOS(EOS_assembled* eos, double nb, double T, double *Y) {

	std::array<double,29> EOS_output;

	const double mb   = eos->GetBaryonMass();
        const double d    = 1.e39*nb*mb*MeV/(c*c);
        const double mu_n = eos->NeutronChemicalPotential(nb, T, Y);
        const double mu_p = eos->ProtonChemicalPotential(nb, T, Y);
        const double mu_e = eos->EOS_leptons<0>::LepChemicalPotential<id_test>(nb, T, Y);
        const double mu_m = eos->EOS_leptons<1>::LepChemicalPotential<id_test>(nb, T, Y);

        const double mu_nue = mu_p - mu_n + mu_e;
        double mu_num = 0.;
        if (with_mu == true) mu_num = mu_p - mu_n + mu_m;

	const double e = MeV*1.e39*(eos->Energy(nb, T, Y) - mb*nb); // + 1.e39*8.265e18*mb/(c*c)*MeV*nb; //erg/cm3
        const double p = MeV*1.e39*eos->Pressure(nb, T, Y);
	const double s = eos->Entropy(nb, T, Y); // kB=1

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

        const double snue  = nu_entropy(nb, T, mu_nue/T);
        const double sanue = anu_entropy(nb, T, mu_nue/T);
        const double snum  = nu_entropy(nb, T, mu_num/T);
        const double sanum = anu_entropy(nb, T, mu_num/T);
        const double snux  = nu_energy(nb, T, 0.);
        
	const double u = e + 1.e39*nb * MeV * (znue + zanue + znum + zanum + 2. * znux);
        const double ptot = p + 1.e39*nb * MeV * (znue + zanue + znum + zanum + 2. * znux) / 3.;
	const double stot = s + (snue + sanue + snum + sanum + 2. * snux);

        EOS_output[0] = d;          // Mass density [g/cm3]
	EOS_output[1] = yle;        // Electron lepton fraction [#/baryon]
	EOS_output[2] = ylm;        // Muon lepton fraction [#/baryon]
	EOS_output[3] = u;          // Internal energy density (including neutrinos) [erg/cm3]
	EOS_output[4] = nb*1.e39;   // Number density [1/cm3]
	EOS_output[5] = T;          // Temperature [MeV]
	EOS_output[6] = ye;         // Electron fraction [#/baryon]
	EOS_output[7] = ym;         // Muon fraction [#/baryon]
	EOS_output[8] = e;          // Internal energy density (including neutrinos) [erg/cm3]
	EOS_output[9] = ptot;       // Pressure (including neutrinos) [erg/cm3]
        EOS_output[10] = stot;      // Entropy per baryon (including neutrinos) [#/baryon]
	EOS_output[11] = yh;        // Heavy nuclei fraction [#/baryon]
	EOS_output[12] = ya;        // Alpha particle fraction [#/baryon]
	EOS_output[13] = yp;        // Proton fraction [#/baryon]
	EOS_output[14] = yn;        // Neutron fraction [#/baryon]
	EOS_output[15] = ynue;      // Electron neutrino fraction [#/baryon]
	EOS_output[16] = yanue;     // Electron antineutrino fraction [#/baryon]
	EOS_output[17] = ynum;      // Muon neutrino fraction [#/baryon]
	EOS_output[18] = yanum;     // Muon antineutrino fraction [#/baryon]
	EOS_output[19] = ynux;      // Tau (anti)neutrino fraction [#/baryon]
	EOS_output[20] = znue;      // Electron neutrino energy per baryon [MeV/baryon]
	EOS_output[21] = zanue;     // Electron antineutrino energy per baryon [MeV/baryon]
	EOS_output[22] = znum;      // Muon neutrino energy per baryon [MeV/baryon]
	EOS_output[23] = zanum;     // Muon antineutrino energy per baryon [MeV/baryon]
	EOS_output[24] = znux;      // Tau (anti)neutrino energy per baryon [MeV/baryon]
	EOS_output[25] = mu_p - mb; // Proton chemical potential wrt baryon (neutron) mass [MeV]
	EOS_output[26] = mu_n - mb; // Neutron chemical potential wrt baryon (neutron) mass [MeV]
	EOS_output[27] = mu_e;      // Electron chemical potential (with rest mass) [MeV]
	EOS_output[28] = mu_m;      // Muon chemical potential (with rest mass) [MeV]

        return EOS_output;
}

int main () {
	std::cout << std::scientific; // scientific format for std output
    	
	/* Print to std output if EOS contains muons contribution or not */
	std::cout << "#####################" << std::endl;
        if (with_mu == true) {
                std::cout << "#   EOS with muons  #" << std::endl;
        } else {
                std::cout << "# EOS without muons #" << std::endl;
        }
        std::cout << "#####################" << std::endl;
        std::cout << std::endl;

	/* Global EOS class */
	EOS_assembled eos;

	/* Read EOS tables (electron and muon tables are read just to initialize the class but they are not used) */
        eos.ReadBarTableFromFile("eos_table/baryons/DD2_bar.h5");
        eos.EOS_leptons<0>::m_lep_active = true;
	if (with_mu == true) eos.EOS_leptons<1>::m_lep_active = true;


	/* Define name of output table */
	std::string table_name;
	if (with_mu == true) {
		table_name = "eos_global_w_mu.txt";
	} else {
		table_name = "eos_global_wo_mu.txt";
	}

	/* Output stream */
	std::ofstream Iout("eos_table/global/"+table_name);
	
	/* Define table input paramters
	
	Input parameters:
	- nb: baryon number density [1/fm3]
	-  T: temperature [MeV]
	- yq: charge fraction [#/baryon]
	- ym: muon fraction [#/baryon]

	where ye = yq - ym (charge neutrality), we discard negative values of the electron fraction ye

	*/

	/* Table dimensions */
	const int n_nb = eos.m_nn; // Number density
	const int n_t  = eos.m_nt; // Temperature
	const int n_yq = eos.m_ny; // Charge fraction
	const int n_ym = 60;       // Muon fraction

	/* Arrays of input variables */
	// Number density
	double const * n_array  = eos.GetRawLogNumberDensity(); // log spaced
	
	// Temperature
	double const * t_array  = eos.GetRawLogTemperature(); // log spaced
    	
	// Charge fraction
	double const * yq_array = eos.GetRawYq(); // lin spaced
	
	// Muon fraction
	std::vector<double> ym_array; // log spaced
	
	const double ymmin = 2.0e-05;
	const double ymmax = 2.0e-01;

	const double log_ymmin = log10(ymmin);
	const double log_ymmax = log10(ymmax);
        
	for (int i=0; i<n_ym; i++) ym_array.push_back(log_ymmin + static_cast<double>(i) * (log_ymmax-log_ymmin) / static_cast<double>(n_ym));


	Iout << std::scientific << std::setprecision(10); // scientific format for output table

	/* Output stream -> legend */
	Iout << "# d [g/cm3], yle [#/baryon], ylm [#/baryon], u [erg/cm3], nb [1/cm3], T [MeV], ye [#/baryon], ";
	Iout <<   "e [erg/cm3], P [erg/cm3], s [#/baryon], yh [#/baryon], ya [#/baryon], yp [#/baryon], yn [#/baryon], ";
	Iout <<   "ynue [#/baryon], yanue [#/baryon], ynum [#/baryon], yanum [#/baryon], ynux [#/baryon], ";
	Iout <<   "znue [MeV/baryon], zanue [MeV/baryon], znum [MeV/baryon], zanum [MeV/baryon], znux [MeV/baryon], ";
	Iout <<   "mup wrt mn [MeV], mun wrt mn [MeV], mue w me [MeV]" << std::endl;

	double nb, T;
	double yq, ye, ym;
	double Y[2] = {0.0};
	
	//Y[1] = 0.0; // set muon fraction equal to zero in case muons are not included

	/* Choose step size for iteration over input variable arrays */
	const int di = 10, dj = 10, dk = 10, dl = 10; //di: nb step, dj: T step, dk: yq step, dl: ym step

	std::array<double,29> eos_array;

	if (with_mu == true) {
		/* Write table WITH MUONS to output file */
		std::cout << "# id_n, id_t, id_ye, id_ym,  nb [1/fm3], T [MeV], Yq [#/baryon], Ym [#/baryon]" << std::endl;
	    	for (int i=0; i<n_nb; i+=di) {
			nb = exp(n_array[i]);
			for (int j=0; j<n_t; j+=dj) {
				T = exp(t_array[j]);
				for (int k=0; k<n_yq; k+=dk) {
					yq = yq_array[k];
					for (int l=0; l<n_ym; l+=dl) {
						ym = pow(10.,ym_array[l]);
						ye = yq - ym;

						if (ye <= 0.) continue; // exclude EOS point if Ye is negative

						Y[0] = ye;
						Y[1] = ym;

						std::cout << i  << "  " << j << "  " << k    << "  " << l    << "  ";
						std::cout << nb << "  " << T << "  " << Y[0] << "  " << Y[1] << std::endl;

						eos_array = compute_EOS(&eos, nb, T, Y);
						for (int i=0; i<29; i++) Iout << eos_array[i] << "  ";
						Iout << std::endl;
					}
				}
			}
		}
	} else {
		/* Write table WITHOUT MUONS to output file */
		std::cout << "# id_n, id_t, id_ye, nb [1/fm3], T [MeV], Yq [#/baryon]" << std::endl;
                for (int i=0; i<n_nb; i+=di) {
                        nb = exp(n_array[i]);
                        for (int j=0; j<n_t; j+=dj) {
                                T = exp(t_array[j]);
                                for (int k=0; k<n_yq; k+=dk) {
                                        Y[0] = yq_array[k];

                                        std::cout << i  << "  " << j << "  " << k    << "  ";
                                        std::cout << nb << "  " << T << "  " << Y[0] << std::endl;

                                        eos_array = compute_EOS(&eos, nb, T, Y);
                                        for (int i=0; i<29; i++) Iout << eos_array[i] << "  ";
                                        Iout << std::endl;
                                }
                        }
                }
	}

	Iout.close(); // close output stream

	return 0;
}
