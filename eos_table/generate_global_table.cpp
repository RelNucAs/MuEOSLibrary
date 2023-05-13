#include <iostream>
#include <fstream>
#include <iomanip>
#include <array>

#include "parameters.hpp"
#include "eos_neutrinos.hpp"
#include "eos_assembled.hpp"

using namespace parameters;

const bool with_mu = true;

std::array<double,28> compute_EOS(EOS_assembled* eos, double nb, double T, double *Y) {

	std::array<double,28> EOS_output;

	const double mb   = eos->GetBaryonMass();
        const double d    = 1.e39*nb*mb*MeV/(c*c);
        const double mu_n = eos->NeutronChemicalPotential(nb, T, Y);
        const double mu_p = eos->ProtonChemicalPotential(nb, T, Y);
        const double mu_e = eos->ElectronChemicalPotential<id_test>(nb, T, Y);
        double mu_m = 0.;
        if (with_mu == true) mu_m = eos->MuonChemicalPotential<id_test>(nb, T, Y);

        const double mu_nue = mu_p - mu_n + mu_e;
        double mu_num = 0.;
        if (with_mu == true) mu_num = mu_p - mu_n + mu_m;

        const double e = MeV*1.e39*(eos->Energy(nb, T, Y) - mb*nb) + 1.e39*8.265e18*mb/(c*c)*MeV*nb; //erg/cm3
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
	EOS_output[2] = ylm;        // Muon lepton fraction [#/baryon]
	EOS_output[3] = u;          // Internal energy density (including neutrinos) [erg/cm3]
	EOS_output[4] = nb*1.e39;   // Number density [1/cm3]
	EOS_output[5] = T;          // Temperature [MeV]
	EOS_output[6] = ye;         // Electron fraction [#/baryon]
	EOS_output[7] = ym;         // Muon fraction [#/baryon]
	EOS_output[8] = e;          // Internal energy density (including neutrinos) [erg/cm3]
	EOS_output[9] = p;          // Pressure (including neutrinos) [erg/cm3]
        //std::cout << "Entropy      = " << eos.Entropy(nb, T, Y)        << std::endl;
	EOS_output[10] = yh;        // Heavy nuclei fraction [#/baryon]
	EOS_output[11] = ya;        // Alpha particle fraction [#/baryon]
	EOS_output[12] = yp;        // Proton fraction [#/baryon]
	EOS_output[13] = yn;        // Neutron fraction [#/baryon]
	EOS_output[14] = ynue;      // Electron neutrino fraction [#/baryon]
	EOS_output[15] = yanue;     // Electron antineutrino fraction [#/baryon]
	EOS_output[16] = ynum;      // Muon neutrino fraction [#/baryon]
	EOS_output[17] = yanum;     // Muon antineutrino fraction [#/baryon]
	EOS_output[18] = ynux;      // Tau (anti)neutrino fraction [#/baryon]
	EOS_output[19] = znue;      // Electron neutrino energy per baryon [MeV/baryon]
	EOS_output[20] = zanue;     // Electron antineutrino energy per baryon [MeV/baryon]
	EOS_output[21] = znum;      // Muon neutrino energy per baryon [MeV/baryon]
	EOS_output[22] = zanum;     // Muon antineutrino energy per baryon [MeV/baryon]
	EOS_output[23] = znux;      // Tau (anti)neutrino energy per baryon [MeV/baryon]
	EOS_output[24] = mu_p - mb; // Proton chemical potential wrt baryon (neutron) mass [MeV]
	EOS_output[25] = mu_n - mb; // Neutron chemical potential wrt baryon (neutron) mass [MeV]
	EOS_output[26] = mu_e;      // Electron chemical potential (with rest mass) [MeV]
	EOS_output[27] = mu_m;      // Muon chemical potential (with rest mass) [MeV]

        return EOS_output;
}

int main () {
	std::cout << std::scientific;
    	
	std::cout << "#####################" << std::endl;
        if (with_mu == true) {
                std::cout << "#   EOS with muons  #" << std::endl;
        } else {
                std::cout << "# EOS without muons #" << std::endl;
        }
        std::cout << "#####################" << std::endl;
        std::cout << std::endl;
        
	EOS_assembled eos;

        eos.ReadBarTableFromFile("/home/leonardo/Desktop/PhD_work/BNS_muons/EOS_module/eos_table/baryons/DD2_bar.h5");
        eos.ReadETableFromFile("/home/leonardo/Desktop/PhD_work/BNS_muons/EOS_module/eos_table/electrons/eos_electrons_primitive_new.txt");
        if (with_mu == true) eos.ReadMTableFromFile("/home/leonardo/Desktop/PhD_work/BNS_muons/EOS_module/eos_table/muons/eos_muons_primitive_new.txt");

	std::string table_name;

	if (with_mu == true) {
		table_name = "eos_global_w_mu.txt";
	} else {
		table_name = "eos_global_wo_mu.txt";
	}

	std::ofstream Iout("/home/leonardo/Desktop/PhD_work/BNS_muons/EOS_module/eos_table/global/"+table_name);
	
	std::vector<double> n_array, t_array, ye_array, ym_array;

	/* Table dimensions */
	const int n_nb = 100; // Number density
	const int n_t  =  10; // Temperature
	const int n_ye =  20; // Electron fraction
	const int n_ym =  10; // Muon fraction

	/* Table limits */
	// Number density [1/fm3]
        const double nmin = 1.00e-12;
        const double nmax = 1.90e+00;
        
	const double log_nmin = log10(nmin);
        const double log_nmax = log10(nmax);

        for (int i=0; i<n_nb; i++) n_array.push_back(log_nmin + static_cast<double>(i) * (log_nmax-log_nmin) / static_cast<double>(n_nb));

	// Temperature [MeV]
	const double tmin = 1.00e-01;
	const double tmax = 1.44e+02;

        const double log_tmin = log10(tmin);
        const double log_tmax = log10(tmax);

        for (int i=0; i<n_t; i++) t_array.push_back(log_tmin + static_cast<double>(i) * (log_tmax-log_tmin) / static_cast<double>(n_t));
	
	// Electron fraction [#/baryon]
	const double yemin = 1.00e-02;
	const double yemax = 6.00e-01;

        for (int i=0; i<n_ye; i++) ye_array.push_back(yemin + static_cast<double>(i) * (yemax-yemin) / static_cast<double>(n_ye));
	
	// Muon fraction [#/baryon]
	const double ymmin = 1.00e-02;
	const double ymmax = 6.00e-01;

        for (int i=0; i<n_ym; i++) ym_array.push_back(ymmin + static_cast<double>(i) * (ymmax-ymmin) / static_cast<double>(n_ym));

	Iout << std::scientific << std::setprecision(10);

	Iout << "# d [g/cm3], yle [#/baryon], ylm [#/baryon], u [erg/cm3], nb [1/cm3], T [MeV], ye [#/baryon], ";
	Iout <<   "e [erg/cm3], P [erg/cm3], yh [#/baryon], ya [#/baryon], yp [#/baryon], yn [#/baryon], ";
	Iout <<   "ynue [#/baryon], yanue [#/baryon], ynum [#/baryon], yanum [#/baryon], ynux [#/baryon], ";
	Iout <<   "znue [MeV/baryon], zanue [MeV/baryon], znum [MeV/baryon], zanum [MeV/baryon], znux [MeV/baryon], ";
	Iout <<   "mup wrt mn [MeV], mun wrt mn [MeV], mue w me [MeV]" << std::endl;

	double nb, T;
	double Y[2];
	Y[1] = 0.0;

	std::array<double,28> eos_array;

	if (with_mu == true) {
		std::cout << "# id_n, id_t, id_ye, id_ym,  nb [1/fm3], T [MeV], Ye [#/baryon], Ym [#/baryon]" << std::endl;
	    	for (int i=0; i<n_nb; i++) {
			nb = pow(10.,n_array[i]);
			for (int j=0; j<n_t; j++) {
				T = pow(10.,t_array[j]);
				for (int k=0; k<n_ye; k++) {
					Y[0] = ye_array[k];
					for (int l=0; l<n_ym; l++) {
						Y[1] = ym_array[l];

						std::cout << i  << "  " << j << "  " << k    << "  " << l    << "  ";
						std::cout << nb << "  " << T << "  " << Y[0] << "  " << Y[1] << std::endl;

						eos_array = compute_EOS(&eos, nb, T, Y);
						for (int i=0; i<28; i++) Iout << eos_array[i] << "  ";
						Iout << std::endl;
					}
				}
			}
		}
	} else {
		std::cout << "# id_n, id_t, id_ye, nb [1/fm3], T [MeV], Ye [#/baryon]" << std::endl;
                for (int i=0; i<n_nb; i++) {
                        nb = pow(10.,n_array[i]);
                        for (int j=0; j<n_t; j++) {
                                T = pow(10.,t_array[j]);
                                for (int k=0; k<n_ye; k++) {
                                        Y[0] = ye_array[k];

                                        std::cout << i  << "  " << j << "  " << k    << "  ";
                                        std::cout << nb << "  " << T << "  " << Y[0] << std::endl;

                                        eos_array = compute_EOS(&eos, nb, T, Y);
                                        for (int i=0; i<28; i++) Iout << eos_array[i] << "  ";
                                        Iout << std::endl;
                                }
                        }
                }
	}

	Iout.close();

	return 0;
}
