#include <stdio.h>
#include <math.h>
#include <iostream> //Input&Output stream model based library
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <time.h>
#include <stdlib.h>
#include <iomanip>

#include "eos_neutrinos.hpp"
#include "eos_assembled.hpp"

using namespace std;

const bool with_mu = true;

void compute_EOS(EOS_assembled *eos, double nb, double T, double *Y) {
	std::cout << "Input quantities: " << std::endl;
	std::cout << "nb = " << nb << " " << "1/fm-3" << std::endl;
	std::cout << "T  = " << T << " " << "MeV" << std::endl;
	std::cout << "Ye = " << Y[0] << " " << "#/baryons" << std::endl;
	std::cout << "Ym = " << Y[1] << " " << "#/baryons" << std::endl;
	std::cout << std::endl;

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

	std::cout << "Output quantities: " << std::endl;
	std::cout << "d            = " << d         << " " << "g/cm3"      << std::endl;	
	std::cout << "yle          = " << yle       << " " << "#/baryon"      << std::endl;	
	std::cout << "ylm          = " << ylm       << " " << "#/baryon"      << std::endl;	
	std::cout << "u            = " << u         << " " << "erg/cm3"      << std::endl;	
	std::cout << "nb           = " << 1.e39*nb  << " " << "1/cm3"      << std::endl;	
	std::cout << "T            = " << T         << " " << "MeV"        << std::endl;	
	std::cout << "Ye           = " << ye        << " " << "#/baryon"   << std::endl;	
	std::cout << "Ym           = " << ym        << " " << "#/baryon"   << std::endl;	
	std::cout << "e            = " << e         << " " << "erg/cm3"    << std::endl;
	std::cout << "P            = " << ptot      << " " << "erg/cm3"    << std::endl;
	//std::cout << "Entropy      = " << eos.Entropy(nb, T, Y)        << std::endl;
	std::cout << "Yh           = " << yh        << " " << "#/baryon"   << std::endl;
	std::cout << "Ya           = " << ya        << " " << "#/baryon"   << std::endl;
	std::cout << "Yp           = " << yp        << " " << "#/baryon"   << std::endl;
	std::cout << "Yn           = " << yn        << " " << "#/baryon"   << std::endl;
	std::cout << "Ynue         = " << ynue      << " " << "#/baryon"   << std::endl;
	std::cout << "Yanue        = " << yanue     << " " << "#/baryon"   << std::endl;
	std::cout << "Ynum         = " << ynum      << " " << "#/baryon"   << std::endl;
	std::cout << "Yanum        = " << yanum     << " " << "#/baryon"   << std::endl;
	std::cout << "Ynux         = " << ynux      << " " << "#/baryon"   << std::endl;
	std::cout << "Znue         = " << znue      << " " << "MeV/baryon" << std::endl;
	std::cout << "Zanue        = " << zanue     << " " << "MeV/baryon" << std::endl;
	std::cout << "Znum         = " << znum      << " " << "MeV/baryon" << std::endl;
	std::cout << "Zanum        = " << zanum     << " " << "MeV/baryon" << std::endl;
	std::cout << "Znux         = " << znux      << " " << "MeV/baryon" << std::endl;
	std::cout << "mup (wrt mn) = " << mu_p - mb << " " << "MeV"        << std::endl;
	std::cout << "mun (wrt mn) = " << mu_n - mb << " " << "MeV"        << std::endl;
	std::cout << "mue (w me)   = " << mu_e      << " " << "MeV"        << std::endl;
	std::cout << "mum (w mm)   = " << mu_m      << " " << "MeV"        << std::endl;

        return;
}

int main() {
	EOS_assembled eos;

	std::cout << std::scientific;

	std::cout << "#####################" << std::endl;
	if (with_mu == true) {
		std::cout << "#   EOS with muons  #" << std::endl;
	} else {
		std::cout << "# EOS without muons #" << std::endl;
	}
	std::cout << "#####################" << std::endl;
	
	std::cout << std::endl;

	eos.ReadBarTableFromFile("/home/leonardo/Desktop/PhD_work/BNS_muons/EOS_module/eos_table/baryons/DD2_bar.h5");
	eos.ReadETableFromFile("/home/leonardo/Desktop/PhD_work/BNS_muons/EOS_module/eos_table/electrons/eos_electrons_primitive_new.txt");
	if (with_mu == true) eos.ReadMTableFromFile("/home/leonardo/Desktop/PhD_work/BNS_muons/EOS_module/eos_table/muons/eos_muons_primitive_new.txt");
	
	double nb = 0.251188644;    // fm-3
        double T  = 19.0546059;     // MeV
        double Y[2] = {0.31,0.01};  // #/baryons
	
	compute_EOS(&eos, nb, T, Y);
	
	return 0;
}
