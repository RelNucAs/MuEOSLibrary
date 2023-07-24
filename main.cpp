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

#include "src/eos_species/eos_assembled.hpp"

using namespace std;

int main() {
	const bool with_mu = true;

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
	eos.EOS_leptons<0>::ReadLepTableFromFile("/home/leonardo/Desktop/PhD_work/BNS_muons/EOS_module/eos_table/electrons/eos_electrons_primitive_new.txt");
	eos.EOS_leptons<1>::ReadLepTableFromFile("/home/leonardo/Desktop/PhD_work/BNS_muons/EOS_module/eos_table/muons/eos_muons_primitive_new.txt");
	eos.EOS_leptons<0>::m_lep_active = true;
	if (with_mu == true) eos.EOS_leptons<1>::m_lep_active = true;


	double nb = 1.20226528E+34*1.e-39;    // fm-3
	double T  = 6.30957362E+00;  //1.31825639E+02;     // MeV
	double Y[2] = {7.00000003E-02,0.01}; //{5.79999983E-01,0.00};  // #/baryons

	std::cout << "Input quantities: " << std::endl;
	std::cout << "nb = " << nb << " " << "1/fm-3" << std::endl;
	std::cout << "T  = " << T << " " << "MeV" << std::endl;
	std::cout << "Ye = " << Y[0] << " " << "#/baryons" << std::endl;
	std::cout << "Ym = " << Y[1] << " " << "#/baryons" << std::endl;
	std::cout << std::endl;

	FullEOSOutput out = eos.compute_full_EOS(nb, T, Y);

	const double mb   = eos.GetBaryonMass();
	const double etot = out.e + 1.0E+39 * nb * MeV * out.nuEOS.Z_tot;
	const double ptot = out.P + 1.0E+39 * nb * MeV * out.nuEOS.Z_tot / 3.;
	const double stot = out.s + out.nuEOS.s_tot;

	std::cout << "Output quantities: " << std::endl;
	std::cout << "d            = " << out.rho                << " " << "g/cm3"      << std::endl;	
	std::cout << "yle          = " << out.Y_part.yle         << " " << "#/baryon"   << std::endl;	
	std::cout << "ylm          = " << out.Y_part.ylm         << " " << "#/baryon"   << std::endl;	
	std::cout << "u            = " << etot                   << " " << "erg/cm3"    << std::endl;	
	std::cout << "nb           = " << 1.0E+39 * nb           << " " << "1/cm3"      << std::endl;	
	std::cout << "T            = " << out.T                  << " " << "MeV"        << std::endl;	
	std::cout << "Ye           = " << out.Y_part.ye          << " " << "#/baryon"   << std::endl;	
	std::cout << "Ym           = " << out.Y_part.ym          << " " << "#/baryon"   << std::endl;	
	std::cout << "e            = " << out.e                  << " " << "erg/cm3"    << std::endl;
	std::cout << "P            = " << ptot                   << " " << "erg/cm3"    << std::endl;
	std::cout << "s            = " << stot                   << " " << "#/baryon"   << std::endl;
	std::cout << "Yh           = " << out.Y_part.yh          << " " << "#/baryon"   << std::endl;
	std::cout << "Ya           = " << out.Y_part.ya          << " " << "#/baryon"   << std::endl;
	std::cout << "Yp           = " << out.Y_part.yp          << " " << "#/baryon"   << std::endl;
	std::cout << "Yn           = " << out.Y_part.yn          << " " << "#/baryon"   << std::endl;
	std::cout << "Ynue         = " << out.Y_part.ynue        << " " << "#/baryon"   << std::endl;
	std::cout << "Yanue        = " << out.Y_part.yanue       << " " << "#/baryon"   << std::endl;
	std::cout << "Ynum         = " << out.Y_part.ynum        << " " << "#/baryon"   << std::endl;
	std::cout << "Yanum        = " << out.Y_part.yanum       << " " << "#/baryon"   << std::endl;
	std::cout << "Ynux         = " << out.Y_part.ynux        << " " << "#/baryon"   << std::endl;
	std::cout << "Znue         = " << out.nuEOS.Z_nue        << " " << "MeV/baryon" << std::endl;
	std::cout << "Zanue        = " << out.nuEOS.Z_anue       << " " << "MeV/baryon" << std::endl;
	std::cout << "Znum         = " << out.nuEOS.Z_num        << " " << "MeV/baryon" << std::endl;
	std::cout << "Zanum        = " << out.nuEOS.Z_anum       << " " << "MeV/baryon" << std::endl;
	std::cout << "Znux         = " << out.nuEOS.Z_nux        << " " << "MeV/baryon" << std::endl;
	std::cout << "mup (wrt mn) = " << out.chem_pot.mu_p - mb << " " << "MeV"        << std::endl;
	std::cout << "mun (wrt mn) = " << out.chem_pot.mu_n - mb << " " << "MeV"        << std::endl;
	std::cout << "mue (w me)   = " << out.chem_pot.mu_e      << " " << "MeV"        << std::endl;
	std::cout << "mum (w mm)   = " << out.chem_pot.mu_m      << " " << "MeV"        << std::endl;
	
	return 0;
}
