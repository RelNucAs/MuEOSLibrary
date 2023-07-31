#include <iomanip>

#include "src/eos_species/eos_assembled.hpp"

int main() {
  /* Boolean variable for including muons */
  const bool with_mu = true;

  /* Name of baryon EOS table */
  std::string BarTableName = "eos_table/baryons/DD2_bar.h5";  // baryon table
  
  /* Initialize global EOS class

  Constructor -> EOS_assembled(const int id_eos, const bool el_bool, const bool mu_bool, std::string BarTableName)

  Inputs:
   - id_EOS: method for EOS computation (1: interpolation, 2: on-the-fly)
   - el_bool: flag for activating electrons
   - mu_bool: flag for activating muons
   - BarTableName: path of baryon EOS table  */
  EOS_assembled eos(2, true, with_mu, BarTableName);

  /* Print header */
  std::cout << std::endl;
  std::cout << "#####################" << std::endl;
  if (with_mu == true) {
    std::cout << "#   EOS with muons  #" << std::endl;
  } else {
    std::cout << "# EOS without muons #" << std::endl;
  }
  std::cout << "#####################" << std::endl;

  std::cout << std::endl;

  
  /* Define input quantities */
  double nb = 1.20226528E-05;  // Baryon number density [fm-3]
  double T  = 6.30957362E+00;	 // Temperature [MeV]
  double Ye = 7.00000003E-02;  // Electron fraction [#/baryons]
  double Ym = 1.00000000E-02;  // Muon fraction [#/baryons] (not needed when "with_mu == false")
  double Y[2] = {Ye, Ym};
  
  /* Print input quantities to standard output */
  std::cout << std::scientific << std::setprecision(8);
  std::cout << "Input quantities: " << std::endl;
  std::cout << "nb = " << nb   << " " << "1/fm-3" << std::endl;
  std::cout << "T  = " << T    << " " << "MeV" << std::endl;
  std::cout << "Ye = " << Y[0] << " " << "#/baryons" << std::endl;
  std::cout << "Ym = " << Y[1] << " " << "#/baryons" << std::endl;
  std::cout << std::endl;

  /* Compute EOS Output */
  FullEOSOutput out = eos.compute_full_EOS(nb, T, Y);

  /*
  FullEOSOutput structure description:

    ! Rest mass contribution is included in the chemical potentials !

   - double rho;                // Mass density [g/cm3]
   - double nb;                 // Baryon number density [1/fm3]
   - double T;                  // Temperature [MeV]
   - double mb;                 // Baryon mass [MeV]
   - double e;                  // Internal energy density [erg/cm3]
   - double P;                  // Pressure [erg/cm3]
   - double s;                  // Entropy per baryon [#/baryon]
   - double chem_pot.mu_n;      // Neutron chemical potential [MeV]
   - double chem_pot.mu_p;      // Proton chemical potential [MeV]
   - double chem_pot.mu_e;      // Electron chemical potential [MeV]
   - double chem_pot.mu_m;      // Muon chemical potential [MeV]
   - double chem_pot.mu_nue;    // Electron neutrino chemical potential [MeV] 
   - double chem_pot.mu_num;    // Muon neutrino chemical potential [MeV]
   - double chem_pot.mu_nux;    // Tau neutrino chemical potential [MeV]
   - double Y_part.yh;          // Heavy nuclei fraction [#/baryon]
   - double Y_part.ya;          // Alpha particle fraction [#/baryon]
   - double Y_part.yn;          // Neutron fraction [#/baryon]
   - double Y_part.yp;          // Proton fraction [#/baryon]
   - double Y_part.ye;          // Electron fraction [#/baryon]
   - double Y_part.ym;          // Muon fraction [#/baryon]
   - double Y_part.ynue;        // Electron neutrino fraction [#/baryon]
   - double Y_part.yanue;       // Electron antineutrino fraction [#/baryon]
   - double Y_part.ynum;        // Muon neutrino fraction [#/baryon]
   - double Y_part.yanum;       // Muon antineutrino fraction [#/baryon]
   - double Y_part.ynux;        // Tau (anti)neutrino fraction [#/baryon]
   - double Y_part.yle;         // Net electronic lepton fraction [#/baryons] (yle = )
   - double Y_part.ylm;         // Net muonic lepton fraction [#/baryon] (ylm = )
   - double nuEOS.z_nue;        // Electron neutrino energy per baryon [MeV/baryon]
   - double nuEOS.z_anue;       // Electron antineutrino energy per baryon [MeV/baryon]
   - double nuEOS.z_num;        // Muon neutrino energy per baryon [MeV/baryon]
   - double nuEOS.z_anum;       // Muon antineutrino energy per baryon [MeV/baryon]
   - double nuEOS.z_nux;        // Tau (anti)neutrino energy per baryon [MeV/baryon]
   - double nuEOS.z_tot;        // Total neutrino energy per baryon [MeV/baryon]
   - double nuEOS.s_nue;        // Electron neutrino entropy per baryon [#/baryon]
   - double nuEOS.s_anue;       // Electron antineutrino entropy per baryon [#/baryon]
   - double nuEOS.s_num;        // Muon neutrino entropy per baryon [#/baryon]
   - double nuEOS.s_anum;       // Muon antineutrino entropy per baryon [#/baryon]
   - double nuEOS.s_nux;        // Tau (anti)neutrino entropy per baryon [#/baryon]
   - double nuEOS.s_tot;        // Total neutrino entropy per baryon [#/baryon]

  */

  /* Add contribution of neutrinos */
  const double etot = out.e + 1.0E+39 * nb * MeV * out.nuEOS.Z_tot;  // Internal energy density [erg/cm3]
  const double ptot = out.P + 1.0E+39 * nb * MeV * out.nuEOS.Z_tot / 3.;  // Pressure [erg/cm3]
  const double stot = out.s + out.nuEOS.s_tot;  // Entropy per baryon [#/baryon]

  /* Print EOS output to standard output */
  std::cout << "Output quantities: " << std::endl;
  std::cout << "d            = " << out.rho                    << " " << "g/cm3"      << std::endl;	
  std::cout << "yle          = " << out.Y_part.yle             << " " << "#/baryon"   << std::endl;	
  std::cout << "ylm          = " << out.Y_part.ylm             << " " << "#/baryon"   << std::endl;	
  std::cout << "u            = " << etot                       << " " << "erg/cm3"    << std::endl;	
  std::cout << "nb           = " << 1.0E+39 * nb               << " " << "1/cm3"      << std::endl;	
  std::cout << "T            = " << out.T                      << " " << "MeV"        << std::endl;	
  std::cout << "Ye           = " << out.Y_part.ye              << " " << "#/baryon"   << std::endl;	
  std::cout << "Ym           = " << out.Y_part.ym              << " " << "#/baryon"   << std::endl;	
  std::cout << "e            = " << out.e                      << " " << "erg/cm3"    << std::endl;
  std::cout << "P            = " << ptot                       << " " << "erg/cm3"    << std::endl;
  std::cout << "s            = " << stot                       << " " << "#/baryon"   << std::endl;
  std::cout << "Yh           = " << out.Y_part.yh              << " " << "#/baryon"   << std::endl;
  std::cout << "Ya           = " << out.Y_part.ya              << " " << "#/baryon"   << std::endl;
  std::cout << "Yp           = " << out.Y_part.yp              << " " << "#/baryon"   << std::endl;
  std::cout << "Yn           = " << out.Y_part.yn              << " " << "#/baryon"   << std::endl;
  std::cout << "Ynue         = " << out.Y_part.ynue            << " " << "#/baryon"   << std::endl;
  std::cout << "Yanue        = " << out.Y_part.yanue           << " " << "#/baryon"   << std::endl;
  std::cout << "Ynum         = " << out.Y_part.ynum            << " " << "#/baryon"   << std::endl;
  std::cout << "Yanum        = " << out.Y_part.yanum           << " " << "#/baryon"   << std::endl;
  std::cout << "Ynux         = " << out.Y_part.ynux            << " " << "#/baryon"   << std::endl;
  std::cout << "Znue         = " << out.nuEOS.Z_nue            << " " << "MeV/baryon" << std::endl;
  std::cout << "Zanue        = " << out.nuEOS.Z_anue           << " " << "MeV/baryon" << std::endl;
  std::cout << "Znum         = " << out.nuEOS.Z_num            << " " << "MeV/baryon" << std::endl;
  std::cout << "Zanum        = " << out.nuEOS.Z_anum           << " " << "MeV/baryon" << std::endl;
  std::cout << "Znux         = " << out.nuEOS.Z_nux            << " " << "MeV/baryon" << std::endl;
  std::cout << "mup (wrt mn) = " << out.chem_pot.mu_p - out.mb << " " << "MeV"        << std::endl;
  std::cout << "mun (wrt mn) = " << out.chem_pot.mu_n - out.mb << " " << "MeV"        << std::endl;
  std::cout << "mue (w me)   = " << out.chem_pot.mu_e          << " " << "MeV"        << std::endl;
  std::cout << "mum (w mm)   = " << out.chem_pot.mu_m          << " " << "MeV"        << std::endl;
	
  return 0;
}
