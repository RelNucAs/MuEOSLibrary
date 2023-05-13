#ifndef PARAMS_H
#define PARAMS_H

//Define parameters

#include <string>

struct eos_table_dim {int nne; int nnm; int nt;};
struct eos_table_limits {const double n_min; const double n_max; const double y_min; const double y_max; const double t_min; const double t_max; const double eta_min; const double eta_max;};

namespace parameters
{
	//numerical integration
	const double alpha = 0.; //power exponent of Gauss-Laguerre weight
	const int ngle = 20;
	const int ngla = 20;

	const double amu_g = 1.66054e-24;

	//fermionic eos
	//electrons
	//const struct eos_table_limits e_tab_lim  = {5.0e+03/amu_g*1.e-39, 2.0e+16/amu_g*1.e-39, 5.0e-03, 5.5e-01, 5.0e-02, 1.05e+02, -2.3e+04, 2.3e+04}; 
	const struct eos_table_limits e_tab_lim  = {1.0000000e-12, 1.9054607e+00, 1.00e-02, 6.00e-01, 1.0000000e-01, 1.4454398e+02, -2.3e+04, 2.3e+04}; 
	//muons
	const struct eos_table_limits mu_tab_lim = {1.00e-12, 1.90e+00, 5.e-08, 5.e-01, 1.00e-01, 1.58e+02, -2.3e+04, 5.0e+04}; //1.3e+04
					   //n_min, n_max, y_min,   y_max,  t_min,    t_max,  eta_min,  eta_max


	constexpr int id_test = 1;

	const bool HR = false;

	const struct eos_table_dim SR_tab = {700 ,  750, 150};
	const struct eos_table_dim HR_tab = {1400, 1500, 300};
	
	const std::string abs_path = "/home/leonardo/Desktop/PhD_work/BNS_muons/EOS_module/"; //to be moved somewhere else
}

#endif
