#ifndef INTERP_H
#define INTERP_H

#include <vector>

struct EOSeta {
        //double d[nne]; //density
        //double t[nt]; //temperature
        //double eta[nne*nt]; //relativistic parameter
	int nn; //number of density points
	int nt; //number of temperature points
	std::vector<double> nL; //density
	std::vector<double> t; //temperature
	std::vector<double> eta; //relativistic parameter
};

struct EOScomplete_old {
	int nn; //number of density points
	int nt; //number of temperature points
	std::vector<double> d; //density
	std::vector<double> t; //temperature
	std::vector<double> mu_part; //chemical potential of particles
	std::vector<double> n_part; //number density of particles 
	std::vector<double> n_antipart; //number density of antiparticles 
	std::vector<double> P_part; //pressure of particles 
	std::vector<double> P_antipart; //pressure of antiparticles 
	std::vector<double> e_part; //internal energy of particles 
	std::vector<double> e_antipart; //internal energy of antiparticles
	std::vector<double> s_part; //entropy of particles
	std::vector<double> s_antipart; //entropy of antiparticles
};

struct EOScomplete {
	int nn; //number of density points
	int nt; //number of temperature points
	std::vector<double> nL; //density
	std::vector<double> t; //temperature
	std::array<std::vector<double>,13> eos_table;
	//i=0: chemical potential of particles
	//i=1: number density of particles 
	//i=2: number density of antiparticles 
	//i=3: pressure of particles 
	//i=4: pressure of antiparticles 
	//i=5: internal energy of particles 
	//i=6: internal energy of antiparticles
	//i=7: entropy of particles
	//i=8: entropy of antiparticles
};

void reset_string(std::stringstream ss);

std::string trim(const std::string& str,
                 const std::string& whitespace); // = " \t");

struct EOSeta read_eta_table(const int sp);

struct EOScomplete read_total_table(const int sp);

void set_idx(const double x, const std::vector<double> &x_arr, int &idx, double &r);

void setd(double d, std::vector<double> d_arr, int* idx, double* r);

double eos_tintep(double d, double t, std::vector<double> &d_arr, std::vector<double> &t_arr, std::vector<double> &y);

std::array<double,13> eos_interp(const double d, const double t, const struct EOScomplete &EOSin);

#endif
