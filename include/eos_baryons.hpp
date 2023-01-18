#pragma once

const int t_entries = 81; //number of temperature entries
const int y_entries = 60; //number of proton fraction entries
const int nb_entries = 326; //303 number of nb entries
const int entries = 19; //number of EOS entries per grid-point

typedef std::vector<std::vector<std::vector<std::vector<double>>>> eos_type;
        
//......define the constants of the t, y, and nb array
const double nbin = 1.0e-12;
const double tin = 0.1e0;
const double yin = 0.0e0;
const double dnb = 0.04e0;
const double dt = 0.04e0;
const double dy = 0.01e0;

std::array<double,t_entries> read_t_baryon();

std::array<double,y_entries> read_y_baryon();

std::array<double,nb_entries> read_nb_baryon();

//This function reads in the ASCII EOS table and calculates the nominal proton fraction and nominal baryon number density arrays
//Input: eos_filename.tab
//std::array<std::array<std::array<std::array<double,entries>,nb_entries>,y_entries>,t_entries> eos_readin(std::string eos_filename) {
eos_type eos_readin();


//.....This function interpolates the table in the temperature
std::array<double,entries> eos_tinterp_bar(double inb, double iy, double it, std::array<double,t_entries> t, std::array<double,y_entries> y, std::array<double,nb_entries> nb, eos_type eos);
