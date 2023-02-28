#include <stdio.h>
#include <math.h>
#include <iostream> //Input&Output stream model based library
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <array>
#include <time.h>
#include <stdlib.h>
#include <chrono>
#include <bits/stdc++.h>
#include <limits>
#include <iomanip>

#include "constants.hpp"
#include "parameters.hpp"
#include "interp.hpp"
#include "find_eta.hpp"
#include "eos_fermions.hpp"
#include "eos_assembled.hpp"
#include "eos_leptons.hpp"

using namespace std;
using namespace std::chrono;
using namespace constants;

typedef std::chrono::time_point<std::chrono::high_resolution_clock> chrono_t; //Chrono type.
typedef std::chrono::duration<long int, std::ratio<1, 1000> > chrono_int;

int main (){
	int niter = 10000;

	const int sp = 1; //1:electrons, 2:muons

	chrono_t begin, stop;
	chrono_int duration;
	
	double r1, r2, r3;
        double nb, t;
	double y[1];

        const double yemin = 5.e-3;
        const double yemax = 5.5e-1;

        const double denmin = 5.e+3;
        const double denmax = 2.e+16;

        const double nb_min = denmin/(1.e39*mu);
        const double nb_max = denmax/(1.e39*mu);

	const double ne_min = yemin*nb_min;
	const double ne_max = yemax*nb_max;

        const double t_min = 5.e-2;
        const double t_max = 1.05e+2;

	const double dn = nb_max-nb_min;
	const double dt = t_max-t_min;
	const double dy = yemax-yemin;

	std::vector<double> nb_arr, t_arr, y_arr;
	std::array<double,9> lep_EOS;

	//Read EOS tables
        struct EOSeta eta_table = read_eta_table(sp);
        struct EOScomplete EOS_table = read_total_table(sp);
	
	EOS_assembled eos;
	eos.ReadETableFromFile(abs_path+"eos_table/electrons/eos_electrons_complete_leo.txt");

	//Initialization of random number generator
	printf("Measuring numerical performance of EOS module, number of iterations: %d\n\n", niter);

	//Initialization of random number generator
	srand(time(NULL));

	//Loop over (ne-t) grid
	for (int i=0; i<niter; i++) {
		r1 = rand() / (double) RAND_MAX; //Returns a pseudo-random number between 0 and 1
		r2 = rand() / (double) RAND_MAX;
		r2 = rand() / (double) RAND_MAX;

		nb_arr.push_back(nb_min + r1*dn);
		t_arr.push_back(t_min  + r2*dt);
		y_arr.push_back(yemin  + r3*dy);
	}


	std::cout << std::setprecision(5);

	printf("Measuring performance: first method\n");
	begin = high_resolution_clock::now();
	//First method: completely on-the-fly (both eta and e.g. energy)
	for (int i=0; i<niter; i++) {
		nb = nb_arr[i];
		t  = t_arr[i];
		y[0]  = y_arr[i];
		lep_EOS = eos_ferm_array<1,sp>(nb*y[0], t, eta_table, EOS_table);

	}
	stop = high_resolution_clock::now();
	//for (int i=0; i<9; i++) printf("EOS_array[%d] = %.10e\n", i, lep_EOS[i]);
	duration = duration_cast<milliseconds>(stop - begin);
	cout << "Elapsed time, first method: " << duration.count()*1.e-3 << " s" << endl << endl;




	printf("Measuring performance: second method\n");
	begin = high_resolution_clock::now();
	//Second method: eta interpolation + e.g. energy on-the-fly
	for (int i=0; i<niter; i++) {
		nb = nb_arr[i];
		t  = t_arr[i];
		y[0]  = y_arr[i];
		lep_EOS = eos_ferm_array<2,sp>(nb*y[0], t, eta_table, EOS_table);
	}
	stop = high_resolution_clock::now();
	//for (int i=0; i<9; i++) printf("EOS_array[%d] = %.10e\n", i, lep_EOS[i]);
	duration = duration_cast<milliseconds>(stop - begin);
	cout << "Elapsed time, second method: " << duration.count()*1.e-3 << " s" << endl << endl;




	printf("Measuring performance: third method\n");
	begin = high_resolution_clock::now();
	//Third method: direct interpolation of e.g. energy
	for (int i=0; i<niter; i++) {
		nb = nb_arr[i];
		t  = t_arr[i];
		y[0]  = y_arr[i];
		//lep_EOS = eos_ferm_array<3,sp>(ne, t, eta_table, EOS_table);
		std::array<double,eos.NVARS> lep_EOS = eos.ComputeFullElectronEOS(nb, t, y);
	}
	stop = high_resolution_clock::now();
	//for (int i=0; i<9; i++) printf("EOS_array[%d] = %.10e\n", i, lep_EOS[i]);
	//duration = duration_cast<microseconds>(stop - begin);
 	std::chrono::duration<double> diff = stop - begin;
	cout << "Elapsed time, third method: " << std::setprecision(2) << diff.count() << " s" << endl << endl;

	return 0;
}

