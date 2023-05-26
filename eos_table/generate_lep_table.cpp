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

#include "parameters.hpp"
#include "find_eta.hpp"
#include "eos_fermions.hpp"
#include "interp.hpp"

using namespace std;

constexpr int id_test = 2;

template<int species>
void make_lep_table() {
        int n1, n2;
        double  nbmin,  nbmax;
        double   ymin,   ymax;
        double   tmin,   tmax;
        double etamin, etamax;
        string sp, res;
        if constexpr(species == 0) {
                nbmin = e_tab_lim.n_min;
                nbmax = e_tab_lim.n_max;
                ymin = e_tab_lim.y_min;
                ymax = e_tab_lim.y_max;
                tmin = e_tab_lim.t_min;
                tmax = e_tab_lim.t_max;
                etamin = e_tab_lim.eta_min;
                etamax = e_tab_lim.eta_max;
		sp = "electrons";
                if (parameters::HR == true){
                        n1 = HR_tab.nne;
                        n2 = HR_tab.nt;
                        res = "_HR";
                } else {
                        n1 = SR_tab.nne;
                        n2 = SR_tab.nt;
                        res = "";
                }

        } else if constexpr(species == 1) {
                nbmin = mu_tab_lim.n_min;
                nbmax = mu_tab_lim.n_max;
                ymin = mu_tab_lim.y_min;
                ymax = mu_tab_lim.y_max;
                tmin = mu_tab_lim.t_min;
                tmax = mu_tab_lim.t_max;
                etamin = mu_tab_lim.eta_min;
                etamax = mu_tab_lim.eta_max;
		sp = "muons";
                if (parameters::HR == true){
                        n1 = HR_tab.nnm;
                        n2 = HR_tab.nt;
                        res = "_HR";
                } else {
                        n1 = SR_tab.nnm;
                        n2 = SR_tab.nt;
                        res = "";
                }

        } else {
                cout << "Make_lep_table: Error in species identifier" << endl;
                exit(EXIT_FAILURE);
        }

        const double nmin = ymin*nbmin; ///(1.e39*mu);
        const double nmax = ymax*nbmax; ///(1.e39*mu);

        const double log_nmin = log10(nmin);
        const double log_nmax = log10(nmax);

        const double log_tmin = log10(tmin);
        const double log_tmax = log10(tmax);

        double nLep, temp, eta;
        double guess;

        std::vector<double>   lne_edges, n_array;
        std::vector<double>    lt_edges, t_array;

	std::vector<std::array<double,9>> eos_out;
        std::vector<std::array<double,4>> der_out;

	string eta_filename   = abs_path + "eos_table/" + sp + "/eos_" + sp + "_eta"+res+".txt";
        string table_filename = abs_path + "eos_table/" + sp + "/eos_" + sp + "_primitive_new_cs2_num_fine_log"+res+".txt";

        cout << "Generating complete table for " << sp << " with :" << endl;
        cout << "# n points = " << n1 << ", # t points = " << n2 << " (HR = " << HR << ")" << endl;
        cout << "n_min = "   <<   nmin << " fm-3, n_max = " <<   nmax << " fm-3" << endl;
        cout << "t_min = "   <<   tmin <<  " MeV, t_max = " <<   tmax <<  " MeV" << endl;
        cout << "eta_min = " << etamin <<    ", eta_max = " << etamax <<            endl;

	ofstream Ieta(eta_filename);
        Ieta << n1 << "\n";
        Ieta << n2 << "\n";
        Ieta << std::scientific << setprecision(8);

	ofstream Iout(table_filename);
        Iout << n1 << "\n";
	Iout << n2 << "\n";
	Iout << std::scientific << setprecision(8);

        //for (int i=0;i<n1+1;i++) lne_edges.push_back(log_nmin + static_cast<double>(i) * (log_nmax-log_nmin) / static_cast<double>(n1));
        //for (int i=0;i<n2+1;i++)  lt_edges.push_back(log_tmin + static_cast<double>(i) * (log_tmax-log_tmin) / static_cast<double>(n2));
        
	//for (int i=0;i<n1;i++) n_array.push_back(0.5 * (lne_edges[i+1]+lne_edges[i]));
        //for (int i=0;i<n2;i++) t_array.push_back(0.5 * ( lt_edges[i+1]+ lt_edges[i]));

        for (int i=0;i<n1;i++) n_array.push_back(log_nmin + static_cast<double>(i) * (log_nmax-log_nmin) / static_cast<double>(n1-1));
        for (int i=0;i<n2;i++) t_array.push_back(log_tmin + static_cast<double>(i) * (log_tmax-log_tmin) / static_cast<double>(n2-1));
	
	cout << std::scientific;

        for (int i=0;i<n1;i++) {
                Ieta << n_array[i] << " ";
		Iout << n_array[i] << " ";
		//cout << pow(10.,n_array[i]) << endl;
        }
        Ieta << "\n";
        Iout << "\n";
	
	cout << endl << endl << endl;
        
	for (int i=0;i<n2;i++) {
                Ieta << t_array[i] << " ";
                Iout << t_array[i] << " ";
		//cout << pow(10.,t_array[i]) << endl;
        }
        Ieta << "\n";
        Iout << "\n";

        cout << "Loop over number density array (from 0 to " << (n1-1) << ")" << endl;
	for (int i=0;i<n1;i++) {
                cout << "i = " << i << endl;
		for (int j=0;j<n2;j++) {
	                nLep = pow(10.,n_array[i]);
			temp = pow(10.,t_array[j]);

			guess = find_guess_eta<species>(1.e39*nLep, temp);
                        eta = rtsafe<species>(1.e39*nLep, temp, guess);
			Ieta << eta << " ";

			// Check species
			eos_out.push_back(eos_ferm_onthefly(eta, temp, species));
                        //der_out.push_back(der_cs2<species>(nLep, temp));
                        der_out.push_back(der_cs2_num<species>(nLep, temp));
		}
		Ieta << "\n";
	}

	Ieta.close();

        // Lepton number density
	for (int i=0;i<n1;i++) {
                for (int j=0;j<n2;j++) {
                        Iout << eos_out[n2*i+j][IL_N] << " ";
                }
                Iout << "\n";
        }

        // Anti-Lepton number density
	for (int i=0;i<n1;i++) {
                for (int j=0;j<n2;j++) {
                        Iout << eos_out[n2*i+j][IA_N] << " ";
                }
                Iout << "\n";
        }

	// Lepton pressure
	for (int i=0;i<n1;i++) {
                for (int j=0;j<n2;j++) {
                        Iout << eos_out[n2*i+j][IL_P] << " ";
                }
                Iout << "\n";
        }

        // Anti-Lepton pressure
	for (int i=0;i<n1;i++) {
                for (int j=0;j<n2;j++) {
                        Iout << eos_out[n2*i+j][IA_P] << " ";
                }
                Iout << "\n";
        }

        // Lepton internal energy density
	for (int i=0;i<n1;i++) {
                for (int j=0;j<n2;j++) {
                        Iout << eos_out[n2*i+j][IL_E] << " ";
                }
                Iout << "\n";
        }

        // Anti-Lepton internal energy density
	for (int i=0;i<n1;i++) {
                for (int j=0;j<n2;j++) {
                        Iout << eos_out[n2*i+j][IA_E] << " ";
                }
                Iout << "\n";
        }

        // Lepton entropy density
	for (int i=0;i<n1;i++) {
                for (int j=0;j<n2;j++) {
                        Iout << eos_out[n2*i+j][IL_S] << " ";
                }
                Iout << "\n";
        }

        // Anti-Lepton entropy density
	for (int i=0;i<n1;i++) {
                for (int j=0;j<n2;j++) {
                        Iout << eos_out[n2*i+j][IA_S] << " ";
                }
                Iout << "\n";
        }

	// Chemical potential
        for (int i=0;i<n1;i++) {
                for (int j=0;j<n2;j++) {
			Iout << eos_out[n2*i+j][IL_MU] << " ";
		}
		Iout << "\n";
	}


	//dP/dn derivative;
        for (int i=0;i<n1;i++) {
                for (int j=0;j<n2;j++) {
                        Iout << der_out[n2*i+j][0] << " ";
                }
                Iout << "\n";
        }

	//ds/dn derivative;
        for (int i=0;i<n1;i++) {
                for (int j=0;j<n2;j++) {
                        Iout << der_out[n2*i+j][1] << " ";
                }
                Iout << "\n";
        }

	//dP/dt derivative;
        for (int i=0;i<n1;i++) {
                for (int j=0;j<n2;j++) {
                        Iout << der_out[n2*i+j][2] << " ";
                }
                Iout << "\n";
        }

	//ds/dt derivative;
        for (int i=0;i<n1;i++) {
                for (int j=0;j<n2;j++) {
                        Iout << der_out[n2*i+j][3] << " ";
                }
                Iout << "\n";
        }


	Iout.close();
        cout << "Done!" << endl << endl;
	return;
}



int main (){

        make_lep_table<0>(); // electrons
        make_lep_table<1>(); // muons

        return 0;

}


