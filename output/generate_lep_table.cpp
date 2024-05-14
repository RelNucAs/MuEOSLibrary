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

#include "../src/helmholtz_eos/helmholtz_eos.hpp"

const double amu_g = 1.66054e-24;

/* Set table limits */
// 0: electrons, 1: muons
const double nb_min[2] = {1.0000000e-12, 1.00e-12}; // {5.0e+03/amu_g*1.e-39, ...}
const double nb_max[2] = {1.9054607e+00, 1.90e+00}; // {2.0e+16/amu_g*1.e-39, ...}
const double  t_min[2] = {1.0000000e-01, 1.00e-01}; // {5.00e-02, ...}
const double  t_max[2] = {1.4454398e+02, 1.58e+02}; // {1.05e+02, ...}  
const double  y_min[2] = {1.00e-02, 5.0e-08}; // {5.0e-03, ...}
const double  y_max[2] = {6.00e-01, 5.0e-01}; // {5.5e-01, ...}

const int n_nl[2] = {700, 750}; // {1400, 1500}
const int  n_t[2] = {150, 150}; // {300, 300}

void make_lep_table(const int id_L) {
  std::string str;
  if (id_L == 0) {
    str = "electrons";
  } else if (id_L == 1) {
    str = "muons"; 
  } else {
    std::cout << "Make_lep_table: Error in species identifier" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  const double nl_min = y_min[id_L] * nb_min[id_L];
  const double nl_max = y_max[id_L] * nb_max[id_L];

  const double log_nmin = log10(nl_min);
  const double log_nmax = log10(nl_max);

  const double log_tmin = log10(t_min[id_L]);
  const double log_tmax = log10(t_max[id_L]);

  const int n1 = n_nl[id_L];
  const int n2 = n_t[id_L];

  std::vector<double>  n_array, t_array;

  for (int i=0; i<n1; i++) n_array.push_back(log_nmin + static_cast<double>(i) * (log_nmax-log_nmin) / static_cast<double>(n1-1));
  for (int i=0; i<n2 ;i++) t_array.push_back(log_tmin + static_cast<double>(i) * (log_tmax-log_tmin) / static_cast<double>(n2-1));

  std::string  eta_table = "eos_table/" + str + "/eos_" + str + "_eta.txt";
  std::string full_table = "eos_table/" + str + "/eos_" + str + "_table.txt";

  std::cout << std::endl;
  std::cout << "Generating complete table for " << str << " with:";
  std::cout << std::endl;
  std::cout << "# n points = " << n1 << ", # t points = " << n2 << std::endl;
  std::cout << std::scientific << std::setprecision(8);
  std::cout << "nl_min = "  << nl_min      << " fm-3, nl_max = " << nl_max      << " fm-3" << std::endl;
  std::cout << "t_min  = "  << t_min[id_L] << " MeV , t_max  = " << t_max[id_L] <<  " MeV" << std::endl;
  std::cout << std::endl;

  /* Output stream for eta table */
  std::ofstream Ieta(eta_table);
  Ieta << n1 << "\n" << n2 << "\n";
  Ieta << std::scientific << std::setprecision(8);
  
  /* Output stream for full table */
  std::ofstream Iout(full_table);
  Iout << n1 << "\n" << n2 << "\n";
  Iout << std::scientific << std::setprecision(8);

  /* Writing to file log10 of lepton number density array */
  for (int i=0; i<n1; i++) {
    Ieta << n_array[i] << " "; // 
    Iout << n_array[i] << " ";
  }
  Ieta << "\n";
  Iout << "\n";

  /* Writing to file log10 of temperature array */
  for (int i=0; i<n2 ;i++) {
     Ieta << t_array[i] << " ";
     Iout << t_array[i] << " ";
  }
  Ieta << "\n";
  Iout << "\n";
  

  double nLep, temp, eta;

  HelmEOSOutput eos_point;
  std::vector<HelmEOSOutput> eos_out;
  std::vector<HelmEOSDer>    eos_der;

  std::cout << "Loop over number density array (from 0 to " << (n1-1) << "):" << std::endl;
  for (int i=0; i<n1; i++) {
    std::cout << "i = " << i << std::endl;
      for (int j=0; j<n2; j++) {
        nLep = pow(10.,n_array[i]);
        temp = pow(10.,t_array[j]);

        eos_point = eos_helm_full(nLep, temp, id_L);

        eos_out.push_back(eos_point);
        eos_der.push_back(der_cs2(nLep, temp, id_L));
        //eos_der.push_back(der_cs2_num(nLep, temp, id_L));

        eta = (eos_point.mul - mL[id_L]) / temp;
        Ieta << eta << " ";
      }
    Ieta << "\n";
  }

  Ieta.close();

  // Lepton number density
  for (int i=0;i<n1;i++) {
    for (int j=0;j<n2;j++) {
      Iout << eos_out[n2*i+j].nl << " ";
    }
    Iout << "\n";
  }

  // Anti-lepton number density
  for (int i=0;i<n1;i++) {
    for (int j=0;j<n2;j++) {
      Iout << eos_out[n2*i+j].a_nl << " ";
    }
    Iout << "\n";
  }
  
  // Lepton pressure
  for (int i=0;i<n1;i++) {
    for (int j=0;j<n2;j++) {
      Iout << eos_out[n2*i+j].pl << " ";
    }
    Iout << "\n";
  }

  // Anti-Lepton pressure
  for (int i=0;i<n1;i++) {
    for (int j=0;j<n2;j++) {
      Iout << eos_out[n2*i+j].a_pl << " ";
    }
    Iout << "\n";
  }

  // Lepton internal energy density
  for (int i=0;i<n1;i++) {
    for (int j=0;j<n2;j++) {
      Iout << eos_out[n2*i+j].el << " ";
    }
    Iout << "\n";
  }

  // Anti-Lepton internal energy density
  for (int i=0;i<n1;i++) {
    for (int j=0;j<n2;j++) {
      Iout << eos_out[n2*i+j].a_el << " ";
    }
    Iout << "\n";
  }

  // Lepton entropy density
  for (int i=0;i<n1;i++) {
    for (int j=0;j<n2;j++) {
      Iout << eos_out[n2*i+j].sl << " ";
    }
    Iout << "\n";
  }

  // Anti-Lepton entropy density
  for (int i=0;i<n1;i++) {
    for (int j=0;j<n2;j++) {
      Iout << eos_out[n2*i+j].a_sl << " ";
    }
    Iout << "\n";
  }

  // Chemical potential
  for (int i=0;i<n1;i++) {
    for (int j=0;j<n2;j++) {
      Iout << eos_out[n2*i+j].mul << " ";
    }
    Iout << "\n";
  }





  //dP/dn derivative;
  for (int i=0;i<n1;i++) {
    for (int j=0;j<n2;j++) {
      Iout << eos_der[n2*i+j].dPdn << " ";
    }
    Iout << "\n";
  }

  //ds/dn derivative;
  for (int i=0;i<n1;i++) {
    for (int j=0;j<n2;j++) {
      Iout << eos_der[n2*i+j].dsdn << " ";
    }
    Iout << "\n";
  }

  //dP/dt derivative;
  for (int i=0;i<n1;i++) {
    for (int j=0;j<n2;j++) {
      Iout << eos_der[n2*i+j].dPdt << " ";
    }
    Iout << "\n";
  }
  
  //ds/dt derivative;
  for (int i=0;i<n1;i++) {
    for (int j=0;j<n2;j++) {
      Iout << eos_der[n2*i+j].dsdt << " ";
    }
    Iout << "\n";
  }

  Iout.close();

  std::cout << "Done!" << std::endl;
  std::cout << std::endl;
  
  return;
}



int main (){
  make_lep_table(0); // electrons
  make_lep_table(1); // muons

  return 0;
}


