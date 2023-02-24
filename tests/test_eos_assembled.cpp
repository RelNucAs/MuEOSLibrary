#include <iostream>
#include <ostream>
#include <cmath>

#include "eos_assembled.hpp"

using namespace std;

int main( ) {
	EOS_assembled eos;
	eos.ReadBarTableFromFile("/home/leonardo/Desktop/PhD_work/BNS_muons/EOS_module/eos_table/baryons/SFHo_noel.h5");
	eos.ReadETableFromFile("/home/leonardo/Desktop/PhD_work/BNS_muons/EOS_module/eos_table/electrons/eos_electrons_complete_leo.txt");
	eos.ReadMTableFromFile("/home/leonardo/Desktop/PhD_work/BNS_muons/EOS_module/eos_table/muons/eos_muons_complete_leo.txt");
	
	cout << eos.IsElectronInitialized() << endl;
	cout << eos.IsMuonInitialized() << endl;
	
	//const double * log_te = eos.GetRawELogTemperature();
	//for (int i=0; i<eos.m_nte; i++) cout << exp(log_te[i]) << endl;
	
	const double * ETable = eos.GetElectronTable();
	const double * MTable = eos.GetMuonTable();
       	
	//for (int i=0; i<eos.NTOT*eos.m_ne*eos.m_nte; i++) cout << ETable[i] << endl;
	//for (int i=0; i<eos.BNVARS*eos.m_nn*eos.m_nt; i++) cout << BTable[i] << endl;

	const double * BTable = eos.GetRawTable();
	const double * log_nb = eos.GetRawLogNumberDensity();
	const double * log_tb = eos.GetRawLogTemperature();
	const double * Yq     = eos.GetRawYq();
        
	//const double * log_nm = eos.GetRawMLogNumberDensity();
	double Ytest[] = {Yq[50], 1.e-1};
       
       	//for (int i=0; i<eos.NTOT*eos.m_ne*eos.m_nte; i++) cout << ETable[i] << endl;
	//for (int i=0; i<eos.BNVARS*eos.m_nn*eos.m_nt; i++) cout << BTable[i] << endl;
	cout << endl << endl << eos.Energy(exp(log_nb[50]),exp(log_tb[50]),Ytest) << endl;
	//cout << endl << endl << eos.MuonEnergy(exp(log_nb[50]),exp(log_tb[50]),Ytest) << endl;
	return 0;
}
