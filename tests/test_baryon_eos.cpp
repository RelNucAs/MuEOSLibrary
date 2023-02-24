#include "parameters.hpp"
#include "eos_baryons.hpp"

int main() {
	std::string h5name = parameters::abs_path + "eos_table/hadrons/SFHo.h5";
	int m_nn, m_nt, m_ny;
	double m_log_nb[];
	double m_log_t[];
	double m_yq[];
	double m_table[];
	ReadTableFromFile(h5name, &m_nn, &m_nt, &m_ny, m_log_nb, m_log_t, m_yq, m_table);
	return 0;
}
