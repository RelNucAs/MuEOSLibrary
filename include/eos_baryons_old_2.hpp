#ifndef EOS_BARYONS_OLD_H
#define EOS_BARYONS_OLD_H

//! \file eos_baryons.hpp
//  \brief Defines EOSTable, which stores information from a tabulated
//         equation of state in CompOSE format.
//
//  Tables should be generated using
//  <a href="https://bitbucket.org/dradice/pycompose">PyCompOSE</a>

///  \warning This code assumes the table to be uniformly spaced in
///           log nb, log t, and yq

#include <string>

const int MAX_SPECIES = 3;

const int ECLOGP  = 0;  //! log (pressure / 1 MeV fm^-3)
const int ECENT   = 1; //! entropy per baryon [kb]
const int ECMUB   = 2;  //! baryon chemical potential [MeV]
const int ECMUQ   = 3;  //! charge chemical potential [MeV]
const int ECMUL   = 4;  //! lepton chemical potential [MeV]
const int ECLOGE  = 5;	//! log (total energy density / 1 MeV fm^-3)
const int ECCS    = 6;	//! sound speed [c]
const int ECNVARS = 7;

// Indexing used to access the data
inline ptrdiff_t index(const int m_nn, const int m_nt, const int m_ny, int iv, int in, int iy, int it) {
	return it + m_nt*(iy + m_ny*(in + m_nn*iv));
}

/// Reads the table file.
void ReadTableFromFile(std::string fname);

//#endif
