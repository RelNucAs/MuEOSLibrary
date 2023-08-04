# MuEOSLibrary

C++ class for implementation of general-purpose Equation of State of dense matter.

Here is the list of the different particle species that the class can account for:
  - **Baryons**<br>
  The EOS of baryons is interpolated on a 3D CompOSE table that includes only the contribution of baryons. The table must be downloaded from the CompOSE website (https://compose.obspm.fr/home) and reprocessed into HDF5 format using **PyCompOSE** (https://bitbucket.org/dradice/pycompose/).<br>
  
  - **Massive leptons (electron and muons)**<br>
  Massive leptons are treated as ideal Fermi gases at finite temperature in thermal equilirium, see *Bludman, S.A., & van Riper, K.A. 1977, ApJ, 212, 859* (https://ui.adsabs.harvard.edu/abs/1977ApJ...212..859B/abstract). To compute the leptonic EOS, one can interpolate on a pre-computed 2D table (that can be generated by the code itself) or perform directly the evaluation of the EOS quantites through the on-the-fly computation of the Generalized Fermi Functions (GFFs).

  - **Photons**<br>
  Photons are treated as an ideal Bose gas in thermal equilibrium with matter and with zero chemical potential. The photon EOS is implemented via the on-the-fly computation of the associated formulae.

  - **Neutrinos**<br>
  Neutrinos are treated as a massless Fermi gas in weak and thermal equilibrium with matter and radiation. The neutrinos EOS is implemented via the on-the-fly computation of the associated formulae.

**Warning:** additive EOS output quantities like the internal energy, pressure and entropy do not include the neutrino contribution by default, it has to be added manually in equilibrium conditions. This choice was made in order to make the EOS call valid also in non-equilibrium conditions, where the user can arbitrarely apply a damping factor to neutrino EOS quantities in order to account for the non equilibrium.

## Preliminary steps
  - Place the HDF5 file containing the baryonic EOS table (must be PyCompOSE output, see description above about baryonic EOS) in the `eos_table/baryons/`  folder.
  - To compute the leptonic EOS through table interpolation, first generate the lepton EOS tables by running **`make lep_table`** in the parent directory.

## EOS class 
The general-purpose EOS of dense matter is implemented through the C++ class named **`EOS_assembled`**, whose definition can be found in `src/eos_species/eos_assembled.cpp`

The class can be initalized by calling the constrctor `EOS_assembled(const int id_eos, const bool el_flag, const bool mu_flag, std::string BarTableName)`, where `id_eos` is an index determining the way in which the leptonic EOS is computed (0: table interpolation, 1: on-the-fly calculation, see description above), `el_flag` (`mu_flag`) ia a boolean variable that can be set to switch on/off the inclusion of electrons (muons) and `BarTableName` is the name of the baryonic EOS table.

The EOS output is computed by calling the class member function `FullEOSOutput EOS_assembled::compute_full_EOS(double nb, double T, double *Y)`.<br>
The input parameters are the following:
```c++
  - double nb               // Baryon number density [fm-3]
  - double T                // Temperature [MeV]
  - double Y[2] = {Ye, Ym}  // Electron and muon fractions [#/baryon]
```

**Warning:** always check that the EOS is called in a thermodynamic point that falls within the boundaries of the input tables!

The EOS output is organized in the `FullEOSOutput` structure which contains the following fields (please note that chemical potentials are including the rest-mass contribution):
```c++
   - double rho                // Mass density [g/cm3]
   - double nb                 // Baryon number density [1/fm3]
   - double T                  // Temperature [MeV]
   - double mb                 // Baryon mass [MeV]
   - double e                  // Internal energy density (without neutrinos) [erg/cm3]
   - double P                  // Pressure (without neutrinos) [erg/cm3]
   - double s                  // Entropy per baryon (without neutrinos) [#/baryon]
   - double chem_pot.mu_n      // Neutron chemical potential [MeV]
   - double chem_pot.mu_p      // Proton chemical potential [MeV]
   - double chem_pot.mu_e      // Electron chemical potential [MeV]
   - double chem_pot.mu_m      // Muon chemical potential [MeV]
   - double chem_pot.mu_nue    // Electron neutrino chemical potential [MeV] 
   - double chem_pot.mu_num    // Muon neutrino chemical potential [MeV]
   - double chem_pot.mu_nux    // Tau neutrino chemical potential [MeV]
   - double Y_part.yh          // Heavy nuclei fraction [#/baryon]
   - double Y_part.ya          // Alpha particle fraction [#/baryon]
   - double Y_part.yn          // Neutron fraction [#/baryon]
   - double Y_part.yp          // Proton fraction [#/baryon]
   - double Y_part.ye          // Electron fraction [#/baryon]
   - double Y_part.ym          // Muon fraction [#/baryon]
   - double Y_part.ynue        // Electron neutrino fraction [#/baryon]
   - double Y_part.yanue       // Electron antineutrino fraction [#/baryon]
   - double Y_part.ynum        // Muon neutrino fraction [#/baryon]
   - double Y_part.yanum       // Muon antineutrino fraction [#/baryon]
   - double Y_part.ynux        // Tau (anti)neutrino fraction [#/baryon]
   - double Y_part.yle         // Net electronic lepton fraction [#/baryons] (yle = ye + ynue - yanue)
   - double Y_part.ylm         // Net muonic lepton fraction [#/baryon] (ylm = ym + ynum - yanum)
   - double nuEOS.z_nue        // Electron neutrino energy per baryon [MeV/baryon]
   - double nuEOS.z_anue       // Electron antineutrino energy per baryon [MeV/baryon]
   - double nuEOS.z_num        // Muon neutrino energy per baryon [MeV/baryon]
   - double nuEOS.z_anum       // Muon antineutrino energy per baryon [MeV/baryon]
   - double nuEOS.z_nux        // Tau (anti)neutrino energy per baryon [MeV/baryon]
   - double nuEOS.z_tot        // Total neutrino energy per baryon [MeV/baryon]
   - double nuEOS.s_nue        // Electron neutrino entropy per baryon [#/baryon]
   - double nuEOS.s_anue       // Electron antineutrino entropy per baryon [#/baryon]
   - double nuEOS.s_num        // Muon neutrino entropy per baryon [#/baryon]
   - double nuEOS.s_anum       // Muon antineutrino entropy per baryon [#/baryon]
   - double nuEOS.s_nux        // Tau (anti)neutrino entropy per baryon [#/baryon]
   - double nuEOS.s_tot        // Total neutrino entropy per baryon [#/baryon]
```

See `main.cpp` for an example on how to initialize the class and compute the full EOS output in a single thermodynamic point.

## Generating an output table: how to
An example of code that generates output EOS tables by combining the contributions from the different species can be found at `output/generate_global_table.cpp`. The code can be compiled and run by running the command **`make output_table`** in the parent directory. The output file is saved in the `output/data/` directory. 

The contribution of muons can be switched on (off) by setting the boolean variable `with_mu` equal to `true` (`false`) at the beginning of the file.

The EOS output is computed by iterating over the following arrays of the input quantities:
```c++
  1. n_array  // Baryon number density [fm-3]
  2. t_array  // Temperature [MeV]
  3. yq_array // Charge fraction [#/baryons] (equal to electron fraction, Ye, if muons are not included)

  If muon contribution is active:
  4. ym_array // Muon fraction [#/baryons] (in this case Ye = Yq - Ym because of charge neutrality,
              //                            if Ye < 0 the code direcly jumps to the following iterartion,
              //                            please note that by doing so the output will not be precisely
              //                            a table but just a collection of output EOS points)
```

The input `n_array`, `t_array`, `yq_array` are read directly from the baryonic EOS table, in order to avoid to introduce uncertainties associated with its interpolation. The first two are uniformly log-spaced, while the third one is uniform in linear space. `ym_array` instead is uniformly log-spaced and can be adjusted by the user by acting on the boundaries, `ymmin` and `ymmax`, and on the number of points, `n_ym`.

If necessary, the user can generate a coarser table by incrementing the variables controlling the step size of the loops over the input arrays (`di` for `n_array`, `dj` for `t_array`, `dk` for `yq_array` and `dl` for `ym_array`). 

## Generating libraries
In order to couple the EOS class to an external code, run `make lib` in the parent directory to generate two static libraries which will be saved under the `lib/` folder. The `libmueos_int.a` (`libmueos_otf.a`) library contains the version of the code where the leptonic EOS is computed via table interpolation (on-the-fly calculation of Fermi-Dirac integrals).
