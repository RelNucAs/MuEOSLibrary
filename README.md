# EOS_module

C++ class for implementation of general-purpose Equation of State of dense matter. Individual species are implemented in the following way:
  - baryons: interpolation on 3D CompOSE EOS table, baryon-contribution only (HDF5 table produced by PyCompOSE tool --> https://bitbucket.org/dradice/pycompose/).<br>
             Input parameters:<br>
              - nb : baryon number density<br>
              - T  : temperature<br>
              - Yq : charge fraction<br>
  - leptons (electron and muons): either interpolation on precomputed 2D lepton table (1) or on-the-fly computation with Fermi integrals (2), as in Timmes EOS.<br>                                    Input parameters:<br>
                                    - nl : lepton number density (nl = nb * Yl)<br>
                                    - T  : temperature<br>
  - photons: on-the-fly evaluation based on ideal photon gas EOS
  - neutrinos (not included in the global class, must be added separately): on-the-fly evaluation based on massless ideal Fermi gas EOS

## Preliminary steps
  - Set the correct absolute path of the repository in "include/parameters.hpp" --> abs_path, this is needed by the code in order to correctly read and write EOS tables.
  - Copy the HDF5 file for the baryonic EOS table (must be PyCompOSE output) in "eos_table/baryons/".
  - To compute the leptonic EOS by table interpolation, generate the lepton EOS tables by running "make lep_table" and "./eos_table/compile/generate_lep_table"
  
# Warning:
Always check that the EOS is called in a thermodynamic point that falls within the boundaries of the input tables!

## Generating an output table: how to
Tables including the contributions to the EOS of all the species can be generated using "eos_table/generate_global_table.cpp". The legend describing the content and the units can be read from the first line of the table.

The contribution of muons can be eventually switched off by setting the boolean variable "with_mu = false" at the beginning of the file.

The input parameters are:
  - nb : number density [fm-3]
  - T  : temperature [MeV]
  - Yq : charge fraction (equal to electron fraction, Ye, if muons are not considered)<br>
  If muon contribution is active:
  - Ym : muon fraction (in this case Ye = Yq-Ym by charge neutrality, please note that if Ye < 0 the EOS is not computed and the code jumps to the following iterartion)

The input nb, T, Yq arrays are directly from the baryonic EOS table, in order to avoid to introduce uncertainties associated with the interpolation of the baryonic table. The first two arrays are uniformly spaced in log space, while the third one is uniform in lin space.
The input Ym array is uniformly log-spaced and can be adjusted by the user by acting on the boundaries, "ymmin" and "ymmax", and on the number of points, "n_ym".

If necessary, the user can generate a coarser table by incrementing the variables controlling the step size of the loop over the input arrays ("di" for "n_array", "dj" for "t_array", "dk" for "yq_array" and "dl" for ym_array"). 

To generate the table:
  - Compile the code by running "make global_table" in the parent directory
  - Execute the code by running "./eos_table/generate_global_table" (the output table will be saved in "eos_table/global/")
