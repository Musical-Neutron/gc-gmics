# Plotting scripts for the E-MOSAICS GM globular cluster paper


**Last reviewed:** 0.3.0

A set of scripts and a repository of reduced data to reproduce the plots in the
genetically modified initial condition E-MOSAICS Globular Cluster (GC)
paper.

This README file is located in the main repository of the plotting scripts.
All plotting scripts should be executed from this directory.

## Installation

The package can be installed by cloning the github repository:

```
    git clone git@github.com:Musical-Neutron/gc-gmics.git
    cd gc-gmics
    pip install .
```

## 1.0 Scripts

There are seven scripts that can be executed independently:

* Fig. 2: [fig_02_key_galaxy_properties.py](/fig_02_key_galaxy_properties.py)
  * Plots the halo mass, <i>M</i><sub>200</sub>,
  stellar mass, <i>M</i><sub>*</sub>, and mass in GCs, <i>M</i><sub>GC</sub>,
  vs. lookback time.
* Fig. 3: [fig_03_gas_bh_properties.py](/fig_03_gas_bh_properties.py)
  * Plots the mass of star-forming gas and the central black hole mass vs.
  lookback time.
* Fig. 4: [fig_04_z0_mass_function.py](/fig_04_z0_mass_function.py)
  * Plots the $z=0$ and birth mass functions of all clusters that survive to
  $z=0$.
* Fig. 5: [fig_05_rates_of_change.py](/fig_05_rates_of_change.py)
  * Plots the star formation rate, GC formation rate, GC destruction rate, and
  net change in GC mass vs. lookback time.
* Figs 6 & 7: [fig_06_07_gc_properties.py](/fig_06_07_gc_properties.py)
  * Plots the stellar and GC birth pressures, M<sub>c,*</sub>, and cluster
  formation efficiency vs. lookback time.
  * Plots the specific frequency, <i>T</i><sub>N</sub> and <i>S</i><sub>M</sub>
  vs. lookback time.
* Fig. 8: [fig_08_fhalo.py](/fig_08_fhalo.py)
  * Plots the fraction of mass in halo stars contributed by the destruction of
  star clusters.
* [print_paper_results.py](/print_paper_results.py)
  * Prints information relevant to the paper to stdout.

There is also a master script, [run_all_scripts.py](/run_all_scripts.py),
that will run all of the above scripts when executed. This produces .svg
and .pdf versions of each figure in the paper.

### Supplementary scripts

* [common_functions.py](/common_functions.py)
  * A set of functions common to more than one of the main scripts.
* [process_data.py](/process_data.py)
  * Contains classes and functions to handle basic processing of data.

## 2.0 Data

The [data](/data) directory that contains all files necessary to reproduce the
figures in the paper. There are XX files:

* [file_one.hdf5](/data/file_one.hdf5)
  * Required for all figures.
<!-- * [17_11_z0_data.hdf5](/data/17_11_z0_data.hdf5)
  * Required for Figs. 2&ndash;4, 6 and A1.
* [ludlow2014_logc_vs_logm200h.csv](/data/ludlow2014_logc_vs_logm200h.csv)
  * Required for Figs. 4 \& A1.
* [halo_positions.hdf5](/data/halo_positions.hdf5)
  * Only required for Fig. 5.
* [gammaldi_2021_data.hdf5](/data/gammaldi_2021_data.hdf5)
  * Only required for Fig. 6. -->

## 3.0 Citations

This code and the accompanying data are freely available.

### If you use this code or derivative work

* [O. Newton et al. (2025)](DOI)
* [O. Newton (2021)](https://doi.org/10.5281/zenodo.4708338)

### If you use these data, a derivative work, or results thereof

* [O. Newton et al. (2022)](https://doi.org/10.1093/mnras/stac1316)

If you have any questions or would like help in using the scripts, please
email:
> onewton 'at' cft.edu.pl
