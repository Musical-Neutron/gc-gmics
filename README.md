# Plotting scripts for the E-MOSAICS GM globular cluster paper

**Last reviewed:** 0.5.0

A set of scripts and a repository of reduced data to reproduce the plots in the
genetically modified initial conditions E-MOSAICS Globular Cluster (GC)
paper.

This README file is located in the main repository of the plotting scripts.
All plotting scripts should be executed from this directory.

## Installation

The package can be installed by cloning the github repository:

```bash
    git clone git@github.com:Musical-Neutron/gc-gmics.git
    cd gc-gmics
    pip install .
```

## Configuration

Key settings in the repository are configured using the
[settings.yaml](/settings.yaml) file. Changes to the file are read at script
runtime.

## 1.0 Scripts

There are seven scripts that can be executed independently:

* Fig. 2: [fig_02_key_galaxy_properties.py](/fig_02_key_galaxy_properties.py)
  * Plots the halo mass, *M*<sub>200</sub>,
  stellar mass, *M*<sub>\*</sub>, and mass in GCs, *M*<sub>GC</sub>,
  vs. lookback time, *t*<sub>lb</sub>.
* Fig. 3: [fig_03_gas_bh_properties.py](/fig_03_gas_bh_properties.py)
  * Plots the mass of star-forming gas, *M*<sub>SF</sub>, and the central black
  hole mass, *M*<sub>SMBH</sub> vs. *t*<sub>lb</sub>.
* Fig. 4: [fig_04_z0_mass_function.py](/fig_04_z0_mass_function.py)
  * Plots the $z=0$ and birth mass functions of all clusters that survive to
  $z=0$. Note that the latter is averaged over time and is *not* equivalent to
  the initial cluster mass function.
* Fig. 5: [fig_05_rates_of_change.py](/fig_05_rates_of_change.py)
  * Plots the star formation rate, GC formation rate, GCFR, GC destruction
  rate, GCDR, and the net change in GC mass, $dM_{\rm GC} / dt$ vs.
  *t*<sub>lb</sub>.
* Figs 6 & 7: [fig_06_07_gc_properties.py](/fig_06_07_gc_properties.py)
  * Plots the stellar and GC birth pressures, M<sub>c,*</sub>, and cluster
  formation efficiency vs. *t*<sub>lb</sub>.
  * Plots the specific frequency, *T*<sub>N</sub>, and the specific mass,
  *S*<sub>M</sub> vs. *t*<sub>lb</sub>.
* Fig. 8: [fig_08_fhalo.py](/fig_08_fhalo.py)
  * Plots the fraction of mass in halo stars contributed by the destruction of
  star clusters.
* [print_paper_results.py](/print_paper_results.py)
  * Prints information relevant to the paper to stdout.

There is also a master script, [run_all_scripts.py](/run_all_scripts.py),
that will run all of the above scripts when executed. This produces .svg
and .pdf versions of each figure in the paper.

### Supplementary scripts

* [common_functions.py](/gcgmics/common_functions.py)
  * A set of functions common to more than one of the main scripts.
* [process_data.py](/gcgmics/process_data.py)
  * Contains classes and functions to handle basic processing of data.
* [settings.py](/gcgmics/settings.py)
  * Contains classes and functions to read in settings from
  [settings.yaml](/settings.yaml).

## 2.0 Data

The [data](/data) directory contains all files necessary to reproduce the
figures in the paper. There are XX files:

### Simulation data

* [z2_od_1p000_mz1p7_0p800_HiRes_gcmass_100000_paper_data.hdf5](/data/z2_od_1p000_mz1p7_0p800_HiRes_gcmass_100000_paper_data.hdf5)
  * Data from simulations produced in this work. Required for all figures except Figs 4 & 8.
* [z2_od_1p000_HiRes_gcmass_100000_paper_data.hdf5](/data/z2_od_1p000_HiRes_gcmass_100000_paper_data.hdf5)
  * Data from simulations produced in this work. Required for all figures except Figs 4 & 8.
* [z2_od_1p000_mz1p7_1p100_HiRes_gcmass_100000_paper_data.hdf5](/data/z2_od_1p000_mz1p7_1p100_HiRes_gcmass_100000_paper_data.hdf5)
  * Data from simulations produced in this work. Required for all figures except Figs 4 & 8.
* [z2_od_1p000_mz1p7_0p800_HiRes_noBH_gcmass_100000_paper_data.hdf5](/data/z2_od_1p000_mz1p7_0p800_HiRes_noBH_gcmass_100000_paper_data.hdf5)
  * Data from simulations produced in this work. Required for Fig. 3.
* [z2_od_1p000_HiRes_noBH_gcmass_100000_paper_data.hdf5](/data/z2_od_1p000_HiRes_noBH_gcmass_100000_paper_data.hdf5)
  * Data from simulations produced in this work. Required for Fig. 3.
* [z2_od_1p000_mz1p7_1p100_HiRes_noBH_gcmass_100000_paper_data.hdf5](/data/z2_od_1p000_mz1p7_1p100_HiRes_noBH_gcmass_100000_paper_data.hdf5)
  * Data from simulations produced in this work. Required for Fig. 3.
* [z2_od_1p000_mz1p7_0p800_HiRes_paper_z0_gc_properties.hdf5](/data/z2_od_1p000_mz1p7_0p800_HiRes_paper_z0_gc_properties.hdf5)
  * Data from simulations produced in this work. Required for Figs 4 & 8.
* [z2_od_1p000_HiRes_paper_z0_gc_properties.hdf5](/data/z2_od_1p000_HiRes_paper_z0_gc_properties.hdf5)
  * Data from simulations produced in this work. Required for Figs 4 & 8.
* [z2_od_1p000_mz1p7_1p100_HiRes_paper_z0_gc_properties.hdf5](/data/z2_od_1p000_mz1p7_1p100_HiRes_paper_z0_gc_properties.hdf5)
  * Data from simulations produced in this work. Required for Figs 4 & 8.

### Observational data

* [caldwell_2009_old_clusters.csv](/data/caldwell_2009_old_clusters.csv)
  * Data from [Caldwell et al. (2009)](https://doi.org/10.1088/0004-6256/137/1/94).
  * Required for Fig. 4.
* [caldwell_2009_intermediate_clusters.csv](/data/caldwell_2009_intermediate_clusters.csv)
  * Data from [Caldwell et al. (2009)](https://doi.org/10.1088/0004-6256/137/1/94).
  * Required for Fig. 4.
* [caldwell_2009_young_clusters.csv](/data/caldwell_2009_young_clusters.csv)
  * Data from [Caldwell et al. (2009)](https://doi.org/10.1088/0004-6256/137/1/94).
  * Required for Fig. 4.
* [caldwell_2011_logMstar_FeH.csv](/data/caldwell_2011_logMstar_FeH.csv)
  * Data from [Caldwell et al. (2011)](https://doi.org/10.1088/0004-6256/141/2/61).
  * Required for Fig. 4.
* [baumgardt_2019_mw_gcs.csv](/data/baumgardt_2019_mw_gcs.csv)
  * Data from [Baumgardt et al. (2019)](https://doi.org/10.1093/mnras/sty2997).
  * Data retrieved from <https://people.smp.uq.edu.au/HolgerBaumgardt/globular/>
  * Required for Fig. 4.
* [hunt_2024_N_M_bins.csv](/data/hunt_2024_N_M_bins.csv)
  * Data extracted from [Hunt et al. (2024, fig. 17)](https://doi.org/10.1051/0004-6361/202348662).
  * Required for Fig. 4.
* [horta_2021_fGC.csv](/data/horta_2021_fGC.csv)
  * Data from [Horta et al. (2021)](https://doi.org/10.1093/mnras/staa3598).
  * Required for Fig. 8.

## 3.0 Citations

This code and the accompanying data are freely available.

### If you use this code or derivative work

* [O. Newton et al. (2025)](https://doi.org/10.1093/mnras/staf1226)
* [O. Newton (2025)](https://doi.org/XX.XXXX/zenodo.XXXXXXX)

### If you use these data, a derivative work, or results thereof

* [O. Newton et al. (2025)](https://doi.org/10.1093/mnras/staf1226)

If you have any questions or would like help in using the scripts, please
email:
> onewton 'at' cft.edu.pl
