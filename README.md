# Plotting scripts for the E-MOSAICS GM globular cluster paper

**Last reviewed:** 0.6.1

A set of scripts and a repository of reduced data to reproduce the plots in the
genetically modified initial conditions E-MOSAICS Globular Cluster (GC)
paper.

This README file is located in the main repository of the plotting scripts.
All plotting scripts should be executed from this directory.

## 1.0 Installation

### 1.1 Pre-installation dependencies

This package uses the [cairo graphics library](https://cairographics.org/).
In particular, as of this release it requires

```
  cairo >= 1.15.10
```

For relevant information about the current cairo dependencies please
review the [Pycairo package documentation](https://pycairo.readthedocs.io/en/latest/).

### 1.2 Package installation

Once you have installed cairo, this package can be installed by cloning it
from the github repository:

```bash
    git clone git@github.com:Musical-Neutron/gc-gmics.git
    cd gc-gmics
    pip install .
```

### 1.3 Tests

If you plan to review the package test suite, run the following from the root
directory

```bash
    pip install .[test]
    pytest
```

## 2.0 Configuration (`settings.yaml`)

Key settings in the repository are configured using the
[settings.yaml](/settings.yaml) file.

### 2.1 Core configuration

- **Data Management**:
  - `data_dir` defines the base directory where relevant data are stored.
  - Relevant observational data sets, set under `Observations`, are also stored
  in `data_dir`.

- **Analysis settings**:
  - `gc_mass` sets the minimum globular cluster mass, in units of solar masses.
  - `cluster_mass_limit` sets the lower mass below which clusters are not
  tracked. This reflects the settings with which the simulations were run.

### 2.2 Output control

Under `Printing`, `print_property_list` specifies which properties will appear
in terminal outputs when called by [print_paper_results.py](/print_paper_results.py).
The full list of available properties is provided in `Property dictionary`

### 2.3 Visualisation

In addition to providing a complete list of available galaxy properties, the
`Property dictionary` controls plot formatting:

- `ylim` sets the y-axis lower and upper limits. `~` is parsed as `None` when
  the settings are read.
- `yscale` sets the scale of the y-axis (linear or log-scale).
- `ylabel` enables the customisation of the y-axis label. LaTeX can be
  provided if the appropriate packages are installed.
- `printlabel` specifies how this property will appear when printed to screen.

Changes to [settings.yaml](/settings.yaml) are read at script runtime.

## 3.0 Scripts

There are seven scripts that can be executed independently:

- Fig. 2: [fig_02_key_galaxy_properties.py](/fig_02_key_galaxy_properties.py)
  - Plots the halo mass, *M*<sub>200</sub>,
  stellar mass, *M*<sub>\*</sub>, and mass in GCs, *M*<sub>GC</sub>,
  vs. lookback time, *t*<sub>lb</sub>.
- Fig. 3: [fig_03_gas_bh_properties.py](/fig_03_gas_bh_properties.py)
  - Plots the mass of star-forming gas, *M*<sub>SF</sub>, and the central black
  hole mass, *M*<sub>SMBH</sub> vs. *t*<sub>lb</sub>.
- Fig. 4: [fig_04_z0_mass_function.py](/fig_04_z0_mass_function.py)
  - Plots the $z=0$ and birth mass functions of all clusters that survive to
  $z=0$. Note that the latter is averaged over time and is *not* equivalent to
  the initial cluster mass function.
- Fig. 5: [fig_05_rates_of_change.py](/fig_05_rates_of_change.py)
  - Plots the star formation rate, GC formation rate, GCFR, GC destruction
  rate, GCDR, and the net change in GC mass, $dM_{\rm GC} / dt$ vs.
  *t*<sub>lb</sub>.
- Figs 6 & 7: [fig_06_07_gc_properties.py](/fig_06_07_gc_properties.py)
  - Plots the stellar and GC birth pressures, M<sub>c,*</sub>, and cluster
  formation efficiency vs. *t*<sub>lb</sub>.
  - Plots the specific frequency, *T*<sub>N</sub>, and the specific mass,
  *S*<sub>M</sub> vs. *t*<sub>lb</sub>.
- Fig. 8: [fig_08_fhalo.py](/fig_08_fhalo.py)
  - Plots the fraction of mass in halo stars contributed by the destruction of
  star clusters.
- [print_paper_results.py](/print_paper_results.py)
  - Prints information relevant to the paper to stdout.

There is also a master script, [run_all_scripts.py](/run_all_scripts.py),
that will run all of the above scripts when executed. This produces .svg
and .pdf versions of each figure in the paper.

### 3.1 Supplementary scripts

- [common_functions.py](/gcgmics/common_functions.py)
  - A set of functions common to more than one of the main scripts.
- [process_data.py](/gcgmics/process_data.py)
  - Contains classes and functions to handle basic processing of data.
- [settings.py](/gcgmics/settings.py)
  - Contains classes and functions to read in settings from
  [settings.yaml](/settings.yaml).

## 4.0 Data

The [data](/data) directory contains all files necessary to reproduce the
figures in the paper. There are nine files containing reduced data from the
simulations produced for this work, and seven files containing relevant
observational data:

### 4.1 Simulation data

The following files contain data from simulations produced in this work.

File name | Figures needed for
:--- | :---
[suppressed_mgc_100000_data.hdf5](/data/suppressed_mgc_100000_data.hdf5) | All except Figs 4 & 8
[organic_mgc_100000_data.hdf5](/data/organic_mgc_100000_data.hdf5) | All except Figs 4 & 8
[enhanced_mgc_100000_data.hdf5](/data/enhanced_mgc_100000_data.hdf5) | All except Figs 4 & 8
[suppressed_noBH_mgc_100000_data.hdf5](/data/suppressed_noBH_mgc_100000_data.hdf5) | Fig. 3
[organic_noBH_mgc_100000_data.hdf5](/data/organic_noBH_mgc_100000_data.hdf5) | Fig. 3
[enhanced_noBH_mgc_100000_data.hdf5](/data/enhanced_noBH_mgc_100000_data.hdf5) | Fig. 3
[z2_od_1p000_mz1p7_0p800_HiRes_paper_z0_gc_properties.hdf5](/data/z2_od_1p000_mz1p7_0p800_HiRes_paper_z0_gc_properties.hdf5) | Figs 4 & 8
[z2_od_1p000_HiRes_paper_z0_gc_properties.hdf5](/data/z2_od_1p000_HiRes_paper_z0_gc_properties.hdf5) | Figs 4 & 8
[z2_od_1p000_mz1p7_1p100_HiRes_paper_z0_gc_properties.hdf5](/data/z2_od_1p000_mz1p7_1p100_HiRes_paper_z0_gc_properties.hdf5) | Figs 4 & 8

### 4.2 Observational data

File name | Figures needed for | Citation | Notes
:--- | :---: | :--- | :---
[caldwell_2009_old_clusters.csv](/data/caldwell_2009_old_clusters.csv) | Fig. 4 | [Caldwell et al. (2009)](https://doi.org/10.1088/0004-6256/137/1/94)
[caldwell_2009_intermediate_clusters.csv](/data/caldwell_2009_intermediate_clusters.csv) | Fig. 4 | [Caldwell et al. (2009)](https://doi.org/10.1088/0004-6256/137/1/94)
[caldwell_2009_young_clusters.csv](/data/caldwell_2009_young_clusters.csv) | Fig. 4 | [Caldwell et al. (2009)](https://doi.org/10.1088/0004-6256/137/1/94)
[caldwell_2011_logMstar_FeH.csv](/data/caldwell_2011_logMstar_FeH.csv) | Fig. 4 | [Caldwell et al. (2011)](https://doi.org/10.1088/0004-6256/141/2/61)
[baumgardt_2019_mw_gcs.csv](/data/baumgardt_2019_mw_gcs.csv) | Fig. 4 | [Baumgardt et al. (2019)](https://doi.org/10.1093/mnras/sty2997) | Data retrieved from <https://people.smp.uq.edu.au/HolgerBaumgardt/globular/>.
[hunt_2024_N_M_bins.csv](/data/hunt_2024_N_M_bins.csv) | Fig. 4 | [Hunt et al. (2024, fig. 17)](https://doi.org/10.1051/0004-6361/202348662)
[horta_2021_fGC.csv](/data/horta_2021_fGC.csv) | Fig. 8 | [Horta et al. (2021)](https://doi.org/10.1093/mnras/staa3598)

## 3.0 Citations

This code and the accompanying data are freely available.

### If you use this code or derivative work

- [O. Newton et al. (2025)](https://doi.org/10.1093/mnras/staf1226)

### If you use these data, a derivative work, or results thereof

- [O. Newton et al. (2025)](https://doi.org/10.1093/mnras/staf1226)

If you have any questions or would like help in using the scripts, please
email:
> onewton 'at' cft.edu.pl
