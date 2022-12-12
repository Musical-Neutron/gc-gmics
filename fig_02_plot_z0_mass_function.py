#!/usr/bin/env python3

# Place import files below
import os

import astro_tools.cosmology as cosmology
# import astro_tools.figure_manipulation as fm
import astro_tools.array_operations as arr_op
# import astropy
# import astro_tools.geometry as gm
# import h5py
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from cycler import cycler
# from matplotlib.ticker import MultipleLocator, LogLocator, Locator
# from plot_galaxy_properties_vs_z import (plot_line_and_spread, tick_function,
#                                          med_spread)
from plot_data import (DataObject, FigureHandler, plot_line_and_spread,
                       plot_main_and_minor_lines)


def main():
    # Plot settings
    try:
        plt.style.use('./paper.mplstyle')
    except OSError:
        pass
    return None


if __name__ == "__main__":
    main()
