#!/usr/bin/env python3

# Place import files below
import fig_01_plot_key_galaxy_properties
import fig_02_plot_z0_mass_function
import fig_03_plot_rates_of_change
import fig_04_05_plot_gc_properties
import fig_06_plot_gas_bh_properties
import fig_07_plot_fhalo
import print_paper_results


def main():
    # Plot Fig. 1
    print("1. Plotting Fig. 1")
    fig_01_plot_key_galaxy_properties.main()

    # Plot Fig. 2
    print("2. Plotting Fig. 2")
    fig_02_plot_z0_mass_function.main()

    # Plot Fig. 3
    print("3. Plotting Fig. 3")
    fig_03_plot_rates_of_change.main()

    # Plot Figs 4 and 5
    print("4. Plotting Figs 4 & 5")
    fig_04_05_plot_gc_properties.main()

    # Plot Fig. 6
    print("5. Plotting Fig. 6")
    fig_06_plot_gas_bh_properties.main()

    # Plot Fig. 7
    print("6. Plotting Fig. 7")
    fig_07_plot_fhalo.main()

    # Print paper results
    print("7. Printing relevant data to screen")
    print_paper_results.main()

    return None


if __name__ == "__main__":
    main()
