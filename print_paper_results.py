#!/usr/bin/env python3

# Place import files below
import numpy as np

from process_data import EvolutionData, return_print_format_lists
from universal_settings import (sim_list, sim_names)


def main():
    # Load relevant data
    property_list = ['M_200', 'M_star', 'M_GC']
    ev_data = [EvolutionData(sim) for sim in sim_list[:3]]
    print_labels = return_print_format_lists(property_list)

    ####################################################################
    # Print z=0 properties of the galaxies, listed in property_list
    ####################################################################
    l_label_pad = 12
    print_fields = "{:.1f}^+{:.1f}_{:.1f} x 10^{:<9.0f}" * len(property_list)
    property_fields = "{:<28}" * len(property_list)

    print("#" * 72)
    print("z = 0")

    # Print first row of output
    print("{0:<{1}}".format("Property", l_label_pad) +
          property_fields.format(*print_labels))

    ####################################################################
    # Remove when no longer needed
    combined_prop_lists = [[] for _ in np.arange(len(property_list))]
    ####################################################################

    # Iterate over simulation data
    for (sim_data, sim_name) in zip(ev_data, sim_names):
        # print()
        z0_value = []
        z0_value_spread = []
        for p_i, property_to_plot in enumerate(property_list):
            med, spread = sim_data.med_spread(property_to_plot)
            z0_value.append(med[0])
            z0_value_spread.append(spread[:, 0][::-1] - med[0])
            ############################################################
            # Remove when no longer needed
            combined_prop_lists[p_i] = np.concatenate(
                (combined_prop_lists[p_i], getattr(sim_data,
                                                   property_to_plot)[0]))
            ############################################################

        # Determine index to which the median value is raised
        value_idx = np.floor(np.log10(z0_value))

        # Prepare data for printing
        print_data = np.column_stack(
            (z0_value / 10**value_idx,
             z0_value_spread / 10**value_idx[:, np.newaxis],
             value_idx)).reshape(len(z0_value) * 4)

        # Print data for this simulation
        print("{0:<{1}}".format(sim_name, l_label_pad) +
              print_fields.format(*print_data))

    ####################################################################
    # Remove when no longer needed
    ####################################################################
    # combined_prop_lists = [
    #     np.concatenate(item) for item in combined_prop_lists
    # ]
    combined_prop_med = np.nanmedian(combined_prop_lists, axis=1)
    combined_prop_spread = (
        np.nanpercentile(combined_prop_lists,
                         (5., 95.), axis=1) - combined_prop_med)[::-1]
    value_idx = np.floor(np.log10(combined_prop_med))

    print_data = np.column_stack(
        (combined_prop_med / 10**value_idx,
         combined_prop_spread.T / 10**value_idx[:, np.newaxis],
         value_idx)).reshape(len(combined_prop_med) * 4)
    print("{0:<{1}}".format('Combined', l_label_pad) +
          print_fields.format(*print_data))
    ####################################################################

    print("#" * 72)

    ####################################################################
    # Remove when no longer needed
    ####################################################################
    tlb_one = 5  # Gyr
    start_idx = find_nearest_idx(sim_data.t_lb, tlb_one)

    m200_ev = np.nanmedian(getattr(sim_data, 'M_200'), axis=1)
    max_m = np.nanmax(m200_ev)
    end_idx = find_nearest_idx(m200_ev, max_m)

    print("Ratio of halo masses: {}".format(max_m / m200_ev[start_idx]))

    ####################################################################

    return None


def find_nearest_idx(array, value):
    idx = (np.abs(array - value)).argmin()
    return idx


if __name__ == "__main__":
    main()
