#!/usr/bin/env python3

# Place import files below
import copy

import matplotlib.pyplot as plt

from common_functions import plot_merger_arrow, save_figures
from process_data import EvolutionData, return_plot_format_lists
from universal_settings import (arrow_length, axis_rescale, figure_handler,
                                mm_arrow_properties, plot_styles, sim_list,
                                sim_names, sim_tlb_major_merger,
                                sim_tlb_target_merger, tm_arrow_properties)


def main():
    # Plot settings
    try:
        plt.style.use('./paper.mplstyle')
    except OSError:
        pass
        # File location
    out_file = 'fig6_mgassf_mbh.pdf'
    out_file_template = '{}_vs_tlb.pdf'

    # Load data for figures
    property_list = ['M_gas,SF', 'M_BH']
    ev_data = [EvolutionData(sim) for sim in sim_list[:3]]
    ylabels, yscales, ylims = return_plot_format_lists(property_list)

    fig, axs = plt.subplots(
        len(property_list),
        1,
        figsize=(8, 8 * (1. + 0.3 * (len(property_list) - 1.))),
        sharex=True,
        gridspec_kw={
            "width_ratios": [1],
            "height_ratios": [1. / len(property_list) for _ in property_list],
            "wspace": 0,
            "hspace": 0
        })
    indiv_figs = [plt.figure(figsize=(8, 8)) for _ in property_list]
    indiv_axs = [fig.add_subplot(111) for fig in indiv_figs]

    # Iterate over each axes object and plot data
    for a_i, (ax, indiv_ax,
              property_to_plot) in enumerate(zip(axs, indiv_axs,
                                                 property_list)):
        for (sim_data, sim, sim_name, tlb_mm,
             tlb_tm) in zip(ev_data, sim_list, sim_names, sim_tlb_major_merger,
                            sim_tlb_target_merger):
            med, spread = sim_data.med_spread(property_to_plot)
            # Plot median
            line, = ax.plot(sim_data.t_lb,
                            med,
                            label=sim_name,
                            **plot_styles[sim])
            indiv_line, = indiv_ax.plot(sim_data.t_lb,
                                        med,
                                        label=sim_name,
                                        **plot_styles[sim])
            if property_list[a_i] == 'M_gas,SF':
                msfgas_r200_med, _ = sim_data.med_spread('M_gas,SF_r200')
                temp_style = copy.deepcopy(plot_styles[sim])
                temp_style.update({'ls': ':'})
                ax.plot(sim_data.t_lb, msfgas_r200_med, **temp_style)
                indiv_ax.plot(sim_data.t_lb, msfgas_r200_med, **temp_style)
            # Plot scatter
            ax.fill_between(sim_data.t_lb,
                            *spread,
                            color=line.get_color(),
                            alpha=0.3)
            indiv_ax.fill_between(sim_data.t_lb,
                                  *spread,
                                  color=indiv_line.get_color(),
                                  alpha=0.3)

            # Plot merger arrows on first and last panels
            # First major merger
            mm_x = 1. - (tlb_mm / axis_rescale)
            mm_arrow_properties.update({'fc': line.get_color()})
            # Target major merger
            tm_arrow_properties.update({'fc': line.get_color()})
            if tlb_tm is not None:
                tm_x = 1. - (tlb_tm / axis_rescale)
            if a_i == 0:
                plot_merger_arrow(ax,
                                  mm_x,
                                  arrow_length,
                                  arrow_properties=mm_arrow_properties,
                                  loc='upper')
                if tlb_tm is not None:
                    plot_merger_arrow(ax,
                                      tm_x,
                                      arrow_length,
                                      arrow_properties=tm_arrow_properties,
                                      loc='upper')
            if a_i == len(axs) - 1:
                plot_merger_arrow(ax,
                                  mm_x,
                                  arrow_length,
                                  arrow_properties=mm_arrow_properties,
                                  loc='lower')
                if tlb_tm is not None:
                    plot_merger_arrow(ax,
                                      tm_x,
                                      arrow_length,
                                      arrow_properties=tm_arrow_properties,
                                      loc='lower')

            plot_merger_arrow(indiv_ax,
                              mm_x,
                              arrow_length,
                              arrow_properties=mm_arrow_properties,
                              loc='upper')
            if tlb_tm is not None:
                plot_merger_arrow(indiv_ax,
                                  tm_x,
                                  arrow_length,
                                  arrow_properties=tm_arrow_properties,
                                  loc='upper')
            plot_merger_arrow(indiv_ax,
                              mm_x,
                              arrow_length,
                              arrow_properties=mm_arrow_properties,
                              loc='lower')
            if tlb_tm is not None:
                plot_merger_arrow(indiv_ax,
                                  tm_x,
                                  arrow_length,
                                  arrow_properties=tm_arrow_properties,
                                  loc='lower')

    # Set common properties of the plots using the figure handler
    figure_handler.set_stacked_figure_properties(axs,
                                                 ylabels,
                                                 yscale=yscales,
                                                 ylims=ylims)
    figure_handler.set_figure_properties(indiv_axs,
                                         ylabels,
                                         yscale=yscales,
                                         ylims=ylims)

    # Save figures
    save_figures(fig, out_file)

    for prop_name, fig in zip(property_list, indiv_figs):
        save_figures(fig, out_file_template.format(prop_name))

    return None


if __name__ == "__main__":
    main()
