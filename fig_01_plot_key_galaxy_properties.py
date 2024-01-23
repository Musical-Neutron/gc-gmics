#!/usr/bin/env python3

# Place import files below
import copy

import matplotlib.pyplot as plt
import numpy as np

from common_functions import plot_merger_arrow, save_figures
from process_data import EvolutionData, return_plot_format_lists
from universal_settings import (
    arrow_length,
    axis_rescale,
    caldwell_2011_m31_mstar_feh_data_file,
    figure_handler,
    gc_mass,
    mm_arrow_properties,
    plot_styles,
    sim_list,
    sim_names,
    sim_tlb_major_merger,
    sim_tlb_target_merger,
    tm_arrow_properties,
)


def main():
    # Plot settings
    try:
        plt.style.use('./paper.mplstyle')
    except OSError:
        pass
        # File location
    out_file = 'fig1_m200_mgc_mstar.pdf'
    out_file_template = '{}_vs_tlb.pdf'
    m31_data = np.genfromtxt(caldwell_2011_m31_mstar_feh_data_file,
                             skip_header=True)

    # Load data for figures
    m31_mass_data = 10.**m31_data[:, 1]  # Msun
    property_list = ['M_200', 'M_star', 'M_GC']
    ev_data = [EvolutionData(sim) for sim in sim_list[:3]]
    ylabels, yscales, ylims = return_plot_format_lists(property_list)
    # (mstar_r200_ylabels, mstar_r200_yscales,
    #  mstar_r200_ylims) = return_plot_format_lists(['M_star_r200'])
    # (mgc_r200_ylabels, mgc_r200_yscales,
    #  mgc_r200_ylims) = return_plot_format_lists(['M_GC_r200'])

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
    _, norm_axs = plt.subplots(
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

    # Iterate over each axis and plot data
    for a_i, (ax, norm_ax, indiv_ax, property_to_plot) in enumerate(
            zip(axs, norm_axs, indiv_axs, property_list)):

        norm_ax.axhline(1., color='k', linestyle=':', zorder=0)
        for (sim_data, sim, sim_name, tlb_mm,
             tlb_tm) in zip(ev_data, sim_list, sim_names, sim_tlb_major_merger,
                            sim_tlb_target_merger):
            med, spread = sim_data.med_spread(property_to_plot)
            ref_tlb = 7.7
            if tlb_tm is not None:
                ref_tlb = tlb_tm
            norm_med = med / med[sim_data.t_lb >= ref_tlb][0]
            norm_spread = spread / med[sim_data.t_lb >= ref_tlb][0]

            # Plot median
            line, = ax.plot(sim_data.t_lb,
                            med,
                            label=sim_name,
                            **plot_styles[sim])
            indiv_line, = indiv_ax.plot(sim_data.t_lb,
                                        med,
                                        label=sim_name,
                                        **plot_styles[sim])
            norm_ax.plot(sim_data.t_lb,
                         norm_med,
                         label=sim_name,
                         **plot_styles[sim])
            # if ('0p800' in sim) or ('1p000_HiRes' in sim):
            #     hr_sim_data = EvolutionData(sim + '_hc')
            #     hr_med, hr_spread = hr_sim_data.med_spread(property_to_plot)
            #     temp_style = copy.deepcopy(plot_styles[sim])
            #     temp_style.update({'ls': '--'})
            #     ax.plot(hr_sim_data.t_lb, hr_med, **temp_style)
            #     indiv_ax.plot(hr_sim_data.t_lb, hr_med, **temp_style)
            if property_list[a_i] == 'M_star':
                star_r200_med, _ = sim_data.med_spread('M_star_r200')
                temp_style = copy.deepcopy(plot_styles[sim])
                temp_style.update({'ls': ':'})
                ax.plot(sim_data.t_lb, star_r200_med, **temp_style)
                indiv_ax.plot(sim_data.t_lb, star_r200_med, **temp_style)
            if property_list[a_i] == 'M_GC':
                gc_r200_med, _ = sim_data.med_spread('M_GC_r200')
                temp_style = copy.deepcopy(plot_styles[sim])
                temp_style.update({'ls': ':'})
                temp_style.update({'ls': ':'})
                ax.plot(sim_data.t_lb, gc_r200_med, **temp_style)
                indiv_ax.plot(sim_data.t_lb, gc_r200_med, **temp_style)
            # Plot scatter
            ax.fill_between(sim_data.t_lb,
                            *spread,
                            color=line.get_color(),
                            alpha=0.3)
            indiv_ax.fill_between(sim_data.t_lb,
                                  *spread,
                                  color=indiv_line.get_color(),
                                  alpha=0.3)
            norm_ax.fill_between(sim_data.t_lb,
                                 *norm_spread,
                                 color=line.get_color(),
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

    mgc_m31 = np.nansum(m31_mass_data[m31_mass_data >= gc_mass])
    axs[-1].scatter(0.15,
                    mgc_m31,
                    color='k',
                    marker='s',
                    label='M31 (observed)')

    # Set common properties of the plots using the figure handler
    figure_handler.set_stacked_figure_properties(axs,
                                                 ylabels,
                                                 yscale=yscales,
                                                 ylims=ylims)
    figure_handler.set_stacked_figure_properties(norm_axs,
                                                 ylabels,
                                                 yscale=['linear'] *
                                                 len(ylabels),
                                                 ylims=[[0, 2], [0, 2.5],
                                                        [0, 1.5]])
    figure_handler.set_figure_properties(indiv_axs,
                                         ylabels,
                                         yscale=yscales,
                                         ylims=ylims)

    plt.show()

    # Save figures
    save_figures(fig, out_file)

    for prop_name, fig in zip(property_list, indiv_figs):
        save_figures(fig, out_file_template.format(prop_name))

    return None


if __name__ == "__main__":
    main()
