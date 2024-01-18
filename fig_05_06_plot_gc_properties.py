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
    figure_handler,
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
    # fig4_out_file = 'fig4_pkstar_pkgc_mcstar_cfe.pdf'
    fig4_out_file = 'fig05_pksfgas_pkgc_mcstar_cfe_newdata.pdf'
    fig5_out_file = 'fig06_TN_SM.pdf'
    out_file_template = '{}_vs_tlb.pdf'

    # Load data for figures
    ev_data = [EvolutionData(sim) for sim in sim_list[:3]]
    # property_list = ['Pk_birth,star', 'Pk_birth,GC', 'Mc_star', 'CFE']
    property_list = ['Pk_SFgas', 'Pk_birth,GC', 'Mc_star', 'CFE']
    ylabels, yscales, ylims = return_plot_format_lists(property_list)

    fig, axs = plt.subplots(
        len(property_list),
        1,
        figsize=(8, 8 * np.ceil(len(property_list) / 3.)),
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
    for a_i, (ax, indiv_ax,
              property_to_plot) in enumerate(zip(axs, indiv_axs,
                                                 property_list)):
        print(property_to_plot)
        for (sim_data, sim, sim_name, tlb_mm,
             tlb_tm) in zip(ev_data, sim_list, sim_names, sim_tlb_major_merger,
                            sim_tlb_target_merger):
            print(sim_data.z)
            med, spread = sim_data.med_spread(property_to_plot)
            print(med)

            # Plot median
            line, = ax.plot(sim_data.t_lb,
                            med,
                            label=sim_name,
                            **plot_styles[sim])
            indiv_line, = indiv_ax.plot(sim_data.t_lb,
                                        med,
                                        label=sim_name,
                                        **plot_styles[sim])
            # if (('0p800' in sim) or
            #     ('1p000_HiRes' in sim)) and ('Pk_SFgas' in property_to_plot):
            #     hr_sim_data = EvolutionData(sim + '_hc')
            #     print('hr')
            #     hr_med, hr_spread = hr_sim_data.med_spread(property_to_plot)
            #     temp_style = copy.deepcopy(plot_styles[sim])
            #     temp_style.update({'ls': '--'})
            #     ax.plot(hr_sim_data.t_lb, hr_med, **temp_style)
            #     indiv_ax.plot(hr_sim_data.t_lb, hr_med, **temp_style)
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
    # axs[0].invert_xaxis()

    # Save figures
    save_figures(fig, fig4_out_file)

    for prop_name, fig in zip(property_list, indiv_figs):
        save_figures(fig, out_file_template.format(prop_name))

    print("TN, SM")

    ################################################################
    # Plot TN and SM
    ################################################################
    new_arrow_length = 0.07
    property_list = ['TN_model0', 'SM_model0']
    ylabels, yscales, ylims = return_plot_format_lists(property_list)

    fig, axs = plt.subplots(
        1,
        len(property_list),
        figsize=(8 * np.ceil(len(property_list)), 8.),
        gridspec_kw={
            "width_ratios": [1. / len(property_list) for _ in property_list],
            "height_ratios": [1],
            "hspace": 0
        })
    indiv_figs = [plt.figure(figsize=(8, 8)) for _ in property_list]
    indiv_axs = [fig.add_subplot(111) for fig in indiv_figs]

    # Iterate over each axis and plot data
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
            # if property_list[a_i] == 'TN_model0':
            #     tn_r200_med, _ = sim_data.med_spread('TN_model0_r200')
            #     temp_style = copy.deepcopy(plot_styles[sim])
            #     temp_style.update({'ls': ':'})
            #     ax.plot(sim_data.t_lb, tn_r200_med, **temp_style)
            #     indiv_ax.plot(sim_data.t_lb, tn_r200_med, **temp_style)
            if property_list[a_i] == 'SM_model0':
                sm_birth_m0_med, _ = sim_data.med_spread('SM_birth_model0')
                # sm_r200_med, _ = sim_data.med_spread('SM_model0_r200')
                temp_style = copy.deepcopy(plot_styles[sim])
                # temp_style.update({'ls': ':'})
                # ax.plot(sim_data.t_lb, sm_r200_med, **temp_style)
                # indiv_ax.plot(sim_data.t_lb, sm_r200_med, **temp_style)
                temp_style.update({'ls': '--'})
                ax.plot(sim_data.t_lb, sm_birth_m0_med, **temp_style)
                indiv_ax.plot(sim_data.t_lb, sm_birth_m0_med, **temp_style)
            # Plot scatter
            ax.fill_between(sim_data.t_lb,
                            *spread,
                            color=line.get_color(),
                            alpha=0.3)
            indiv_ax.fill_between(sim_data.t_lb,
                                  *spread,
                                  color=indiv_line.get_color(),
                                  alpha=0.3)

            # Plot merger arrows
            # First major merger
            mm_x = 1. - (tlb_mm / axis_rescale)
            mm_arrow_properties.update({'fc': line.get_color()})
            # Target major merger
            tm_arrow_properties.update({'fc': line.get_color()})
            if tlb_tm is not None:
                tm_x = 1. - (tlb_tm / axis_rescale)
            plot_merger_arrow(ax,
                              mm_x,
                              new_arrow_length,
                              arrow_properties=mm_arrow_properties,
                              loc='lower')
            if tlb_tm is not None:
                plot_merger_arrow(ax,
                                  tm_x,
                                  new_arrow_length,
                                  arrow_properties=tm_arrow_properties,
                                  loc='lower')
            plot_merger_arrow(indiv_ax,
                              mm_x,
                              new_arrow_length,
                              arrow_properties=mm_arrow_properties,
                              loc='lower')
            if tlb_tm is not None:
                plot_merger_arrow(indiv_ax,
                                  tm_x,
                                  new_arrow_length,
                                  arrow_properties=tm_arrow_properties,
                                  loc='lower')

    # Set common properties of the plots using the figure handler
    figure_handler.set_stacked_figure_properties(axs,
                                                 ylabels,
                                                 yscale=yscales,
                                                 ylims=ylims,
                                                 direction='horizontal')
    figure_handler.set_figure_properties(indiv_axs,
                                         ylabels,
                                         yscale=yscales,
                                         ylims=ylims)

    orig_legend = figure_handler.legend

    legend_markers = [
        plt.Line2D([0, 1], [0, 0], color='k', linestyle='--'),
        plt.Line2D([0, 1], [0, 0], color='k', linestyle='-')
    ]
    legend_labels = ['Without dynamical evolution', 'With dynamical evolution']

    axs[1].legend(legend_markers, legend_labels, loc='upper center')
    axs[1].add_artist(orig_legend)

    # Save figures
    save_figures(fig, fig5_out_file)

    for prop_name, fig in zip(property_list, indiv_figs):
        save_figures(fig, out_file_template.format(prop_name))

    return None


if __name__ == "__main__":
    main()
