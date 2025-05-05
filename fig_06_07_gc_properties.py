#!/usr/bin/env python3

# Place import files below
import matplotlib.pyplot as plt
import numpy as np

from common_functions import (
    get_scaled_arrow_properties,
    plot_merger_arrow,
    save_figures,
)
from process_data import EvolutionData, return_plot_format_lists
from universal_settings import (
    arrow_length,
    axis_rescale,
    baumgardt_2019_mw_cluster_file,
    caldwell_2011_m31_mstar_feh_data_file,
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
    fig4_out_file = 'fig06_pksfgas_pkgc_mcstar_cfe_newdata.pdf'
    fig5_out_file = 'fig07_TN_SM.pdf'
    out_file_template = '{}_vs_tlb.pdf'

    mw_cl_data = np.genfromtxt(baumgardt_2019_mw_cluster_file,
                               skip_header=True)
    m31_data = np.genfromtxt(caldwell_2011_m31_mstar_feh_data_file,
                             skip_header=True)

    # Load data for figures
    mw_cl_mass_data = mw_cl_data[:, 1]  # Msun
    m31_mass_data = 10.**m31_data[:, 1]  # Msun
    ev_data = [EvolutionData(sim) for sim in sim_list[:3]]
    property_list = ['Pk_SFgas', 'Pk_birth,GC', 'Mc_star', 'CFE']
    ylabels, yscales, ylims = return_plot_format_lists(property_list)

    ####################################################################
    # Create figures
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

    ####################################################################
    # Iterate over each axis and plot data
    for a_i, (ax, indiv_ax,
              property_to_plot) in enumerate(zip(axs, indiv_axs,
                                                 property_list)):
        for (sim_data, sim, sim_name, tlb_mm,
             tlb_tm) in zip(ev_data, sim_list, sim_names, sim_tlb_major_merger,
                            sim_tlb_target_merger):
            med, spread = sim_data.med_spread(property_to_plot)

            ############################################################
            # Plot lines and fills
            line, = ax.plot(sim_data.t_lb,
                            med,
                            label=sim_name,
                            **plot_styles[sim])
            ax.fill_between(sim_data.t_lb,
                            *spread,
                            lw=0,
                            color=line.get_color(),
                            alpha=0.3)

            indiv_line, = indiv_ax.plot(sim_data.t_lb,
                                        med,
                                        label=sim_name,
                                        **plot_styles[sim])
            indiv_ax.fill_between(sim_data.t_lb,
                                  *spread,
                                  lw=0,
                                  color=indiv_line.get_color(),
                                  alpha=0.3)

            ############################################################
            # Plot merger arrows on first and last panels
            # First major merger
            mm_x = 1. - (tlb_mm / axis_rescale)
            mm_props = {**mm_arrow_properties, 'fc': line.get_color()}
            (new_arrow_length, mm_props) = get_scaled_arrow_properties(
                arrow_length, mm_props,
                ax.get_gridspec()._row_height_ratios[a_i])

            # Target major merger
            tm_props = {**tm_arrow_properties, 'fc': line.get_color()}
            _, tm_props = get_scaled_arrow_properties(
                arrow_length, tm_props,
                ax.get_gridspec()._row_height_ratios[a_i])
            if tlb_tm is not None:
                tm_x = 1. - (tlb_tm / axis_rescale)

            # Main axis arrows
            if a_i == 0:
                plot_merger_arrow(ax, mm_x, new_arrow_length, mm_props,
                                  'upper')
                if tlb_tm is not None:
                    plot_merger_arrow(ax, tm_x, new_arrow_length, tm_props,
                                      'upper')
            elif a_i == len(axs) - 1:
                plot_merger_arrow(ax, mm_x, new_arrow_length, mm_props,
                                  'lower')
                if tlb_tm is not None:
                    plot_merger_arrow(ax, tm_x, new_arrow_length, tm_props,
                                      'lower')

            # Individual axes and arrows
            for loc in ['upper', 'lower']:
                plot_merger_arrow(indiv_ax, mm_x, new_arrow_length, mm_props,
                                  loc)
                if tlb_tm is not None:
                    plot_merger_arrow(indiv_ax, tm_x, new_arrow_length,
                                      tm_props, loc)
            ############################################################

    ####################################################################
    # Set common properties of the plots using the figure handler
    figure_handler.set_stacked_figure_properties(axs,
                                                 ylabels,
                                                 yscale=yscales,
                                                 ylims=ylims)
    figure_handler.set_figure_properties(indiv_axs,
                                         ylabels,
                                         yscale=yscales,
                                         ylims=ylims)

    ####################################################################
    # Save figures
    save_figures(fig, fig4_out_file)

    for prop_name, fig in zip(property_list, indiv_figs):
        save_figures(fig, out_file_template.format(prop_name))

    ####################################################################
    # Plot TN and SM
    ####################################################################
    new_arrow_length = 0.07
    property_list = ['TN_model0', 'SM_model0']
    ylabels, yscales, ylims = return_plot_format_lists(property_list)
    birth_map = {'SM_model0': 'SM_birth_model0'}

    ####################################################################
    # Create figures
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

    ####################################################################
    # Iterate over each axis and plot data
    for a_i, (ax, indiv_ax,
              property_to_plot) in enumerate(zip(axs, indiv_axs,
                                                 property_list)):

        for (sim_data, sim, sim_name, tlb_mm,
             tlb_tm) in zip(ev_data, sim_list, sim_names, sim_tlb_major_merger,
                            sim_tlb_target_merger):
            med, spread = sim_data.med_spread(property_to_plot)

            ############################################################
            # Plot lines and fills
            line, = ax.plot(sim_data.t_lb,
                            med,
                            label=sim_name,
                            **plot_styles[sim])
            ax.fill_between(sim_data.t_lb,
                            *spread,
                            lw=0,
                            color=line.get_color(),
                            alpha=0.3)

            indiv_line, = indiv_ax.plot(sim_data.t_lb,
                                        med,
                                        label=sim_name,
                                        **plot_styles[sim])
            indiv_ax.fill_between(sim_data.t_lb,
                                  *spread,
                                  lw=0,
                                  color=indiv_line.get_color(),
                                  alpha=0.3)

            ############################################################
            # R200 lines
            if property_to_plot in birth_map:
                birth_med, _ = sim_data.med_spread(birth_map[property_to_plot])
                style = {**plot_styles[sim], 'ls': '--'}
                ax.plot(sim_data.t_lb, birth_med, **style)
                indiv_ax.plot(sim_data.t_lb, birth_med, **style)

            ############################################################
            # Plot merger arrows
            # First major merger
            mm_x = 1. - (tlb_mm / axis_rescale)
            mm_props = {**mm_arrow_properties, 'fc': line.get_color()}
            (new_arrow_length, mm_props) = get_scaled_arrow_properties(
                arrow_length, mm_props,
                ax.get_gridspec()._row_height_ratios[a_i])

            # Target major merger
            tm_props = {**tm_arrow_properties, 'fc': line.get_color()}
            _, tm_props = get_scaled_arrow_properties(
                arrow_length, tm_props,
                ax.get_gridspec()._row_height_ratios[a_i])
            if tlb_tm is not None:
                tm_x = 1. - (tlb_tm / axis_rescale)

            plot_merger_arrow(ax, mm_x, new_arrow_length, mm_props, 'lower')
            plot_merger_arrow(indiv_ax, mm_x, new_arrow_length, mm_props,
                              'lower')
            if tlb_tm is not None:
                plot_merger_arrow(ax, tm_x, new_arrow_length, tm_props,
                                  'lower')
                plot_merger_arrow(indiv_ax, tm_x, new_arrow_length, tm_props,
                                  'lower')

    ####################################################################
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

    ####################################################################
    # Legend
    orig_legend = figure_handler.legend

    legend_markers = [
        plt.Line2D([0, 1], [0, 0], color='k', linestyle='--'),
        plt.Line2D([0, 1], [0, 0], color='k', linestyle='-')
    ]
    legend_labels = ['Without dynamical evolution', 'With dynamical evolution']

    axs[1].legend(legend_markers,
                  legend_labels,
                  loc='upper right',
                  alignment='right')
    axs[0].add_artist(orig_legend)

    ####################################################################
    # Save figures
    save_figures(fig, fig5_out_file)

    for prop_name, fig in zip(property_list, indiv_figs):
        save_figures(fig, out_file_template.format(prop_name))

    return None


if __name__ == "__main__":
    main()
