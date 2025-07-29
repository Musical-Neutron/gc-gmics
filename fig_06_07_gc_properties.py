#!/usr/bin/env python3
# Place import files below
import matplotlib.pyplot as plt
import numpy as np

from gcgmics.common_functions import (
    get_scaled_arrow_properties,
    plot_merger_arrow,
    save_figures,
)
from gcgmics.process_data import EvolutionData, return_plot_format_lists
from gcgmics.settings import Plotting, Simulations


def main():
    # Plot settings
    try:
        plt.style.use("./paper.mplstyle")
    except OSError:
        pass

    # File locations
    fig6_out_file = "fig06_pksfgas_pkgc_mcstar_cfe.pdf"
    fig7_out_file = "fig07_TN_SM.pdf"
    out_file_template = "{}_vs_tlb.pdf"

    # Load data for figures
    ev_data = [EvolutionData(sim) for sim in Simulations["sim_list"][:3]]
    property_list = ["Pk_SFgas", "Pk_birth,GC", "Mc_star", "CFE"]
    ylabels, yscales, ylims = return_plot_format_lists(property_list)

    ####################################################################
    # Plot Pk, Mc*, CFE
    ####################################################################
    # Create figures
    fig, axs = plt.subplots(
        len(property_list),
        1,
        figsize=(8, 8 * np.ceil(len(property_list) / 3.0)),
        sharex=True,
        gridspec_kw={
            "width_ratios": [1],
            "height_ratios": [1.0 / len(property_list) for _ in property_list],
            "wspace": 0,
            "hspace": 0,
        },
    )
    indiv_figs = [plt.figure(figsize=(8, 8)) for _ in property_list]
    indiv_axs = [fig.add_subplot(111) for fig in indiv_figs]

    ####################################################################
    # Iterate over each axis and plot data
    for a_i, (ax, indiv_ax, property_to_plot) in enumerate(
        zip(axs, indiv_axs, property_list)
    ):
        for sim_data, sim, sim_name, tlb_mm, tlb_tm in zip(
            ev_data,
            Simulations["sim_list"],
            Simulations["sim_names"],
            Simulations["sim_tlb_major_merger"],
            Simulations["sim_tlb_target_merger"],
        ):
            med, spread = sim_data.med_spread(property_to_plot)

            ############################################################
            # Plot lines and fills
            (line,) = ax.plot(
                sim_data.t_lb, med, label=sim_name, **Plotting["plot_styles"][sim]
            )
            ax.fill_between(
                sim_data.t_lb, *spread, lw=0, color=line.get_color(), alpha=0.3
            )

            (indiv_line,) = indiv_ax.plot(
                sim_data.t_lb, med, label=sim_name, **Plotting["plot_styles"][sim]
            )
            indiv_ax.fill_between(
                sim_data.t_lb, *spread, lw=0, color=indiv_line.get_color(), alpha=0.3
            )

            ############################################################
            # Plot merger arrows on first and last panels
            # First major merger
            mm_x = 1.0 - (tlb_mm / Plotting["axis_rescale"])
            mm_props = {
                **Plotting["mm_arrow_properties"],
                "fc": line.get_color(),
            }
            (new_arrow_length, stack_mm_props) = get_scaled_arrow_properties(
                Plotting["arrow_length"],
                mm_props,
                ax.get_gridspec()._row_height_ratios[a_i] * len(property_list) / 3.0,
            )

            # Target major merger
            tm_props = {
                **Plotting["tm_arrow_properties"],
                "fc": line.get_color(),
            }
            _, stack_tm_props = get_scaled_arrow_properties(
                Plotting["arrow_length"],
                tm_props,
                ax.get_gridspec()._row_height_ratios[a_i] * len(property_list) / 3.0,
            )
            if tlb_tm is not None:
                tm_x = 1.0 - (tlb_tm / Plotting["axis_rescale"])

            # Main axis arrows
            if a_i == 0:
                plot_merger_arrow(ax, mm_x, new_arrow_length, stack_mm_props, "upper")
                if tlb_tm is not None:
                    plot_merger_arrow(
                        ax, tm_x, new_arrow_length, stack_tm_props, "upper"
                    )
            elif a_i == len(axs) - 1:
                plot_merger_arrow(ax, mm_x, new_arrow_length, stack_mm_props, "lower")
                if tlb_tm is not None:
                    plot_merger_arrow(
                        ax, tm_x, new_arrow_length, stack_tm_props, "lower"
                    )

            ############################################################
            # Individual figures
            # First major merger
            (new_arrow_length, indiv_mm_props) = get_scaled_arrow_properties(
                Plotting["arrow_length"],
                mm_props,
                indiv_ax.get_gridspec()._row_height_ratios[0] / 3.0,
            )

            # Target major merger
            _, indiv_tm_props = get_scaled_arrow_properties(
                Plotting["arrow_length"],
                tm_props,
                indiv_ax.get_gridspec()._row_height_ratios[0] / 3.0,
            )

            # Individual axes and arrows
            for loc in ["upper", "lower"]:
                plot_merger_arrow(indiv_ax, mm_x, new_arrow_length, indiv_mm_props, loc)
                if tlb_tm is not None:
                    plot_merger_arrow(
                        indiv_ax, tm_x, new_arrow_length, indiv_tm_props, loc=loc
                    )
            ############################################################

    ####################################################################
    # Set common properties of the plots using the figure handler
    Plotting["figure_handler"].set_stacked_figure_properties(
        axs, ylabels, yscale=yscales, ylims=ylims
    )
    Plotting["figure_handler"].set_figure_properties(
        indiv_axs, ylabels, yscale=yscales, ylims=ylims
    )

    ####################################################################
    # Save figures
    save_figures(fig, fig6_out_file)

    for prop_name, fig in zip(property_list, indiv_figs):
        save_figures(fig, out_file_template.format(prop_name))

    ####################################################################
    # Plot TN and SM
    ####################################################################
    property_list = ["TN_model0", "SM_model0"]
    ylabels, yscales, ylims = return_plot_format_lists(property_list)
    birth_map = {"SM_model0": "SM_birth_model0"}

    ####################################################################
    # Create figures
    fig, axs = plt.subplots(
        1,
        len(property_list),
        figsize=(8 * np.ceil(len(property_list)), 8.0),
        gridspec_kw={
            "width_ratios": [1.0 / len(property_list) for _ in property_list],
            "height_ratios": [1],
            "hspace": 0,
        },
    )
    indiv_figs = [plt.figure(figsize=(8, 8)) for _ in property_list]
    indiv_axs = [fig.add_subplot(111) for fig in indiv_figs]

    ####################################################################
    # Iterate over each axis and plot data
    for a_i, (ax, indiv_ax, property_to_plot) in enumerate(
        zip(axs, indiv_axs, property_list)
    ):

        for sim_data, sim, sim_name, tlb_mm, tlb_tm in zip(
            ev_data,
            Simulations["sim_list"],
            Simulations["sim_names"],
            Simulations["sim_tlb_major_merger"],
            Simulations["sim_tlb_target_merger"],
        ):
            med, spread = sim_data.med_spread(property_to_plot)

            ############################################################
            # Plot lines and fills
            (line,) = ax.plot(
                sim_data.t_lb, med, label=sim_name, **Plotting["plot_styles"][sim]
            )
            ax.fill_between(
                sim_data.t_lb, *spread, lw=0, color=line.get_color(), alpha=0.3
            )

            (indiv_line,) = indiv_ax.plot(
                sim_data.t_lb, med, label=sim_name, **Plotting["plot_styles"][sim]
            )
            indiv_ax.fill_between(
                sim_data.t_lb, *spread, lw=0, color=indiv_line.get_color(), alpha=0.3
            )

            ############################################################
            # R200 lines
            if property_to_plot in birth_map:
                birth_med, _ = sim_data.med_spread(birth_map[property_to_plot])
                style = {**Plotting["plot_styles"][sim], "ls": "--"}
                ax.plot(sim_data.t_lb, birth_med, **style)
                indiv_ax.plot(sim_data.t_lb, birth_med, **style)

            ############################################################
            # Plot merger arrows on first and last panels
            # First major merger
            mm_x = 1.0 - (tlb_mm / Plotting["axis_rescale"])
            mm_props = {
                **Plotting["mm_arrow_properties"],
                "fc": line.get_color(),
            }
            (new_arrow_length, stack_mm_props) = get_scaled_arrow_properties(
                Plotting["arrow_length"],
                mm_props,
                ax.get_gridspec()._col_width_ratios[a_i],
            )

            # Target major merger
            tm_props = {
                **Plotting["tm_arrow_properties"],
                "fc": line.get_color(),
            }
            _, stack_tm_props = get_scaled_arrow_properties(
                Plotting["arrow_length"],
                tm_props,
                ax.get_gridspec()._col_width_ratios[a_i],
            )
            if tlb_tm is not None:
                tm_x = 1.0 - (tlb_tm / Plotting["axis_rescale"])

            # Main axis arrows
            plot_merger_arrow(ax, mm_x, new_arrow_length, stack_mm_props, "upper")
            plot_merger_arrow(ax, mm_x, new_arrow_length, stack_mm_props, "lower")
            if tlb_tm is not None:
                plot_merger_arrow(ax, tm_x, new_arrow_length, stack_tm_props, "upper")
                plot_merger_arrow(ax, tm_x, new_arrow_length, stack_tm_props, "lower")

            ############################################################
            # Individual figures
            # First major merger
            (new_arrow_length, indiv_mm_props) = get_scaled_arrow_properties(
                Plotting["arrow_length"],
                mm_props,
                indiv_ax.get_gridspec()._row_height_ratios[0] / 2.0,
            )

            # Target major merger
            _, indiv_tm_props = get_scaled_arrow_properties(
                Plotting["arrow_length"],
                tm_props,
                indiv_ax.get_gridspec()._row_height_ratios[0] / 2.0,
            )

            # Individual axes and arrows
            for loc in ["upper", "lower"]:
                plot_merger_arrow(indiv_ax, mm_x, new_arrow_length, indiv_mm_props, loc)
                if tlb_tm is not None:
                    plot_merger_arrow(
                        indiv_ax, tm_x, new_arrow_length, indiv_tm_props, loc=loc
                    )
            ############################################################

    ####################################################################
    # Set common properties of the plots using the figure handler
    Plotting["figure_handler"].set_stacked_figure_properties(
        axs, ylabels, yscale=yscales, ylims=ylims, direction="horizontal"
    )
    Plotting["figure_handler"].set_figure_properties(
        indiv_axs, ylabels, yscale=yscales, ylims=ylims
    )

    ####################################################################
    # Legend
    orig_legend = Plotting["figure_handler"].legend

    legend_markers = [
        plt.Line2D([0, 1], [0, 0], color="k", linestyle="--"),
        plt.Line2D([0, 1], [0, 0], color="k", linestyle="-"),
    ]
    legend_labels = ["Without dynamical evolution", "With dynamical evolution"]

    axs[1].legend(legend_markers, legend_labels, loc="upper right", alignment="right")
    axs[0].add_artist(orig_legend)

    ####################################################################
    # Save figures
    save_figures(fig, fig7_out_file)

    for prop_name, fig in zip(property_list, indiv_figs):
        save_figures(fig, out_file_template.format(prop_name))

    return None


if __name__ == "__main__":
    main()
