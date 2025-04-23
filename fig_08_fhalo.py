#!/usr/bin/env python3

# Place import files below
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patches

from common_functions import save_figures
from process_data import Z0Data
from universal_settings import horta_2021_zetagc_file, plot_styles, sim_list, sim_names


def main():
    # Plot settings
    try:
        plt.style.use('./paper.mplstyle')
    except OSError:
        pass
    # File location
    fig7_out_file = 'fig08_mcfield_mfield.pdf'
    out_file_template = '{}_vs_tlb.pdf'
    horta_mw_data = np.genfromtxt(horta_2021_zetagc_file, skip_header=True)

    # Analysis settings
    outer_radius = 50.  # kpc
    gc_mass = 1.e5  # Msun
    n_bins = 16
    CL = [16., 84.]

    # Load data for figures
    all_z0_data = [Z0Data(sim) for sim in sim_list[:3]]

    property_list = ['fhalo_cl', 'fhalo_gc']

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

    plt.rc('hatch', color='k', linewidth=0.5)

    # Select only Enhanced and Suppressed simulations for plotting
    selected_sim_list = np.asarray(sim_list)[[0, 2]]
    selected_sim_names = np.asarray(sim_names)[[0, 2]]
    selected_z0_data = np.asarray(all_z0_data)[[0, 2]]
    # selected_sim_list = np.asarray(sim_list)[[0, 1, 2]]
    # selected_sim_names = np.asarray(sim_names)[[0, 1, 2]]
    # selected_z0_data = np.asarray(all_z0_data)[[0, 1, 2]]

    # for sim, sim_name, z0_data in zip(sim_list, sim_names, all_z0_data):
    for sim, sim_name, z0_data in zip(selected_sim_list, selected_sim_names,
                                      selected_z0_data):
        # Reina-Campos definition
        rc_cl, rc_gcs, rc_init_cl, rc_init_gcs = z0_data.reina_campos_f_halo(
            outer_radius=outer_radius, m_gc=gc_mass)

        # Radial calculation
        (radius_mid_bins, f_halo_cl, f_halo_gcs, f_halo_init_cl,
         f_halo_init_gcs) = z0_data.f_halo_vs_distance(
             n_bins=n_bins, outer_radius=outer_radius, m_gc=gc_mass)
        max_rhalf = z0_data.max_rhalf

        # FHalo in radial bins
        med_fcl, spread_fcl = med_spread(f_halo_cl, confidence=CL)
        med_fgc, spread_fgc = med_spread(f_halo_gcs, confidence=CL)
        (med_fcl_init, spread_fcl_init) = med_spread(f_halo_init_cl,
                                                     confidence=CL)
        (med_fgc_init, spread_fgc_init) = med_spread(f_halo_init_gcs,
                                                     confidence=CL)

        # FHalo (Reina-Campos definition) data
        # Radius at which to plot Reina-Campos-style data
        rc_x = (max_rhalf + outer_radius) / 2.
        # y-values
        med_rc_cl, spread_rc_cl = med_spread(rc_cl, confidence=CL)
        med_rc_gc, spread_rc_gc = med_spread(rc_gcs, confidence=CL)
        med_rc_init_cl, spread_rc_init_cl = med_spread(rc_init_cl,
                                                       confidence=CL)
        med_rc_init_gc, spread_rc_init_gc = med_spread(rc_init_gcs,
                                                       confidence=CL)

        # Plot median
        line, = axs[0].plot(radius_mid_bins,
                            med_fcl,
                            label=sim_name,
                            **plot_styles[sim])
        indiv_line, = indiv_axs[0].plot(radius_mid_bins,
                                        med_fcl,
                                        label=sim_name,
                                        **plot_styles[sim])
        # Plot scatter
        axs[0].fill_between(radius_mid_bins,
                            *spread_fcl,
                            hatch='xxx',
                            color='None',
                            edgecolor=indiv_line.get_color(),
                            linewidth=0.5,
                            alpha=0.3)
        indiv_axs[0].fill_between(radius_mid_bins,
                                  *spread_fcl,
                                  hatch='xxx',
                                  color='None',
                                  edgecolor=indiv_line.get_color(),
                                  linewidth=0.5,
                                  alpha=0.3)

        # Plot evolved initial mass
        # Median
        axs[0].plot(radius_mid_bins,
                    med_fcl_init,
                    ls='--',
                    color=line.get_color())
        indiv_axs[0].plot(radius_mid_bins,
                          med_fcl_init,
                          ls='--',
                          color=indiv_line.get_color())
        # # Scatter
        # axs[0].fill_between(radius_mid_bins,
        #                     *spread_fcl_init,
        #                     color=line.get_color(),
        #                     alpha=0.3)

        # Plot FHalo (Reina-Campos definition)

        print(sim)
        print("   Clusters")
        print("   ", med_rc_cl)
        axs[0].axhline(med_rc_cl, ls=':', color=line.get_color(), zorder=1)
        axs[0].axhspan(*spread_rc_cl,
                       color=line.get_color(),
                       zorder=0,
                       ec=None,
                       alpha=0.2)
        indiv_axs[0].axhline(med_rc_cl,
                             ls=':',
                             color=indiv_line.get_color(),
                             zorder=1)
        indiv_axs[0].axhspan(*spread_rc_cl,
                             color=indiv_line.get_color(),
                             zorder=0,
                             ec=None,
                             alpha=0.2)
        # axs[0].errorbar(
        #     rc_x,
        #     med_rc_init_cl,
        #     yerr=[[err] for err in np.abs(spread_rc_init_cl - med_rc_init_cl)],
        #     marker='.',
        #     color=line.get_color(),
        #     markersize=5,
        #     lw=None,
        #     elinewidth=1,
        #     capsize=3,
        #     capthick=1,
        #     label=r'$F^{\rm halo}$',
        #     zorder=0)
        # axs[0].errorbar(rc_x,
        #                 med_rc_cl,
        #                 yerr=[[err]
        #                       for err in np.abs(spread_rc_cl - med_rc_cl)],
        #                 marker='.',
        #                 color=line.get_color(),
        #                 markersize=5,
        #                 lw=None,
        #                 elinewidth=1,
        #                 capsize=3,
        #                 capthick=1,
        #                 label=r'$F^{\rm halo}$',
        #                 zorder=0)

        # Plot median
        line, = axs[1].plot(radius_mid_bins,
                            med_fgc,
                            label=sim_name,
                            **plot_styles[sim])
        indiv_line, = indiv_axs[1].plot(radius_mid_bins,
                                        med_fgc,
                                        label=sim_name,
                                        **plot_styles[sim])
        # Plot scatter
        axs[1].fill_between(radius_mid_bins,
                            *spread_fgc,
                            hatch='xxx',
                            color='None',
                            edgecolor=indiv_line.get_color(),
                            linewidth=0.5,
                            alpha=0.3)
        indiv_axs[1].fill_between(radius_mid_bins,
                                  *spread_fgc,
                                  hatch='xxx',
                                  color='None',
                                  edgecolor=indiv_line.get_color(),
                                  linewidth=0.1,
                                  alpha=0.5)

        # Plot evolved initial mass
        # Median
        axs[1].plot(radius_mid_bins,
                    med_fgc_init,
                    ls='--',
                    color=line.get_color())
        indiv_axs[1].plot(radius_mid_bins,
                          med_fgc_init,
                          ls='--',
                          color=indiv_line.get_color())
        # # Scatter
        # axs[1].fill_between(radius_mid_bins,
        #                     *spread_fgc_init,
        #                     color=line.get_color(),
        #                     alpha=0.3)

        # Plot FHalo (Reina-Campos definition)
        print("   Globular Clusters")
        print("   ", med_rc_gc)
        axs[1].axhline(med_rc_gc, ls=':', color=line.get_color(), zorder=1)
        axs[1].axhspan(*spread_rc_gc,
                       color=line.get_color(),
                       zorder=0,
                       ec=None,
                       alpha=0.2)
        indiv_axs[1].axhline(med_rc_gc,
                             ls=':',
                             color=indiv_line.get_color(),
                             zorder=1)
        indiv_axs[1].axhspan(*spread_rc_gc,
                             color=indiv_line.get_color(),
                             zorder=0,
                             ec=None,
                             alpha=0.2)
        # axs[1].errorbar(
        #     rc_x,
        #     med_rc_init_gc,
        #     yerr=[[err] for err in np.abs(spread_rc_init_gc - med_rc_init_gc)],
        #     marker='.',
        #     color=line.get_color(),
        #     markersize=5,
        #     lw=None,
        #     elinewidth=1,
        #     capsize=3,
        #     capthick=1,
        #     label=r'$F^{\rm halo}$',
        #     zorder=0)
        # axs[1].errorbar(rc_x,
        #                 med_rc_gc,
        #                 yerr=[[err]
        #                       for err in np.abs(spread_rc_gc - med_rc_gc)],
        #                 marker='.',
        #                 color=line.get_color(),
        #                 markersize=5,
        #                 lw=None,
        #                 elinewidth=1,
        #                 capsize=3,
        #                 capthick=1,
        #                 label=r'$F^{\rm halo}$',
        #                 zorder=0)

    ####################################################################
    # Plot Horta+2021 MW data
    ####################################################################
    for ax, indiv_ax in zip(axs, indiv_axs):
        # Plot median
        line, = ax.plot(horta_mw_data[:, 0],
                        1.5 * horta_mw_data[:, 1] * 0.01,
                        color='grey',
                        marker='*',
                        markerfacecolor='yellow',
                        markeredgecolor='yellow',
                        markevery=8)
        indiv_line, = indiv_ax.plot(horta_mw_data[:, 0],
                                    1.5 * horta_mw_data[:, 1] * 0.01,
                                    color='grey',
                                    marker='*',
                                    markerfacecolor='yellow',
                                    markeredgecolor='yellow',
                                    markevery=8)
        # Plot scatter
        ax.fill_between(horta_mw_data[:, 0],
                        *horta_mw_data[:, 2:].T * 1.5 * 0.01,
                        color=line.get_color(),
                        edgecolor=line.get_color(),
                        linewidth=0.5,
                        alpha=0.3)
        indiv_ax.fill_between(horta_mw_data[:, 0],
                              *horta_mw_data[:, 2:].T * 1.5 * 0.01,
                              color=indiv_line.get_color(),
                              edgecolor=indiv_line.get_color(),
                              linewidth=0.5,
                              alpha=0.3)
    ####################################################################

    legend_markers = [
        plt.Line2D([0, 1], [0, 0], color='k', linestyle='--'),
        plt.Line2D([0, 1], [0, 0], color='k', linestyle='-'),
        plt.Line2D([0, 1], [0, 0], color='k', linestyle=':'),
        (patches.Patch(color='k', alpha=0.3, linewidth=0),
         plt.Line2D([0, 1], [0, 0],
                    color='grey',
                    linestyle='-',
                    marker='*',
                    markerfacecolor='yellow',
                    markeredgecolor='yellow')),
    ]
    legend_labels = [
        'Fully disrupted',
        'Partially disrupted',
        r'$F^{\rm halo}$',
        'Horta+(2021b)',
    ]

    for ax in axs:
        ax.minorticks_on()
        ax.set(yscale='log')
    axs[0].set(
        ylabel=
        r'$\zeta_{\rm CL} = M^{\rm halo}_{\rm CL}\, /\, M^{\rm halo}_{\rm tot}$',
        ylim=[None, 2.5e-1])
    axs[1].set(
        xlabel=r'$r\, \left[{\rm kpc}\right]$',
        ylabel=
        r'$\zeta_{\rm GC} = M^{\rm halo}_{\rm GC}\, /\, M^{\rm halo}_{\rm tot}$',
        ylim=[None, 8.e-2])
    orig_legend = axs[0].legend(markerfirst=False,
                                loc='upper right',
                                handlelength=0.,
                                handletextpad=0.,
                                labelspacing=0.)
    for t_item, line in zip(orig_legend.get_texts(), orig_legend.get_lines()):
        t_item.set_color(line.get_color())
    axs[0].legend(legend_markers,
                  legend_labels,
                  bbox_to_anchor=(0., 1.02, 1., .102),
                  loc='lower left',
                  ncols=2,
                  mode="expand",
                  borderaxespad=0.,
                  frameon=True,
                  fancybox=True,
                  framealpha=0.75)
    # .legend(legend_markers, legend_labels, loc='upper right', markerfirst=False)
    axs[0].add_artist(orig_legend)

    for indiv_ax in indiv_axs:
        indiv_ax.minorticks_on()
        ax0_orig_legend = indiv_ax.legend(markerfirst=True,
                                          loc='upper right',
                                          handlelength=0.,
                                          handletextpad=0.,
                                          labelspacing=0.)
        for t_item, line in zip(ax0_orig_legend.get_texts(),
                                ax0_orig_legend.get_lines()):
            t_item.set_color(line.get_color())
        indiv_ax.legend(legend_markers,
                        legend_labels,
                        bbox_to_anchor=(0., 1.02, 1., .102),
                        loc='lower left',
                        ncols=2,
                        mode="expand",
                        borderaxespad=0.,
                        frameon=True,
                        fancybox=True,
                        framealpha=0.75)
        # .legend(legend_markers,legend_labels,loc='upper right',markerfirst=False)
        indiv_ax.add_artist(ax0_orig_legend)
        indiv_ax.set(xlabel=r'$r\, \left[{\rm kpc}\right]$', yscale='log')
    indiv_axs[0].set(
        ylabel=
        r'$\zeta_{\rm CL} = M_{\rm field,\, CL}\, /\, M_{\rm field,\, tot}$',
        ylim=[None, 4.e-1])
    indiv_axs[1].set(
        ylabel=
        r'$\zeta_{\rm GC} = M_{\rm field,\, GC}\, /\, M_{\rm field,\, tot}$',
        ylim=[None, 2.e-1])

    # plt.show()

    # Save figures
    save_figures(fig, fig7_out_file)

    for prop_name, fig in zip(property_list, indiv_figs):
        save_figures(fig, out_file_template.format(prop_name))

    return None


def med_spread(data, confidence=[16., 84.], axis=None):
    if axis is not None:
        axis = axis
    else:
        axis = len(data.shape) - 1

    med = np.nanmedian(data, axis=axis)
    percentiles = np.nanpercentile(data, confidence, axis=axis)
    return med, percentiles


if __name__ == "__main__":

    main()
