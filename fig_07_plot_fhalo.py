#!/usr/bin/env python3

# Place import files below
import matplotlib.pyplot as plt
import numpy as np

from common_functions import save_figures
from process_data import Z0Data
from universal_settings import (plot_styles, sim_list, sim_names)


def main():
    # Plot settings
    try:
        plt.style.use('./paper.mplstyle')
    except OSError:
        pass
    # File location
    fig7_out_file = 'fig7_mcfield_mfield.pdf'
    out_file_template = '{}_vs_tlb.pdf'

    # Analysis settings
    outer_radius = 50.  # kpc
    gc_mass = 1.e5  # Msun
    n_radial_bins = 21
    CL = [16., 84.]

    # Load data for figures
    all_z0_data = [Z0Data(sim) for sim in sim_list[:3]]

    fig, axs = plt.subplots(2,
                            1,
                            figsize=(8, 8 * (1. + 0.3 * (2 - 1.))),
                            sharex=True,
                            gridspec_kw={
                                "width_ratios": [1],
                                "height_ratios": [1. / 2] * 2,
                                "wspace": 0,
                                "hspace": 0
                            })

    property_list = ['fhalo_cl', 'fhalo_gc']
    indiv_figs = [plt.figure(figsize=(8, 8)) for _ in property_list]
    indiv_axs = [fig.add_subplot(111) for fig in indiv_figs]

    for sim, sim_name, z0_data in zip(sim_list, sim_names, all_z0_data):
        # Reina-Campos definition
        rc_cl, rc_gcs, rc_init_cl, rc_init_gcs = z0_data.reina_campos_f_halo(
            outer_radius=outer_radius, m_gc=gc_mass)

        # Radial calculation
        (radius_mid_bins, f_halo_cl, f_halo_gcs, f_halo_init_cl,
         f_halo_init_gcs) = z0_data.f_halo_vs_distance(
             n_bins=n_radial_bins, outer_radius=outer_radius, m_gc=gc_mass)
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
                            color=line.get_color(),
                            alpha=0.3)
        indiv_axs[0].fill_between(radius_mid_bins,
                                  *spread_fcl,
                                  color=indiv_line.get_color(),
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
        print("Clusters")
        print(med_rc_cl)
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
                            color=line.get_color(),
                            alpha=0.3)
        indiv_axs[1].fill_between(radius_mid_bins,
                                  *spread_fgc,
                                  color=indiv_line.get_color(),
                                  alpha=0.3)

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
        print("Globular Clusters")
        print(med_rc_gc)
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

    axs[0].minorticks_on()
    axs[1].minorticks_on()
    axs[0].legend()
    axs[0].set(ylabel=r'$M_{\rm field,\, CL}\, /\, M_{\rm field,\, tot}$',
               yscale='log',
               ylim=[None, 4.e-1])
    axs[1].set(xlabel=r'$r\, \left[{\rm kpc}\right]$',
               ylabel=r'$M_{\rm field,\, GC}\, /\, M_{\rm field,\, tot}$',
               yscale='log',
               ylim=[None, 8.e-2])
    indiv_axs[0].minorticks_on()
    indiv_axs[1].minorticks_on()
    indiv_axs[0].legend()
    indiv_axs[0].set(
        xlabel=r'$r\, \left[{\rm kpc}\right]$',
        ylabel=r'$M_{\rm field,\, CL}\, /\, M_{\rm field,\, tot}$',
        yscale='log',
        ylim=[None, 4.e-1])
    indiv_axs[1].set(
        xlabel=r'$r\, \left[{\rm kpc}\right]$',
        ylabel=r'$M_{\rm field,\, GC}\, /\, M_{\rm field,\, tot}$',
        yscale='log',
        ylim=[None, 8.e-2])

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
