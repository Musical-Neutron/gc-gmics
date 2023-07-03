#!/usr/bin/env python3

# Place import files below
import matplotlib.pyplot as plt
import numpy as np

from common_functions import save_figures
from process_data import Z0Data
from universal_settings import (
    caldwell_2011_m31_mstar_feh_data_file, johnson_2017_m31_dNdlogM_data_file,
    caldwell_2009_m31_old_data_file, caldwell_2009_m31_int_data_file,
    caldwell_2009_m31_young_data_file, cosmology_parameters, plot_styles,
    sim_list, sim_names)


def main():
    # Plot settings
    try:
        plt.style.use('./paper.mplstyle')
    except OSError:
        pass
    # File location
    fig2_out_file = 'fig2_dn_dlogM.pdf'
    caldwell_m31_data = np.genfromtxt(caldwell_2011_m31_mstar_feh_data_file,
                                      skip_header=True)
    johnson_m31_data = np.genfromtxt(johnson_2017_m31_dNdlogM_data_file,
                                     skip_header=True)
    caldwell_m31_old_data = np.genfromtxt(caldwell_2009_m31_old_data_file,
                                          skip_header=True)
    caldwell_m31_int_data = np.genfromtxt(caldwell_2009_m31_int_data_file,
                                          skip_header=True)
    caldwell_m31_yng_data = np.genfromtxt(caldwell_2009_m31_young_data_file,
                                          skip_header=True)
    johnson_dlogM = 0.1
    caldwell_dlogM = 0.2

    # Analysis settings
    n_bins = 30
    m_bins = np.logspace(2, 8, n_bins)  # Msun / h
    dlogM = np.log10(m_bins[1:]) - np.log10(m_bins[:-1])
    mid_logmbins = 10.**(np.log10(m_bins[:-1]) + dlogM / 2.)
    data_m_bins = np.logspace(2, 8, int(n_bins / 2))  # Msun
    data_dlogM = np.log10(data_m_bins[1:]) - np.log10(data_m_bins[:-1])
    mid_data_logmbins = 10.**(np.log10(data_m_bins[:-1]) + data_dlogM / 2.)
    # aperture = 30. * cosmology_parameters['h']  # kpc / h
    # aperture_volume = (4. * np.pi * aperture**3.) / 3.  # h^-3 kpc^3
    CL = [16., 84.]

    # Load data for figures
    all_z0_data = [Z0Data(sim) for sim in sim_list[:3]]
    m31_mass_data = 10.**caldwell_m31_data[:, 1]  # Msun
    dN_dlogM_m31 = np.histogram(m31_mass_data[~np.isnan(m31_mass_data)],
                                data_m_bins)[0] / data_dlogM
    dN_dlogM_m31_old = caldwell_m31_old_data[:, 1] / caldwell_dlogM
    dN_dlogM_m31_int = caldwell_m31_int_data[:, 1] / caldwell_dlogM
    dN_dlogM_m31_yng = caldwell_m31_yng_data[:, 1] / caldwell_dlogM
    dN_dlogM_m31_all = (caldwell_m31_old_data[:, 1] +
                        caldwell_m31_int_data[:, 1] +
                        caldwell_m31_yng_data[:, 1]) / caldwell_dlogM
    # dn_dlogM_m31 = np.histogram(m31_mass_data[~np.isnan(m31_mass_data)],
    #                             m_bins)[0] / dlogM  / aperture_volume

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)

    ################################################################
    # Plot surviving cluster mass function in M31
    ################################################################
    # Plot median
    ax.plot(10.**caldwell_m31_old_data[:, 0],
            dN_dlogM_m31_all,
            label=r'Caldwell et al. (2009)',
            color='c')
    ax.plot(mid_data_logmbins,
            dN_dlogM_m31,
            label='Caldwell et al. (2011)',
            color='k')
    ax.plot(10.**johnson_m31_data[:, 0],
            4 * johnson_m31_data[:, 1] / johnson_dlogM,
            label=r'$4 \times$ YCs (Johnson et al., 2017)',
            color='m')

    current_med = []
    birth_med = []
    all_birth_med = []

    for sim, sim_name, z0_data in zip(sim_list[:3], sim_names[:3],
                                      all_z0_data):
        temp_fig = plt.figure(figsize=(8, 8))
        temp_ax = temp_fig.add_subplot(111)

        mcl_current = z0_data.current_cluster_m_current  # Msun
        mcl_birth = z0_data.current_cluster_m_birth  # Msun
        all_mcl_birth = z0_data.all_cluster_m_birth  # Msun

        # prefixes = [
        #     'current_cluster_',
        #     'disrupted_',
        #     'allcluster_',
        # ]
        # all_mcl_birth = np.row_stack(
        #     [getattr(z0_data, prefix + 'm_birth') for prefix in prefixes])

        dN_dlogM_current = np.column_stack([
            np.histogram(mass_data, m_bins)[0] / dlogM
            for mass_data in mcl_current.T
        ])
        # dn_dlogM_current = np.column_stack([
        #     np.histogram(mass_data, m_bins)[0] / dlogM / aperture_volume
        #     for mass_data in mcl_current.T
        # ])
        dN_dlogM_birth = np.column_stack([
            np.histogram(mass_data, m_bins)[0] / dlogM
            for mass_data in mcl_birth.T
        ])
        dallN_dlogM_birth = np.column_stack([
            np.histogram(mass_data, m_bins)[0] / dlogM
            for mass_data in all_mcl_birth.T
        ])
        # dn_dlogM_birth = np.column_stack([
        #     np.histogram(mass_data, m_bins)[0] / dlogM / aperture_volume
        #     for mass_data in mcl_birth.T
        # ])

        # z0_mass_pdfs = [
        #     np.histogram(mass_data[selection],
        #                  int(2. * selection.sum()**(1. / 3.)),
        #                  density=True)
        #     for mass_data, selection in zip(mcl_current.T, (
        #         mcl_birth <= 2.e4).T)
        # ]
        for current_mass_data, birth_mass_data, selection in zip(
                mcl_current.T, mcl_birth.T, (mcl_birth <= 2.e4).T):
            temp_ax.hist(
                current_mass_data[selection] / birth_mass_data[selection],
                np.logspace(start=np.log10(0.005), stop=0, num=15),
                # np.linspace(start=0, stop=1, num=15),
                density=True,
                histtype='step',
                #  log=True
            )

        temp_ax.set(
            xlabel=r'$M_{\rm cl,\, current}\, /\, M_{\rm cl,\, birth}$',
            ylabel=r'PDF',
            xscale='log')
        temp_ax.tick_params(axis='x', which='major', pad=7)
        temp_ax.minorticks_on()
        temp_ax.set(title=sim_name +
                    r' $(M_{\rm cl,\, birth} \leq 2\times10^4)$')

        med_dNlogM_current, spread_dNlogM_current = med_spread(
            dN_dlogM_current, confidence=CL)
        # med_dnlogM_current, spread_dnlogM_current = med_spread(
        #     dN_dlogM_current, confidence=CL)
        med_dNlogM_birth, spread_dNlogM_birth = med_spread(dN_dlogM_birth,
                                                           confidence=CL)
        # med_dnlogM_birth, spread_dnlogM_birth = med_spread(dN_dlogM_birth,
        #                                                    confidence=CL)
        med_dallNlogM_birth, spread_dallNlogM_birth = med_spread(
            dallN_dlogM_birth, confidence=CL)

        med_dNlogM_current[med_dNlogM_current == 0] = np.nan
        med_dNlogM_birth[med_dNlogM_birth == 0] = np.nan
        med_dallNlogM_birth[med_dallNlogM_birth == 0] = np.nan

        current_med.append(med_dNlogM_current)
        birth_med.append(med_dNlogM_birth)
        all_birth_med.append(med_dallNlogM_birth)

        ################################################################
        # Plot surviving cluster mass function at z=0
        ################################################################
        # Plot median
        line, = ax.plot(mid_logmbins,
                        med_dNlogM_current,
                        label=sim_name,
                        **plot_styles[sim])
        # Plot scatter
        ax.fill_between(mid_logmbins,
                        *spread_dNlogM_current,
                        color=line.get_color(),
                        alpha=0.3)

        # ################################################################
        # # Plot birth mass function of surviving clusters
        # ################################################################
        # err_line_dict = {'color': line.get_color(), 'ls': ':'}

        # ax.errorbar(mid_logmbins,
        #             med_dNlogM_birth,
        #             yerr=np.abs(spread_dNlogM_birth - med_dNlogM_birth),
        #             **err_line_dict)
        # ################################################################

        ################################################################
        # Plot birth mass function of all clusters
        ################################################################
        err_line_dict = {'color': line.get_color(), 'ls': '--'}

        ax.errorbar(mid_logmbins,
                    med_dallNlogM_birth,
                    yerr=np.abs(spread_dallNlogM_birth - med_dallNlogM_birth),
                    **err_line_dict)
        ################################################################

    ax.set(xlabel=r'$M_{\rm cl}\, \left[{\rm M_\odot}\right]$',
           ylabel=r'$dN\, /\, d \log M_{\rm cl}$',
           xscale='log',
           yscale='log',
           xlim=np.array([10.**(3 * 0.95), 5.e7]))
    legend = ax.legend(
        markerfirst=False,
        loc='upper right',
        handlelength=0.,
        handletextpad=0.,
        labelspacing=0.,
    )
    for t_item, line in zip(legend.get_texts(), legend.get_lines()):
        t_item.set_color(line.get_color())
    ax.tick_params(axis='x', which='major', pad=7)
    ax.minorticks_on()

    # Save figures
    save_figures(fig, fig2_out_file)

    ####################################################################
    # DELETE WHEN DONE
    ####################################################################
    print("#" * 75)
    print("Current cluster mass")
    print("#" * 75)
    print("Enhanced / Organic")
    print(current_med[0] / current_med[1])
    print("Suppressed / Organic")
    print(current_med[2] / current_med[1])

    print()
    print("#" * 75)
    print("Birth cluster mass")
    print("#" * 75)
    print("Enhanced / Organic")
    print(all_birth_med[0] / all_birth_med[1])
    print("Suppressed / Organic")
    print(all_birth_med[2] / all_birth_med[1])
    # print("Enhanced / Organic")
    # print(birth_med[0] / birth_med[1])
    # print("Suppressed / Organic")
    # print(birth_med[2] / birth_med[1])
    ####################################################################

    plt.show()

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
