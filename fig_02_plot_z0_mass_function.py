#!/usr/bin/env python3

# Place import files below
import matplotlib.pyplot as plt
import numpy as np

from common_functions import save_figures
from process_data import Z0Data
from universal_settings import (cosmology_parameters, plot_styles, sim_list,
                                sim_names)


def main():
    # Plot settings
    try:
        plt.style.use('./paper.mplstyle')
    except OSError:
        pass
    # File location
    fig2_out_file = 'fig2_dn_dlogM.pdf'

    # Analysis settings
    m_bins = np.logspace(2, 8)  # Msun / h
    dlogM = np.log10(m_bins[1:]) - np.log10(m_bins[:-1])
    mid_logmbins = 10.**(np.log10(m_bins[:-1]) + dlogM / 2.)
    aperture = 30. * cosmology_parameters['h']  # kpc / h
    aperture_volume = (4. * np.pi * aperture**3.) / 3.  # h^-3 kpc^3
    CL = [16., 84.]

    # Load data for figures
    all_z0_data = [Z0Data(sim) for sim in sim_list[:3]]

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)

    current_med = []
    birth_med = []

    for sim, sim_name, z0_data in zip(sim_list, sim_names, all_z0_data):
        temp_fig = plt.figure(figsize=(8, 8))
        temp_ax = temp_fig.add_subplot(111)

        mcl_current = z0_data.current_cluster_m_current  # Msun
        mcl_birth = z0_data.current_cluster_m_birth  # Msun

        dn_dlogM_current = np.column_stack([
            np.histogram(mass_data, m_bins)[0] / dlogM / aperture_volume
            for mass_data in mcl_current.T
        ])
        dn_dlogM_birth = np.column_stack([
            np.histogram(mass_data, m_bins)[0] / dlogM / aperture_volume
            for mass_data in mcl_birth.T
        ])

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

        med_dnlogM_current, spread_dnlogM_current = med_spread(
            dn_dlogM_current, confidence=CL)
        med_dnlogM_birth, spread_dnlogM_birth = med_spread(dn_dlogM_birth,
                                                           confidence=CL)

        med_dnlogM_current[med_dnlogM_current == 0] = np.nan
        med_dnlogM_birth[med_dnlogM_birth == 0] = np.nan

        current_med.append(med_dnlogM_current)
        birth_med.append(med_dnlogM_birth)

        ################################################################
        # Plot surviving cluster mass function at z=0
        ################################################################
        # Plot median
        line, = ax.plot(mid_logmbins,
                        med_dnlogM_current,
                        label=sim_name,
                        **plot_styles[sim])
        # Plot scatter
        ax.fill_between(mid_logmbins,
                        *spread_dnlogM_current,
                        color=line.get_color(),
                        alpha=0.3)

        ################################################################
        # Plot birth mass function of surviving clusters
        ################################################################
        err_line_dict = {'color': line.get_color(), 'ls': ':'}

        ax.errorbar(mid_logmbins,
                    med_dnlogM_birth,
                    yerr=np.abs(spread_dnlogM_birth - med_dnlogM_birth),
                    **err_line_dict)
        ################################################################

    ax.set(
        xlabel=r'$M_{\rm cl}\, \left[{\rm M_\odot}\right]$',
        ylabel=
        r'$dn\, /\, d \log M_{\rm cl}\, \left[h^{-3} {\rm kpc^{-3}}\right]$',
        xscale='log',
        yscale='log')
    ax.legend()
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
    print(birth_med[0] / birth_med[1])
    print("Suppressed / Organic")
    print(birth_med[2] / birth_med[1])
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
