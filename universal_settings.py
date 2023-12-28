#!/usr/bin/env python3
import copy
import os

import numpy as np

from common_functions import Cosmology, FigureHandler

########################################################################
# File locations
data_dir = "data"
data_file_template = os.path.join(data_dir,
                                  "{}_gcmass_{:g}_paper_data_new.hdf5")
z0_data_file = os.path.join(data_dir, '{}_paper_z0_gc_properties.hdf5')

########################################################################
# Analysis settings
aperture = 30.  # ckpc
gc_mass = 1.e5  # Msun
cosmology_parameters = {
    'omega_M': 0.3156,
    'omega_lambda': 0.6844,
    'h': 0.6727,
}  # PLANCK 2015
cosmology_object = Cosmology(cosmology_parameters)

########################################################################
# Observations
caldwell_2011_m31_mstar_feh_data_file = os.path.join(
    data_dir, 'caldwell_2011_logMstar_FeH.csv')
johnson_2017_m31_dNdlogM_data_file = os.path.join(
    data_dir, 'johnson_2017_dNdlogM_M.csv')
caldwell_2009_m31_old_data_file = os.path.join(
    data_dir, 'caldwell_2009_old_clusters.csv')
caldwell_2009_m31_int_data_file = os.path.join(
    data_dir, 'caldwell_2009_intermediate_clusters.csv')
caldwell_2009_m31_young_data_file = os.path.join(
    data_dir, 'caldwell_2009_young_clusters.csv')
baumgardt_2019_mw_cluster_file = os.path.join(data_dir,
                                              'baumgardt_2019_mw_gcs.csv')

########################################################################
# Simulations
bastian_2020_sm_mstar_data_file = os.path.join(data_dir,
                                               'bastian2020_SM_Mstar_z0.csv')
bastian_2020_sminit_mstar_data_file = os.path.join(
    data_dir, 'bastian2020_SMinit_Mstar_z0.csv')
sim_list = [
    'z2_od_1p000_mz1p7_1p100_HiRes',
    'z2_od_1p000_HiRes',
    'z2_od_1p000_mz1p7_0p800_HiRes',
    'z2_od_1p000_mz1p7_1p100_HiRes_noBH',
    'z2_od_1p000_HiRes_noBH',
    'z2_od_1p000_mz1p7_0p800_HiRes_noBH',
]
# Onset of merger
sim_z_major_merger = [1.63, 1.7, 1.68, 1.63, 1.7, 1.68]
sim_tlb_major_merger = [9.85, 10.00, 9.95, 9.85, 10.00, 9.95]  # Gyr
sim_dtlb_mm = [0.325, 0.259, 0.129, 0.325, 0.259, 0.129]  # Gyr
sim_z_target_merger = [0.89, 0.77, None, 0.89, 0.77, None]
sim_tlb_target_merger = [7.49, 6.89, None, 7.49, 6.89, None]  # Gyr
sim_dtlb_tm = [0.486, 0.107, None, 0.486, 0.107, None]  # Gyr
sim_names = [
    r'$\textsc{enhanced}$',
    r'$\textsc{organic}$',
    r'$\textsc{suppressed}$',
    r'$\textsc{enhanced - noAGN}$',
    r'$\textsc{organic - noAGN}$',
    r'$\textsc{suppressed - noAGN}$',
]

########################################################################
# Plotting settings
tlb_lim = [0, 12.]
z_lim = [0., 8.]
fine_spacing = 0.2
z_cut = 6.5
coarse_spacing = 1.
all_redshifts = np.arange(*z_lim, fine_spacing)
major_redshifts = all_redshifts[~(all_redshifts % 1).astype(bool)]
minor_redshifts = all_redshifts[(all_redshifts % 1).astype(bool)]
minor_redshifts = np.concatenate((minor_redshifts[minor_redshifts < z_cut],
                                  np.arange(z_cut, z_lim[-1], coarse_spacing)))
arrow_length = 0.1
axis_rescale = np.abs(np.diff(tlb_lim)[0])
common_arrow_properties = {
    'head_width': arrow_length * 0.2,
    'head_length': arrow_length * 0.3,
    'width': arrow_length * 0.075,
    'shape': 'full',
    'length_includes_head': True,
    'zorder': 99
}
mm_arrow_properties = copy.deepcopy(common_arrow_properties)
mm_arrow_properties.update({'ls': ':', 'hatch': '/' * 5})
tm_arrow_properties = copy.deepcopy(common_arrow_properties)
tm_arrow_properties.update({'ls': '-'})

figure_handler = FigureHandler(major_redshifts=major_redshifts,
                               minor_redshifts=minor_redshifts,
                               tlb_lim=tlb_lim,
                               cos_obj=cosmology_object)
plot_styles = {
    'z2_od_1p000_mz1p7_1p100_HiRes': {
        'color': '#7E317B'
    },
    'z2_od_1p000_HiRes': {
        'color': '#34B83B'
    },
    'z2_od_1p000_mz1p7_0p800_HiRes': {
        'color': '#2F3BD6'
    },
    'z2_od_1p000_mz1p7_1p100_HiRes_noBH': {
        'color': '#7E317B',
        'marker': '.',
        'ls': '--'
    },
    'z2_od_1p000_HiRes_noBH': {
        'color': '#34B83B',
        'marker': '.',
        'ls': '--'
    },
    'z2_od_1p000_mz1p7_0p800_HiRes_noBH': {
        'color': '#2F3BD6',
        'marker': '.',
        'ls': '--'
    },
}
evo_property_dict = {
    'CFE': {
        'ylim': [0., 0.8],
        'yscale': 'linear',
        'ylabel': r'$\Gamma$',
        'printlabel': 'CFE(r < {} ckpc)'.format(aperture)
    },
    'GCDR': {
        'ylim': [-1.45, 0.2],
        'yscale': 'linear',
        'ylabel': r'${\rm GCDR\, \left[M_\odot\, yr^{-1}\right]}$',
        'printlabel': 'GCDR(r < {} ckpc) [Msun / yr]'.format(aperture)
    },
    'GCFR': {
        'ylim': [-0.2, 4.6],
        'yscale': 'linear',
        'ylabel': r'${\rm GCFR\, \left[M_\odot\, yr^{-1}\right]}$',
        'printlabel': 'GCFR(r < {} ckpc) [Msun / yr]'.format(aperture)
    },
    'M_200': {
        'ylim': [1.e11, None],
        'yscale': 'log',
        'ylabel': r'$M_{200}\, \left[{\rm M_\odot}\right]$',
        'printlabel': 'M_200 [Msun]'
    },
    'M_BH': {
        'ylim': [1.e5, 2.e8],
        'yscale': 'log',
        'ylabel': r'$M_{\rm SMBH}\, \left[{\rm M_\odot}\right]$',
        'printlabel': 'M_AGN [Msun]'
    },
    'M_GC': {
        'ylim': [5.e7, None],
        'yscale': 'log',
        'ylabel': r'$M_{\rm GC}\, \left[{\rm M_\odot}\right]$',
        'printlabel': 'M_GC(r < {} ckpc) [Msun]'.format(aperture)
    },
    'M_GC_r200': {
        'ylim': [5.e7, None],
        'yscale': 'log',
        'ylabel': r'$M_{\rm GC}\!\left(r < R_{200}\right)\, ' +
        r'\left[{\rm M_\odot}\right]$',
        'printlabel': 'M_GC(r < R_200) [Msun]'
    },
    'M_gas,SF': {
        'ylim': [8.e8, None],
        'yscale': 'log',
        'ylabel': r'$M_{\rm gas,\, SF}\, \left[{\rm M_\odot}\right]$',
        'printlabel': 'M_gas,SF(r < {} ckpc) [Msun]'.format(aperture)
    },
    'M_gas,SF_r200': {
        'ylim': [8.e8, None],
        'yscale': 'log',
        'ylabel': r'$M_{\rm gas,\, SF}\, \left[{\rm M_\odot}\right]$',
        'printlabel': 'M_gas,SF(r < R_200) [Msun]'
    },
    'M_star': {
        'ylim': [5.e9, 1.e11],
        'yscale': 'log',
        'ylabel': r'$M_\ast\, \left[{\rm M_\odot}\right]$',
        'printlabel': 'M_star(r < {} ckpc) [Msun]'.format(aperture)
    },
    'M_star_r200': {
        'ylim': [5.e9, 1.e12],
        'yscale': 'log',
        'ylabel': r'$M_\ast\!\left(r < R_{200}\right)\, ' +
        r'\left[{\rm M_\odot}\right]$',
        'printlabel': 'M_star(r < R_200) [Msun]'
    },
    'Mc_star': {
        'ylim': None,
        'yscale': 'log',
        'ylabel': r'$M_{\rm c,\, \ast} \left[{\rm M_\odot}\right]$',
        'printlabel': 'M_c,star [Msun]'
    },
    'Net_dM_GCdt': {
        'ylim': [-1., 3.5],
        'yscale': 'linear',
        'ylabel':
        r'$dM_{\rm GC}\, /\, dt\, \left[{\rm M_\odot\, yr^{-1}}\right]$',
        'printlabel': 'Net dM_GC/dt [Msun / yr]'
    },
    'Pk_birth,GC': {
        'ylim': [1.1e3, 8.e8],
        'yscale':
        'log',
        'ylabel':
        r'$\left(P\, /\, k_{\rm B}\right)_{\rm birth,\, GC}\!' +
        r'\left(r < {:.0f}\, {{\rm ckpc}}\right)\, '.format(aperture) +
        r'\left[{\rm K\, cm^{-3}}\right]$',
        'printlabel':
        'P_k,birth,GC [K cm^-3]'
    },
    'Pk_birth,star': {
        'ylim': [1.1e3, 8.e8],
        'yscale':
        'log',
        'ylabel':
        r'$\left(P\, /\, k_{\rm B}\right)_{\rm birth,\, \ast}\!' +
        r'\left(r < {:.0f}\, {{\rm ckpc}}\right)\, '.format(aperture) +
        r'\left[{\rm K\, cm^{-3}}\right]$',
        'printlabel':
        'P_k,birth,star [K cm^-3]'
    },
    'Pk_SFgas': {
        'ylim': [1.1e3, 8.e8],
        'yscale':
        'log',
        'ylabel':
        r'$\left(P\, /\, k_{\rm B}\right)_{\rm SF\, gas}\!' +
        r'\left(r < {:.0f}\, {{\rm ckpc}}\right)\, '.format(aperture) +
        r'\left[{\rm K\, cm^{-3}}\right]$',
        'printlabel':
        'P_k,SFgas [K cm^-3]'
    },
    'SFR': {
        'ylim': [-1.4, None],
        'yscale':
        'linear',
        'ylabel':
        r'${{\rm SFR}}\!\left(r < {:.0f}\, {{\rm ckpc}}\right)\, '.format(
            aperture) + r'\left[{\rm M_\odot\, yr^{-1}}\right]$',
        'printlabel':
        'SFR(r < {} ckpc) [Msun / yr]'.format(aperture)
    },
    'SM_model0': {
        'ylim': [0., 12.5],
        'yscale': 'linear',
        'ylabel': r'$S_{\rm M} = 100\, M_{\rm GC}\, /\, M_\ast$',
        'printlabel': 'S_M,0(r < {} ckpc)'.format(aperture)
    },
    'SM_model0_r200': {
        'ylim': [0., 12.5],
        'yscale': 'linear',
        'ylabel': r'$S_{\rm M} = 100\, M_{\rm GC}\, /\, M_\ast$',
        'printlabel': 'S_M,0(r < R_200)'
    },
    'SM_model3': {
        'ylim': [0., 12.5],
        'yscale': 'linear',
        'ylabel':
        r'$S_{\rm M}_{\rm no\, phys} = 100\, M_{\rm GC}\, /\, M_\ast$',
        'printlabel': 'S_M,3(r < {} ckpc)'.format(aperture)
    },
    'SM_birth_model0': {
        'ylim': [0., 13.],
        'yscale': 'linear',
        'ylabel':
        r'$S_{\rm M,\, birth} = 100\, M_{\rm GC,\, birth}\, /\, M_\ast$',
        'printlabel': 'S_M_birth,0(r < {} ckpc)'.format(aperture)
    },
    'TN_model0': {
        'ylim': [0., 2.25e2],
        'yscale': 'linear',
        'ylabel': r'$T_{\rm N} = N_{\rm GC}\, /\, M_\ast\, ' +
        r'\left[\left({\rm 10^9\, M_\odot}\right)^{-1}\right]$',
        'printlabel': 'T_N,0(r < {} ckpc) [10^-9 Msun^-1]'.format(aperture)
    },
    'TN_model0_r200': {
        'ylim': [0., 2.25e2],
        'yscale': 'linear',
        'ylabel': r'$T_{\rm N} = N_{\rm GC}\, /\, M_\ast\, ' +
        r'\left[\left({\rm 10^9\, M_\odot}\right)^{-1}\right]$',
        'printlabel': 'T_N,0(r < R_200) [10^-9 Msun^-1]'
    },
    'TN_model3': {
        'ylim': [0., 2.25e2],
        'yscale': 'linear',
        'ylabel': r'$T_{\rm N\, no\, phys} = N_{\rm GC}\, /\, M_\ast\, ' +
        r'\left[\left({\rm 10^9\, M_\odot}\right)^{-1}\right]$',
        'printlabel': 'T_N,3(r < {} ckpc) [10^-9 Msun^-1]'.format(aperture)
    },
    'SFR_MgasSF': {
        'ylim': np.asarray([0.3, 3]) * 1.e-9,
        'yscale': 'log',
        'ylabel': r'${\rm SFR}\, /\, M_{\rm gas,\, SF}\, ' +
        r'\left[{\rm yr^{-1}}\right]$',
        'printlabel': 'SFR/MgasSF(r < {} ckpc) [yr^-1]'.format(aperture)
    },
    'dMgc_dt_MgasSF': {
        'ylim': np.asarray([-2.5, 5]) * 1.e-10,
        'yscale': 'linear',
        'ylabel': r'$dM_{\rm GC}\, /\, dt\, /\, M_{\rm gas,\, SF}\, ' +
        r'\left[{\rm yr^{-1}}\right]$',
        'printlabel': 'dMgc_dt/MgasSF(r < {} ckpc) [yr^-1]'.format(aperture)
    },
}
