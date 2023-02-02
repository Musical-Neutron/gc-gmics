#!/usr/bin/env python3
import copy
import os

import numpy as np

from common_functions import Cosmology, FigureHandler

########################################################################
# File locations
data_dir = "data"
data_file_template = os.path.join(data_dir, "{}_gcmass_{:g}_paper_data.hdf5")
z0_data_file = os.path.join(data_dir, '{}_paper_z0_gc_properties.hdf5')

########################################################################
# Analysis settings
aperture = 30.  # kpc
gc_mass = 1.e5  # Msun
cosmology_parameters = {
    'omega_M': 0.3156,
    'omega_lambda': 0.6844,
    'h': 0.6727,
}  # PLANCK 2015
cosmology_object = Cosmology(cosmology_parameters)

########################################################################
# Simulations
sim_list = [
    'z2_od_1p000_mz1p7_1p100_HiRes',
    'z2_od_1p000_HiRes',
    'z2_od_1p000_mz1p7_0p800_HiRes',
    'z2_od_1p000_mz1p7_1p100_HiRes_noBH',
    'z2_od_1p000_HiRes_noBH',
    'z2_od_1p000_mz1p7_0p800_HiRes_noBH',
]
sim_tlb_major_merger = [9.3, 10., 9.6, 9.3, 10., 9.6]  # Gyr
sim_tlb_target_merger = [8.07, 6.71, None, 8.07, 6.71, None]  # Gyr
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
all_redshifts = np.arange(0, 8.2, 0.2)
major_redshifts = all_redshifts[~(all_redshifts % 1).astype(bool)]
minor_redshifts = all_redshifts[(all_redshifts % 1).astype(bool)]
minor_redshifts = np.concatenate((minor_redshifts[:24], np.arange(6.5, 8.,
                                                                  1.)))
arrow_length = 0.1
axis_rescale = np.abs(np.diff(tlb_lim)[0])
common_arrow_properties = {
    'head_width': arrow_length * 0.2,
    'head_length': arrow_length * 0.3,
    'shape': 'full',
    'length_includes_head': True,
    'zorder': 99
}
mm_arrow_properties = copy.deepcopy(common_arrow_properties)
mm_arrow_properties.update({'ls': ':'})
tm_arrow_properties = copy.deepcopy(common_arrow_properties)
tm_arrow_properties.update({'ls': '-'})

figure_handler = FigureHandler(major_redshifts=major_redshifts,
                               minor_redshifts=minor_redshifts,
                               tlb_lim=tlb_lim,
                               cos_obj=cosmology_object)
plot_styles = {
    'z2_od_1p000_mz1p7_1p100_HiRes': {
        'color': '#2F3BD6'
    },
    'z2_od_1p000_HiRes': {
        'color': '#34B83B'
    },
    'z2_od_1p000_mz1p7_0p800_HiRes': {
        'color': '#7E317B'
    },
    'z2_od_1p000_mz1p7_1p100_HiRes_noBH': {
        'color': '#2F3BD6',
        'marker': '.',
        'ls': '--'
    },
    'z2_od_1p000_HiRes_noBH': {
        'color': '#34B83B',
        'marker': '.',
        'ls': '--'
    },
    'z2_od_1p000_mz1p7_0p800_HiRes_noBH': {
        'color': '#7E317B',
        'marker': '.',
        'ls': '--'
    },
}
evo_property_dict = {
    'CFE': {
        'ylim': [0., 0.75],
        'yscale': 'linear',
        'ylabel': r'$\Gamma\!\left(r < 30\, {\rm kpc}\right)$',
        'printlabel': 'CFE(r < 30 kpc)'
    },
    'GCDR': {
        'ylim': [-1.45, 0.3],
        'yscale': 'linear',
        'ylabel': r'${\rm GCDR\, \left[M_\odot\, yr^{-1}\right]}$',
        'printlabel': 'GCDR(r < 30 kpc) [Msun / yr]'
    },
    'GCFR': {
        'ylim': [-0.2, 4.6],
        'yscale': 'linear',
        'ylabel': r'${\rm GCFR\, \left[M_\odot\, yr^{-1}\right]}$',
        'printlabel': 'GCFR(r < 30 kpc) [Msun / yr]'
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
        'ylabel': r'$M_{\rm GC}\!\left(r < 30\, {\rm kpc}\right)\, ' +
        r'\left[{\rm M_\odot}\right]$',
        'printlabel': 'M_GC(r < 30 kpc) [Msun]'
    },
    'M_gas,SF': {
        'ylim': [8.e8, None],
        'yscale': 'log',
        'ylabel': r'$M_{\rm gas,\, SF}\, \left[{\rm M_\odot}\right]$',
        'printlabel': 'M_gas,SF(r < 30 kpc) [Msun]'
    },
    'M_star': {
        'ylim': [5.e9, 1.e11],
        'yscale': 'log',
        'ylabel': r'$M_\ast\!\left(r < 30\, {\rm kpc}\right)\, ' +
        r'\left[{\rm M_\odot}\right]$',
        'printlabel': 'M_star(r < 30 kpc) [Msun]'
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
        'ylim':
        None,
        'yscale':
        'log',
        'ylabel':
        r'$\left(P\, /\, k_{\rm B}\right)_{\rm birth,\, GC}\!' +
        r'\left(r < 30\, {\rm kpc}\right)\, ' +
        r'\left[{\rm K\, cm^{-3}}\right]$',
        'printlabel':
        'P_k,birth,GC [K cm^-3]'
    },
    'Pk_birth,star': {
        'ylim':
        None,
        'yscale':
        'log',
        'ylabel':
        r'$\left(P\, /\, k_{\rm B}\right)_{\rm birth,\, \ast}\!' +
        r'\left(r < 30\, {\rm kpc}\right)\, ' +
        r'\left[{\rm K\, cm^{-3}}\right]$',
        'printlabel':
        'P_k,birth,star [K cm^-3]'
    },
    'SFR': {
        'ylim': [-1.4, None],
        'yscale': 'linear',
        'ylabel': r'${\rm SFR}\!\left(r < 30\, {\rm kpc}\right)\, ' +
        r'\left[{\rm M_\odot\, yr^{-1}}\right]$',
        'printlabel': 'SFR(r < 30 kpc) [Msun / yr]'
    },
    'SM_model0': {
        'ylim': [0., 10.],
        'yscale': 'linear',
        'ylabel': r'$S_{\rm M}\!\left(r < 30\, {\rm kpc}\right)' +
        r'= 100\, M_{\rm GC}\, /\, M_\ast$',
        'printlabel': 'S_M,0(r < 30 kpc)'
    },
    'SM_model3': {
        'ylim': [0., 10.],
        'yscale': 'linear',
        'ylabel':
        r'$S_{\rm M}_{\rm no\, phys}\!\left(r < 30\, {\rm kpc}\right)' +
        r'= 100\, M_{\rm GC}\, /\, M_\ast$',
        'printlabel': 'S_M,3(r < 30 kpc)'
    },
    'TN_model0': {
        'ylim': [0., 2.e2],
        'yscale':
        'linear',
        'ylabel':
        r'$T_{\rm N}\!\left(r < 30\, {\rm kpc}\right)' +
        r'= N_{\rm GC}\, /\, M_\ast\, ' +
        r'\left[\left({\rm 10^9\, M_\odot}\right)^{-1}\right]$',
        'printlabel':
        'T_N,0(r < 30 kpc) [10^-9 Msun^-1]'
    },
    'TN_model3': {
        'ylim': [0., 2.e2],
        'yscale':
        'linear',
        'ylabel':
        r'$T_{\rm N\, no\, phys}\!\left(r < 30\, {\rm kpc}\right)' +
        r'= N_{\rm GC}\, /\, M_\ast\, ' +
        r'\left[\left({\rm 10^9\, M_\odot}\right)^{-1}\right]$',
        'printlabel':
        'T_N,3(r < 30 kpc) [10^-9 Msun^-1]'
    },
}
