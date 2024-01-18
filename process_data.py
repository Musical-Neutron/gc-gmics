#!/usr/bin/env python3

# Place import files below
import h5py
import numpy as np

from universal_settings import (
    data_file_template,
    evo_property_dict,
    gc_mass,
    z0_data_file,
)


def return_plot_format_lists(properties_to_plot,
                             property_dict=evo_property_dict):
    """Returns plotting properties associated with a given attribute.

    Args:
        properties_to_plot (list): List of str designating the desired
            properties to plot.
        property_dict (dict, optional): Dictionary of properties
            associated with each attribute.
            Defaults to evo_property_dict.

    Returns:
        tuple:  [0]: y-axis labels.
                [1]: y-axis scale (usually 'linear' or 'log').
                [2]: y-axis limits.
    """

    ylabels = []
    yscales = []
    ylims = []
    for sub in properties_to_plot:
        ylabels.append(property_dict[sub]['ylabel'])
        yscales.append(property_dict[sub]['yscale'])
        ylims.append(property_dict[sub]['ylim'])
    return ylabels, yscales, ylims


def return_print_format_lists(properties_to_print,
                              property_dict=evo_property_dict):
    """Returns stdout print properties associated with a given attribute

    Args:
        properties_to_print (list): List of str designating the desired
            properties to print.
        property_dict (dict, optional): Dictionary of properties
            associated with each attribute.
            Defaults to evo_property_dict.

    Returns:
        list: Print labels.
    """

    labels = []
    for sub in properties_to_print:
        labels.append(property_dict[sub]['printlabel'])
    return labels


def quantity_in_shells(quantity_to_sum, distances, shell_r):
    """Calculates the amount of quantity_to_sum in radial shells.

    Args:
        quantity_to_sum (len(N) arr): Array of values to sum.
        distances (len(N) arr): Distances associated with each value in
            quantity_to_sum
        shell_r (len(M) arr): Array of bin edges.

    Returns:
        len(M)-1 arr: Total amount of quantity_to_sum in each radial
            shell.
    """
    q_enc = np.asarray([
        np.nansum(quantity_to_sum * (distances <= r), axis=0) for r in shell_r
    ])
    q_shell = q_enc[1:] - q_enc[:-1]
    return q_shell


class Z0Data(object):

    def __init__(self, sim) -> None:
        self.sim = sim
        self.process_data()
        pass

    def _iter_read_file_data(self, f, group_name, prefix):
        for key in f[group_name].keys():
            data = f[group_name][key][()]
            str_key = prefix + key.strip().lower().replace(' ', '_')
            setattr(self, str_key, data)
        return None

    def process_data(self):
        with h5py.File(z0_data_file.format(self.sim), 'r') as f:
            groups = [
                'All clusters', 'Current clusters', 'Fully disrupted',
                'Partially disrupted', 'Field', 'Undisrupted'
            ]
            group_prefixes = [
                'all_cluster_', 'current_cluster_', 'disrupted_',
                'part_disrupted_', 'field_', 'undisrupted_'
            ]
            # Read in current cluster data
            for group, prefix in zip(groups, group_prefixes):
                self._iter_read_file_data(f, group, prefix)

            self.r_half = f['R_half'][()]
        return None

    def f_halo_vs_distance(self,
                           n_bins=40,
                           outer_radius=50.,
                           m_gc=1.e5,
                           metallicity_range=None):
        # Identify largest stellar half-mass radius
        self.max_rhalf = np.nanmax(self.r_half)

        # Construct radial bins
        radius_bins = np.linspace(self.max_rhalf, outer_radius, n_bins)
        radius_mid_bins = (radius_bins[1:] + radius_bins[:-1]) / 2.

        # Select objects in the given metallicity range
        if metallicity_range is not None:
            select_field_metallicity = (
                (self.field_fe_h >= np.nanmin(metallicity_range)) *
                (self.field_fe_h < np.nanmax(metallicity_range)))
        else:
            select_field_metallicity = np.ones(self.field_distance.shape,
                                               dtype=bool)

        # Select objects between the max stellar half-mass radius and
        # outer_radius
        # Field
        select_field_radius = ((self.field_distance >= self.max_rhalf) *
                               (self.field_distance <= outer_radius))
        field_m_shell = quantity_in_shells(
            self.field_m_current * select_field_radius *
            select_field_metallicity, self.field_distance, radius_bins)

        ################################################################
        # Disrupted objects
        # Select objects in the given metallicity range
        if metallicity_range is not None:
            select_disrupted_metallicity = (
                (self.disrupted_fe_h >= np.nanmin(metallicity_range)) *
                (self.disrupted_fe_h < np.nanmax(metallicity_range)))
        else:
            select_disrupted_metallicity = np.ones(
                self.disrupted_distance.shape, dtype=bool)

        select_disrupted_clusters = (
            (self.disrupted_distance >= self.max_rhalf) *
            (self.disrupted_distance <= outer_radius) *
            select_disrupted_metallicity)
        select_disrupted_gcs = (select_disrupted_clusters *
                                (self.disrupted_m_birth >= m_gc))

        # Mass from fully disrupted objects
        disrupted_cl_m_shell = quantity_in_shells(
            self.disrupted_m_disrupted * select_disrupted_clusters,
            self.disrupted_distance, radius_bins)
        disrupted_gcs_m_shell = quantity_in_shells(
            self.disrupted_m_disrupted * select_disrupted_gcs,
            self.disrupted_distance, radius_bins)
        disrupted_cl_m_init_shell = quantity_in_shells(
            self.disrupted_m_birth *
            self.disrupted_stellar_evolution_fraction *
            select_disrupted_clusters, self.disrupted_distance, radius_bins)
        disrupted_gcs_m_init_shell = quantity_in_shells(
            self.disrupted_m_birth *
            self.disrupted_stellar_evolution_fraction * select_disrupted_gcs,
            self.disrupted_distance, radius_bins)

        ################################################################
        # Partially disrupted objects
        # Select objects in the given metallicity range
        if metallicity_range is not None:
            select_part_disrupted_metallicity = (
                (self.part_disrupted_fe_h >= np.nanmin(metallicity_range)) *
                (self.part_disrupted_fe_h < np.nanmax(metallicity_range)))
        else:
            select_part_disrupted_metallicity = np.ones(
                self.part_disrupted_distance.shape, dtype=bool)
        select_part_disrupted_clusters = (
            (self.part_disrupted_distance >= self.max_rhalf) *
            (self.part_disrupted_distance <= outer_radius) *
            select_part_disrupted_metallicity)
        select_part_disrupted_gcs = (select_part_disrupted_clusters *
                                     (self.part_disrupted_m_birth >= m_gc))

        # Mass sloughed by partially disrupted objects
        part_disrupted_cl_m_shell = quantity_in_shells(
            self.part_disrupted_m_disrupted * select_part_disrupted_clusters,
            self.part_disrupted_distance, radius_bins)
        part_disrupted_gcs_m_shell = quantity_in_shells(
            self.part_disrupted_m_disrupted * select_part_disrupted_gcs,
            self.part_disrupted_distance, radius_bins)
        part_disrupted_cl_m_init_shell = quantity_in_shells(
            self.part_disrupted_m_birth *
            self.part_disrupted_stellar_evolution_fraction *
            select_part_disrupted_clusters, self.part_disrupted_distance,
            radius_bins)
        part_disrupted_gcs_m_init_shell = quantity_in_shells(
            self.part_disrupted_m_birth *
            self.part_disrupted_stellar_evolution_fraction *
            select_part_disrupted_gcs, self.part_disrupted_distance,
            radius_bins)

        f_halo_cl = (disrupted_cl_m_shell +
                     part_disrupted_cl_m_shell) / field_m_shell
        f_halo_gcs = (disrupted_gcs_m_shell +
                      part_disrupted_gcs_m_shell) / field_m_shell
        f_halo_init_cl = (disrupted_cl_m_init_shell +
                          part_disrupted_cl_m_init_shell) / field_m_shell
        f_halo_init_gcs = (disrupted_gcs_m_init_shell +
                           part_disrupted_gcs_m_init_shell) / field_m_shell

        return (radius_mid_bins, f_halo_cl, f_halo_gcs, f_halo_init_cl,
                f_halo_init_gcs)

    def reina_campos_f_halo(self,
                            outer_radius=50.,
                            m_gc=1.e5,
                            metallicity_range=None):
        # Select objects between the stellar half-mass radius and
        # outer_radius
        select_field_radius = ((self.field_distance >= self.r_half) *
                               (self.field_distance <= outer_radius))

        # Select objects in the given metallicity range
        if metallicity_range is not None:
            select_field_metallicity = (
                (self.field_fe_h >= np.nanmin(metallicity_range)) *
                (self.field_fe_h < np.nanmax(metallicity_range)))
            select_undisrupted_metallicity = (
                (self.undisrupted_fe_h >= np.nanmin(metallicity_range)) *
                (self.undisrupted_fe_h < np.nanmax(metallicity_range)))
            select_disrupted_metallicity = (
                (self.disrupted_fe_h >= np.nanmin(metallicity_range)) *
                (self.disrupted_fe_h < np.nanmax(metallicity_range)))
            select_part_disrupted_metallicity = (
                (self.part_disrupted_fe_h >= np.nanmin(metallicity_range)) *
                (self.part_disrupted_fe_h < np.nanmax(metallicity_range)))
        else:
            select_field_metallicity = np.ones(self.field_distance.shape,
                                               dtype=bool)
            select_undisrupted_metallicity = np.ones(
                self.undisrupted_distance.shape, dtype=bool)
            select_disrupted_metallicity = np.ones(
                self.disrupted_distance.shape, dtype=bool)
            select_part_disrupted_metallicity = np.ones(
                self.part_disrupted_distance.shape, dtype=bool)

        select_disrupted_clusters = (
            (self.disrupted_distance >= self.r_half) *
            (self.disrupted_distance <= outer_radius) *
            select_disrupted_metallicity)
        select_disrupted_gcs = (select_disrupted_clusters *
                                (self.disrupted_m_birth >= m_gc))

        select_part_disrupted_clusters = (
            (self.part_disrupted_distance >= self.r_half) *
            (self.part_disrupted_distance <= outer_radius) *
            select_part_disrupted_metallicity)
        select_part_disrupted_gcs = (select_part_disrupted_clusters *
                                     (self.part_disrupted_m_birth >= m_gc))

        select_undisrupted_clusters = (
            (self.undisrupted_distance >= self.r_half) *
            (self.undisrupted_distance <= outer_radius) *
            select_undisrupted_metallicity)
        select_undisrupted_gcs = (select_undisrupted_clusters *
                                  (self.undisrupted_m_birth >= m_gc))

        prefixes = ['disrupted_', 'part_disrupted_', 'undisrupted_']
        prefix_selections = [
            [select_disrupted_clusters, select_disrupted_gcs],
            [select_part_disrupted_clusters, select_part_disrupted_gcs],
            [select_undisrupted_clusters, select_undisrupted_gcs]
        ]

        m_cl = []
        m_gcs = []
        m_init_cl = []
        m_init_gcs = []

        for prefix, selections in zip(prefixes, prefix_selections):
            cl_selection, gc_selection = selections

            # Mass lost from clusters
            m_cl.append(
                np.nansum(getattr(self, prefix + 'm_disrupted') * cl_selection,
                          axis=0))
            # Mass lost from GCs
            m_gcs.append(
                np.nansum(getattr(self, prefix + 'm_disrupted') * gc_selection,
                          axis=0))

            m_init_cl.append(
                np.nansum(
                    getattr(self, prefix + 'm_birth') *
                    getattr(self, prefix + 'stellar_evolution_fraction') *
                    cl_selection,
                    axis=0))
            m_init_gcs.append(
                np.nansum(
                    getattr(self, prefix + 'm_birth') *
                    getattr(self, prefix + 'stellar_evolution_fraction') *
                    gc_selection,
                    axis=0))

        # Sum masses of disrupted objects
        m_disrupted_cl = np.nansum(self.disrupted_m_disrupted *
                                   select_disrupted_clusters,
                                   axis=0)
        m_disrupted_gc = np.nansum(self.disrupted_m_disrupted *
                                   select_disrupted_gcs,
                                   axis=0)
        m_init_disrupted_cl = np.nansum(
            self.disrupted_m_birth *
            self.disrupted_stellar_evolution_fraction *
            select_disrupted_clusters,
            axis=0)
        m_init_disrupted_gc = np.nansum(
            self.disrupted_m_birth *
            self.disrupted_stellar_evolution_fraction * select_disrupted_gcs,
            axis=0)

        # Sum mass sloughed by partially disrupted objects (i.e. zero)
        m_part_disrupted_cl = np.nansum(self.part_disrupted_m_disrupted *
                                        select_part_disrupted_clusters,
                                        axis=0)
        m_part_disrupted_gc = np.nansum(self.part_disrupted_m_disrupted *
                                        select_part_disrupted_clusters *
                                        select_part_disrupted_gcs,
                                        axis=0)
        m_init_part_disrupted_cl = np.nansum(
            self.part_disrupted_m_birth *
            self.part_disrupted_stellar_evolution_fraction *
            select_part_disrupted_clusters,
            axis=0)
        m_init_part_disrupted_gc = np.nansum(
            self.part_disrupted_m_birth *
            self.part_disrupted_stellar_evolution_fraction *
            select_part_disrupted_gcs,
            axis=0)

        # Sum mass sloughed by undisrupted objects
        m_undisrupted_cl = np.nansum(self.undisrupted_m_disrupted *
                                     select_undisrupted_clusters,
                                     axis=0)
        m_undisrupted_gc = np.nansum(self.undisrupted_m_disrupted *
                                     select_undisrupted_clusters *
                                     select_undisrupted_gcs,
                                     axis=0)
        m_init_undisrupted_cl = np.nansum(
            self.undisrupted_m_birth *
            self.undisrupted_stellar_evolution_fraction *
            select_undisrupted_clusters,
            axis=0)
        m_init_undisrupted_gc = np.nansum(
            self.undisrupted_m_birth *
            self.undisrupted_stellar_evolution_fraction *
            select_undisrupted_gcs,
            axis=0)

        # Sum mass of halo stars
        m_halo = np.nansum(self.field_m_current * select_field_radius *
                           select_field_metallicity,
                           axis=0)

        # self.f_halo_all_clusters = (m_disrupted_cl + m_part_disrupted_cl +
        #                             m_undisrupted_cl) / m_halo
        self.f_halo_all_clusters = np.nansum(m_cl, axis=0) / m_halo
        # self.f_halo_formed_as_gcs = (m_disrupted_gc + m_part_disrupted_gc +
        #                              m_undisrupted_gc) / m_halo
        self.f_halo_formed_as_gcs = np.nansum(m_gcs, axis=0) / m_halo
        # self.f_halo_init_all_clusters = (m_init_disrupted_cl +
        #                                  m_init_part_disrupted_cl +
        #                                  m_init_undisrupted_cl) / m_halo
        self.f_halo_init_all_clusters = np.nansum(m_init_cl, axis=0) / m_halo
        # self.f_halo_init_formed_as_gcs = (m_init_disrupted_gc +
        #                                   m_init_part_disrupted_gc +
        #                                   m_init_undisrupted_gc) / m_halo
        self.f_halo_init_formed_as_gcs = np.nansum(m_init_gcs, axis=0) / m_halo

        return (self.f_halo_all_clusters, self.f_halo_formed_as_gcs,
                self.f_halo_init_all_clusters, self.f_halo_init_formed_as_gcs)


class EvolutionData(object):

    def __init__(self, sim) -> None:
        self.sim = sim
        self.gc_mass = gc_mass
        self.process_data()
        self.make_specific_gas_dmgc_sfr()

        return None

    def process_data(self):
        with h5py.File(data_file_template.format(self.sim, gc_mass), 'r') as f:
            for key in f['z evolution'].keys():
                data = f['z evolution'][key][()]
                str_key = key.strip().replace(' ', '_')
                setattr(self, str_key, data)
        return None

    def make_specific_gas_dmgc_sfr(self):
        sfr = getattr(self, 'SFR')
        net_dmgc = getattr(self, 'Net_dM_GCdt')
        mgas_sf = getattr(self, 'M_gas,SF')

        sfr_mgas_sf = sfr / mgas_sf
        dmgc_mgas_sf = net_dmgc / mgas_sf

        sfr_mgas_sf[np.isinf(sfr_mgas_sf)] = np.nan
        dmgc_mgas_sf[np.isinf(dmgc_mgas_sf)] = np.nan

        setattr(self, 'SFR_MgasSF', sfr_mgas_sf)
        setattr(self, 'dMgc_dt_MgasSF', dmgc_mgas_sf)

        return None

    def med_spread(self, attr_name, confidence=[16., 84.]):
        """Returns the median and percentiles of a given property.

        Args:
            attr_name (str): One of the attributes set by
                self.process_data.
            confidence (arr/list, optional): Lower and upper percentiles
                to return. Defaults to [16., 84.].

        Returns:
            tuple:  [0]: Median.
                    [1]: Lower and upper percentiles.
        """
        # print(np.all(~np.isnan(getattr(self, attr_name)), axis=0).sum())
        print(np.sum(~np.isnan(getattr(self, attr_name)), axis=1))
        med = np.nanmedian(getattr(self, attr_name), axis=1)
        spread = np.nanpercentile(getattr(self, attr_name), confidence, axis=1)
        return (med, spread)
