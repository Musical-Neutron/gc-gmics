#!/usr/bin/env python3

import h5py
import numpy as np
import pytest

from gcgmics.common_functions import (
    Cosmology,
    FigureHandler,
    MergerArrow,
    MinorSymLogLocator,
    SecondXAxis,
    embed_symbols,
    get_scaled_arrow_properties,
    plot_merger_arrow,
    save_figures,
    tick_function,
)


class TestCosmology:

    def test_default_initialisation(self):
        cosmo = Cosmology()
        assert cosmo.omega_M == 0.25
        assert cosmo.omega_lambda == 0.75
        assert cosmo.omega_k == 0.0
        assert cosmo.omega_r == 0.0
        assert cosmo.h == 0.72
        return None

    def test_custom_initialisation(self):
        args = {"omega_M": 0.3, "omega_lambda": 0.7, "h": 0.7}
        cosmo = Cosmology(args)
        assert cosmo.omega_M == 0.3
        assert cosmo.omega_lambda == 0.7
        assert cosmo.omega_k == 0.0
        assert cosmo.omega_r == 0.0
        assert cosmo.h == 0.7
        return None

    @pytest.mark.parametrize("z", ["s", 2])
    def test_compute_lookback_time_invalid(self, z):
        """Test output of Cosmology.compute_lookback_time"""
        # Should raise error if z not float, list or np.ndarray
        cosmo = Cosmology()
        with pytest.raises(TypeError):
            cosmo.compute_lookback_time(z)

    @pytest.mark.parametrize(
        "z, expected",
        [
            (0.3, 3.35871824),
            (0.1, 1.269832863),
            ([0.3, 0.1], [3.35871824, 1.269832863]),
        ],
    )
    def test_compute_lookback_time_valid(self, z, expected):
        """Test output of Cosmology.compute_lookback_time"""
        cosmo = Cosmology()
        assert cosmo.compute_lookback_time(z) == pytest.approx(expected)

    @pytest.mark.parametrize(
        "z1, z2",
        [
            ("s", 0.1),
            (0.2, "s"),
            ("s", "s"),
            (2, 0.1),
            (2.0, 1),
        ],
    )
    def test_compute_time_interval_invalid(self, z1, z2):
        cosmo = Cosmology()
        with pytest.raises(TypeError):
            cosmo.compute_time_interval(z1, z2)

        return None

    @pytest.mark.parametrize(
        "z1, z2",
        [
            ([0.3, 0.4, 0.5], [0.1, 0.2]),
            ([0.4, 0.5], [0.1, 0.2, 0.3]),
        ],
    )
    def test_compute_time_interval_mismatched_inputs(self, z1, z2):
        cosmo = Cosmology()
        with pytest.raises(IOError):
            cosmo.compute_time_interval(z1, z2)

        return None

    @pytest.mark.parametrize(
        "z1, z2",
        [
            (0.1, 0.2),
            ([0.1, 0.2], [0.2, 0.3]),
            ([0.1, 0.2], [0.2, 0.1]),
            ([0.2, 0.1], [0.1, 0.2]),
        ],
    )
    def test_compute_time_interval_invalid_z_order(self, z1, z2):
        cosmo = Cosmology()
        with pytest.raises(IOError):
            cosmo.compute_time_interval(z1, z2)

        return None

    @pytest.mark.parametrize(
        "z1, z2",
        [
            (np.nan, 0.1),
            (0.1, np.nan),
            (np.nan, np.nan),
        ],
    )
    def test_compute_time_interval_nan_scalar(self, z1, z2):
        cosmo = Cosmology()

        assert np.isnan(cosmo.compute_time_interval(z1, z2))

        return None

    @pytest.mark.parametrize(
        "z1, z2",
        [
            ([0.2, 0.3], [0.1, np.nan]),
            ([0.2, 0.3], [np.nan, 0.1]),
            ([0.2, np.nan], [0.1, 0.2]),
            ([np.nan, 0.3], [0.1, 0.2]),
        ],
    )
    def test_compute_time_interval_nan_array(self, z1, z2):
        cosmo = Cosmology()

        assert any(np.isnan(cosmo.compute_time_interval(z1, z2)))

        return None

    @pytest.mark.parametrize(
        "z1, z2, expected",
        [
            (0.3, 0.1, 2.088885),
            ([0.3], 0.1, 2.088885),
        ],
    )
    def test_compute_time_interval_valid(self, z1, z2, expected):
        """Test output of Cosmology.compute_time_interval"""
        cosmo = Cosmology()
        assert cosmo.compute_time_interval(z1, z2) == pytest.approx(expected)

    def test_E_z(self):
        """Test output of Cosmology.E_z"""
        cosmo = Cosmology()
        z = 1.0

        expected = 0.25 * 2**3 + 0.75
        assert np.isclose(cosmo.E_z(z), expected)

        z = np.asarray([1.0, 2.0, 3.0])
        expected = [0.25 * 2**3 + 0.75, 0.25 * 3**3 + 0.75, 0.25 * 4**3 + 0.75]
        assert np.array_equal(cosmo.E_z(z), expected)

        return None

    def test_hubble_parameter_in_sec(self):
        """Test output of Cosmology.hubble_parameter_in_sec"""
        cosmo = Cosmology()
        assert cosmo.hubble_parameter_in_sec() == 0.72 / (10.0 * 3.085677581491367e16)
