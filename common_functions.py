#!/usr/bin/env python3
# Place import files below
import os
import string

import numpy as np
import scipy.integrate as integrate
from astropy import constants as const
from matplotlib.ticker import Locator

panel_labels = [it + ")" for it in list(string.ascii_lowercase)]


def embed_symbols(pdf_file):
    """Embeds symobls in pdf files.

    Args:
        pdf_file (str): Filepath to the file

    Returns:
        None
    """
    os.system('gs -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress '
              '-dEmbedAllFonts=true -sOutputFile={} -f {}'.format(
                  pdf_file.replace('.pdf', '_embedded.pdf'), pdf_file))
    return None


def plot_merger_arrow(
        ax,
        x,
        length,
        #   ls='-',
        #   fc='k',
        arrow_properties={},
        loc='upper'):
    if loc == 'upper':
        y = 1. - length
        dy = length
    elif loc == 'lower':
        y = length
        dy = -length
    else:
        y = 1.
        dy = 1.

    ax.arrow(x, y, 0., dy, transform=ax.transAxes, **arrow_properties)
    return None


def save_figures(fig, location, embed=False):
    """Saves svg and pdf versions of figures.

    Args:
        fig (Matplotlib figure object): The figure to save
        location (str): Filepath to the save file
        embed (bool): If True, embeds the symbols in the pdf file.
            Default: False.

    Returns:
        None
    """
    if '.pdf' in location:
        pdf_file = location
        svg_file = location.replace('.pdf', '.svg')
    else:
        pdf_file = location + '.pdf'
        svg_file = location + '.svg'

    fig.savefig(pdf_file, dpi=600, format='pdf', transparent=False)
    fig.savefig(svg_file, dpi=600, format='svg', transparent=False)

    if embed:
        embed_symbols(pdf_file)

    return None


def tick_function(redshifts, cos_obj):
    return cos_obj.compute_lookback_time(redshifts)


class Cosmology(object):

    def __init__(self, args={}) -> None:
        self.omega_M = 0.25
        self.omega_lambda = 0.75
        self.omega_k = 0.0
        self.omega_r = 0.0
        self.h = 0.72

        for key in args:
            setattr(self, key, args[key])

        return None

    def compute_lookback_time(self, z):
        """Computes lookback time via numerical integration over z.

        Args:
            z (fl/list): Redshifts to calculate lookback time to.

        Returns:
            list: Time interval(s)
                [0]: Lookback time [s]
                [1]: dict containing 'Units'
        """
        # Initial type-checking
        permitted_types = (float, list, np.ndarray)
        if not isinstance(z, permitted_types):
            raise TypeError('z must be float, list or np.ndarray')

        # Convert lists to np.ndarray
        if isinstance(z, list):
            z = np.asarray(z)

        return self.compute_time_interval(z, 0.)

    def compute_time_interval(self, z1, z2):
        """Computes time interval via numerical integration over z

        Args:
            z1 (fl/list): Starting redshift (high)
            z2 (fl/list): End redshift (low)

        Returns:
            list: Time interval(s)
                [0]: Lookback time [s]
                [1]: dict containing 'Units'
        """
        array_flag = False

        # Initial type-checking
        permitted_types = (float, list, np.ndarray)
        if not isinstance(z1, permitted_types):
            raise TypeError('z1 must be float, list or np.ndarray')
        if not isinstance(z2, permitted_types):
            raise TypeError('z2 must be float, list or np.ndarray')

        # Convert lists to np.ndarray
        if isinstance(z1, list):
            z1 = np.asarray(z1)
        if isinstance(z2, list):
            z2 = np.asarray(z2)
        if isinstance(z1, float) and isinstance(z2, np.ndarray):
            z1 = np.repeat(np.asarray([z1]), len(z2))
        if isinstance(z2, float) and isinstance(z1, np.ndarray):
            z2 = np.repeat(np.asarray([z2]), len(z1))

        if isinstance(z1, np.ndarray):
            if len(z1) != len(z2):
                raise IOError(
                    'z1 (len: {}) and z2 (len: {}) '.format(len(z1), len(z2)) +
                    'arrays must be the same length')
            if (z1 < z2).any():
                raise IOError(
                    'All z1 must be greater than corresponding z2 values')
            array_flag = True
        else:
            if z1 < z2:
                raise IOError('z1 must be greater than z2 values')

        def integrand(z):
            return 1. / ((1. + z) * np.sqrt(self.E_z(z)))

        if array_flag:
            integrate_out = np.empty(len(z1))
            for z_i, (z1_item, z2_item) in enumerate(zip(z1, z2)):
                integrate_out[z_i] = integrate.quad(integrand, z1_item,
                                                    z2_item)[0]
        else:
            integrate_out = integrate.quad(integrand, z1, z2)[0]

        delta_t = integrate_out / -self.hubble_parameter_in_sec()

        return delta_t / (1.e9 * 365.25 * 24. * 3600.)  # Gyr

    def E_z(self, z):
        """Calculates the square of :math:`H(z) / H_0`

            Args:
                z (float/list): Redshifts at which to calculate :math:`E(z)`

            Returns:
                float/list: Calculated :math:`E(z)` values.
            """

        output = (self.omega_r * (1. + z)**4. + self.omega_M * (1. + z)**3. +
                  self.omega_k * (1. + z)**2. + self.omega_lambda)
        return output

    def hubble_parameter_in_sec(self):
        """Returns the Hubble parameter in units of [1 / s]

        Returns:
            list:
                [0]: Value of H_0
                [1]: Dict of attributes, esp. units
        """
        return self.h / (10. * const.pc.value)  # [1 / s]


class FigureHandler(object):

    def __init__(
            self,
            major_redshifts,
            minor_redshifts,
            tlb_lim,
            cos_obj,
            xlabel=r'$t_{\rm lookback}\, \left[{\rm Gyr}\right]$') -> None:
        self.major_redshifts = major_redshifts
        self.minor_redshifts = minor_redshifts
        self.tlb_lim = tlb_lim
        self.cos_obj = cos_obj
        self.x_label = xlabel
        return None

    def set_figure_properties(
        self,
        axs,
        ylabels,
        yscale='log',
        ylims=None,
        #   leg_loc='lower right',
        leg_loc='best',
        linthresh=1.e-3,
        no_legend=False,
    ):
        if yscale is None:
            yscale_prop_list = [None for _ in axs]
        elif isinstance(yscale, str):
            yscale_prop_list = [yscale for _ in axs]
        else:
            yscale_prop_list = yscale
        sym_scale = False
        # yscale_prop = yscale
        # if yscale == 'symlog':
        #     yscale_prop = None
        if ylims is None:
            ylim_list = [None for _ in axs]
        else:
            ylim_list = ylims

        for i, (ax, ylabel, ylim, yscale_prop) in enumerate(
                zip(axs, ylabels, ylim_list, yscale_prop_list)):
            new_x_ax = SecondXAxis(ax, tick_function)
            new_x_ax.set_major_x_ticks(self.major_redshifts,
                                       self.cos_obj,
                                       label_override=[
                                           '{:d}'.format(int(num))
                                           for num in self.major_redshifts
                                       ])
            new_x_ax.set_minor_x_ticks(self.minor_redshifts, self.cos_obj)
            new_x_ax.set_axis_limits(self.tlb_lim)
            new_x_ax.set_xlabel(r'$z$')
            if yscale_prop == 'symlog':
                yscale_prop = None
                sym_scale = True
            ax.set(xlabel=self.x_label,
                   ylabel=ylabel,
                   yscale=yscale_prop,
                   ylim=ylim)
            ax.minorticks_on()
            if sym_scale:
                ax.set_yscale('symlog', linthresh=linthresh)
                ax.yaxis.set_minor_locator(MinorSymLogLocator(linthresh))

            new_x_ax.invert_axis()
            if not no_legend:
                ax.legend(fontsize='small', loc=leg_loc)

            sym_scale = False

        return None

    def set_figure_with_ratio_properties(
            self,
            axs,
            ylabels,
            yscale='log',
            ylims=None,
            #  leg_loc='lower right',
            leg_loc='best',
            linthresh=1.e-3):
        yscale_prop = yscale
        if yscale == 'symlog':
            yscale_prop = None
        if ylims is None:
            ylim_list = [[None, None] for _ in axs]
        else:
            ylim_list = ylims

        for i, (ax, ylabel,
                (ylim_0, ylim_1)) in enumerate(zip(axs, ylabels, ylim_list)):
            ax[0].axhline(1, color='k', linestyle=':')
            new_x_ax = SecondXAxis(ax[0], tick_function)
            new_x_ax.set_major_x_ticks(self.major_redshifts,
                                       self.cos_obj,
                                       label_override=[
                                           '{:d}'.format(int(num))
                                           for num in self.major_redshifts
                                       ])
            new_x_ax.set_minor_x_ticks(self.minor_redshifts, self.cos_obj)
            new_x_ax.set_axis_limits(self.tlb_lim)
            new_x_ax.set_xlabel(r'$z$')
            ax[0].set(ylabel=ylabel[0], yscale='linear', ylim=ylim_0)
            ax[1].set(xlabel=self.x_label,
                      ylabel=ylabel[1],
                      yscale=yscale_prop,
                      ylim=ylim_1)
            ax[0].minorticks_on()
            ax[1].minorticks_on()
            if yscale == 'symlog':
                ax[1].set_yscale('symlog', linthresh=linthresh)
                ax[1].yaxis.set_minor_locator(MinorSymLogLocator(linthresh))
            new_x_ax.invert_axis()
            ax[1].legend(fontsize='small', loc=leg_loc)
        return None

    def set_stacked_figure_properties(
            self,
            axs,
            ylabels,
            yscale='log',
            ylims=None,
            #   leg_loc='lower right',
            leg_loc='best',
            linthresh=1.e-3,
            text_labels=None,
            text_label_args=None,
            direction='vertical'):
        """Sets the figure properties for stacked axes.

        Args:
            axs (list): List of Axes objects.
            ylabels (list): List of y-axis labels to be applied to each
                axis.
            yscale (str, optional): Scale of y-axis. Defaults to 'log'.
            ylims (list, optional): Lower and upper y-axis values.
                Defaults to None.
            leg_loc (str, optional): Location of legend. Defaults to
                'best'.
            linthresh (fl, optional): Threshold at which symlog axes
                switch from linear to logarithmic. Defaults to 1.e-3.
            text_labels (list, optional): List of labels to be applied
                to each axis. Defaults to None.
            text_label_args (dict, optional): args to set the properties
                of each text label. Defaults to None.

        Returns:
            None
        """
        if yscale is None:
            yscale_prop_list = [None for _ in axs]
        else:
            yscale_prop_list = yscale
        # if yscale == 'symlog':
        #     yscale_prop = None
        if ylims is None:
            ylim_list = [[None, None] for _ in axs]
        else:
            ylim_list = ylims
        if text_labels is None:
            text_label_list = [None for _ in axs]
        else:
            text_label_list = text_labels
        if text_label_args is None:
            text_label_arg_list = [None for _ in axs]
        else:
            text_label_arg_list = text_label_args

        for i, (ax, ax_label, ylabel, ylim, yscale, t_label,
                t_label_args) in enumerate(
                    zip(axs, panel_labels, ylabels, ylim_list,
                        yscale_prop_list, text_label_list,
                        text_label_arg_list)):
            # Add redshift x-axis on the topmost panel only
            if (((direction == 'vertical') and (i == 0))
                    or (direction == 'horizontal')):
                new_x_ax = SecondXAxis(ax, tick_function)
                new_x_ax.set_major_x_ticks(self.major_redshifts,
                                           self.cos_obj,
                                           label_override=[
                                               '{:d}'.format(int(num))
                                               for num in self.major_redshifts
                                           ])
                new_x_ax.set_minor_x_ticks(self.minor_redshifts, self.cos_obj)
                new_x_ax.set_axis_limits(self.tlb_lim)
                new_x_ax.set_xlabel(r'$z$')

            if yscale != 'symlog':
                ax.set(yscale=yscale)
            if (ylim is not None):
                if (ylim[0] is not None) and (ylim[0] < 0.):
                    ax.axhline(0., color='k', linestyle=':', zorder=0)
            ax.set(ylim=ylim)
            ax.set_ylabel(ylabel=ylabel, fontsize='small')

            # Set lower x-label on the last axis object in vertical stack
            if (((direction == 'vertical') and (i == len(axs) - 1))
                    or (direction == 'horizontal')):
                ax.set(xlabel=self.x_label)
            ax.minorticks_on()
            if yscale == 'symlog':
                ax.set_yscale('symlog', linthresh=linthresh)
                ax.yaxis.set_minor_locator(MinorSymLogLocator(linthresh))

            # Sets legend in the first panel
            if i == 0:
                ax.legend(fontsize='small', loc=leg_loc)

            if t_label is not None:
                ax.text(s=t_label, transform=ax.transAxes, **t_label_args)

            if direction == 'horizontal':
                new_x_ax.invert_axis()

            if len(axs) > 2:
                ax.text(0.02,
                        0.96,
                        ax_label,
                        color='k',
                        verticalalignment='top',
                        transform=ax.transAxes)

        if direction == 'vertical':
            new_x_ax.invert_axis()

        return None

    def save_figures(figs, file_paths):
        for fig, file_path in zip(figs, file_paths):
            save_figures(fig, file_path, embed=True)
        return None


class MinorSymLogLocator(Locator):
    """
    Dynamically find minor tick positions based on the positions of
    major ticks for a symlog scaling.

    Courtesy of stackoverflow users 'Matt Pitkin' and 'Peruz':
    https://stackoverflow.com/questions/20470892/how-to-place-minor-ticks-on-symlog-scale
    """

    def __init__(self, linthresh, nints=10):
        """
        Ticks will be placed between the major ticks.
        The placement is linear for x between -linthresh and linthresh,
        otherwise its logarithmically. nints gives the number of
        intervals that will be bounded by the minor ticks.
        """
        self.linthresh = linthresh
        self.nintervals = nints

    def __call__(self):
        # Return the locations of the ticks
        majorlocs = self.axis.get_majorticklocs()

        first_major = majorlocs[0]
        if first_major == 0:
            outrange_first = -self.linthresh
        else:
            outrange_first = first_major * float(10)**(-np.sign(first_major))
        # top of the axis (high values)
        last_major = majorlocs[-1]
        if last_major == 0:
            outrange_last = self.linthresh
        else:
            outrange_last = last_major * float(10)**(np.sign(last_major))
        majorlocs = np.concatenate(
            ([outrange_first], majorlocs, [outrange_last]))

        if len(majorlocs) == 1:
            return self.raise_if_exceeds(np.array([]))

        # add temporary major tick locs at either end of the current range
        # to fill in minor tick gaps
        dmlower = majorlocs[1] - majorlocs[
            0]  # major tick difference at lower end
        dmupper = majorlocs[-1] - majorlocs[
            -2]  # major tick difference at upper end

        # add temporary major tick location at the lower end
        if majorlocs[0] != 0. and (((majorlocs[0] != self.linthresh) and
                                    (dmlower > self.linthresh)) or
                                   ((dmlower == self.linthresh) and
                                    (majorlocs[0] < 0))):
            majorlocs = np.insert(majorlocs, 0, majorlocs[0] * 10.)
        else:
            majorlocs = np.insert(majorlocs, 0, majorlocs[0] - self.linthresh)

        # add temporary major tick location at the upper end
        if majorlocs[-1] != 0. and ((np.abs(
            (majorlocs[-1]) != self.linthresh) and
                                     (dmupper > self.linthresh)) or
                                    ((dmupper == self.linthresh) and
                                     (majorlocs[-1] > 0))):
            majorlocs = np.append(majorlocs, majorlocs[-1] * 10.)
        else:
            majorlocs = np.append(majorlocs, majorlocs[-1] + self.linthresh)

        # iterate through minor locs
        minorlocs = []

        # handle the lowest part
        for i in np.arange(1, len(majorlocs)):
            majorstep = majorlocs[i] - majorlocs[i - 1]
            if abs(majorlocs[i - 1] + majorstep / 2) < self.linthresh:
                ndivs = self.nintervals
            else:
                ndivs = self.nintervals - 1.

            minorstep = majorstep / ndivs
            locs = np.arange(majorlocs[i - 1], majorlocs[i], minorstep)[1:]
            minorlocs.extend(locs)

        return self.raise_if_exceeds(np.array(minorlocs))

    def tick_values(self, vmin, vmax):
        raise NotImplementedError('Cannot get tick locations for a '
                                  '{} type.'.format(type(self)))


class SecondXAxis(object):

    def __init__(self, ax, conversion_func) -> None:
        """Add a second x-axis to a matplotlib figure instance.

        Args:
            ax (Axis object): Matplotlib axis object
            conversion_func (function): Function to convert values from
                the second x-axis to the main x-axis.

        Returns:
            None
        """
        self.ax1 = ax
        self.ax2 = ax.twiny()
        self.conversion_func = conversion_func
        return None

    def set_axis_limits(self, values):
        """Sets the axis limits of *both* axes.

        Args:
            values (arr): Upper and lower axis limits for the *main*
                axis.

        Returns:
            None
        """
        self.ax1.set_xlim(values)
        self.ax2.set_xlim(values)
        return None

    def invert_axis(self):
        """Inverts *both* axes.

        Args:
            None

        Returns:
            None
        """
        self.ax1.invert_xaxis()
        self.ax2.invert_xaxis()
        return None

    def set_major_x_ticks(self, values, *args, **kwargs):
        """Sets major x-ticks and labels on the new x-axis

        Args:
            values (fl/arr): Array of values to be set on the second
                x-axis.

        Returns:
            None
        """
        self.ax2.set_xticks(self.conversion_func(values, *args))
        if 'label_override' in kwargs:
            self.ax2.set_xticklabels(kwargs['label_override'])
        else:
            self.ax2.set_xticklabels(values)

        return None

    def set_minor_x_ticks(self, values, *args):
        """Sets minor x-ticks on the new x-axis

        Args:
            values (fl/arr): Array of values to be set on the second
                x-axis.

        Returns:
            None
        """
        self.ax2.set_xticks(self.conversion_func(values, *args), minor=True)
        return None

    def set_xlabel(self, label):
        """Sets x-axis label on the new x-axis

        Args:
            label (str): Label to set on the new x-axis.

        Returns:
            None
        """
        self.ax2.set(xlabel=label)
        return None
