# SPDX-FileCopyrightText: BioPAL <biopal@esa.int>
# SPDX-License-Identifier: MIT

import warnings
import numpy as np
from matplotlib import pyplot as plt

__all__ = ["plot", "plot_db", "plot_abs", "plot_angle", "plot_rad2deg"]


def __plot(rasterObj, data_to_plot, ax=None, title=None, clims=None, xlims=None, ylims=None):

    if ax is None:
        _, ax = plt.subplots()

    extent = [
        rasterObj.x_axis[0],
        rasterObj.x_axis[-1],
        rasterObj.y_axis[-1],
        rasterObj.y_axis[0],
    ]
    if rasterObj._data_type == "SlantRange_Azimuth":
        cmap = "gray"
    elif rasterObj._data_type == "ground":
        cmap = "YlGn"
    else:
        cmap = None
        warnings('Cannot recognaze input data object data_type "{}"'.format(rasterObj._data_type))

    image = ax.imshow(data_to_plot, interpolation="none", origin="upper", cmap=cmap, extent=extent,)
    if clims:
        image.set_clim(clims)
    plt.colorbar(image, ax=ax)
    ax.set_xlabel(rasterObj.x_axis_description)
    ax.set_ylabel(rasterObj.y_axis_description)
    if title:
        ax.set_title(title)
    if xlims:
        ax.set_xlim(xlims)
    if ylims:
        ax.set_ylim(ylims)

    return ax


def plot(rasterObj, ax=None, title=None, clims=None, xlims=None, ylims=None):
    """Plotting function for BiomassL1cRaster and BiomassL2Raster data objects
    
    Parameters
    ----------
    rasterObj : BiomassL1cRaster or BiomassL2Raster object
        initialized object of the data to plot 
    ax : pyplot.subplot axes, optional
        for custom settings
    title : str, optional
    clims : list of numbers, optional
        plot saturation limits [minvalue, maxvalue]
        unit is [km], [deg] or [km], depending on the rasterObj.data content (respectively L1, intermediate or L2)
    xlims : list of numbers, optional
        plot x axis limits [minvalue, maxvalue]
        slant range [km], longitude [deg] or east [km], depending on rasterObj.data content (respectively L1, intermediate or L2)
    ylims : list of numbers, optional
        plot y axis limits [minvalue, maxvalue]
        azimuth [km], latitude [deg] or north [km], depending on rasterObj.data content (respectively L1, intermediate or L2)


    See Also
    --------
    plot_db : plot the data in dB
    plot_abs : plot the absolute value of data
    plot_angle : plot the angle (degrees) of complex data
    plot_rad2deg : plot the data after conversion from radiants to degrees

    Notes
    -----
    plot, plots the data as it is:

    .. math:: rasterObj.data

    Examples
    --------
    >>> rasterObj = BiomassL1cRaster(folder_path)
    >>> plot(rasterObj)
    >>> plt.show()

    >>> rasterObj = BiomassL2Raster(tif_path)
    >>> fig, ax = plt.subplots()
    >>> plot(rasterObj, ax=ax)
    >>> ax.set_xlabel("my xlabel")
    >>> fig.show()
    """
    return __plot(rasterObj, rasterObj.data, ax, title, clims, xlims, ylims,)


def plot_db(rasterObj, ax=None, title=None, clims=None, xlims=None, ylims=None):
    """Plotting function (dB) for BiomassL1cRaster and BiomassL2Raster data objects
    
    Parameters
    ----------
    rasterObj : BiomassL1cRaster or BiomassL2Raster object
        initialized object of the data to plot 
    ax : pyplot.subplot axes, optional
        for custom settings
    title : str, optional
    clims : list of numbers, optional
        plot saturation limits [minvalue, maxvalue]
        unit is [km], [deg] or [km], depending on the rasterObj.data content (respectively L1, intermediate or L2)
    xlims : list of numbers, optional
        plot x axis limits [minvalue, maxvalue]
        slant range [km], longitude [deg] or east [km], depending on rasterObj.data content (respectively L1, intermediate or L2)
    ylims : list of numbers, optional
        plot y axis limits [minvalue, maxvalue]
        azimuth [km], latitude [deg] or north [km], depending on rasterObj.data content (respectively L1, intermediate or L2)

    See Also 
    --------
    plot : plot the data as it is
    plot_abs : plot the absolute value of data
    plot_angle : plot the angle (degrees) of complex data
    plot_rad2deg : plot the data after conversion from radiants to degrees

    Notes
    -----
    plot_db, plots the data in dB:

    .. math:: 10*log10( abs( rasterObj.data )**2 )

    Examples
    --------
    >>> rasterObj = BiomassL1cRaster(folder_path)
    >>> plo_db(rasterObj)
    >>> plt.show()

    >>> rasterObj = BiomassL2Raster(tif_path)
    >>> fig, ax = plt.subplots()
    >>> plot_db(rasterObj, ax=ax)
    >>> ax.set_xlabel("my xlabel")
    >>> fig.show()
    """
    if title is None:
        title = "[dB]"
    return __plot(rasterObj, 10 * np.log10(np.abs(rasterObj.data) ** 2), ax, title, clims, xlims, ylims,)


def plot_abs(rasterObj, ax=None, title=None, clims=None, xlims=None, ylims=None):
    """Plotting function (ABS) for BiomassL1cRaster and BiomassL2Raster data objects

    Parameters
    ----------
    rasterObj : BiomassL1cRaster or BiomassL2Raster object
        initialized object of the data to plot 
    ax : pyplot.subplot axes, optional
        for custom settings
    title : str, optional
    clims : list of numbers, optional
        plot saturation limits [minvalue, maxvalue]
        unit is [km], [deg] or [km], depending on the rasterObj.data content (respectively L1, intermediate or L2)
    xlims : list of numbers, optional
        plot x axis limits [minvalue, maxvalue]
        slant range [km], longitude [deg] or east [km], depending on rasterObj.data content (respectively L1, intermediate or L2)
    ylims : list of numbers, optional
        plot y axis limits [minvalue, maxvalue]
        azimuth [km], latitude [deg] or north [km], depending on rasterObj.data content (respectively L1, intermediate or L2)


    See Also
    --------
    plot : plot the data as it is
    plot_db : plot the data in dB
    plot_angle : plot the angle (degrees) of complex data
    plot_rad2deg : plot the data after conversion from radiants to degrees

    Notes
    -----
    plot_abs, plots the absolute value of data:

    .. math:: abs( rasterObj.data )

    Examples
    --------
    >>> rasterObj = BiomassL1cRaster(folder_path)
    >>> plot_abs(rasterObj)
    >>> plt.show()

    >>> rasterObj = BiomassL2Raster(tif_path)
    >>> fig, ax = plt.subplots()
    >>> plot_abs(rasterObj, ax=ax)
    >>> ax.set_xlabel("my xlabel")
    >>> fig.show()
    """
    if title is None:
        title = "ABS"
    return __plot(rasterObj, np.abs(rasterObj.data), ax, title, clims, xlims, ylims,)


def plot_angle(rasterObj, ax=None, title=None, clims=None, xlims=None, ylims=None):
    """Plotting function (angle) for BiomassL1cRaster and BiomassL2Raster data objects
    
    Parameters
    ----------
    rasterObj : BiomassL1cRaster or BiomassL2Raster object
        initialized object of the data to plot 
    ax : pyplot.subplot axes, optional
        for custom settings
    title : str, optional
    clims : list of numbers, optional
        plot saturation limits [minvalue, maxvalue]
        unit is [km], [deg] or [km], depending on the rasterObj.data content (respectively L1, intermediate or L2)
    xlims : list of numbers, optional
        plot x axis limits [minvalue, maxvalue]
        slant range [km], longitude [deg] or east [km], depending on rasterObj.data content (respectively L1, intermediate or L2)
    ylims : list of numbers, optional
        plot y axis limits [minvalue, maxvalue]
        azimuth [km], latitude [deg] or north [km], depending on rasterObj.data content (respectively L1, intermediate or L2)


    See Also
    --------
    plot : plot the data as it is
    plot_db : plot the data in dB
    plot_abs : plot the absolute value of data
    plot_rad2deg : plot the data after conversion from radiants to degrees

    Notes
    ----- 
    plot_angle, plots the angle (degrees) of complex data:

    .. math:: angle ( rasterObj.data )

    Examples
    --------
    >>> rasterObj = BiomassL1cRaster(folder_path)
    >>> plot_angle(rasterObj)
    >>> plt.show()

    >>> rasterObj = BiomassL2Raster(tif_path)
    >>> fig, ax = plt.subplots()
    >>> plot_angle(rasterObj, ax=ax)
    >>> ax.set_xlabel("my xlabel")
    >>> fig.show()
    """
    if title is None:
        title = "Angle [deg]"
    return __plot(rasterObj, np.rad2deg(np.angle(rasterObj.data)), ax, title, clims, xlims, ylims,)


def plot_rad2deg(rasterObj, ax=None, title=None, clims=None, xlims=None, ylims=None):

    """Plotting function (rad2deg) for BiomassL1cRaster and BiomassL2Raster data objects
    
    Parameters
    ----------
    rasterObj : BiomassL1cRaster or BiomassL2Raster object
        initialized object of the data to plot 
    ax : pyplot.subplot axes, optional
        for custom settings
    title : str, optional
    clims : list of numbers, optional
        plot saturation limits [minvalue, maxvalue]
        unit is [km], [deg] or [km], depending on the rasterObj.data content (respectively L1, intermediate or L2)
    xlims : list of numbers, optional
        plot x axis limits [minvalue, maxvalue]
        slant range [km], longitude [deg] or east [km], depending on rasterObj.data content (respectively L1, intermediate or L2)
    ylims : list of numbers, optional
        plot y axis limits [minvalue, maxvalue]
        azimuth [km], latitude [deg] or north [km], depending on rasterObj.data content (respectively L1, intermediate or L2)


    Examples
    --------
    >>> rasterObj = BiomassL1cRaster(folder_path)
    >>> plot_rad2deg(rasterObj)
    >>> plt.show()

    >>> rasterObj = BiomassL2Raster(tif_path)
    >>> fig, ax = plt.subplots()
    >>> plot_rad2deg(rasterObj, ax=ax)
    >>> ax.set_xlabel("my xlabel")
    >>> fig.show()

    See Also
    --------
    plot : plot the data as it is
    plot_db : plot the data in dB
    plot_abs : plot the absolute value of data
    plot_angle : plot the angle (degrees) of complex data

    Notes
    -----
    plot_rad2deg, plots the data after conversion from radiants to degrees:

    .. math:: rasterObj.data*180/pi
    """

    if title is None:
        title = "[deg]"
    return __plot(rasterObj, np.rad2deg(rasterObj.data), ax, title, clims, xlims, ylims,)
