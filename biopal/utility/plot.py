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
    
    plot: plot the data as it is
    
    See also: plot_db, plot_abs, plot_angle, plot_rad2deg

    Parameters
    ----------
        rasterObj: object of the data to plot, see BiomassL1cRaster, BiomassL2Raster classes
        ax:    (optional) pyplot axis, for custom settings
        title: (optional) string
        clims: (optional) list of plot saturation limits (unit depends on the type of plot) [minvalue, maxvalue]
        xlims: (optional) list of plot range axis limits in km [minvalue, maxvalue]
        ylims: (optional) list of plot azimuth axis limits in km [minvalue, maxvalue]

    Examples
    ----------
        rasterObj = BiomassL1cRaster(data_path)
        
        plot(rasterObj)
        plt.show()

        fig, ax = plt.subplots()
        plot(rasterObj, ax=ax)
        ax.set_xlabel("my xlabel")
        fig.show()
    """
    return __plot(rasterObj, rasterObj.data, ax, title, clims, xlims, ylims,)


def plot_db(rasterObj, ax=None, title=None, clims=None, xlims=None, ylims=None):

    """Plotting function (dB) for BiomassL1cRaster and BiomassL2Raster data objects
        
    plot_db: plot the data in dB (10*log10( abs(data)**2 )

    See also: plot, plot_abs, plot_angle, plot_rad2deg

    Parameters
    ----------
        rasterObj: object of the data to plot, see BiomassL1cRaster, BiomassL2Raster classes
        ax:    (optional) pyplot axis, for custom settings
        title: (optional) string
        clims: (optional) list of plot saturation limits (unit depends on the type of plot) [minvalue, maxvalue]
        xlims: (optional) list of plot range axis limits in km [minvalue, maxvalue]
        ylims: (optional) list of plot azimuth axis limits in km [minvalue, maxvalue]

    Examples
    ----------
        rasterObj = BiomassL1cRaster(data_path)

        plot_db(rasterObj)
        plt.show()
        
        fig, ax = plt.subplots()
        plot_db(rasterObj, ax=ax)
        ax.set_xlabel("my xlabel")
        fig.show()
    """
    if title is None:
        title = "[dB]"
    return __plot(rasterObj, 10 * np.log10(np.abs(rasterObj.data) ** 2), ax, title, clims, xlims, ylims,)


def plot_abs(rasterObj, ax=None, title=None, clims=None, xlims=None, ylims=None):

    """Plotting function (ABS) for BiomassL1cRaster and BiomassL2Raster data objects
    
    plot_abs: plot the absolute value of data ( abs(data) )

    See also: plot, plot_db, plot_angle, plot_rad2deg
    
    Parameters
    ----------
        rasterObj: object of the data to plot, see BiomassL1cRaster, BiomassL2Raster classes
        ax:    (optional) pyplot axis, for custom settings
        title: (optional) string
        clims: (optional) list of plot saturation limits (unit depends on the type of plot) [minvalue, maxvalue]
        xlims: (optional) list of plot range axis limits in km [minvalue, maxvalue]
        ylims: (optional) list of plot azimuth axis limits in km [minvalue, maxvalue]

    Examples
    ----------
        rasterObj = BiomassL1cRaster(data_path)

        plot_abs(rasterObj)
        plt.show()
        
        fig, ax = plt.subplots()
        plot_abs(rasterObj, ax=ax)
        ax.set_xlabel("my xlabel")
        fig.show()
    """
    if title is None:
        title = "ABS"
    return __plot(rasterObj, np.abs(rasterObj.data), ax, title, clims, xlims, ylims,)


def plot_angle(rasterObj, ax=None, title=None, clims=None, xlims=None, ylims=None):
    """Plotting function (angle) for BiomassL1cRaster and BiomassL2Raster data objects
    
    plot_angle: plot the angle (degrees) of complex data ( angle(data) )
    
    See also: plot, plot_db, plot_abs, plot_rad2deg

    Parameters
    ----------
        rasterObj: object of the data to plot, see BiomassL1cRaster, BiomassL2Raster classes
        ax:    (optional) pyplot axis, for custom settings
        title: (optional) string
        clims: (optional) list of plot saturation limits (unit depends on the type of plot) [minvalue, maxvalue]
        xlims: (optional) list of plot range axis limits in km [minvalue, maxvalue]
        ylims: (optional) list of plot azimuth axis limits in km [minvalue, maxvalue]

    Examples
    ----------
        rasterObj = BiomassL1cRaster(data_path)

        plot_angle(rasterObj)
        plt.show()
        
        fig, ax = plt.subplots()
        plot_angle(rasterObj, ax=ax)
        ax.set_xlabel("my xlabel")
        fig.show()
    """
    if title is None:
        title = "Angle [deg]"
    return __plot(rasterObj, np.rad2deg(np.angle(rasterObj.data)), ax, title, clims, xlims, ylims,)


def plot_rad2deg(rasterObj, ax=None, title=None, clims=None, xlims=None, ylims=None):

    """Plotting function (rad2deg) for BiomassL1cRaster and BiomassL2Raster data objects
    
    plot_rad2deg: plot the data after conversion from radiants to degrees ( rad2deg(data) )

    See also: plot, plot_db, plot_abs, plot_angle

    Parameters
    ----------
        rasterObj: object of the data to plot, see BiomassL1cRaster, BiomassL2Raster classes
        ax:    (optional) pyplot axis, for custom settings
        title: (optional) string
        clims: (optional) list of plot saturation limits (unit depends on the type of plot) [minvalue, maxvalue]
        xlims: (optional) list of plot range axis limits in km [minvalue, maxvalue]
        ylims: (optional) list of plot azimuth axis limits in km [minvalue, maxvalue]

    Examples
    ----------
        rasterObj = BiomassL1cRaster(data_path)

        plot_rad2deg(rasterObj)
        plt.show()
        
        fig, ax = plt.subplots()
        plot_rad2deg(rasterObj, ax=ax)
        ax.set_xlabel("my xlabel")
        fig.show()
    """
    if title is None:
        title = "[deg]"
    return __plot(rasterObj, np.rad2deg(rasterObj.data), ax, title, clims, xlims, ylims,)

