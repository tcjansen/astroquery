import os

import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.table import Table
from astropy.utils.data import download_file
from astropy.modeling.tabular import Tabular1D

__all__ = ['Filter', 'FilterObject']

data_path = os.path.join(os.path.dirname(__file__), 'data', 'filters.fits')


class Filter(object):
    """
    Astronomical filter object generator.
    """

    def __init__(self, path=data_path):
        """
        Parameters
        ----------
        path : str (optional)
            Path to `filters.fits` file
        """
        self.path = path
        self.table = Table(fits.getdata(path))
        self.table.add_index('col0')

    def available_filters(self):
        """
        Return a list of the available filters in the archive
        """
        return self.table['col0'].data

    def download(self, filter_id):
        """
        Query the SVO service for a given filter,
        return the transmittance curve.

        Parameters
        ----------
        filter_id : str
            Name of the filter. To see available filter IDs, run
            `~tynt.Filter.available_filters()`

        Returns
        -------
        filt : `~tynt.Filter`
            Astronomical filter object.
        """
        filter_path = download_file('http://svo2.cab.inta-csic.es/theory/'
                                    'fps3/fps.php?ID={0}'.format(filter_id))

        transmittance = Table.read(filter_path, format='votable')

        filt = FilterObject(transmittance['Wavelength'].data.data * u.Angstrom,
                            transmittance['Transmission'].data.data)
        return filt


class FilterObject(object):
    """
    Astronomical filter object.
    """
    def __init__(self, wavelength, transmittance, model=None):
        """
        Parameters
        ----------
        wavelength : `~numpy.ndarray`
            Wavelength array
        transmittance : `~numpy.ndarray`
            Transmittance array
        model : `~astropy.modeling.Model`
            Astropy model for the transmittance curve
        """
        self.wavelength = wavelength
        self.transmittance = transmittance
        self.model = model

    @property
    def table(self):
        return Tabular1D(points=self.wavelength,
                         lookup_table=self.transmittance)
