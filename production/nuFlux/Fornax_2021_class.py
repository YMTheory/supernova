import numpy as np
import logging
import os
import sys
import tarfile

import h5py
from astropy import units as u
from astropy.io import ascii, fits
from astropy.table import Table
from scipy.special import gamma, lpmv
from neutrino import Flavor
import re
from SupernovaModel import SupernovaModel

class Fornax_2021(SupernovaModel):
    """Model based on axisymmetric simulations from A. Burrows and D.  Vartanyan, Nature 589:29, 2021. Data available at https://www.astro.princeton.edu/~burrows/nu-emissions.2d/.
    """
    def __init__(self, filename):
        """
        Parameters
        ----------
        filename : str
            Absolute or relative path to FITS file with model data.
        """
        # Set up model metadata. 
        self.progenitor_mass = float(filename.split('/')[-1].split('_')[2][:-1]) * u.Msun
        self.metadata = {
                        'Progenitor mass':self.progenitor_mass,
                        }
        # Read time and flux information from HDF5 data file.
        _h5file = h5py.File(filename, 'r')

        self.time = _h5file['nu0'].attrs['time'] * u.s

        self.luminosity = {}
        self._E = {}
        self._dLdE = {}
        for flavor in Flavor:
            # Convert flavor to key name in the model HDF5 file
            key = {Flavor.NU_E: 'nu0',
                   Flavor.NU_E_BAR: 'nu1',
                   Flavor.NU_X: 'nu2',
                   Flavor.NU_X_BAR: 'nu2'}[flavor]

            self._E[flavor] = np.asarray(_h5file[key]['egroup'])
            self._dLdE[flavor] = {f"g{i}": np.asarray(_h5file[key][f'g{i}']) for i in range(12)}

            # Compute luminosity by integrating over model energy bins.
            dE = np.asarray(_h5file[key]['degroup'])
            n = len(dE[0])
            dLdE = np.zeros((len(self.time), n), dtype=float)
            for i in range(n):
                dLdE[:, i] = self._dLdE[flavor][f"g{i}"]

            # Note factor of 0.25 in nu_x and nu_x_bar.
            factor = 1. if flavor.is_electron else 0.25
            self.luminosity[flavor] = np.sum(dLdE*dE, axis=1) * factor * 1e50 * u.erg/u.s

    def get_initial_spectra(self, t, E, flavors=Flavor, interpolation='linear'):
        """Get neutrino spectra/luminosity curves after oscillation.
        Parameters
        ----------
        t : astropy.Quantity
            Time to evaluate initial spectra.
        E : astropy.Quantity or ndarray of astropy.Quantity
            Energies to evaluate the initial spectra.
        flavors: iterable of snewpy.neutrino.Flavor
            Return spectra for these flavors only (default: all)
        interpolation : str
            Scheme to interpolate in spectra ('nearest', 'linear').
        Returns
        -------
        initialspectra : dict
            Dictionary of model spectra, keyed by neutrino flavor.
        """
        initialspectra = {}

        # Avoid "division by zero" in retrieval of the spectrum.
        E[E == 0] = np.finfo(float).eps * E.unit
        logE = np.log10(E.to_value('MeV'))

        # Make sure the input time uses the same units as the model time grid.
        # Convert input time to a time index.
        t = t.to(self.time.unit)
        j = (np.abs(t - self.time)).argmin()

        for flavor in flavors:
            # Energy bin centers (in MeV)
            _E = self._E[flavor][j]
            _logE = np.log10(_E)
            _dlogE = np.diff(_logE)

            # Model flavors (internally) are nu_e, nu_e_bar, and nu_x, which stands
            # for nu_mu(_bar) and nu_tau(_bar), making the flux 4x higher than nu_e and nu_e_bar.
            factor = 1. if flavor.is_electron else 0.25

            # Linear interpolation in flux.
            if interpolation.lower() == 'linear':
                # Pad log(E) array with values where flux is fixed to zero.
                _logEbins = np.insert(_logE, 0, np.log10(np.finfo(float).eps))
                _logEbins = np.append(_logEbins, _logE[-1] + _dlogE[-1])

                # Luminosity spectrum _dLdE is in units of 1e50 erg/s/MeV.
                # Pad with values where flux is fixed to zero, then divide by E to get number luminosity
                _dNLdE = np.asarray([0.] + [self._dLdE[flavor]['g{}'.format(i)][j] for i in range(12)] + [0.])
                initialspectra[flavor] = (np.interp(logE, _logEbins, _dNLdE) / E * factor * 1e50 * u.erg/u.s/u.MeV).to('1 / (MeV s)')

            elif interpolation.lower() == 'nearest':
                # Find edges of energy bins and identify which energy bin (each entry of) E falls into
                _logEbinEdges = _logE - _dlogE[0] / 2
                _logEbinEdges = np.concatenate((_logEbinEdges, [_logE[-1] + _dlogE[-1] / 2]))
                _EbinEdges = 10**_logEbinEdges
                idx = np.searchsorted(_EbinEdges, E) - 1
                select = (idx > 0) & (idx < len(_E))

                # Divide luminosity spectrum by energy at bin center to get number luminosity spectrum
                _dNLdE = np.zeros(len(E))
                _dNLdE[np.where(select)] = np.asarray([self._dLdE[flavor]['g{}'.format(i)][j] / _E[i] for i in idx[select]])
                initialspectra[flavor] = ((_dNLdE << 1/u.MeV) * factor * 1e50 * u.erg/u.s/u.MeV).to('1 / (erg s)')

            else:
                raise ValueError('Unrecognized interpolation type "{}"'.format(interpolation))

        return initialspectra

    def get_averageE(self, t, E, flavors=Flavor, interpolation='linear'):
        initialspectra = {}
        averageE = {}

        # Avoid "division by zero" in retrieval of the spectrum.
        E[E == 0] = np.finfo(float).eps * E.unit
        logE = np.log10(E.to_value('MeV'))

        # Make sure the input time uses the same units as the model time grid.
        # Convert input time to a time index.
        t = t.to(self.time.unit)
        j = (np.abs(t - self.time)).argmin()

        for flavor in flavors:
            # Energy bin centers (in MeV)
            _E = self._E[flavor][j]
            _logE = np.log10(_E)
            _dlogE = np.diff(_logE)

            # Model flavors (internally) are nu_e, nu_e_bar, and nu_x, which stands
            # for nu_mu(_bar) and nu_tau(_bar), making the flux 4x higher than nu_e and nu_e_bar.
            factor = 1. if flavor.is_electron else 0.25

            # Linear interpolation in flux.
            if interpolation.lower() == 'linear':
                # Pad log(E) array with values where flux is fixed to zero.
                _logEbins = np.insert(_logE, 0, np.log10(np.finfo(float).eps))
                _logEbins = np.append(_logEbins, _logE[-1] + _dlogE[-1])

                # Luminosity spectrum _dLdE is in units of 1e50 erg/s/MeV.
                # Pad with values where flux is fixed to zero, then divide by E to get number luminosity
                _dNLdE = np.asarray([0.] + [self._dLdE[flavor]['g{}'.format(i)][j] for i in range(12)] + [0.])
                initialspectra[flavor] = (np.interp(logE, _logEbins, _dNLdE) / E * factor * 1e50 * u.erg/u.s/u.MeV).to('1 / (erg s)')
                averageE[flavor] = (np.sum(np.interp(logE, _logEbins, _dNLdE))) / np.sum((np.interp(logE, _logEbins, _dNLdE) / E)) 

            else:
                raise ValueError('Unrecognized interpolation type "{}"'.format(interpolation))

        return averageE
        
