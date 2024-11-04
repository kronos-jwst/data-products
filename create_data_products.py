import numpy as np
from astropy.io import fits
from astropy.table import Table
from time import gmtime, strftime

import utils

__all__ = ['set_header', 'stellar_spec', 'white_light_curve',
           'spec_light_curves', 'transmission_spec']

# Helper functions to ensure data product uniformity for the kronos program.
# These will be in the standard format required by MAST to be hosted as HLSPs

# Standard filename convention:
#    hlsp_proj-id_observatory_instrument_target_opt-elem_version_product-type.extension

class CreateDataProducts(object):

    def __init__(self, creator, creator_email, target, instrument,
                 element, pipeline, stage2_file, version=1.0):
        """
        Initializes the CreateDataProducts object.

        Parameters
        ----------
        creator : str
           The name of the person who created the data product.
        creator_email : str
           The email of the person who created the data product.
        target : str
           The name of the planet. Spaces in the target name should be replaced
           with hyphens (e.g., V1298-Tau-c or TOI-451-d).
        instrument : str
           The name of the instrument. Should be either NIRISS or NIRSpec.
        element : str
           The grating used during the observations. Should be either SOSS or G395H.
        pipeline : str
           The name of the pipeline used to reduce the data.
        stage2_file : str
           The name of a FITS file created during Stage 2 of the JWST reduction
           pipeline. This will be used to copy the header information over.
        version : float, optional
           The version of the data product. Default is 1.0.
        """
        self.creator = creator
        self.creator_email = creator_email

        self.target = target.replace(" ", "-") # ensures no spaces in the target name

        if instrument != 'NIRISS' and instrument != 'NIRSpec':
            return('instrument should be "NIRISS" or "NIRSpec".')
        else:
            self.instrument = instrument

        if element != 'SOSS' and element != 'G395H':
            return('element should be "SOSS" or "G395H".')
        else:
            self.element = element

        self.pipeline = pipeline
        self.stage2_file = fits.open(stage2_file)

        self.version  = float(version)

        return


    def set_main_header(self):
        """
        Creates the main header for the FITS file.

        Parameters
        ----------


        Returns
        -------
        hdr : fits.HDU.Header
        """
        keys, setval = utils.required_keys_all(self.stage2_file[0].header)

        hdr = fits.Header()

        for i in range(len(keys)):
            if setval[i] == False:
                hdr.set(keys[i], self.stage2_file[0].header[keys[i]])
            else:
                if keys[i] == 'CONTACT':
                    hdr.set(keys[i], self.creator_email, 'Conact email for the person who produced the product.')
                if keys[i] == 'CREATOR':
                    hdr.set(keys[i], self.creator, 'Name of the person who produced the product.')
                if keys[i] == 'DATEMADE':
                    hdr.set(keys[i], strftime("%Y-%m-%d %H:%M:%S", gmtime()), 'Date the product was created.')
                if keys[i] == 'DOI':
                    hdr.set(keys[i], 'XXX.XXX')
                if keys[i] == 'EQUINOX':
                    hdr.set(keys[i], 0)
                if keys[i] == 'HLSPID':
                    hdr.set(keys[i], 0)
                if keys[i] == 'HLSPNAME':
                    hdr.set(keys[i], 'KRONOS: Keys to Revealing the Origin and Nature Of sub-neptune Systems')
                if keys[i] == 'HLSPLEAD':
                    hdr.set(keys[i], 'Adina Feinstein')
                if keys[i] == 'HLSPTARG':
                    hdr.set(keys[i], self.target, 'Target for this data product.')
                if keys[i] == 'HLSPVER':
                    hdr.set(keys[i], self.version, 'Data product version.')
                if keys[i] == 'LICENSE':
                    hdr.set(keys[i], 'MIT CC BY 4.0')
                if keys[i] == 'LICENURL':
                    hdr.set(keys[i], 'https://creativecommons.org/licenses/by/4.0/')
                if keys[i] == 'MJD-BEG':
                    hdr.set(keys[i], 0)
                if keys[i] == 'MJD-END':
                    hdr.set(keys[i], 0)
                if keys[i] == 'MJD-MID':
                    hdr.set(keys[i], 0)
                if keys[i] == 'OBSERVAT':
                    hdr.set(keys[i], 'JWST')
                if keys[i] == 'PLANET':
                    hdr.set(keys[i], self.target, 'Planet for this data product.')
                if keys[i] == 'PIPELINE':
                    hdr.set(keys[i], self.pipeline, 'Reduction pipeline used.')
                if keys[i] == 'PROPOSID':
                    hdr.set(keys[i], 5959)
                if keys[i] == 'XPOSURE':
                    hdr.set(keys[i], 0)

        return hdr


    def set_spectra_hdr(self, method, aperture=None):
        """
        Creates the header for the stellar spectra extension.

        Parameters
        ----------
        method : str
           How the spectra were extracted. Should be either 'box' or 'optimal'.
        aperture : float, optional
           The aperture of the box used to extract the stellar spectra. This is
           not applicable if the spectra were extracted via 'box' method. Default
           is None.
        """
        spec_hdr = fits.Header()

        """
        spec_hdr['METHOD']   =
        spec_hdr['APERTURE'] = aperture
        spec_hdr['XPOSURE']  =
        spec_hdr['BUNIT']    =
        spec_hdr['CDi_j']    =
        spec_hdr['CDELTi']   =
        spec_hdr['CRPIXj']   =
        spec_hdr['CRVALi']   =
        spec_hdr['CTYPEi']   =
        spec_hdr['PCi_j']    =
        spec_hdr['RADESYS']  =
        spec_hdr['WCSAXES']  =
        """
        return spec_hdr


    def set_whitelc_hdr(self, wave_start, wave_end):
        """
        Creates the header for the stellar spectra extension.

        Parameters
        ----------
        cenwave : float
           The central wavelength value.
        cenwave_err : float
           The error on the bin size of the central wavelength value.
        channel : int, optional
           The channel number for the spectroscopic light curve. Default is 0.
        """
        whitelc_hdr = fits.Header()

        whitelc_hdr['XTENSION'] = 1

        whitelc_hdr['FILT']   = 'WHITE'
        whitelc_hdr['CHANNEL'] = 0
        whitelc_hdr['WAVE-START'] = wave_start
        whitelc_hdr['WAVE-END'] = wave_end

        return speclc_hdr


    def set_speclc_hdr(self, cenwave, cenwave_err, channel=0):
        """
        Creates the header for the stellar spectra extension.

        Parameters
        ----------
        cenwave : float
           The central wavelength value.
        cenwave_err : float
           The error on the bin size of the central wavelength value.
        channel : int, optional
           The channel number for the spectroscopic light curve. Default is 0.
        """
        speclc_hdr = fits.Header()

        speclc_hdr['XTENSION'] = channel + 1

        speclc_hdr['FILT']   = 'MULTI'
        speclc_hdr['CHANNEL'] = channel

        speclc_hdr.set('CENWAVE', cenwave, 'microns')
        speclc_hdr.set('WAVEERR', cenwave_err, 'microns')
        speclc_hdr.set('LC-UNIT', 'normalized')

        return speclc_hdr


    def stellar_spec(self):
        """
        Creates an h5py file for the stellar spectra from a given pipeline.

        Parameters
        ----------
        """
        filename = 'hlsp_kronos_jwst_{0}_{1}_{2}_{3}_v{4}_stellarspec.fits'.format(self.instrument,
                                                                                   self.target,
                                                                                   self.element,
                                                                                   self.pipeline,
                                                                                   self.version)
        hdr = set_main_header()
        return


    def white_light_curve(self, pipeline):
        """
        Creates an h5py file for the white light curve from a given pipeline.

        Parameters
        ----------
        """
        filename = 'hlsp_kronos_jwst_{0}_{1}_{2}_{3}_v{4}_whitelc.fits'.format(self.instrument,
                                                                               self.target,
                                                                               self.element,
                                                                               self.pipeline,
                                                                               self.version)
        hdr = self.set_main_header()
        return


    def spec_light_curves(self, cenwave, cenwave_err, lc, lc_err, models, R):
        """
        Creates an h5py file for the spectroscopic light curves from a given pipeline.

        Parameters
        ----------
        cenwave : np.array
           Array of the central wavelengths for each spectroscopic channel.
        cenwave_err : np.array
           Array of half of the bin width for each spectroscopic channel.
        lc : np.ndarray
           2D array of the light curves for each spectroscopic channel.
        lc_err : np.ndarray
           2D array of the light curve errors for each spectroscopic channel.
        models : np.ndarray
           2D array of the best-fit light curve model for each spectroscopic
           channel.
        """
        filename = 'hlsp_kronos_jwst_{0}_{1}_{2}_{3}_v{4}_R{5}_speclc.fits'.format(self.instrument,
                                                                                   self.target,
                                                                                   self.element,
                                                                                   self.pipeline,
                                                                                   self.version,
                                                                                   R)
        dset = fits.HDUList()
        hdr = self.set_main_header()

        if (lc.shape != lc_err.shape) and (lc.shape != models.shape):
            return('lc, lc_err, and models should all be the same shape.')
        if len(cenwave) != len(lc):
            return('cenwave and lc should be the same length.')

        hdulist = [hdr]

        for i in range(len(lc)):
            shdr = self.set_speclc_hdr(cenwave[i], cenwave_err[i], i)

            tab = Table()
            tab['LC'] = lc[i]
            tab['LC-ERR'] = lc_err[i]
            tab['MODEL'] = models[i]

            tab_hdu = fits.BinTableHDU(tab, header=shdr)

            hdulist.append(tab_hdu)

        print(hdulist[40].header)
        return


    def transmission_spec(self, wave, wave_err, depth, depth_err, R):
        """
        Creates a text file for the transmission spectrum.

        Parameters
        ----------
        wave : np.array
           Array of the central wavelength for the transmission spectrum, in units
           of microns.
        wave_err : np.array
           Array of the half width of each wavelength bin for the transmission
           spectrum, in units of microns.
        depth : np.array
           Array of the measured transit depth as a function of wavelength, in units
           of parts-per-million (ppm).
        depth_err : np.array
           Array of the error of the measured transit depth as a function of
           wavelength, in units of parts-per-million (ppm).
        R : int
           The resolution of the transmission spectrum.
        pipeline : str
           The name of the pipeline used to reduce the data.
        """
        if np.nanmedian(depth) < 100:
            return('Transit depth should be given in units of ppm.')

        tab = Table()
        tab['wave'] = wave
        tab['wave_err'] = wave_err
        tab['depth'] = depth
        tab['depth_err'] = depth_err

        # Writes a CSV for the modelers
        tab.write('{0}_R{1}_spectrum.txt'.format(pipeline, R),
                  format='ascii')

        # Writes a FITS file for MAST
        filename = 'hlsp_kronos_jwst_{0}_{1}_{2}_{3}_v{4}_transmission.fits'.format(self.instrument,
                                                                                    self.target,
                                                                                    self.element,
                                                                                    self.pipeline,
                                                                                    self.version)

        return
