import numpy as np
from astropy.io import fits
from astropy.table import Table
from time import gmtime, strftime

import utils

__all__ = ['CreateDataProducts']

# Helper functions to ensure data product uniformity for the kronos program.
# These will be in the standard format required by MAST to be hosted as HLSPs

# Standard filename convention:
#    hlsp_proj-id_observatory_instrument_target_opt-elem_version_product-type.extension

class CreateDataProducts(object):

    def __init__(self, creator, creator_email, target, instrument,
                 element, pipeline, stage2_file, version=1):
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

        self.version  = int(version)

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

        primary = fits.PrimaryHDU()
        hdr = primary.header

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
                    hdr.set(keys[i], self.target, 'Target in this data product.')
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
                    hdr.set(keys[i], self.target, 'Planet in this data product.')
                if keys[i] == 'PIPELINE':
                    hdr.set(keys[i], self.pipeline, 'Reduction pipeline used.')
                if keys[i] == 'PROPOSID':
                    hdr.set(keys[i], 5959)
                if keys[i] == 'XPOSURE':
                    hdr.set(keys[i], 0)

        return primary

    def stellar_spectra_niriss(self, filename, flux_calibrated=False, overwrite=False):
        """
        Creates an h5py file for the stellar spectra from a given pipeline.

        Parameters
        ----------
        filename : str
           Name of the file to load data from.
        """
        if flux_calibrated == False:
            extname = 'uncalstellarspec'
        else:
            extname = 'fluxcalstellarspec'

        outputname = 'hlsp_kronos_jwst_{0}_{1}_{2}_v{3}_{4}.fits'.format(self.instrument.lower(),
                                                                     self.target.lower(),
                                                                     self.pipeline.lower(),
                                                                     self.version,
                                                                     extname)
        hdu = fits.open(filename)

        ext_list = [self.set_main_header()]

        for i in range(1, len(hdu)):
            datatype = hdu[i].header['EXTNAME'].split(' ')

            if (datatype[0] == 'Wave') and (len(datatype) == 2):
                extname = 'WAVELENGTH'
            elif (datatype[0] == 'Wave') and (len(datatype) > 2):
                extname = 'WAVELEGTH_ERROR'
            elif (datatype[0] == 'Flux') and (len(datatype) == 2):
                extname = 'FLUX'
            elif (datatype[0] == 'Flux') and (len(datatype) > 2):
                extname = 'FLUX_ERROR'
            elif datatype[0] == 'Time':
                extname = 'TIME'

            ext = fits.ImageHDU(data=hdu[i].data, name=extname)

            if i == 0:
                ext.header['METHOD']  = hdu[0].header['METHOD']
                ext.header['APWIDTH'] = hdu[0].header['WIDTH']

            if datatype[-1] == 'O1':
                ext.header['ORDER'] = 1
            elif datatype[-1] == 'O2':
                ext.header['ORDER'] = 2
            else:
                pass

            ext.header['UNITS'] = hdu[i].header['UNITS']
            ext.header['XTENSION'] = hdu[i].header['XTENSION']
            ext.header['NAXIS'] = hdu[i].data.shape[0]

            try:
                ext.header['NAXIS2'] = hdu[i].data.shape[1]
            except IndexError:
                pass

            ext_list.append(ext)

        hdulist = fits.HDUList(ext_list)
        hdulist.writeto(outputname, overwrite=overwrite)
        return


    def transmission_spec(self, files, resolutions, cases, overwrite=True):
        """
        Creates a text file for the transmission spectrum.

        Parameters
        ----------
        files : list, np.array
           List of files to open and include in the single FITS file. Files
           should be .txt files with four columns: wavelength, wavelength error,
           flux, and flux error. It is assumed the file is in this order.
        resolutions : list, np.array
           List of resolutions of the transmission spectra. Should be in the
           same order as the files.
        cases : list, array
           List of test cases of the transmission spectra. Should be in the same
           order as the files array.
        """
        files = np.array(files)
        cases = np.array(cases)
        resolutions = np.array(resolutions)

        ext_list = [self.set_main_header()]


        # Writes a FITS file for MAST
        outputname = 'hlsp_kronos_jwst_{0}_{1}_{2}_v{3}_transmission-spectra.fits'.format(self.instrument.lower(),
                                                                                          self.target.lower(),
                                                                                          self.pipeline.lower(),
                                                                                          self.version)

        # Sort files by Case #
        argsort = np.argsort(cases)

        for i in range(len(files)):
            dat = np.loadtxt(files[argsort][i])

            tab = Table()
            tab['wave'] = dat[:,0]
            tab['wave_err'] = dat[:,1]

            if np.nanmedian(dat[:,2]) < 1:
                depth = dat[:,2]
                depth_err = dat[:,3]
            else:
                print('Converting ppm to transit depth for consistent units.')
                depth = (dat[:,2] / 1e6)**2.0
                depth_err = (dat[:,3] / 1e6)**2.0

            tab['depth'] = depth
            tab['depth_err'] = depth_err

            ext = fits.BinTableHDU(data=tab, name='Transmission spectrum')

            ext.header['WAVEUNITS'] = 'micron'
            ext.header['DEPTH_UNITS'] = '(rp/rstar)^2'
            ext.header['RES'] = resolutions[argsort][i]
            ext.header['CASE'] = cases[argsort][i]

            if cases[argsort][i] == 0:
                ext.header['CASE_EXP'] = 'Spots are ignored'
            elif cases[argsort][i] == 1:
                ext.header['CASE_EXP'] = 'Spots are masked'
            elif cases[argsort][i] == 2:
                ext.header['CASE_EXP'] = 'Spots modeled as Gaussians'
            elif cases[argsort][i] == 3:
                ext.header['CASE_EXP'] = 'Spots modeled with fleck'
            elif cases[argsort][i] == 4:
                ext.header['CASE_EXP'] = 'Full Gaussian process fit'
            else:
                return('Case numbers must be between [0,4].')

            ext_list.append(ext)

        hdulist = fits.HDUList(ext_list)
        hdulist.writeto(outputname, overwrite=overwrite)

        return
