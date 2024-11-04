import os
import sys
import numpy as np


def required_keys_all(header):
    """
    Determines which required keys for all FITS files exist or need to be set.

    Parameter
    ---------
    header : hdu.fits.Header

    Returns
    -------
    keys : np.array
       Array of keyword that need to be set in the header.
    set_key : np.array
       Boolean array of which keys need to be set by the user. True = must be set
       by the user. False = the keyword is already in the Stage 2 data product.
    """

    keys = np.sort(['DATE-BEG', 'DATE-END', 'DOI', 'EQUINOX', 'HLSPID', 'HLSPLEAD',
                     'HLSPNAME', 'HLSPTARG', 'HLSPVER', 'INSTRUME', 'LICENSE', 'LICENURL',
                     'MJD-BEG', 'MJD-END', 'MJD-MID', 'OBSERVAT', 'PROPOSID', 'TELESCOP',
                     'TIMESYS', 'XPOSURE', 'PLANET', 'PIPELINE', 'CREATOR', 'CONTACT',
                     'DATEMADE', 'NGROUPS', 'NFRAMES', 'NINTS', 'TARG_RA', 'TARG_DEC'
                     ])

    set_key = np.zeros(len(keys), dtype=bool)

    for i, k in enumerate(keys):
        try:
            header[k]
            set_key[i] = False
        except KeyError:
            set_key[i] = True

    return keys, set_key
