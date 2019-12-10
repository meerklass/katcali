
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.coordinates import *
from astropy import units as u

class CalSource(object):
    
    def __init__(self, calname, ra, dec, frame='icrs', 
                 flux_408=None, flux_1410=None, alpha=None):
        """
        Parameters
        ----------
        calname : str
            Name of calibrator source.
        
        ra, dec : float
            RA and Dec of calibrator source, in degrees.
        
        frame : str, optional
            Frame to use for coordinates. Default: 'icrs'.
        
        flux_408, flux_1410 : float, optional
            Flux of source at 408 MHz and 1410 MHz respectively, in Jy.
            
            If alpha is not specified, both must be specified. Otherwise, only 
            flux_1410 must be specified to give an overall flux scale. 
            Default: None, None.
        
        alpha : float, optional
            Spectral index of frequency spectrum. Default: None.
        """
        # Store arguments
        self.calname = calname
        self.ra = ra
        self.de = dec
        
        # Build SkyCoord object for source position
        self.coords = SkyCoord(ra*u.deg,  dec*u.deg, frame=frame)
        
        # Construct SED parameters
        self.set_params(alpha=alpha, flux_408=flux_408, flux_1410=flux_1410)
        
    
    def set_params(self, alpha=None, flux_408=None, flux_1410=None):
        """
        Build a consistent set of SED parameters for this calibration source.
        
        Parameters
        ----------
        flux_408, flux_1410 : float, optional
            Flux of source at 408 MHz and 1410 MHz respectively, in Jy.
            
            If alpha is not specified, both must be specified. Otherwise, only 
            one must be specified to give an overall flux scale. 
            Default: None, None.
        
        alpha : float, optional
            Spectral index of frequency spectrum. Default: None.
        """
        assert flux_408 is not None or flux_1410 is not None, \
            "Either flux_410 or flux_1410 (or both) must be set."
        
        # Spectral index
        if alpha is not None:
            assert flux_408 is None, \
                "If alpha is specified, flux_408 must not be specified."
        else:
            alpha = - np.log10(flux_1410 / flux_408) / np.log10(1410. / 408.)
        self.alpha = alpha
        
        # Reference frequency and flux
        self.nu_ref = None
        self.flux_ref = None
        if flux_408 is not None and flux_1410 is None:
            self.nu_ref = 0.408
            self.flux_ref = flux_408
        if flux_408 is None and flux_1410 is not None:
            self.nu_ref = 1.410
            self.flux_ref = flux_1410
        if flux_408 is not None and flux_1410 is not None:
            self.nu_ref = 1.410
            self.flux_ref = flux_1410
        
    
    def flux_model(self, freq):
        """
        Return model for SED of calibration source as a function of frequency.
        
        Parameters
        ----------
        freq : array_like
            Array of frequencies, in GHz.
        
        Returns
        -------
        flux : array_like
            Flux of source as a function of frequency, in Jy.
        """
        return self.flux_ref * (freq / self.nu_ref)**(-alpha)
        
