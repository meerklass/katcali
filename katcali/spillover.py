import numpy as np
import scikits.fitting as fit

class Spillover(object):
    
    def __init__(self, filename=None, mode='spline'):
        """
        Spillover model.
        
        Parameters
        ----------
        filename : str, optional
            Filename of the spillover model. Only used if mode = 'spline'. The 
            model file must have 3 columns:
                - theta, degrees (0 at zenith)
                - temperature, K
                - frequency, MHz
        
        mode : str, optional
            What kind of spillover model to build. Options are:
                - 'spline': A 2D spline using scikits.fitting.Spline2DGridFit
                - 'zero':   An all-zero model.
            Default: 'spline'.
        """
        # TODO Need to sort out better frequency interpolation & example
        if mode == 'zero':
            self.zero_spline()
        elif mode == 'spline':
            self.filename = filename
            self.spline_from_file(filename)
        else:
            raise ValueError("Unrecognised mode '%s'" % mode)
    
    
    def temp(self, el, freq, pol):
        """
        Evaluate spillover model for the HH polarisation, as a function of 
        elevation and frequency. Returns values in Kelvin.
        
        Parameters
        ----------
        el : array_like
            Array of elevation values (in degrees).
        
        freq : array_like
            Array of frequencies, in MHz.
        
        pol : str
            Which polarisation to evaluate the model for. Can be 'HH' or 'VV'.
        
        Returns
        -------
        spill : array_like
            Spillover temperature evaluated at each elevation/frequency, in K.
        """
        if pol.upper() == 'HH':
            return self.T_H(el, freq)
        elif pol.upper() == 'VV':
            return self.T_V(el, freq)
        else:
            ValueError("Unrecognised pol '%s'" % pol)
    
    
    def spline_from_file(self, filename):
        """
        Construct spillover model from values in a file.
        
        Parameters
        ----------
        filename : str
            Name of file to load data from. The model file must have 3 columns:
                - theta, degrees (0 at zenith)
                - temperature, K
                - frequency, MHz
        """
        # Models use convention (theta=0) == (el=90).
        # Load data from file
        datafile = np.loadtxt(filename)
        elevation = 90. - datafile[1:,0]
        numfreqs = (datafile.shape[1] - 1) // 2
        freqs = datafile[0,1::2]
        
        # HH and VV data
        spill_hh_data = datafile[1:, 1::2]
        spill_vv_data = datafile[1:, 2::2]
        
        # Fit splines to H and V polarisations
        self.T_H = fit.Spline2DGridFit(degree=(3,3))
        self.T_H.fit((elevation, freqs), spill_hh_data)
        
        self.T_V = fit.Spline2DGridFit(degree=(3,3))
        self.T_V.fit((elevation, freqs), spill_vv_data)
    
    
    def zero_spline(self):
        """
        Construct spillover model that is zero for all frequencies/elevations.
        """
        # Build array of elevations, spillover temperatures, and frequencies
        spillover = np.array([ [0.,90.,0.,90.],
                               [0.,0.,0.,0.],
                               [900.,900.,2000.,2000.] ])
        spillover[0] = 90. - spillover[0]
        
        # Construct splines on zero data
        self.T_H = fit.Delaunay2DScatterFit()
        self.T_H.fit(spillover[[0,2],:], spillover[1,:])
        self.T_V = self.T_H


