import numpy as np
from numpy.polynomial.legendre import legval
import os

KBOLTZ = 1.38e-23 # Boltzmann constant
CLIGHT = 2.99792485e8 # Speed of light, m/s

def file_size(filename):
    """
    Calculate the size of a file and return in human-readable units.
    
    Parameters
    ----------
    filename : str
        Path to file.
    
    Returns
    -------
    fsize : float
        Size of file in some units.
    
    units : str
        Units (e.g. MB, GB).
    """
    # Get size of file
    fsize = os.path.getsize(filename)
    
    # Convert to human-readable units
    for u in ['bytes', 'kB', 'MB', 'GB', 'TB', 'PB']:
        if fsize < 1024.:
            return fsize, u
        fsize /= 1024.
        

def log_normal(x, mu, sigma):
    return -((x - mu)**2. /(2* sigma**2.)) \
           - np.log(np.sqrt(2 * np.pi) * np.abs(sigma))


def legendre_timeseries(t, coeffs):
    """
    Return a Legendre polynomial as a function of time, where the time series 
    is rescaled to the domain [-1, 1]. 
    
    Parameters
    ----------
    t : array_like
        Array of times.
    
    coeffs : array_like
        Coefficients of Legendre polynomial up to some order.
    
    Returns
    -------
    poly : array_like
        Legendre polynomial evaluated at each time.
    """
    x = (t - t[0]) / (t[-1] - t[0]) * 2 - 1  # rescale time to [-1,1]
    return legval(x, coeffs)
    

def valid_filename(filename):
    """
    Return a standardised, valid form of the filename, or raise an error if it 
    is invalid. This function is needed because filenames are often integers, 
    but need to be passed in as strings. 
    
    Parameters
    ----------
    filename : str or int
        Filename to be validated and standardised.
    
    Returns
    -------
    filename : str
        Standardised filename as string.
    """
    # Check for valid filename types
    if isinstance(filename, str): return filename
    if isinstance(filename, (int, float, np.integer, np.float)):
        return "%s" % filename
    
    # Raise error if invalid type found
    raise TypeError("filename '%s' must be a string or integer" % filename)
