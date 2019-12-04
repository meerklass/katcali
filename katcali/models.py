import numpy as np
import pysm
from pysm.nominal import models
import healpy as hp
from astropy.coordinates import SkyCoord
from astropy.coordinates import *
from astropy import units as u
from scipy.interpolate import Rbf
import matplotlib.pylab as plt
from .calsource import CalSource
from .spillover import Spillover
from . import atmosphere

#-------------------------------------------------------------------------------
# Calibration point sources
#-------------------------------------------------------------------------------

calsrc_3C273 = CalSource("3C273",
                         ra=187.2779154, dec=2.0523883, frame='icrs', 
                         flux_408=55.1, flux_1410=42.)

calsrc_3C237 = CalSource("3C237",
                         ra=152.000125, dec=7.504541, frame='icrs', 
                         flux_408=15.4, flux_1410=6.6)

calsrc_PicA  = CalSource("PictorA",
                         ra=79.9571708, dec=-45.7788278, frame='icrs', 
                         flux_408=0.85, flux_1410=66.)

cal_sources = {
    '3C273':    calsrc_3C273,
    '3C237':    calsrc_3C237,
    'PictorA':  calsrc_PicA
}

#-------------------------------------------------------------------------------
# Spillover model
#-------------------------------------------------------------------------------

def cal_Tspill(el, pol, freqs, ch, version=2, 
               filename="/users/jywang/MeerKAT/model_test/mkat_model/spillover-models/mkat/MK_L_Tspill_AsBuilt_atm_mask.dat"):
    """
    Compatibility function for spillover model (single frequency).
    """
    print("cal_Tspill() is for single channel only! cal_Tspill_func has higher "
          "efficiency for multi channel calibration")
    
    # Fail if unsupported version requested
    if version==1:
        raise NotImplementedError("Spillover model version 1 has been removed")
    
    # Backwards compatibility
    if pol == 'h': pol = 'HH'
    if pol == 'v': pol = 'VV'
    
    # Load spillover model and evaluate
    spill = Spillover(filename=filename, mode='spline')
    return spill.temp(el, freqs[ch]/1e6, pol)[:,0]


def cal_Tspill_func(el, pol, freqs, filename='/users/jywang/MeerKAT/model_test/mkat_model/spillover-models/mkat/MK_L_Tspill_AsBuilt_atm_mask.dat'):
    """
    Compatibility function for spillover model.
    """
    # Backwards compatibility
    if pol == 'h': pol = 'HH'
    if pol == 'v': pol = 'VV'
    
    # Load spillover model
    spill = Spillover(filename=filename, mode='spline')
    
    # Wrap interpolation function for specific polarisation and return 
    # wrapper function
    return lambda _el, _freq: spill.temp(_el, _freq, pol)


#-------------------------------------------------------------------------------
# Atmosphere model
#-------------------------------------------------------------------------------

def calc_atmosphere_model_1ch(autodata, ch):
    """
    ???
    
    Parameters
    ----------
    autodata : ?
        ???
    
    ch : ?
        ???
    
    Returns
    -------
    atmos : array_like
        Atmopsheric temperature (in K) for each timestamp in autodata.timestamps.
    """
    # Function based on KAT-7 equivalent
    
    surface_temperature = autodata.temperature
    # This is some equation, I can't remember where
    T_atm = 1.12 * (273.15 + surface_temperature) - 50.0
    
    air_relative_humidity = autodata.humidity / 100.  # Percentage in katdal
    pressure = autodata.pressure
    
    assert len(autodata.ants) == 1, "More than one antenna found in autodata"
    
    # MeerKAT, height in kilometers above sea level
    height = autodata.ants[0].observer.elevation / 1000.
    frequency = autodata.freqs[ch] / 1e9  # frequency in GHz.
    
    # Build array of atmospheric temperature
    atmo_table = np.zeros_like(autodata.timestamps)
    for i in range(0, atmo_table.shape[0]):
        tau = atmosphere.opacity(surface_temperature[i], 
                                 air_relative_humidity[i], 
                                 pressure[i], 
                                 height, 
                                 frequency)
        atmo_v = T_atm[i] * (1 - np.exp(-tau / np.sin(np.radians(autodata.el[i]))))
        atmo_table[i] = atmo_v
    return atmo_table


#-------------------------------------------------------------------------------
# Galactic synchrotron model
#-------------------------------------------------------------------------------

def cal_Gal_model_np(vis, freqs, ra, dec, ch_min, ch_max, nside):
    sky_config = {
        'synchrotron': models("s1", nside),
    }
    sky = pysm.Sky(sky_config)
    c = SkyCoord(ra=ra * u.degree, dec=dec * u.degree, frame='icrs')
    theta = 90 - (c.galactic.b / u.degree).value
    phi = (c.galactic.l / u.degree).value
    result = np.zeros_like(vis)
    
    for i in range(ch_min, ch_max):
        syn = sky.synchrotron(nu=freqs[i]/1e9) / 1e6 # K
        I = hp.pixelfunc.get_interp_val(syn[0, :], theta/180*np.pi, phi/180*np.pi)
        result[:, i] = I / 2.
        
    return result


#-------------------------------------------------------------------------------
# Noise diode model
#-------------------------------------------------------------------------------

def call_Tnd(data, ant, pol, freqs, ch, plot_key=False):
    """
    
    Parameters
    ----------
    data
    
    ant
    
    pol
    
    freqs
    
    ch : 
    
    plot_key : bool, optional
        
    
    Returns
    -------
    Tnd_std : float
        Empirical standard deviation of noise diode data.
    
    Tnd_ref : 
        
    
    noise : 
        
    
    Tnd_spl : func
        Spline for noise diode temperature as a function of frequency [UNITS?].
    """
    print ("#cal_Tnd is for single channel only! Tnd_spl has higher efficiency for multi channel calibration")
    noise_id = data.receivers[ant]
    
    noise = {}
    for pol_i in ['h','v'] :  
        print(noise_id, pol_i)
    
        filename = '/users/jywang/MeerKAT/model_test/mkat_model/noise-diode-models/mkat/rx.%s.%s.csv'%(noise_id,pol_i)
        
        noise[pol_i] = np.loadtxt(filename,delimiter=',')
    
    x = noise[pol][:,0]/1e9
    y = noise[pol][:,1]
    Tnd_std = y.std()
    #print(Tnd_std)
    
    # Construct spline
    Tnd_spl = Rbf(x, y)
    Tnd_ref = Tnd_spl(freqs[ch] / 1e9)
    
    # Plot fits
    if plot_key == 1 or plot_key == True:
        plt.figure(figsize=(6,3))
        plt.plot(noise['h'][:,0]/1e9,noise['h'][:,1],'b.', ms=5)
        plt.plot(noise['v'][:,0]/1e9,noise['v'][:,1],'g.', ms=5)
        plt.plot(freqs/1e9,Tnd_spl(freqs/1e9),'m.',ms=.5) #Rbf line
        plt.plot(freqs[ch]/1e9,Tnd_ref,'ro',ms=6)
        plt.xlabel('freq (GHz)')
        plt.ylabel('Tnd (K)')
        plt.legend(['HH','VV'])
        plt.title(str(ant)+', noise id '+str(noise_id))
        #plt.savefig(str(fname)+'_nd_model_'+str(recv1)+'.pdf',bbox_inches='tight')
        plt.show()

    return Tnd_std, Tnd_ref, noise, Tnd_spl



def Tnd_spl(data, ant, pol, template='/users/jywang/MeerKAT/model_test/mkat_model/noise-diode-models/mkat/rx.%s.%s.csv'):
    """
    
    Parameters
    ----------
    data : 
        x
    
    ant : str
        Name of antenna to fetch noise diode data for.
    
    pol : str
        Polarisation ('h' or 'v').
    
    template : str, optional
        Filename template for noise diode model. Must have the form:
        "path/filename*{receiver}*{pol}*"
    
    Returns
    -------
    Tnd_std : float
        Empirical standard deviation of noise diode data.
    
    Tnd_spl : func
        Spline for noise diode temperature as a function of frequency [UNITS?].
    """
    pol = str(pol).lower()
    
    # 
    noise_id = data.receivers[ant]
    
    noise = {}
    for pol_i in ['h', 'v'] :  
        #print(noise_id, pol_i)
        filename = template % (noise_id, pol_i)
        noise[pol_i] = np.loadtxt(filename, delimiter=',')
    
    # Calculate standard deviation and spline fit to noise diode data
    x = noise[pol][:,0] / 1e9
    y = noise[pol][:,1]
    Tnd_std = y.std()
    Tnd_spl = Rbf(x, y)
    return Tnd_std, Tnd_spl

            
