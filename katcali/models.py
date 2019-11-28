import numpy as np
import pysm
from pysm.nominal import models
import healpy as hp
from astropy.coordinates import SkyCoord
from astropy.coordinates import *
from astropy import units as u
from scipy.interpolate import Rbf
import matplotlib.pylab as plt
from calsource import CalSource
from spillover import Spillover

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
    
def calc_atmospheric_opacity(T, RH, P, h, f):
    """
    Calculates zenith opacity according to ITU-R P.676-9. This is for 
    elevations > 10 deg. Use as "Tsky * (1 - exp(-opacity/sin(el)))" for 
    elevation dependence. (N.B. Same as used for KAT-7.)
    
    Parameters
    ----------
    temp : float
        Temperature in deg C.
        
    rel_hum : float
        Relative humidity, 0 < RH < 1.
    
    pressure : float
        Dry air pressure in hPa (equiv. mbar).
        
    height : float
        Height above sea level in km.
    
    freq : float
        Frequency in GHz (must be < 55 GHz).
    
    Returns
    -------
    opac : float
        Approximate atmospheric opacity at zenith [Nepers].
    """
    # [hPa] from A. L. Buck research manual 1996
    es = 6.1121*np.exp((18.678-T/234.5) * T / (257.14+T))
    
    # [g/m^3] from A. L. Buck research manual 1996
    # (ITU-R omitted the factor "RH" - a mistake)
    rho = RH*es*216.7/(T+273.15)

    # The following is taken directly from ITU-R P.676-9
    p_tot = P + es # from eq 3

    rho = rho*np.exp(h/2) # Adjust to sea level as per eq 32

    # eq 22
    r_t = 288./(273.+T)
    r_p = p_tot/1013.
    phi = lambda a, b, c, d: r_p**a*r_t**b*np.exp(c*(1-r_p)+d*(1-r_t))
    E_1 = phi(0.0717,-1.8132,0.0156,-1.6515)
    E_2 = phi(0.5146,-4.6368,-0.1921,-5.7416)
    E_3 = phi(0.3414,-6.5851,0.2130,-8.5854)
    
    # Following is valid only for f <= 54 GHz
    yo = ( 7.2*r_t**2.8 / (f**2+0.34*r_p**2*r_t**1.6) + 0.62*E_3 / ((54-f)**(1.16*E_1)+0.83*E_2) ) * f**2 * r_p**2 *1e-3
    
    # eq 23
    n_1 = 0.955*r_p*r_t**0.68 + 0.006*rho
    n_2 = 0.735*r_p*r_t**0.5 + 0.0353*r_t**4*rho
    g = lambda f, f_i: 1+(f-f_i)**2/(f+f_i)**2
    yw = (  3.98*n_1*np.exp(2.23*(1-r_t))/((f-22.235)**2+9.42*n_1**2)*g(f,22) + 11.96*n_1*np.exp(0.7*(1-r_t))/((f-183.31)**2+11.14*n_1**2)
          + 0.081*n_1*np.exp(6.44*(1-r_t))/((f-321.226)**2+6.29*n_1**2) + 3.66*n_1*np.exp(1.6*(1-r_t))/((f-325.153)**2+9.22*n_1**2)
          + 25.37*n_1*np.exp(1.09*(1-r_t))/(f-380)**2 + 17.4*n_1*np.exp(1.46*(1-r_t))/(f-448)**2
          + 844.6*n_1*np.exp(0.17*(1-r_t))/(f-557)**2*g(f,557) + 290*n_1*np.exp(0.41*(1-r_t))/(f-752)**2*g(f,752)
          + 8.3328e4*n_2*np.exp(0.99*(1-r_t))/(f-1780)**2*g(f,1780)
          ) * f**2*r_t**2.5*rho*1e-4

    # eq 25
    t_1 = 4.64/(1+0.066*r_p**-2.3) * np.exp(-((f-59.7)/(2.87+12.4*np.exp(-7.9*r_p)))**2)
    t_2 = 0.14*np.exp(2.12*r_p) / ((f-118.75)**2+0.031*np.exp(2.2*r_p))
    t_3 = 0.0114/(1+0.14*r_p**-2.6) * f * (-0.0247+0.0001*f+1.61e-6*f**2) / (1-0.0169*f+4.1e-5*f**2+3.2e-7*f**3)
    ho = 6.1/(1+0.17*r_p**-1.1)*(1+t_1+t_2+t_3)

    # eq 26
    sigma_w = 1.013/(1+np.exp(-8.6*(r_p-0.57)))
    hw = 1.66*( 1 + 1.39*sigma_w/((f-22.235)**2+2.56*sigma_w) + 3.37*sigma_w/((f-183.31)**2+4.69*sigma_w) + 1.58*sigma_w/((f-325.1)**2+2.89*sigma_w) )

    # Attenuation from dry & wet atmosphere relative to a point outside of the 
    # atmosphere
    A = yo*ho*np.exp(-h/ho) + yw*hw*np.exp(-h/hw) # [dB] from equations 27, 30 & 31

    return A*np.log(10)/10.0 # Convert dB to Nepers


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
        tau = calc_atmospheric_opacity(surface_temperature[i], 
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

def call_Tnd(data, ant,pol,freqs,ch,plot_key):
    print ("#cal_Tnd is for single channel only! Tnd_spl has higher efficiency for multi channel calibration")
    noise_id=data.receivers[ant]
    
    noise={}
    for pol_i in ['h','v'] :  
        print(noise_id, pol_i)
    
        filename = '/users/jywang/MeerKAT/model_test/mkat_model/noise-diode-models/mkat/rx.%s.%s.csv'%(noise_id,pol_i)
        
        noise[pol_i] = np.loadtxt(filename,delimiter=',')
    x=noise[pol][:,0]/1e9
    y=noise[pol][:,1]
    Tnd_std=y.std()
    print(Tnd_std)
    Tnd_spl = Rbf(x, y)
    Tnd_ref=Tnd_spl(freqs[ch]/1e9)
    if plot_key==1:
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

    return Tnd_std,Tnd_ref,noise,Tnd_spl



def Tnd_spl(data, ant,pol):
    noise_id=data.receivers[ant]
    
    noise={}
    for pol_i in ['h','v'] :  
        print(noise_id, pol_i)
    
        filename = '/users/jywang/MeerKAT/model_test/mkat_model/noise-diode-models/mkat/rx.%s.%s.csv'%(noise_id,pol_i)
        
        noise[pol_i] = np.loadtxt(filename,delimiter=',')
    x=noise[pol][:,0]/1e9
    y=noise[pol][:,1]
    Tnd_std=y.std()
    #print Tnd_std
    Tnd_spl = Rbf(x, y)
    return Tnd_std,Tnd_spl

            
