import numpy as np
import pickle
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.coordinates import EarthLocation
import astropy.coordinates as ac
from astropy.coordinates.angles import Angle
import astropy.time as at
import glob
import astropy.io.fits as pyfits
from .utils import KBOLTZ, CLIGHT
#from . import label_dump as kl


def beam_pattern_ptr(freq, flux, Aeffmax, Pn):
    """
    
    """
    flux_freq = flux(freq/1e9) * 1e-26
    Aeff = Aeffmax * Pn
    TA = 0.5 * flux_freq * Aeff/ KBOLTZ
    return TA
    
def gaussian_beam_ptr(freq, flux, Aeff, sep):
    """
    
    """
    flux_freq = flux(freq/1e9) * 1e-26
    TA = 0.5 * flux_freq * Aeff(freq, sep) / KBOLTZ # half
    return TA 

# BM
def calc_Omega_A_from_sigma(sigma_deg):
    return 2 * np.pi * np.radians(sigma_deg) ** 2

def calc_Aeff_max_from_sigma(sigma_deg, freq_Hz):
    lbd = CLIGHT / freq_Hz
    return lbd ** 2 / calc_Omega_A_from_sigma(sigma_deg)
    
def cal_Pn(pattern,timestamps,dp_ca,dp_cb,x_pix,y_pix):
    Pn=np.zeros_like(timestamps)
    for i in range(len(timestamps)):
        if i in dp_ca or i in dp_cb:
            Pn[i]=pattern[x_pix[i],y_pix[i]]
        else:
            Pn[i]=np.NaN
    return Pn

############ BM-I: beam from equation ####################
def cal_BMI(freqs, ch, flux_model, ang_deg, Ddish=14., eta=0.83):
    """
    
    Parameters
    ----------
    freqs : array_like
        xxx
    
    ch : int
        xxx
    
    flux_model : ???
        xxx
    
    ang_deg : ???
        xxx
    
    Ddish : float, optional
        Dish diameter, in m. Default: 14 (MeerKAT-like)
    
    eta : float, optional
        Aperture efficiency. Default: 0.83 (mean value from katconfig)
    
    Returns
    -------
    xxx
    """
    lbd = CLIGHT / freqs[ch] # wavelength
    Bsigma0 = np.degrees(0.45 * lbd / Ddish)
    Aeffmax0 = eta* calc_Aeff_max_from_sigma(Bsigma0, freqs[ch])
    
    Aeff0 = lambda freq, r: Aeffmax0*np.exp(-r**2./(2.*Bsigma0**2.)) #freq useless here
    T_ptr0 = gaussian_beam_ptr(freqs[ch], flux_model, Aeff0, ang_deg)
    return T_ptr0
    

"""
def load_Bdata(ch, beam_select):
    print ("#load_Bdata is for single channel only! load_Bdata_fband has higher efficiency for multi channel calibration")
    if ch < 1024:
        beam_file='p513_d5_ch4096/p1'
        ch_local=ch
    if ch >= 1024 and ch < 2048:
        beam_file='p513_d5_ch4096/p2'
        ch_local=ch-1024
    if ch >= 2048 and ch < 3072:
        beam_file='p513_d5_ch4096/p3'
        ch_local=ch-2048
    if ch >= 3072:
        beam_file='p513_d5_ch4096/p4'
        ch_local=ch-3072
    print beam_file
    Bdata=pickle.load(open('/users/jywang/MeerKAT/model_test/beam_model/eidos_sim/'+beam_file+'/beam_fit_'+beam_select+'_data','rb'))
    return Bdata,ch_local,beam_file
"""

def load_Bdata_fband(beam_select, template='/users/jywang/MeerKAT/model_test/beam_model/eidos_sim/p513_d5_ch4096/beam_fit_%s_p513_ch4096_d5_data'):
    """
    
    Parameters
    ----------
    
    Returns
    -------
    """
    Bdata = pickle.load(open(template % beam_select,'rb'))
    return Bdata


# BM-II: beam from gaussian fitting on Pattern
def cal_BMII(freqs, ch, pol, flux_model, ang_deg, beam_select):
    """
    
    Parameters
    ----------
    
    Returns
    -------
    """
    # FIXME: Single channel?
    #Bdata, ch_local, beam_file = load_Bdata(ch, beam_select)
    Bdata, ch_local, beam_file = load_Bdata_fband(beam_select)  
    grid_freq = Bdata['grid_freq']
    sigma_HH = np.array(Bdata['sigma_HH'])
    sigma_VV = np.array(Bdata['sigma_VV'])
    assert(grid_freq[ch_local]==freqs[ch]/1e6)
    
    if pol=='h':
        Bsigma1=sigma_HH[ch_local]
    elif pol=='v':
        Bsigma1=sigma_VV[ch_local]
    else:
        raise ValueError("Unknown polarisation '%s'" % pol)
        
    Aeffmax1 = calc_Aeff_max_from_sigma(Bsigma1, freqs[ch])
    Aeff1 = lambda freq, r: Aeffmax1 * np.exp(-r**2./(2.*Bsigma1**2.))
    T_ptr1 = gaussian_beam_ptr(freqs[ch], flux_model, Aeff1, ang_deg)
    return T_ptr1


"""
############BM-III: beam from pattern#######################################
def cal_BMIII(fname,data,ch,ant,pol,flux_model,c0, dp_ca,dp_cb,ang_deg,beam_select):
    print ("#cal_BMIII is for single channel only! cal_BMIII_1ch has higher efficiency for multi channel calibration")
    timestamps=data.timestamps
    freqs=data.freqs
    az=data.az[:,0]
    el=data.el[:,0]
    print np.shape(az),np.shape(el)
    Bdata,ch_local,beam_file= load_Bdata(ch,beam_select)

    grid_freq=Bdata['grid_freq']
    Aeff_max_HH=Bdata['Aeff_max_HH']
    Aeff_max_VV=Bdata['Aeff_max_VV']
    Aeff_max_HH=np.array(Aeff_max_HH)
    Aeff_max_VV=np.array(Aeff_max_VV)

    if pol=='h':
        Aeffmax2=Aeff_max_HH[ch_local]
    if pol=='v':
        Aeffmax2=Aeff_max_VV[ch_local]
    print Aeffmax2

    lon=Angle(data.ants[0].observer.lon,unit='rad')
    lat=Angle(data.ants[0].observer.lat, unit='rad')
    height=data.ants[0].observer.elevation
    ant_location=EarthLocation(lon=lon,lat=lat,height=height)

    radec_ptr = ac.SkyCoord(ra=c0.ra, dec=c0.dec) 
    altaz_frame = ac.AltAz(obstime=at.Time(timestamps * u.second, format='unix'),
                            location=ant_location)  # alt-az Coordinate System based on ant position
    altaz = radec_ptr.transform_to(altaz_frame)
    #print altaz
    ptr_el = (altaz.alt / u.deg).value
    ptr_az = (altaz.az / u.deg).value
    for i in range(len(ptr_az)):
        if ptr_az[i]>180:
            ptr_az[i]=ptr_az[i]-360

    y_sep=el-ptr_el#ptr_el-el####test
    azel=SkyCoord(az*u.deg, el*u.deg, frame='altaz')
    azel_ptr=SkyCoord(ptr_az*u.deg, ptr_el*u.deg, frame='altaz')
    r=azel.separation(azel_ptr).degree
    theta=np.arcsin(y_sep/r) 
    x_sep=r*np.cos(theta)
    assert(np.all([i > 0 for i in x_sep]))

    for i in range(len(ptr_az)):
        if ptr_az[i]-az[i]< 0: #only +/-, value is meanless!
            x_sep[i]=-1*x_sep[i]
        
    print beam_file
    Ddeg=5.0
    print Ddeg

    file1=glob.glob('/users/jywang/MeerKAT/model_test/beam_model/eidos_sim/'+beam_file+'/*'+beam_select+'*re*.fits')
    assert(len(file1)==1)
    file1=file1[0]
    print file1
    file2=glob.glob('/users/jywang/MeerKAT/model_test/beam_model/eidos_sim/'+beam_file+'/*'+beam_select+'*im*.fits')
    assert(len(file2)==1)
    file2=file2[0]
    print file2

    imgs1=pyfits.open(file1)[0].data
    imgs2=pyfits.open(file2)[0].data

    print np.shape(imgs1)
    Npix=np.shape(imgs1)[-1]
    print Npix

    grid_i=ch_local
    print grid_i
    img_HH=imgs1[grid_i,0,0,:,:]**2+imgs2[grid_i,0,0,:,:]**2 #squared!!!
    img_VV=imgs1[grid_i,1,1,:,:]**2+imgs2[grid_i,1,1,:,:]**2
    print img_HH.max(),img_VV.max()
    img_HH=img_HH/img_HH.max()
    img_VV=img_VV/img_VV.max()
    
    #in Khan's paper, smaller pix number (top in plots) is higher elevation
    img_HH=np.flip(img_HH,axis=0) #updown flip to make sure smaller pix number is lower elevation
    img_VV=np.flip(img_VV,axis=0)
    #img_HH=np.flip(img_HH,axis=1)#leftright
    #img_VV=np.flip(img_VV,axis=1)
    
    print img_HH.max(),img_VV.max()
    if pol=='h':
        pattern=img_HH 
    if pol=='v':
        pattern=img_VV

    x_pc=(Npix-1)/2.
    y_pc=(Npix-1)/2.
    print x_pc, y_pc

    print np.where(pattern==pattern.max())
    x_pix_max=np.where(pattern==pattern.max())[0][0]
    y_pix_max=np.where(pattern==pattern.max())[1][0]
    print x_pix_max,y_pix_max

    x_pix=x_pc+x_sep/Ddeg*Npix #relative on beam pattern
    
    y_pix=y_pc+y_sep/Ddeg*Npix
    
    assert(len(x_pix)==len(y_pix))
    for i in range(len(x_pix)):
        x_pix[i]=round(x_pix[i])
        y_pix[i]=round(y_pix[i])
    
    x_pix=x_pix.astype('int')
    y_pix=y_pix.astype('int')
        
    Pn=cal_Pn(pattern,timestamps,dp_ca,dp_cb,x_pix,y_pix)
    T_ptr2=beam_pattern_ptr(freqs[ch], flux_model, Aeffmax2, Pn)
    pix_label=x_pix,y_pix,x_pix_max,y_pix_max
    return T_ptr2,pattern,pix_label
"""


####split cal_BMIII to avoid repeat calucluation ################################################
def cal_pix_params(data,c0,Npix,Ddeg):
    timestamps=data.timestamps
    az=data.az[:,0]
    el=data.el[:,0]
    print(np.shape(az),np.shape(el))
    lon=Angle(data.ants[0].observer.lon,unit='rad')
    lat=Angle(data.ants[0].observer.lat, unit='rad')
    height=data.ants[0].observer.elevation
    ant_location=EarthLocation(lon=lon,lat=lat,height=height)

    radec_ptr = ac.SkyCoord(ra=c0.ra, dec=c0.dec) 
    altaz_frame = ac.AltAz(obstime=at.Time(timestamps * u.second, format='unix'),
                            location=ant_location)  # alt-az Coordinate System based on ant position
    altaz = radec_ptr.transform_to(altaz_frame)
    #print altaz
    ptr_el = (altaz.alt / u.deg).value
    ptr_az = (altaz.az / u.deg).value
    for i in range(len(ptr_az)):
        if ptr_az[i]>180:
            ptr_az[i]=ptr_az[i]-360

    y_sep=el-ptr_el#ptr_el-el####test
    azel=SkyCoord(az*u.deg, el*u.deg, frame='altaz')
    azel_ptr=SkyCoord(ptr_az*u.deg, ptr_el*u.deg, frame='altaz')
    r=azel.separation(azel_ptr).degree
    theta=np.arcsin(y_sep/r) 
    x_sep=r*np.cos(theta)
    assert(np.all([i > 0 for i in x_sep]))

    for i in range(len(ptr_az)):
        if ptr_az[i]-az[i]< 0: #only +/-, value is meanless!
            x_sep[i]=-1*x_sep[i]
        
    x_pc=(Npix-1)/2.
    y_pc=(Npix-1)/2.
    print(x_pc, y_pc)
    x_pix=x_pc+x_sep/Ddeg*Npix #relative on beam pattern
    
    y_pix=y_pc+y_sep/Ddeg*Npix
    
    assert(len(x_pix)==len(y_pix))
    for i in range(len(x_pix)):
        x_pix[i]=round(x_pix[i])
        y_pix[i]=round(y_pix[i])
    
    x_pix=x_pix.astype('int')
    y_pix=y_pix.astype('int')

    return x_pix,y_pix

def load_Aeff_max_fband(beam_select, pol):
    Bdata= load_Bdata_fband(beam_select)
    if pol=='h':
        Aeff_max_HH=Bdata['Aeff_max_HH']
        Aeff_max=np.array(Aeff_max_HH)
    if pol=='v':
        Aeff_max_VV=Bdata['Aeff_max_VV']
        Aeff_max=np.array(Aeff_max_VV)
    return  Aeff_max

def load_pattern_fband(beam_select, pol):
    if pol=='h':
        fname = glob.glob('/users/jywang/MeerKAT/model_test/beam_model/eidos_sim/p513_d5_ch4096/primary_beam_%s_p513_ch4096_d5_HH.fits' % beam_select)
    if pol=='v':
        fname = glob.glob('/users/jywang/MeerKAT/model_test/beam_model/eidos_sim/p513_d5_ch4096/primary_beam_%s_p513_ch4096_d5_VV.fits' % beam_select)
    
    assert len(fname)==1, "More than one file was found that matches the pattern."
    fname = fname[0]
    imgs = pyfits.open(fname)[0].data
    return imgs

def cal_BMIII_1ch(data, ch, flux_model, dp_ca, dp_cb, pattern_fband, 
                  x_pix, y_pix, Aeff_max_fband):
    timestamps = data.timestamps
    freqs = data.freqs
       
    Aeffmax2 = Aeff_max_fband[ch] #add in notebook
            
    pattern = pattern_fband[ch,:,:]
    pattern = pattern/pattern.max()
    
    #in Khan's paper, smaller pix number (top in plots) is higher elevation
    pattern = np.flip(pattern, axis=0) #updown flip to make sure smaller pix number is lower elevation
      
    print(np.where(pattern==pattern.max()))
    x_pix_max = np.where(pattern==pattern.max())[0][0]
    y_pix_max = np.where(pattern==pattern.max())[1][0]
    print(x_pix_max,y_pix_max)
      
    Pn = cal_Pn(pattern, timestamps, dp_ca, dp_cb, x_pix, y_pix)
    T_ptr2 = beam_pattern_ptr(freqs[ch], flux_model, Aeffmax2, Pn)
    return T_ptr2, pattern, x_pix_max, y_pix_max


