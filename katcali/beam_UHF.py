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
#from . import label_dump as kl
def beam_pattern_ptr(freq, flux, Aeffmax, Pn):
    flux_freq = flux(freq/1e9) * 1e-26
    k_B = 1.38E-23
      
    Aeff=Aeffmax*Pn
    TA = flux_freq / 2 * Aeff/ k_B #half
    return TA 
def gaussian_beam_ptr(freq, flux, Aeff, sep):
    flux_freq = flux(freq/1e9) * 1e-26
    k_B = 1.38E-23
    TA = flux_freq / 2 * Aeff(freq, sep) / k_B #half
    return TA 

#############BM###########################
def calc_Omega_A_from_sigma(sigma_deg):
    return 2 * np.pi * np.radians(sigma_deg) ** 2

def calc_Aeff_max_from_sigma(sigma_deg, freq_Hz):
    light_speed = 2.99792485e8
    lbd = light_speed / freq_Hz
    return lbd ** 2 / calc_Omega_A_from_sigma(sigma_deg)

def cal_Pn(pattern,timestamps,dp_ca,dp_cb,x_pix,y_pix):
    assert(pattern.ndim==2)
    Pn=np.zeros_like(timestamps)
    for i in range(len(timestamps)):
        if i in dp_ca or i in dp_cb:
            #Pn[i]=pattern[x_pix[i],y_pix[i]]
            Pn[i]=pattern[y_pix[i],x_pix[i]] #pattern is (el,az) so (y,x)
        else:
            Pn[i]=np.NaN
    return Pn

############ BM-I: beam from equation ####################
def cal_BMI(freqs,ch,flux_model,ang_deg):
    light_speed = 2.99792485e8
    lbd=light_speed/freqs[ch]
    #D_MeerKAT=14.0 #suggested by Khan
    D_MeerKAT=13.5 #for UHF
    Bsigma0=np.degrees(0.45*lbd/D_MeerKAT)
    #print Bsigma0
    eta=0.83 #mean value from katconfig 
    print ('*** eta=0.83 is used, not sure also works for UHF')
    Aeffmax0 = eta* calc_Aeff_max_from_sigma(Bsigma0, freqs[ch])
    print (Aeffmax0)
    Aeff0=lambda freq,r: Aeffmax0*np.exp(-r**2/(2*Bsigma0**2)) #freq useless here
    T_ptr0=gaussian_beam_ptr(freqs[ch], flux_model, Aeff0, ang_deg)
    return T_ptr0
###########################################################

def load_Bdata_fband(beam_select):
    if beam_select=='HiRes':
        Bdata=pickle.load(open('/users/jywang/MeerKAT/model_test/beam_model_UHF/beam_HiRes_spline_data','rb'),encoding='latin-1') 
        print ('# beam_HiRes_spline_data loaded')
    return Bdata

############BM-II: beam from gaussian fitting on Pattern ##############
def cal_BMII(freqs,ch,pol,flux_model,ang_deg,beam_select):
    
    #Bdata,ch_local,beam_file= load_Bdata(ch,beam_select)
    Bdata = load_Bdata_fband(beam_select)
    ch_local=ch
    grid_freqs=Bdata['freqs']
    sigma_HH=Bdata['sigma_HH']
    sigma_VV=Bdata['sigma_VV']
    sigma_HH=np.array(sigma_HH)
    sigma_VV=np.array(sigma_VV)
    #print np.shape(sigma_HH),np.shape(sigma_VV)
    #print grid_freq[ch_local],freqs[ch]/1e6
    assert(grid_freqs[ch]==freqs[ch]/1e6)
    
    if pol=='h':
        Bsigma1=sigma_HH[ch_local]
    if pol=='v':
        Bsigma1=sigma_VV[ch_local]
       
    #print Bsigma1
    Aeffmax1 = calc_Aeff_max_from_sigma(Bsigma1, freqs[ch])
    #print Aeffmax1
    Aeff1=lambda freq,r: Aeffmax1*np.exp(-r**2/(2*Bsigma1**2)) #freq useless here
    T_ptr1=gaussian_beam_ptr(freqs[ch], flux_model, Aeff1, ang_deg)
    return T_ptr1



############BM-III: beam from pattern#######################################

####split cal_BMIII to avoid repeat calucluation ################################################
def cal_pix_params(data,c0,Npix,Ddeg,key=0):
    timestamps=data.timestamps
    if key==0:
        az=data.az[:,0]
        el=data.el[:,0]
        print (np.shape(az),np.shape(el))
        lon=Angle(data.ants[0].observer.lon,unit='rad')
        lat=Angle(data.ants[0].observer.lat, unit='rad')
        height=data.ants[0].observer.elevation
    if key==-1:
        az=data.az
        el=data.el
        print (np.shape(az),np.shape(el))
        lon=Angle(data.ant_observer_lon,unit='rad')
        lat=Angle(data.ant_observer_lat, unit='rad')
        height=data.ant_observer_elevation
        
    ant_location=EarthLocation(lon=lon,lat=lat,height=height)

    radec_ptr = ac.SkyCoord(ra=c0.ra, dec=c0.dec) 
    altaz_frame = ac.AltAz(obstime=at.Time(timestamps * u.second, format='unix'),
                            location=ant_location)  # alt-az Coordinate System based on ant position
    altaz = radec_ptr.transform_to(altaz_frame)
    #print altaz
    ptr_el = (altaz.alt / u.deg).value
    ptr_az = (altaz.az / u.deg).value
    #print ('ptr az,el: ', ptr_az, ptr_el)
    ptr_az[ptr_az>180]-=360
    ####to avoid some issue in UHF data###
    dp_i=np.where(el>90)[0]
    for i in dp_i:
        print ('WARNNING: [dump'+str(i)+'] az and el are reset:', az[i], el[i],'-->',end=' ')
        az[i]+=180
        el[i]=180-el[i]
        print (az[i],el[i])
    ######################################
    az[az>180]-=360
    
    y_sep=ptr_el-el
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
    print (x_pc, y_pc)
    x_pix=x_pc+x_sep/Ddeg*Npix #relative on beam pattern
    
    y_pix=y_pc+y_sep/Ddeg*Npix
    
    assert(len(x_pix)==len(y_pix))
    for i in range(len(x_pix)):
        x_pix[i]=round(x_pix[i])
        y_pix[i]=round(y_pix[i])
    
    x_pix=x_pix.astype('int')
    y_pix=y_pix.astype('int')

    return x_pix,y_pix

def load_Aeff_max_fband(beam_select,pol):
    Bdata= load_Bdata_fband(beam_select)
    if pol=='h':
        Aeff_max_HH=Bdata['Aeff_max_HH']
        Aeff_max=np.array(Aeff_max_HH)
    if pol=='v':
        Aeff_max_VV=Bdata['Aeff_max_VV']
        Aeff_max=np.array(Aeff_max_VV)
    return  Aeff_max

def load_pattern_fband(beam_select,pol):
    if pol=='h':
        d=np.load('/users/jywang/MeerKAT/model_test/beam_model_UHF/data_saved/HH_pattern_spline_data.npz')
        imgs=d['HH_pattern_fband']
        print ('# HH_pattern_spline_data.npz loaded')
    if pol=='v':
        d=np.load('/users/jywang/MeerKAT/model_test/beam_model_UHF/data_saved/VV_pattern_spline_data.npz')
        imgs=d['VV_pattern_fband']
        print ('# VV_pattern_spline_data.npz loaded')

    print (np.shape(imgs))
    for i in range(np.shape(imgs)[0]):
        imgs[i,:,:]/=imgs[i,:,:].max()
        
    return imgs

def cal_BMIII_1ch(data,ch,flux_model, dp_ca,dp_cb,pattern_fband,x_pix,y_pix,Aeff_max_fband):
    assert(np.array(ch).ndim==0)
    timestamps=data.timestamps
    freqs=data.freqs
       
    Aeffmax2=Aeff_max_fband[ch] #add in notebook
    print (Aeffmax2)
            
    pattern=pattern_fband[ch,:,:]
    pattern=pattern/pattern.max()
    
    #in Khan's paper, smaller pix number (top in plots) is higher elevation
    ### pattern=np.flip(pattern,axis=0) #updown flip to make sure smaller pix number is lower elevation @modified at 2025-01-13 for UHF-band
    print ('### no pattern flip for UHF-band')  
    print (np.where(pattern==pattern.max()))
    #x_pix_max=np.where(pattern==pattern.max())[0][0]
    #y_pix_max=np.where(pattern==pattern.max())[1][0]
    x_pix_max=np.where(pattern==pattern.max())[1][0] #pattern is [el,az]
    y_pix_max=np.where(pattern==pattern.max())[0][0]
    print (x_pix_max,y_pix_max)
      
    Pn=cal_Pn(pattern,timestamps,dp_ca,dp_cb,x_pix,y_pix)
    T_ptr2=beam_pattern_ptr(freqs[ch], flux_model, Aeffmax2, Pn)
    return T_ptr2,pattern,x_pix_max,y_pix_max

############END of split###################################################

