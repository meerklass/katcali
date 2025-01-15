#!/usr/bin/env python3

import numpy as np
import scipy.optimize as opt
import sys
import astropy.io.fits as pyfits
import math
import matplotlib.pylab as plt
import pickle
import time

def gauss_model(xi, yi, center_x, center_y, sigma):
    r2=(xi-center_x)**2.0+(yi-center_y)**2.0
    result= 1*np.exp(-r2/(2*sigma**2.0))
    return result

def func_obj(param, *args):
    sigma = param[0]
    center_x=param[1]
    center_y=param[2]
    xi=args[0]
    yi=args[1]
    data=args[2]
    
    model_value=gauss_model(xi, yi, center_x, center_y, sigma)
    result= np.sum((model_value-data)**2)
    #print(sigma, result)
    return result

print ('start @ ' + time.asctime(time.localtime(time.time())) +'#')

#d=np.load('MeerKAT_U_band_primary_beam.npz')
d=np.load('MeerKAT_U_band_primary_beam_aa_highres.npz')

print (d['pols'])
print (d['antnames'])
freqs_MHz=(d['freq_MHz'])
print (len(freqs_MHz))
print (np.shape(d['beam']))
#beam_ave=d['beam'][:,-1,:,:,:] #use the array_average
beam_ave=d['beam'][:,0,:,:,:] 
print ('use the array_average')
print (np.shape(beam_ave))
Npix=np.shape(beam_ave)[-1]
print (Npix)

Ddeg=12. #UHF
D_MeerKAT=13.5 #14.0 Khan preferred 14.0

r_cut=6.
pix_cut=r_cut/Ddeg*Npix
print (r_cut,pix_cut)

grid_freq=[]
sigma_HH=[]
sigma_VV=[]
Aeff_max_HH=[]
Aeff_max_VV=[]

for grid_i in range(len(freqs_MHz)):
      
    freq_local=freqs_MHz[grid_i]
    print (freq_local)
    
    img_HH=abs(beam_ave[0,grid_i,:,:])**2 #squared!!!
    img_VV=abs(beam_ave[3,grid_i,:,:])**2
    print (img_HH.max(),img_VV.max())
    img_HH=img_HH/img_HH.max()
    img_VV=img_VV/img_VV.max()
    print (img_HH.max(),img_VV.max())
    
    freq0=freqs_MHz[grid_i]*1e6 #freq*1e6
    sigma0=0.45*(3e8/freq0)/D_MeerKAT/np.pi*180.
    sigma0_pix=sigma0/Ddeg*Npix
    print (sigma0,sigma0_pix)
    
    for pol in ['h','v']:

        if pol=='h':
            img=img_HH
        if pol=='v':
            img=img_VV
        print (pol)
        
        center_x0=np.where(img[:,:]==img[:,:].max())[1][0]
        center_y0=np.where(img[:,:]==img[:,:].max())[0][0]
        print (center_x0,center_y0)
        #print (np.std(img_VV))
        
        xlist=[]
        ylist=[]
        zlist=[]
        imgsize=img.shape[0]
        img_cut=np.zeros([imgsize,imgsize])
        
        for i in range(0, imgsize):
            for j in range(0, imgsize):
                r=np.sqrt((j-center_x0)**2+(i-center_y0)**2)
                x=float(j)
                y=float(i)
                if r<=pix_cut:
                    img_cut[i,j]=img[i,j]
                    
                    xlist.append(x)
                    ylist.append(y)
                    zlist.append(img[i,j])
                   
        
        z=np.array(zlist)
        assert(z.max()==1)
        Omega_A=(Ddeg/Npix/180.*np.pi)**2*z.sum()
        light_speed = 2.99792485e8
        lbd = light_speed / (freqs_MHz[grid_i]*1e6)
        Aeff_max=lbd**2/Omega_A
        print (Aeff_max)
        
        fit_result=opt.fmin_powell(func=func_obj, x0=[sigma0_pix,center_x0,center_y0], args=(np.array(xlist), np.array(ylist), np.array(zlist)),xtol=1e-10,ftol=1e-10)
        sigma=fit_result[0]
        sigma_deg=sigma/Npix*Ddeg
        print (sigma_deg, sigma)
                
        sig_factor=sigma0/sigma_deg
        print ('# ', end='')
        print (grid_i, pol, Aeff_max, sigma_deg, sig_factor)
        
        if pol=='h':
            sigma_HH.append(sigma_deg)
            Aeff_max_HH.append(Aeff_max)
        if pol=='v':
            sigma_VV.append(sigma_deg)
            Aeff_max_VV.append(Aeff_max)
    grid_freq.append(freq_local)

data_p={}
data_p['grid_freq']=grid_freq
data_p['sigma_HH']=sigma_HH
data_p['sigma_VV']=sigma_VV
data_p['Aeff_max_HH']=Aeff_max_HH
data_p['Aeff_max_VV']=Aeff_max_VV
fs=open('beam_HiRes_fit_data','wb')
pickle.dump(data_p,fs,protocol=2)
fs.close()

print ('end @ ' + time.asctime(time.localtime(time.time())) +'#')
print ('Hakuna Matata')
################################################
#### Jingying Wang <astro.jywang@gmail.com> ####
################################################
