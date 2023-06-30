#!/usr/bin/env python3
import matplotlib
matplotlib.use('AGG')

#imports
import katdal
import numpy as np
import matplotlib.pylab as plt
import astropy.coordinates as ac
import functools
import healpy as hp
import optparse
import warnings
from matplotlib.backends.backend_pdf import PdfPages
import healpy as hp
from astropy import units as u
from matplotlib.offsetbox import AnchoredText
import time
import pickle
import sys
import katcali
import katcali.visualizer as kv
import katcali.models as km
import katcali.rfi as kr
import katcali.solver as ks
import katcali.io as kio
import katcali.label_dump as kl
import katcali.diode as kd
import katcali.filter as kf
import katcali.beam as kb
from astropy.coordinates import SkyCoord
from astropy import units as u
from scipy.sparse import coo_matrix
import astropy.wcs
from astropy.wcs.utils import skycoord_to_pixel
from astropy.wcs.utils import pixel_to_skycoord
from astropy.io import fits
from astropy.wcs import WCS
from scipy.interpolate import UnivariateSpline
import matplotlib.colors as colors
Tcmb=2.725
k_B = 1.38E-23

print ('start @ ' + time.asctime(time.localtime(time.time())) +'...')

pix_deg=0.3 #0.15 #0.3
std_sigma=4. #the smaller number will make more data be delated 
total_count=961
i_iter=2
print (pix_deg, std_sigma, i_iter)
ch_ref=800
ch_plot=800
c_lim=40

output_file='../level3/'

#round1
#Fits='/idia/projects/hi_im/raw_vis/MeerKLASS2021/level6/ALL'+str(total_count)+'/sigma_'+str(int(std_sigma))+'/data/Nscan'+str(total_count)+'_Tsky_cube_p'+str(pix_deg)+'d_sigma'+str(std_sigma)+'_iter'+str(i_iter)+'.fits'
#Fits2='/idia/projects/hi_im/raw_vis/MeerKLASS2021/level6/ALL'+str(total_count)+'/sigma_'+str(int(std_sigma))+'/data/Nscan'+str(total_count)+'_Npix_count_cube_p'+str(pix_deg)+'d_sigma'+str(std_sigma)+'_iter'+str(i_iter)+'.fits'

#Fits='/idia/projects/hi_im/raw_vis/MeerKLASS2021/level6/'+str(pix_deg)+'/sigma_'+str(int(std_sigma))+'/data/Nscan'+str(total_count)+'_Tsky_cube_p'+str(pix_deg)+'d_sigma'+str(std_sigma)+'_iter'+str(i_iter)+'.fits'
#Fits2='/idia/projects/hi_im/raw_vis/MeerKLASS2021/level6/'+str(pix_deg)+'/sigma_'+str(int(std_sigma))+'/data/Nscan'+str(total_count)+'_Npix_count_cube_p'+str(pix_deg)+'d_sigma'+str(std_sigma)+'_iter'+str(i_iter)+'.fits'

#round2-3

Fits='/users/jywang/MeerKAT/calibration2021/level6/re_cali1/level6/'+str(pix_deg)+'/sigma'+str(int(std_sigma))+'_count'+str(c_lim)+'/re_cali1_round4/Nscan'+str(total_count)+'_Tsky_cube_p'+str(pix_deg)+'d_sigma'+str(std_sigma)+'_iter'+str(i_iter)+'.fits'

Fits2='/users/jywang/MeerKAT/calibration2021/level6/re_cali1/level6/'+str(pix_deg)+'/sigma'+str(int(std_sigma))+'_count'+str(c_lim)+'/re_cali1_round4/Nscan'+str(total_count)+'_Npix_count_cube_p'+str(pix_deg)+'d_sigma'+str(std_sigma)+'_iter'+str(i_iter)+'.fits'

print (Fits)
print (Fits2)

print ('# data cube loaded')

Tsky_cube0 = fits.open(Fits)[0].data
print (np.shape(Tsky_cube0))
count_cube0 = fits.open(Fits2)[0].data

Tsky_cube=np.ma.array(Tsky_cube0,mask=True)
Tsky_cube.mask[~np.isnan(Tsky_cube0)]=False
Tsky_cube.mask[count_cube0<c_lim]=True 

Npix_x= np.shape(Tsky_cube)[0]
Npix_y= np.shape(Tsky_cube)[1]

w=WCS(Fits).dropaxis(-1)
print (w)

pix_ra_map=np.zeros([Npix_x,Npix_y])
pix_dec_map=np.zeros([Npix_x,Npix_y])
print (np.shape(pix_ra_map),np.shape(pix_dec_map))

for pix_x in range(Npix_x):
    for pix_y in range(Npix_y):
        pix_radec=pixel_to_skycoord(pix_x,pix_y,w)
        pix_ra,pix_dec=(pix_radec.ra/u.deg).value,(pix_radec.dec/u.deg).value
        pix_ra_map[pix_x,pix_y]=pix_ra
        pix_dec_map[pix_x,pix_y]=pix_dec

pix_radec_map = SkyCoord(pix_ra_map*u.deg,  pix_dec_map*u.deg, frame='icrs')
print (np.shape(pix_radec_map)) 

fname=sys.argv[1]
ant=sys.argv[2]
pol='h' #no need to change
recv=ant+pol

data=kio.load_data(fname)
target,c0,bad_ants,flux_model=kio.check_ants(fname)
vis,flags= kio.call_vis(fname,recv)

print (ch_plot)

data.select(ants=ant,pol=pol)
ra,dec,az,el=kio.load_coordinates(data)
timestamps,freqs=kio.load_tf(data)
#dp_tt,dp_ss,dp_f,dp_t,dp_s=kl.label_dump_1ch(data,ant,pol,flags,ch_plot)
ang_deg=kio.load_ang_deg(ra,dec,c0)
dp_tt,dp_ss,dp_f,dp_w, dp_t,dp_s,dp_slew,dp_stop=kl.cal_dp_label(data,flags,ant,pol,ch_ref,ang_deg)
dp_sb=dp_ss[0]
dp_se=dp_ss[-1]

Tsky_block=np.ma.array(np.zeros_like(vis),mask=True)
print (np.shape(Tsky_block))

print (pix_deg)
for i in range(dp_sb,dp_se+1):
    p=SkyCoord(ra[i]*u.deg, dec[i]*u.deg, frame='icrs')
    
    p_ang=(p.separation(pix_radec_map)/u.deg)
    l=np.where(p_ang==p_ang.min())
    x,y=l[0][0],l[1][0]
    
    assert(p_ang[x,y]<pix_deg)
    Tsky_block[i,:]=Tsky_cube[x,y,:]

d={}
d['fname']=fname
d['ant']=ant
d['ra']=ra
d['dec']=dec
d['Tsky_block']=Tsky_block
fs=open(output_file+str(fname)+'_'+str(ant)+'_Tsky_block_data','wb')
pickle.dump(d,fs,protocol=2)
fs.close()

print (fname, recv, np.ma.mean(Tsky_block),np.ma.std(Tsky_block))

print ('end @ ' + time.asctime(time.localtime(time.time())) +'#')
print ('Hakuna Matata')
################################################
#### Jingying Wang <astro.jywang@gmail.com> ####
################################################
