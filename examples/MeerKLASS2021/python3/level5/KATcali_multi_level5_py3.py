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
Tcmb=2.725
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

print ('start @ ' + time.asctime(time.localtime(time.time())) +'...')

print  (plt.rcParams['font.size'], plt.rcParams[u'axes.linewidth'],plt.rcParams['lines.linewidth'])
plt.rcParams['font.size'], plt.rcParams[u'axes.linewidth'],plt.rcParams['lines.linewidth'] = 14, 1.5, 1.5
#plt.rcParams['font.size'], plt.rcParams[u'axes.linewidth'],plt.rcParams['lines.linewidth'] = 10.0, 0.8, 1.5
print  (plt.rcParams['font.size'], plt.rcParams[u'axes.linewidth'],plt.rcParams['lines.linewidth'])

p_radec=np.loadtxt('radio_source2021.txt')

ch_plot=800

#input_file='/idia/projects/hi_im/raw_vis/MeerKLASS2021/level4/data/'
input_file='/scratch3/users/jywang/MeerKLASS2021/level4/data/'
output_file='/scratch3/users/jywang/MeerKLASS2021/level5/'

fname=sys.argv[1] #'1634252028'
ant=sys.argv[2] #'m001'

print (fname, ant)

d1=pickle.load(open(input_file+fname+'_'+ant+'_level4_data','rb'), encoding='latin-1')
print (d1.keys())

ra=d1['ra']
dec=d1['dec']
Tresi_map=d1['Tresi_map']
Tsky_map=d1['Tsky_map']
nd_s0=d1['nd_s0_h']

freqs=kio.cal_freqs(range(4096))

plot_gsize=90

#set the sky area to be pixelized 
x_cen=-18 #deg #RA
x_half=20 #deg

y_cen=-34#deg #DEC
y_half=11#deg

pix_deg=0.3

N_half_x=int(x_half/pix_deg)
N_half_y=int(y_half/pix_deg)
Npix_x=2*N_half_x+1
Npix_y=2*N_half_y+1

print (Npix_x,Npix_y)

p_choice='ZEA'

w = astropy.wcs.WCS(naxis=2)
w.wcs.crval = [x_cen-x_half, y_cen-y_half] # reference pointing of the image #deg
w.wcs.crpix = [1.0, 1.0] # pixel index corresponding to the reference pointing (try either 1 or 0 to see if the behaviour agrees to your expectation!)
w.wcs.cdelt = np.array([pix_deg, pix_deg]) # resolution 
#w.wcs.ctype = ['RA---ZEA', 'DEC--ZEA'] #projection 
#w.wcs.ctype = ['RA---AIT', 'DEC--AIT'] #ref_p0 can't go back to zero
w.wcs.ctype = ['RA---'+p_choice, 'DEC--'+p_choice] 

#other choice: https://docs.astropy.org/en/stable/wcs/
#skyview comment: https://skyview.gsfc.nasa.gov/current/help/fields.html#Projection

print (w)

##check the (min ra, min dec) of sky area will fall into pix (0,0)
p0=ac.SkyCoord(ra=(x_cen-x_half)*u.deg, dec=(y_cen-y_half)*u.deg) 
ref_p=skycoord_to_pixel(p0, w)
print (ref_p[0],ref_p[1])
assert(ref_p[0]<1e-12) #should be zero
assert(ref_p[1]<1e-12) #shoule be zero

p_list=ac.SkyCoord(ra=ra*u.deg, dec=dec*u.deg) #pointings in observation
x_pix_list,y_pix_list=skycoord_to_pixel(p_list,w) #observation (ra,dec) to pix

x_pix_list=np.round(x_pix_list).astype(int)
y_pix_list=np.round(y_pix_list).astype(int)

#out range due to track data, will filter when add to the map cube
print (np.min(x_pix_list),np.max(x_pix_list),Npix_x)
print (np.min(y_pix_list),np.max(y_pix_list),Npix_y)


fits_temp=np.zeros([Npix_x,Npix_y,4096])
plt.figure(figsize=(5,1.5))
plt.imshow(fits_temp[:,:,ch_plot].T,aspect='auto')
plt.colorbar()
plt.show()

Sum_Tsky_xy=fits_temp.copy()
Sum_Tresi_xy=fits_temp.copy()
Npix_xy_count1=fits_temp.copy()
Npix_xy_count2=fits_temp.copy()


assert((Tsky_map.mask==Tresi_map.mask).all()==True)

for i in range(len(ra)):
    x_pix,y_pix=x_pix_list[i],y_pix_list[i]
    
    mask1=Tsky_map.mask[i,:]
    if (mask1==True).all()==False:
        Sum_Tsky_xy[x_pix,y_pix,~mask1]+=Tsky_map[i,~mask1]
        Npix_xy_count1[x_pix,y_pix,~mask1]+=1
    
    mask2=Tresi_map.mask[i,:]
    if (mask2==True).all()==False:
        Sum_Tresi_xy[x_pix,y_pix,~mask2]+=Tresi_map[i,~mask2]
        Npix_xy_count2[x_pix,y_pix,~mask2]+=1

assert((Npix_xy_count1==Npix_xy_count2).all())

Tsky_xy=Sum_Tsky_xy/Npix_xy_count1
Tresi_xy=Sum_Tresi_xy/Npix_xy_count2

ptr_list=ac.SkyCoord(ra=p_radec[:,0]*u.deg, dec=p_radec[:,1]*u.deg)
ptr_ra_pix,ptr_dec_pix=skycoord_to_pixel(ptr_list,w)

plt.figure(figsize=(18,12))
plt.subplot(211,projection=w)
plt.imshow(Tsky_xy[:,:,ch_plot].T,cmap=kv.cmap1(),aspect='auto')
#plt.gca().invert_yaxis() #only works when no coordinates
plt.colorbar(label='Kelvin')
plt.ylabel('Dec',fontsize=16)
plt.xlim([0,Npix_x])
plt.ylim([0,Npix_y])
plt.grid()
plt.title('Pixelized ('+p_choice+') Tsky of '+fname+', '+ant+' (0.3deg/pixel)')
plt.subplot(212,projection=w)
plt.imshow(Tresi_xy[:,:,ch_plot].T,cmap=kv.cmap2(),aspect='auto')
#plt.gca().invert_yaxis() #only works when no coordinates
plt.plot(ptr_ra_pix,ptr_dec_pix,'mo')
plt.colorbar(label='Kelvin')
plt.xlabel('R. A.',fontsize=16)
plt.ylabel('Dec',fontsize=16)
plt.xlim([0,Npix_x])
plt.ylim([0,Npix_y])
plt.grid()
plt.title('Pixelized ('+p_choice+') Tres of '+fname+', '+ant+' (0.3deg/pixel)')
plt.savefig(output_file+'F_'+fname+'_'+ant+'_pix'+str(pix_deg)+'d_map.png', bbox_inches='tight')
#plt.show()

assert((Npix_xy_count1==Npix_xy_count2).all())


print (w)

list=[Sum_Tsky_xy,Sum_Tresi_xy,Npix_xy_count1,Tsky_xy,Tresi_xy]
list_str=['Sum_Tsky_xy','Sum_Tresi_xy','Npix_xy_count','Tsky_xy','Tresi_xy']
for i in range(len(list)):
    hdu=w.to_fits()
    hdu[0].data=list[i]
    hdu.writeto(output_file+fname+'_'+ant+'_'+list_str[i]+'_p'+str(pix_deg)+'d.fits', overwrite=True)
    
print ('end @ ' + time.asctime(time.localtime(time.time())) +'#')
print ('Hakuna Matata') 
################################################
#### Jingying Wang <astro.jywang@gmail.com> ####
################################################
