#!/usr/bin/env python2
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
from matplotlib.gridspec import GridSpec

print 'start @ ' + time.asctime(time.localtime(time.time())) +'#'

print katcali.__version__

print  plt.rcParams['font.size'], plt.rcParams[u'axes.linewidth'],plt.rcParams['lines.linewidth']
plt.rcParams['font.size'], plt.rcParams[u'axes.linewidth'],plt.rcParams['lines.linewidth'] = 14, 1.5, 1.5
#plt.rcParams['font.size'], plt.rcParams[u'axes.linewidth'],plt.rcParams['lines.linewidth'] = 10.0, 0.8, 1.5
print  plt.rcParams['font.size'], plt.rcParams[u'axes.linewidth'],plt.rcParams['lines.linewidth']

input_file='../level3/level3/'
output_file='./level4/'

def cal_map_I(map_h, map_v):    
    assert(np.shape(map_h)==np.shape(map_v))
    map=(map_h+map_v)/2.   
    print 'min value h,v:',
    print np.ma.min(map_h),np.ma.min(map_v)
    assert(np.ma.min(map_h)>0)
    assert(np.ma.min(map_v)>0)
    diff_ratio=np.ma.exp(abs(np.ma.log(map_h)-np.ma.log(map_v))) #difference betwwen map h,v; 
                                                        #function to make sure maps are exchangable
    return map, diff_ratio



fname=sys.argv[1]
ant=sys.argv[2]

print fname,ant

try:
    dh=pickle.load(open(input_file+fname+'_'+ant+'h_level3_data'))
    dv=pickle.load(open(input_file+fname+'_'+ant+'v_level3_data'))

    assert((dh['ra']==dv['ra']).all()==True)
    assert((dh['dec']==dv['dec']).all()==True)
    assert((dh['timestamps']==dv['timestamps']).all()==True)
    assert((dh['Tgal_map']==dv['Tgal_map']).all()==True)
    ra=dh['ra']
    dec=dh['dec']
    timestamps=dh['timestamps']
    Tgal_map=dh['Tgal_map']

    gain_map_h=dh['gain_map']
    T_map_h=dh['T_map']
    Tresi_map_h=dh['Tresi_map']
    Tsm_map_h=dh['Tsm_map']
    nd_s0_h=dh['nd_s0'] 
    Tel_map_h=dh['Tel_map']
    Tsky_map_h=T_map_h-Tel_map_h-Tsm_map_h

    gain_map_v=dv['gain_map']
    T_map_v=dv['T_map']
    Tresi_map_v=dv['Tresi_map']
    Tsm_map_v=dv['Tsm_map']
    nd_s0_v=dv['nd_s0'] 
    Tel_map_v=dv['Tel_map']
    Tsky_map_v=T_map_v-Tel_map_v-Tsm_map_v


    freqs=kio.cal_freqs(range(4096))

    #ch_e=1051

    assert((T_map_h.mask==Tresi_map_h.mask).all()==True)
    assert((T_map_v.mask==Tresi_map_v.mask).all()==True)
    div=T_map_h/T_map_v
    mask_clean=div.mask.copy()
    mask_clean_backup=mask_clean.copy()

    l=[T_map_h,T_map_v,Tresi_map_h,Tresi_map_v,Tsm_map_h,Tsm_map_v,gain_map_h,gain_map_v]
    l_str=[r'$T_{cal}(\nu)$'+', pol=HH',r'$T_{cal}(\nu)$'+', pol=VV',r'$T_{res}(\nu)$'+', pol=HH',r'$T_{res}(\nu)$'+', pol=VV',
           r'$T_{rec}(\nu)$'+', pol=HH',r'$T_{rec}(\nu)$'+', pol=VV','gain'+r'$(\nu)$'+', pol=HH','gain'+r'$(\nu)$'+', pol=VV']

    del_point0,del_point1,del_point2,del_point3,del_point4,del_point5,del_point6,del_point7=[],[],[],[],[],[],[],[]


    #####filter######
    niter=6

    for i in range(niter):
        print '> The iteration '+str(i+1)+' is in progress...'
        mask_clean2=mask_clean.copy()

        del_point=[]

        for l_i in range(len(l)):
            a=l[l_i]
            a_ch=np.ma.mean(np.ma.array(a,mask=mask_clean),axis=0)

            ax1=kf.curve_filter_ma_array(range(550,1051),a_ch[550:1051],sigma=4.)#550-1050 filter
            mask_clean2[:,ax1]=True ###apply to the mask
            del_point+=list(ax1)
            '''
            ax2=kf.curve_filter_ma_array(range(2150,3101),a_ch[2150:3101],sigma=4.)#2150-3100 filter
            mask_clean2[:,ax2]=True ###apply to the mask
            del_point+=list(ax2)

            ax3=kf.curve_filter_ma_array(range(550,3101),a_ch[550:3101],sigma=4.,k=1)
            mask_clean2[:,ax3]=True ###apply to the mask
            del_point+=list(ax3)
            '''
            if l_i==0:
                del_point0+=list(ax1)#+list(ax2)+list(ax3)
            if l_i==1:
                del_point1+=list(ax1)#+list(ax2)+list(ax3)
            if l_i==2:
                del_point2+=list(ax1)#+list(ax2)+list(ax3)
            if l_i==3:
                del_point3+=list(ax1)#+list(ax2)+list(ax3)
            if l_i==4:
                del_point4+=list(ax1)#+list(ax2)+list(ax3)
            if l_i==5:
                del_point5+=list(ax1)#+list(ax2)+list(ax3)
            if l_i==6:
                del_point6+=list(ax1)#+list(ax2)+list(ax3)
            if l_i==7:
                del_point7+=list(ax1)#+list(ax2)+list(ax3)

        print 'masks before and after filter are same:',
        print (mask_clean==mask_clean2).all()
        mask_clean=mask_clean2.copy()

        print str(len(set(del_point))) + ' channels masked by '+str(len(del_point))+' times'
        if len(del_point)==0:
            print '***filter done after iter time '+str(i+1) 
            break


    print del_point0
    print del_point1
    print del_point2
    print del_point3
    print del_point4
    print del_point5
    print del_point6
    print del_point7

    del_point_all=del_point0+del_point1+del_point2+del_point3+del_point4+del_point5+del_point6+del_point7
    del_point_all=list(set(del_point_all))
    print len(del_point_all)

    plt.figure(figsize=(14,80))
    plt.subplots_adjust (wspace=0.2, hspace=0.3) 

    for l_i in range(len(l)):
        a=l[l_i]
        plt.subplot(len(l)/2*niter,2,i*len(l)+l_i+1)


        if l_i==0:
            del_point_plot=del_point0
        if l_i==1:
            del_point_plot=del_point1
        if l_i==2:
            del_point_plot=del_point2
        if l_i==3:
            del_point_plot=del_point3
        if l_i==4:
            del_point_plot=del_point4
        if l_i==5:
            del_point_plot=del_point5
        if l_i==6:
            del_point_plot=del_point6
        if l_i==7:
            del_point_plot=del_point7

        a_ch=np.ma.mean(a,axis=0) #raw
        plt.plot(freqs[del_point_all]/1e6,a_ch[del_point_all],'.',color='gold',ms=6) #all deleted data
        plt.plot(freqs[del_point_plot]/1e6,a_ch[del_point_plot],'.',color='r',ms=6) #deleted data

        a_ch1=np.ma.mean(np.ma.array(a,mask=mask_clean),axis=0) #clean 
        plt.plot(freqs/1e6,a_ch1,'.', ms=4) #clean data

        if l_i !=2 and l_i!=3:
            plt.ylim(a_ch.min()*0.95,a_ch.max()*1.05)

        plt.title('time mean '+ l_str[l_i])

        if (l_i+2) % 2 == 0 and (l_i+1) % 7 !=0:
            plt.ylabel('Temperature (K)')

        if (l_i+1) % 7 ==0:
            plt.ylabel('gain')    

        if (l_i+2) % 8 == 0 or (l_i+1) % 8 == 0:
            plt.xlabel('Frequency (MHz)')

        if l_i==1:
            plt.legend(['deleted data due to other parameters', 'deleted data due to this parameter','clean data'], fontsize=11.5)
    plt.savefig(output_file+'F_'+fname+'_'+ant+'_ch_filter_final.png', bbox_inches='tight')
    #plt.show()

    Tsky_map_h=np.ma.array(Tsky_map_h,mask=mask_clean)
    Tsky_map_v=np.ma.array(Tsky_map_v,mask=mask_clean)

    Tresi_map_h=np.ma.array(Tresi_map_h,mask=mask_clean)
    Tresi_map_v=np.ma.array(Tresi_map_v,mask=mask_clean)


    Tsky_map,Tsky_ratio=cal_map_I(Tsky_map_h, Tsky_map_v)
    Tresi_map=(Tresi_map_h+Tresi_map_v)/2.
    assert((Tsky_map.mask==Tresi_map.mask).all()==True)

    assert((Tsky_ratio.mask==mask_clean).all()==True)
    print np.ma.max(Tsky_ratio)

    assert((Tsky_map.mask==mask_clean).all()==True)
    assert((Tresi_map.mask==mask_clean).all()==True)
    assert((Tsky_ratio.mask==mask_clean).all()==True)

    d={}
    d['Tsky_map']=Tsky_map
    d['Tresi_map']=Tresi_map
    d['Tsky_ratio']=Tsky_ratio

    d['timestamps']=timestamps
    d['ra']=ra
    d['dec']=dec
    assert((nd_s0_h==nd_s0_v).all())
    d['nd_s0']=nd_s0_h
    fs=open(output_file+str(fname)+'_'+str(ant)+'_level4_data','wb')
    pickle.dump(d,fs,protocol=2)
    fs.close()

    d1={}

    d1['Inten_mask']=mask_clean
    fs=open(output_file+str(fname)+'_'+str(ant)+'_level4_mask','wb')
    pickle.dump(d1,fs,protocol=2)
    fs.close()

except IOError:
    print '#### no output.'


print 'end @ ' + time.asctime(time.localtime(time.time())) +'#'
print 'Hakuna Matata' 
################################################
#### Jingying Wang <astro.jywang@gmail.com> ####
################################################