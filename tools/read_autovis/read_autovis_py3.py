#!/usr/bin/env python3

import katdal
import numpy as np
import pickle
import sys

fname=sys.argv[1]

print ('read data in...')

#if fname in ['1689090392']:
#    data = katdal.open('/idia/projects/hi_im/temp_raw/SCI-20220822-MS-01/'+fname+'/'+fname+'/'+fname+'_sdp_l0.full.rdb')
#else:
if fname in ['1684087370','1678726283','1678734987','1679615321','1666370606']:
    data = katdal.open('/idia/projects/hi_im/SCI-20220822-MS-01/'+fname+'/'+fname+'/'+fname+'_sdp_l0.full.rdb')

'''
if fname=='1666111882':
    data = katdal.open('https://archive-gw-1.kat.ac.za/1666111882/1666111882_sdp_l0.full.rdb?token=eyJ0eXAiOiJKV1QiLCJhbGciOiJFUzI1NiJ9.eyJpc3MiOiJrYXQtYXJjaGl2ZS5rYXQuYWMuemEiLCJhdWQiOiJhcmNoaXZlLWd3LTEua2F0LmFjLnphIiwiaWF0IjoxNjY2NTk3MTMxLCJwcmVmaXgiOlsiMTY2NjExMTg4MiJdLCJleHAiOjE2NjcyMDE5MzEsInN1YiI6ImFzdHJvLmp5d2FuZ0BnbWFpbC5jb20iLCJzY29wZXMiOlsicmVhZCJdfQ.Ph4JtFzFLWHMd4W6n8A1evFCGiEH4JdO3hGi5kXHfupW4uWCmcinExCltTEFqY0yFr3uP_ct8DLNxskmRMCbvA')
'''
#data = katdal.open('/idia/projects/hi_im/SCI-20210212-MS-01/'+fname+'/'+fname+'/'+fname+'_sdp_l0.full.rdb')
    
print (data)

nt,nch,ncorr=data.vis.shape
print ('nt, nch, ncorr = '+ str(nt)+', '+str(nch)+', '+str(ncorr))

autocorr_ids=np.where(data.corr_products[:,0]==data.corr_products[:,1])[0]
print ('autocorr_ids = '+ str(autocorr_ids))

recvs=data.corr_products[autocorr_ids,0]
print ('recvs: '+ str(recvs))

auto_vis=np.zeros([nt,nch,len(autocorr_ids)],dtype=float)
auto_flags=np.zeros([nt,nch,len(autocorr_ids)],dtype=bool)

iter_vis=iter(data.vis) #for nt
iter_flags=iter(data.flags) #for nt

print ('split vis data...')
print (str(nt)+' in total')

for nt_i in range(nt):
    print (nt_i)
    auto_vis[nt_i,:,:]=next(iter_vis)[:,autocorr_ids].real
    auto_flags[nt_i,:,:]=next(iter_flags)[:,autocorr_ids]
    
print ('saving data...')

for i in autocorr_ids:
    print (i)
    d={}
    d['recv_pair']=data.corr_products[i]
    d['vis']=auto_vis[:,:,i]
    d['flags']=auto_flags[:,:,i]
    
    fs=open(fname+'_'+str(recvs[i])+'_vis_data','wb')
    pickle.dump(d,fs,protocol=2)
    fs.close()
    
print ('Hakuna Matata')

