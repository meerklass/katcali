import katdal
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
import pickle
import yaml
from pkg_resources import resource_string
from . import models as km 


def load_data(fname):
    if fname in ['1551037708','1551055211', '1553966342','1554156377']:
        data = katdal.open('/idia/projects/hi_im/MEERKLASS-1/SCI-20180330-MS-01/'+fname+'/'+fname+'/'+fname+'_sdp_l0.full.rdb')
    if fname in ['1555775533','1555793534', '1555861810', '1556034219', '1556052116', '1556120503', '1556138397','1555879611','1561650779']:
        data = katdal.open('/idia/projects/hi_im/MEERKLASS-1/SCI-20190418-MS-01/'+fname+'/'+fname+'/'+fname+'_sdp_l0.full.rdb')
    if fname in ['1558464584','1558472940']:
        data = katdal.open('/idia/projects/hi_im/MEERKLASS-1/COM-20190418-MS-01/'+fname+'/'+fname+'/'+fname+'_sdp_l0.full.rdb')
    if fname=='1562857793':
        data = katdal.open('/idia/projects/hi_im//1562857793/1562857793/1562857793_sdp_l0.full.rdb')
    if fname in ['1630519596','1631552188','1631667564','1631810671','1631990463','1632184922','1633365980','1634402485','1637346562','1637699408',
                '1638301944','1638647186','1639331184','1640712986','1631379874','1631559762','1631724508','1631818149','1632069690','1632505883',
                '1633970780','1634748682','1637354605','1638130295','1638386189','1638898468','1639935088','1640799689','1631387336','1631659886',
                '1631732038','1631982988','1632077222','1632760885','1634252028','1634835083','1637691677','1638294319','1638639082','1639157507',
                '1640540184']:
        data = katdal.open('/idia/projects/hi_im/MEERKLASS-1/SCI-20210212-MS-01/'+fname+'/'+fname+'/'+fname+'_sdp_l0.full.rdb')
    ###UHF####        
    if fname in ['1665679673', '1666205857',  '1675021905', '1675816512', '1677174749', '1678295187', '1678899080', '1679615321', 
                 '1681920680',  '1688399183', '1665938490', '1666284501', '1675106912', '1676313206', '1677183387', '1678381591', 
                 '1679247986', '1680626188', '1682448988', '1689003684', '1666024766', '1666293386', '1675210948', '1676657789',
                 '1677195529', '1678467685', '1679333668', '1680644082', '1683492604', '1689090392', '1666032415', '1666370606', 
                 '1675623808', '1677002481', '1677777992', '1678726283', '1679419886', '1680798562', '1684087370',  '1689176790', 
                 '1666111882', '1666543907', '1675632179', '1677011008', '1677795989', '1678734987', '1679592842', '1681143685', '1684781618',
        '1666198047', '1666630286', '1675643846', '1677020482', '1678122565', '1678743988', '1679605292', '1681229848', '1685641589']:
        data= katdal.open('/idia/projects/hi_im/MEERKLASS-1/SCI-20220822-MS-01/'+fname+'/'+fname+'/'+fname+'_sdp_l0.full.rdb')
    if fname in ['1709659886', '1715108064', '1715012489', '1715098888', '1715032586', '1715208401','1709832691', '1710869782', '1712685146', '1710282085', '1710888742', '1712791459']:
        data = katdal.open('/idia/projects/hi_im/MEERKLASS-1/SCI-20230907-MS-01/'+fname+'/'+fname+'/'+fname+'_sdp_l0.full.rdb')
    return data


def check_ants(fname):
    
    ant_param = resource_string(__name__, 'config.yaml')
    config = yaml.safe_load(ant_param)

        
    bad_ants=[]
    target_list=[]
    target = None

    entry = config.get('entries', {}).get(str(fname), {})
    if entry:
        if 'target' in entry:
            target = entry['target']
        if 'target_list' in entry:
            target_list = entry['target_list']
        bad_ants = entry.get('bad_ants', [])
        
    ''' 
    #below have been moved to config.yaml, these codes will be deleted soon if YANL works well   
    #####UHF data#######

    if fname in ['1684087370','1675210948','1675643846','1679615321', '1680644082','1710888742']:
        target='3C273'
    if fname in ['1675632179', '1679247986', '1679592842']:
        target='PictorA'
    if fname in ['1715032586', '1715208401','1712791459']:
        target = '3C353'
        
    if fname in ['1709659886','1709832691', '1710869782']:
        target_list=['HydraA','PictorA']
    if fname in ['1712685146']:
        target_list=['PictorA', 'HydraA']
    if fname in ['1715012489', '1715098888']:
        target_list = ['PictorA', '3C273']
    if fname in ['1715108064','1710282085']:
        target_list=['', '3C273']
        
    if fname in ['1684087370']:
        bad_ants=['m050']
    if fname in ['1675632179']:
        bad_ants=['m017','m020','m036']
    if fname in ['1679247986']:
        bad_ants = ['m009', 'm016', 'm017']
    if fname in ['1679592842']:
        bad_ants = ['m017', 'm018']
    if fname in ['1679615321']:
        bad_ants = ['m017', 'm018', 'm046']
    if fname in ['1680644082']:
        bad_ants = ['m017', 'm030', 'm033', 'm034', 'm035']
    if fname in ['1715012489', '1715098888', '1715032586', '1715208401']:
        bad_ants = ['m044', 'm045', 'm046', 'm058', 'm061']
    if fname in ['1710282085']:
        bad_ants=['m004', 'm013', 'm024', 'm025', 'm026']
    if fname in ['1710869782', '1710888742']:
        bad_ants=['m027', 'm028', 'm029', 'm035', 'm038']
    if fname in ['1712685146']:
        bad_ants=['m036', 'm037', 'm038', 'm042']
    if fname in ['1712791459']:
        bad_ants=['m042']
    ###########2021  data ############

    if fname in ['1634835083', '1634748682', '1634402485', '1633970780', '1633365980', '1632760885', '1632505883', '1632077222', '1632069690', '1631990463', '1631982988', '1631818149', '1631810671', '1631732038', '1631724508', '1631559762', '1631552188', '1631387336', '1631379874', '1630519596']:

        target='PKS1934-638'

    if fname in ['1640799689', '1640712986', '1640540184', '1639935088', '1639331184', '1639157507', '1638898468', '1638647186', '1638639082', '1638386189', '1638301944', '1638294319', '1638130295', '1637699408', '1637691677', '1637354605', '1637346562', '1634252028', '1632184922', '1631667564', '1631659886']:

        target='PictorA'

    

    #######from Marta######
    if fname=='1634835083': 
        bad_ants=['m012','m044','m047']
    if fname=='1634748682': 
        bad_ants=['m034','m041']
    if fname=='1634402485': 
        bad_ants=['m003','m006','m035']
    if fname=='1633970780': 
        bad_ants=['m027']
    if fname=='1632505883': 
        bad_ants=['m063']
    if fname=='1632077222': 
        bad_ants=['m033']
    if fname=='1632069690': 
        bad_ants=['m033']
    if fname=='1631990463': 
        bad_ants=['m031']
    if fname=='1631818149': 
        bad_ants=['m020','m052']
    if fname=='1631732038': 
        bad_ants=['m013','m017']
    if fname=='1631724508': 
        bad_ants=['m013','m017']
    if fname=='1631552188': 
        bad_ants=['m013','m014','m032']
    if fname=='1639157507': 
        bad_ants=['m022','m034','m054']
    if fname=='1638898468': 
        bad_ants=['m036']
    if fname=='1638639082': 
        bad_ants=['m009','m011']
    if fname=='1638386189': 
        bad_ants=['m063']
    if fname=='1638294319': 
        bad_ants=['m025','m031']
    if fname=='1637699408': 
        bad_ants=['m052','m060']
    if fname=='1637691677': 
        bad_ants=['m022','m026','m052','m063']
    if fname=='1637346562': 
        bad_ants=['m037','m049']
    if fname=='1634835083': 
        bad_ants=['m012','m044','m047']
    if fname=='1634748682': 
        bad_ants=['m034','m041']
    if fname=='1634402485': 
        bad_ants=['m003','m006','m035']
    if fname=='1633970780': 
        bad_ants=['m027']
    if fname=='1632505883': 
        bad_ants=['m063']
    if fname=='1632077222': 
        bad_ants=['m033']
    if fname=='1632069690': 
        bad_ants=['m033']
    if fname=='1631990463': 
        bad_ants=['m031']
    if fname=='1631818149': 
        bad_ants=['m020','m052']
    if fname=='1631732038': 
        bad_ants=['m013','m017']
    if fname=='1631724508': 
        bad_ants=['m013','m017']
    if fname=='1631552188': 
        bad_ants=['m013','m014','m032']
    #######from Marta#######


################################################################################################################################    
    if fname=='1551037708':
        bad_ants=['m001', 'm007', 'm008', 'm018', 'm023', 'm025', 'm032', 'm036', 'm038', 'm040', 'm041', 'm059']
        target='3C237'
        
    if fname=='1551055211':
        bad_ants=['m018', 'm025', 'm032', 'm036', 'm041']
        target='3C273'
        
    if fname=='1553966342':
        bad_ants=['m004', 'm017', 'm024', 'm032', 'm036', 'm041']
        target='PictorA'
    
    if fname=='1554156377':
        bad_ants=['m014', 'm017', 'm032', 'm036', 'm041', 'm054']
        target='3C273'
    
    if fname=='1556034219':
        #bad_ants=['m000', 'm001', 'm002', 'm003', 'm004', 'm005', 'm006', 'm007', 'm008', 'm009', 'm010', 'm011', 'm012', 'm013', 'm014', 'm015', 'm016', 'm017', 'm018', 'm019', 'm020', 'm021', 'm022', 'm023', 'm024', 'm025', 'm026', 'm027', 'm028', 'm029', 'm030', 'm031', 'm032', 'm033', 'm034', 'm035', 'm036', 'm037', 'm038', 'm039', 'm040', 'm041', 'm042', 'm043', 'm044', 'm045', 'm046', 'm047', 'm048', 'm049', 'm050', 'm051', 'm052', 'm053', 'm054', 'm055', 'm056', 'm057', 'm058', 'm059', 'm060', 'm061', 'm062', 'm063']
        bad_ants=[] #for test
        target='PictorA'
    
    if fname=='1556052116':
        #bad_ants=['m000', 'm001', 'm002', 'm003', 'm004', 'm005', 'm006', 'm007', 'm008', 'm009', 'm010', 'm011', 'm012', 'm013', 'm014', 'm015', 'm016', 'm017', 'm018', 'm019', 'm020', 'm021', 'm022', 'm023', 'm024', 'm025', 'm026', 'm027', 'm028', 'm029', 'm030', 'm031', 'm032', 'm033', 'm034', 'm035', 'm036', 'm037', 'm038', 'm039', 'm040', 'm041', 'm042', 'm043', 'm044', 'm045', 'm046', 'm047', 'm048', 'm049', 'm050', 'm051', 'm052', 'm053', 'm054', 'm055', 'm056', 'm057', 'm058', 'm059', 'm060', 'm061', 'm062', 'm063']
        bad_ants=[]
        target='3C273'
    
    if fname=='1556120503':
        bad_ants=['m014', 'm017', 'm024', 'm032', 'm035', 'm041', 'm042', 'm057', 'm059']
        target='PictorA'
     
    if fname=='1556138397':
        bad_ants=['m014', 'm015', 'm016', 'm017', 'm024', 'm032', 'm033', 'm041', 'm042', 'm057', 'm059']
        target='3C273'
        
    if fname=='1555775533':
        #bad_ants=['m000', 'm001', 'm002', 'm003', 'm004', 'm005', 'm006', 'm007', 'm008', 'm009', 'm010', 'm011', 'm012', 'm013', 'm014', 'm015', 'm016', 'm017', 'm018', 'm019', 'm020', 'm021', 'm022', 'm023', 'm024', 'm025', 'm026', 'm027', 'm028', 'm029', 'm030', 'm031', 'm032', 'm033', 'm034', 'm035', 'm036', 'm037', 'm038', 'm039', 'm040', 'm041', 'm042', 'm043', 'm044', 'm045', 'm046', 'm047', 'm048', 'm049', 'm050', 'm051', 'm052', 'm053', 'm054', 'm055', 'm056', 'm057', 'm058', 'm059', 'm060', 'm061', 'm062', 'm063']
        bad_ants=[] #for test
        target='PictorA'

    if fname=='1555793534':
        #bad_ants=['m000', 'm001', 'm002', 'm003', 'm004', 'm005', 'm006', 'm007', 'm008', 'm009', 'm010', 'm011', 'm012', 'm013', 'm014', 'm015', 'm016', 'm017', 'm018', 'm019', 'm020', 'm021', 'm022', 'm023', 'm024', 'm025', 'm026', 'm027', 'm028', 'm029', 'm030', 'm031', 'm032', 'm033', 'm034', 'm035', 'm036', 'm037', 'm038', 'm039', 'm040', 'm041', 'm042', 'm043', 'm044', 'm045', 'm046', 'm047', 'm048', 'm049', 'm050', 'm051', 'm052', 'm053', 'm054', 'm055', 'm056', 'm057', 'm058', 'm059', 'm060', 'm061', 'm062', 'm063']
        bad_ants=[] #for test
        target='3C273'
        
    if fname=='1555861810':
        bad_ants=['m000', 'm001', 'm002', 'm003', 'm004', 'm005', 'm006', 'm007', 'm008', 'm009', 'm010', 'm011', 'm012', 'm013', 'm014', 'm015', 'm016', 'm017', 'm018', 'm019', 'm020', 'm021', 'm022', 'm023', 'm024', 'm025', 'm026', 'm027', 'm028', 'm029', 'm030', 'm031', 'm032', 'm033', 'm034', 'm035', 'm036', 'm037', 'm038', 'm039', 'm040', 'm041', 'm042', 'm043', 'm044', 'm045', 'm046', 'm047', 'm048', 'm049', 'm050', 'm051', 'm052', 'm053', 'm054', 'm055', 'm056', 'm057', 'm058', 'm059', 'm060', 'm061', 'm062', 'm063']
        target='PictorA'
    
    if fname=='1555879611':
        bad_ants=['m000', 'm001', 'm002', 'm003', 'm004', 'm005', 'm006', 'm007', 'm008', 'm009', 'm010', 'm011', 'm012', 'm013', 'm014', 'm015', 'm016', 'm017', 'm018', 'm019', 'm020', 'm021', 'm022', 'm023', 'm024', 'm025', 'm026', 'm027', 'm028', 'm029', 'm030', 'm031', 'm032', 'm033', 'm034', 'm035', 'm036', 'm037', 'm038', 'm039', 'm040', 'm041', 'm042', 'm043', 'm044', 'm045', 'm046', 'm047', 'm048', 'm049', 'm050', 'm051', 'm052', 'm053', 'm054', 'm055', 'm056', 'm057', 'm058', 'm059', 'm060', 'm061', 'm062', 'm063']
        target='3C273'
    
    if fname=='1561650779':
        bad_ants=['m001', 'm006', 'm039', 'm054', 'm056']
        target='3C273'

    if fname=='1558464584':
        bad_ants=['m006', 'm017', 'm024', 'm032', 'm036', 'm051', 'm052', 'm060', 'm061']
        target='3C273'

    if fname=='1558472940':
        bad_ants=['m006', 'm017', 'm024', 'm032', 'm036', 'm051', 'm052', 'm059', 'm061', 'm062']
        target='3C273'
    
    if fname=='1562857793':
        bad_ants=['m000', 'm001', 'm002', 'm003', 'm035', 'm039', 'm043', 'm045', 'm048', 'm056', 'm058', 'm059', 'm060']
        target='3C273'
                  
    if fname=='1580260015':#test only
        bad_ants=[]
        target='3C273'
        
    if fname=='1579725085':#test only
        bad_ants=[]
        target='PictorA'
    #above have been moved to config.yaml, these codes will be deleted soon if YANL works well
    ''' 
    
    print ('bad_ants: '+ str(bad_ants))

    if len(target_list) == 0 and target is not None:    
        print ('calibrator number = 1')
        flux_model,cc=target_para(target)
        print ('calibrator: '+str(target)+', ra,dec= '+str(cc.ra)+', '+str(cc.dec))
    elif len(target_list) == 2:
        print ('calibrator number = 2')
        flux_model,cc=[],[]
        for target_local in target_list:
            flux_model_i,cc_i=target_para(target_local)
            flux_model.append(flux_model_i)
            cc.append(cc_i)
            if cc_i != None:
                print ('calibrator: '+str(target_local)+', ra,dec= '+str(cc_i.ra)+', '+str(cc_i.dec))
            else:
                print ('calibrator: EMPTY')
        target=target_list
    else:
        flux_model, cc = None, None
        print ('#failed to load flux_model')
    return target,cc,bad_ants,flux_model


def target_para(target):
    flux_model,cc= None, None
    ################set calibrator coordinates####################################
    if target=='HydraA':
        flux_model=km.flux_HydraA
        cc = SkyCoord(139.523617*u.deg, -12.095502*u.deg, frame='icrs') #HydraA
    if target=='3C353':
        flux_model=km.flux_3C353
        cc = SkyCoord(260.117381*u.deg, -0.979642*u.deg, frame='icrs') #3C353
    
    if target=='3C273':
        flux_model=km.flux_3C273
        cc = SkyCoord(187.2779154*u.deg,  2.0523883*u.deg, frame='icrs') #3C273
    if target=='3C237':
        flux_model=km.flux_3C237
        cc = SkyCoord(152.000125*u.deg,  7.504541*u.deg, frame='icrs') #3C237
    if target=='PictorA':
        flux_model=km.flux_PictorA
        cc = SkyCoord(79.9571708*u.deg,  -45.7788278*u.deg, frame='icrs') #Pictor A
    if target=='PKS1934-638':
        flux_model=km.flux_1934_sd
        cc = SkyCoord(294.854275*u.deg, -63.712674*u.deg, frame='icrs') #PKS 1934-63 (1934-638)

    return flux_model,cc


def ant_list(data):
    ant_list=[]
    for corr_i in range(len(data.ants)):
        assert(data.corr_products[corr_i][0]==data.corr_products[corr_i][1]) #check auto-corr
        ant_list.append(data.corr_products[corr_i][0][0:4])
    return np.array(ant_list)

def call_vis(fname,recv):
    ####UHF band####
    #if fname in ['1684087370','1675632179','1675210948','1675643846']:
    if fname in ['1665679673', '1666205857',  '1675021905', '1675816512', '1677174749', '1678295187', '1678899080', '1679615321', '1681920680', 
                 '1688399183', '1665938490', '1666284501', '1675106912',   '1676313206', '1677183387', '1678381591', '1679247986', '1680626188', 
                 '1682448988', '1689003684', '1666024766', '1666293386', '1675210948', '1676657789', '1677195529', '1678467685', '1679333668', 
                 '1680644082', '1683492604', '1689090392', '1666032415', '1666370606', '1675623808', '1677002481', '1677777992', '1678726283', 
                 '1679419886', '1680798562', '1684087370',  '1689176790', '1666111882', '1666543907', '1675632179', '1677011008', '1677795989', 
                 '1678734987', '1679592842', '1681143685', '1684781618', '1666198047', '1666630286', '1675643846', '1677020482', '1678122565',
                 '1678743988', '1679605292', '1681229848', '1685641589', '1709659886','1715108064']:
        data1 = pickle.load(open('/idia/projects/hi_im/raw_vis/MeerKLASS2023/'+str(fname)+'/'+str(fname)+'_'+str(recv)+'_vis_data','rb'))

    if fname in ['1715012489', '1715098888', '1715032586', '1715208401','1709832691', '1710869782', '1712685146', '1710282085', '1710888742', '1712791459']:
        data1 = pickle.load(open('/idia/projects/hi_im/raw_vis/MeerKLASS2024/'+str(fname)+'/'+str(fname)+'_'+str(recv)+'_vis_data','rb'))
        
    #####L-band####
    if fname in ['1551037708','1551055211', '1553966342','1554156377']:
        data1 = pickle.load(open('/idia/projects/hi_im/raw_vis/SCI-20180330-MS-01/'+str(fname)+'/'+str(fname)+'_'+str(recv)+'_vis_data','rb'))
    if fname in ['1555775533','1555793534', '1555861810', '1556034219', '1556052116', '1556120503', '1556138397','1555879611','1561650779', '1562857793', '1579725085', '1580260015']:
        data1 = pickle.load(open('/idia/projects/hi_im/raw_vis/SCI-20190418-MS-01/'+str(fname)+'/'+str(fname)+'_'+str(recv)+'_vis_data','rb'))
    if fname in['1558464584','1558472940']:
        data1 = pickle.load(open('/idia/projects/hi_im/raw_vis/COM-20190418-MS-01/'+str(fname)+'/'+str(fname)+'_'+str(recv)+'_vis_data','rb'))
    if fname in ['1630519596','1631552188','1631667564','1631810671','1631990463','1632184922','1633365980','1634402485','1637346562','1637699408',
                '1638301944','1638647186','1639331184','1640712986','1631379874','1631559762','1631724508','1631818149','1632069690','1632505883',
                '1633970780','1634748682','1637354605','1638130295','1638386189','1638898468','1639935088','1640799689','1631387336','1631659886',
                '1631732038','1631982988','1632077222','1632760885','1634252028','1634835083','1637691677','1638294319','1638639082','1639157507',
                '1640540184']:
        try:
            data1 = pickle.load(open('/idia/projects/hi_im/raw_vis/SCI-20210212-MS-01/'+str(fname)+'/'+str(fname)+'_'+str(recv)+'_vis_data','rb'))
        except(Exception):
            data1 = pickle.load(open('/idia/projects/hi_im/raw_vis/SCI-20210212-MS-01/'+str(fname)+'/'+str(fname)+'_'+str(recv)+'_vis_data','rb'), encoding='latin-1')
            print ('# loaded data was saved in python2')
    print (data1['recv_pair'].astype(np.str_))
    recv1=data1['recv_pair'].astype(np.str_)[0]
    assert(recv1==recv)

    vis=data1['vis']
    flags=data1['flags']
    return vis,flags

def cal_corr_id(data,recv):
    for i in range(len(data.ants)*2):
        a=data.corr_products[i]
        if (a[0]==recv and a[1]==recv):
            corr_id=i
            break
    return corr_id

def load_coordinates(data):
    ra=data.ra[:,0]
    dec=data.dec[:,0]    
    az=data.az[:,0]
    el=data.el[:,0]

    return ra,dec,az,el

def load_ndparam(fname, data):
    nd_set=float(fname)
    nd_time=1.799235 
    nd_cycle=19.9915424299 #only for stable diode noise pattern
    nd_ratio=1.8/data.dump_period #1.8???test
    return nd_set,nd_time,nd_cycle,nd_ratio

def load_tf(data):
    freqs = data.freqs
    timestamps=data.timestamps
    return timestamps,freqs

def load_ang_deg(ra,dec,cc):
    c_obs = SkyCoord(ra*u.deg,  dec*u.deg, frame='icrs')
    ang_deg=cc.separation(c_obs)/u.deg
    return ang_deg

def load_ang_deg2(ra,dec,cc):
    c_obs = SkyCoord(ra*u.deg,  dec*u.deg, frame='icrs')
    print ('# calibrator number: '+str(len(cc)))
    ang_deg=[]
    for i in range(len(cc)):
        if cc[i] != None:
            ang_deg_i=cc[i].separation(c_obs)/u.deg
            ang_deg.append(ang_deg_i)
        else:
            print ('* no track here')
            ang_deg.append([np.NaN]*len(ra))
    return np.array(ang_deg)

def cal_freq(ch):
    v_min=856.0 
    v_max=1712.0 
    dv=0.208984375
    assert((v_max-v_min)/dv==4096)
    freq_MHz=ch*dv+v_min
    freq=freq_MHz*1e6
    return freq

def cal_freqs(ch_list):
    freq_list=[]
    v_min=856.0 
    v_max=1712.0 
    dv=0.208984375
    assert((v_max-v_min)/dv==4096)
    for ch in ch_list:
        freq_MHz=ch*dv+v_min
        freq=freq_MHz*1e6
        freq_list.append(freq)
    freq_list=np.array(freq_list)
    return freq_list

def cal_freq_UHF(ch):
    f_min=5.44000000e+08
    freq_gap=132812.5
    return f_min+ch*freq_gap

def cal_freqs_UHF(ch_list):
    freq_list=[]
    f_min=5.44000000e+08
    freq_gap=132812.5
    for ch in ch_list:
        freq_list.append(f_min+ch*freq_gap)
    return freq_list
    
################ for data has no rdb link #####################
#tools to save the parameters
def is_dict(x):
    return issubclass(dict, type(x))

class RecDictWrapper(object):
    def __init__(self, data):
        self.keys=[]
        for k,v in data.items():
            self.keys.append(k)
            if is_dict(v):
                self.__setattr__(k,RecDictWrapper(v))
            else:
                self.__setattr__(k, v)
###############################################################

def mask_None(list_input):
    return np.ma.masked_invalid(np.array(list_input,dtype=float), copy=False).astype(float)
