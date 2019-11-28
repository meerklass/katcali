

def legacy_metadata(fname):
    """
    Aggregate all file-specific metadata into a single function, so it can be 
    converted into a more uniform external data format.
    """
    
    # File paths (from io.load_data)
    data = None
    if fname in ['1551037708','1551055211', '1553966342','1554156377']:
        data = '/idia/projects/hi_im/SCI-20180330-MS-01/'
    if fname in ['1555775533','1555793534', '1555861810', '1556034219', 
                 '1556052116', '1556120503', '1556138397','1555879611', 
                 '1561650779']:
        data = '/idia/projects/hi_im/SCI-20190418-MS-01/'
    if fname in['1558464584','1558472940']:
        data = '/idia/projects/hi_im/COM-20190418-MS-01/'
    if fname=='1562857793':
        data = '/idia/projects/hi_im//'
    
    # Antennas (from io.check_ants)
    if fname=='1551037708':
        bad_ants=['m001', 'm007', 'm008', 'm018', 'm023', 'm025', 'm032', 
                  'm036', 'm038', 'm040', 'm041', 'm059']
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
        bad_ants=['m000', 'm001', 'm002', 'm003', 'm004', 'm005', 'm006', 'm007', 'm008', 'm009', 'm010', 'm011', 'm012', 'm013', 'm014', 'm015', 'm016', 'm017', 'm018', 'm019', 'm020', 'm021', 'm022', 'm023', 'm024', 'm025', 'm026', 'm027', 'm028', 'm029', 'm030', 'm031', 'm032', 'm033', 'm034', 'm035', 'm036', 'm037', 'm038', 'm039', 'm040', 'm041', 'm042', 'm043', 'm044', 'm045', 'm046', 'm047', 'm048', 'm049', 'm050', 'm051', 'm052', 'm053', 'm054', 'm055', 'm056', 'm057', 'm058', 'm059', 'm060', 'm061', 'm062', 'm063']
        target='PictorA'
    
    if fname=='1556052116':
        bad_ants=['m000', 'm001', 'm002', 'm003', 'm004', 'm005', 'm006', 'm007', 'm008', 'm009', 'm010', 'm011', 'm012', 'm013', 'm014', 'm015', 'm016', 'm017', 'm018', 'm019', 'm020', 'm021', 'm022', 'm023', 'm024', 'm025', 'm026', 'm027', 'm028', 'm029', 'm030', 'm031', 'm032', 'm033', 'm034', 'm035', 'm036', 'm037', 'm038', 'm039', 'm040', 'm041', 'm042', 'm043', 'm044', 'm045', 'm046', 'm047', 'm048', 'm049', 'm050', 'm051', 'm052', 'm053', 'm054', 'm055', 'm056', 'm057', 'm058', 'm059', 'm060', 'm061', 'm062', 'm063']
        target='3C273'
    
    if fname=='1556120503':
        bad_ants=['m014', 'm017', 'm024', 'm032', 'm035', 'm041', 'm042', 'm057', 'm059']
        target='PictorA'
     
    if fname=='1556138397':
        bad_ants=['m014', 'm015', 'm016', 'm017', 'm024', 'm032', 'm033', 'm041', 'm042', 'm057', 'm059']
        target='3C273'
        
    if fname=='1555775533':
        bad_ants=['m000', 'm001', 'm002', 'm003', 'm004', 'm005', 'm006', 'm007', 'm008', 'm009', 'm010', 'm011', 'm012', 'm013', 'm014', 'm015', 'm016', 'm017', 'm018', 'm019', 'm020', 'm021', 'm022', 'm023', 'm024', 'm025', 'm026', 'm027', 'm028', 'm029', 'm030', 'm031', 'm032', 'm033', 'm034', 'm035', 'm036', 'm037', 'm038', 'm039', 'm040', 'm041', 'm042', 'm043', 'm044', 'm045', 'm046', 'm047', 'm048', 'm049', 'm050', 'm051', 'm052', 'm053', 'm054', 'm055', 'm056', 'm057', 'm058', 'm059', 'm060', 'm061', 'm062', 'm063']
        target='PictorA'

    if fname=='1555793534':
        bad_ants=['m000', 'm001', 'm002', 'm003', 'm004', 'm005', 'm006', 'm007', 'm008', 'm009', 'm010', 'm011', 'm012', 'm013', 'm014', 'm015', 'm016', 'm017', 'm018', 'm019', 'm020', 'm021', 'm022', 'm023', 'm024', 'm025', 'm026', 'm027', 'm028', 'm029', 'm030', 'm031', 'm032', 'm033', 'm034', 'm035', 'm036', 'm037', 'm038', 'm039', 'm040', 'm041', 'm042', 'm043', 'm044', 'm045', 'm046', 'm047', 'm048', 'm049', 'm050', 'm051', 'm052', 'm053', 'm054', 'm055', 'm056', 'm057', 'm058', 'm059', 'm060', 'm061', 'm062', 'm063']
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
    
    # Jump limit (from label_nd_injection) 
    f=10.
    if fname in ['1555793534','1551055211','1551037708','1555775533']:
        f=10
    if fname=='1556120503':
        f=2
    if fname=='1556052116':
        f=15.
    
    # Noise diode initial step (diode.call_nd_1a_param)
    nd_1a_gap=10
    nd_1a0=-999
    if fname=='1551055211':
        nd_1a0=1
    if fname=='1555793534':
        nd_1a0 =8
    if fname=='1551037708':
        nd_1a0 =4
    if fname=='1553966342':
        nd_1a0 =4
    if fname=='1554156377':
        nd_1a0 =5
    if fname=='1555775533':
        nd_1a0 =0
    if fname=='1556034219':
        nd_1a0 =-1 #to keep nd_1b in the list
    if fname=='1556120503':
        nd_1a0 =0
    if fname=='1556138397':
        nd_1a0 =2
    if fname=='1556052116':
        nd_1a0 =1
    if fname=='1561650779':
        nd_1a0 =7
    if fname=='1562857793':
        nd_1a0 =7
    
    # Diode lead time (from diode.cal_t_line)
    trange = None
    if fname=='1551037708':
        trange = (5, None)
    if fname=='1551055211':
        trange = (6, None)
    if fname=='1553966342':
        trange = (2, None)
    if fname=='1556034219':
        trange = (2, None)
    if fname=='1554156377':
        trange = (1, None)
    if fname=='1556138397':
        trange = (2,-1)
    if fname=='1556052116':
        trange = (2, None)
    if fname=='1555775533':
        trange = (2, None)
    if fname=='1556120503':
        trange = (2, None)
    if fname=='1555793534':
        trange = (2, None)
    if fname=='1561650779':
        trange = (3,-1)
    if fname=='1562857793':
        trange = (3, None)
    
    return target, bad_ants, f, data, nd_1a0, trange


if __name__ == '__main__':
    # Process filenames
    
    # Prepare DataSet object
    from dataset import DataSet
    ds = DataSet("meerkat2019.json")
    
    # Filenames
    fnames = [
        '1551037708',
        '1551055211',
        '1553966342',
        '1554156377',
        '1556034219',
        '1556052116',
        '1556120503',
        '1556138397',
        '1555775533',
        '1555793534',
        '1555861810',
        '1555879611',
        '1561650779',
        '1558464584',
        '1558472940',
        '1562857793'
    ]
    
    # Loop over filenames
    for f in fnames:
        target, bad_ants, fjump, data, nd_1a0, trange = legacy_metadata(f)
        
        root = data
        
        diode = {
            'offset':       nd_1a0,
            'jump_limit':   fjump,
            'time_range':   trange,
        }
        
        print("%s" % f)
        print("  %s" % target)
        print("  %s" % fjump)
        print("  %s" % data)
        print("  %s" % nd_1a0)
        print("  %s" % str(trange))
        print("  %s" % len(bad_ants))
        print("")
    
        # Load this file into DataSet
        ds.add_file(f, target=target, root=root, diode=diode, bad_ants=bad_ants)
    
    # Save dataset
    ds.save()
    
    # Try loading it as a new dataset
    ds2 = DataSet("meerkat2019.json")
    print(list(ds2.files.keys()))
    
    
    
