import numpy as np
from scipy import stats
from . import filter as kf
from .utils import valid_filename


def select_track(data, ant, pol):
    """

    Parameters
    ----------
    data : metadata object
        x

    ant : str
        Antenna, e.g. 'm006'.

    pol : str
        Polarisation, e.g. 'h'.

    Returns
    -------
    idxs_track : ?
        ?
    """
    data.select(ants=ant, pol=pol, scans='track')
    scans_t = []
    scans_tw = []
    for s in data.scans():
        if data.shape[0] > 50:
            scans_t.append(data.scan_indices[0])
        else:
            scans_tw.append(data.scan_indices[0])
    data.select(ants=ant, pol=pol, scans=scans_t)
    dp_t = data.dumps
    data.select()  # recover after select!!!
    data.select(ants=ant, pol=pol)
    return dp_t, scans_tw


def select_scan(data, ant, pol):
    """

    Parameters
    ----------
    data : metadata object
        x

    ant : str
        Antenna, e.g. 'm006'.

    pol : str
        Polarisation, e.g. 'h'.

    Returns
    -------
    idxs_waste : ?
        ?
    """
    data.select(ants=ant, pol=pol, scans='scan')
    scans_s = []
    scans_sw = []
    for s in data.scans():
        if data.shape[0] > 50:
            scans_s.append(data.scan_indices[0])
        else:
            scans_sw.append(data.scan_indices[0])
    data.select(ants=ant, pol=pol, scans=scans_s)
    dp_s = data.dumps
    # dp_sb=data.dumps[0]
    # dp_se=data.dumps[-1]
    data.select()  # recover after select!!!
    data.select(ants=ant, pol=pol)
    return dp_s, scans_sw


def select_waste(data, ant, pol):
    """

    Parameters
    ----------
    data : metadata object
        x

    ant : str
        Antenna, e.g. 'm006'.

    pol : str
        Polarisation, e.g. 'h'.

    Returns
    -------
    idxs_waste : ?
        ?
    """
    data.select(ants=ant, pol=pol, scans=('slew', 'stop'))
    dp_w1 = data.dumps

    dp_t, scans_tw = select_track(data, ant, pol)
    dp_s, scans_sw = select_scan(data, ant, pol)

    data.select(ants=ant, pol=pol, scans=scans_tw+scans_sw)
    dp_w2 = data.dumps
    dp_w = list(dp_w1)+list(dp_w2)
    dp_w.sort()
    dp_w = np.array(dp_w)
    data.select()  # recover after select!!!
    data.select(ants=ant, pol=pol)
    return dp_w


def label_dump_1ch(data, ant, pol, flags, ch):
    """

    Parameters
    ----------

    Returns
    -------

    """
    dp_t, scans_tw = select_track(data, ant, pol)
    dp_s, scans_sw = select_scan(data, ant, pol)
    dp_w = select_waste(data, ant, pol)

    flags_1ch = flags[:, ch]

    dp_tt = []
    dp_ss = []
    dp_f = []

    for i in dp_t:
        if flags_1ch[i] == False:
            dp_tt.append(i)
        else:
            dp_f.append(i)

    for i in dp_s:
        if flags_1ch[i] == False:
            dp_ss.append(i)
        else:
            dp_f.append(i)
    return dp_tt, dp_ss, dp_f, dp_t, dp_s


def cal_dp_c(fname, data, ant, pol, flags, ch, dp_tt, dp_ss, ang_deg):
    """

    Parameters
    ----------

    Returns
    -------

    """
    sigma_level = 10
    n_iter = 3
    flags_1ch = flags[:, ch]
    dp_sb = dp_ss[0]
    dp_se = dp_ss[-1]
    data.select(ants=ant, pol=pol, scans='track', targets=0)
    dp_c0 = data.dumps
    for i in dp_c0:
        if i not in dp_tt or flags_1ch[i] == True:
            dp_c0 = list(dp_c0)
            dp_c0.remove(i)
            #print 'rm '+str(i)+' from dp_c0'
    dp_c0 = kf.deg_filter(dp_c0, ang_deg, sigma_level, n_iter)
    dp_c0 = np.array(dp_c0)
    dp_c0a = dp_c0[dp_c0 < dp_sb]
    dp_c0b = dp_c0[dp_c0 > dp_se]

    data.select(ants=ant, pol=pol, scans='track', targets=1)
    dp_c1 = data.dumps
    for i in dp_c1:
        if i not in dp_tt or flags_1ch[i] == True:
            dp_c1 = list(dp_c1)
            dp_c1.remove(i)
            #print 'rm '+str(i)+' from dp_c1'
    dp_c1 = kf.deg_filter(dp_c1, ang_deg, sigma_level, n_iter)
    dp_c1 = np.array(dp_c1)
    dp_c1a = dp_c1[dp_c1 < dp_sb]
    dp_c1b = dp_c1[dp_c1 > dp_se]

    # +list(dp_c2a)+list(dp_c3a)+list(dp_c4a)
    dp_ca = list(dp_c0a)+list(dp_c1a)
    # +list(dp_c2b)+list(dp_c3b)+list(dp_c4b)
    dp_cb = list(dp_c0b)+list(dp_c1b)

    result = dp_ca, dp_cb, dp_c0a, dp_c1a, dp_c0b, dp_c1b

    #######################all have above#######################

    if fname in ['1551055211', '1551037708']:
        data.select(ants=ant, pol=pol, scans='track', targets=2)
        dp_c2 = data.dumps
        for i in dp_c2:
            if i not in dp_tt or flags_1ch[i] == True:
                dp_c2 = list(dp_c2)
                dp_c2.remove(i)
                #print 'rm '+str(i)+' from dp_c2'
        dp_c2 = kf.deg_filter(dp_c2, ang_deg, sigma_level, n_iter)
        dp_c2 = np.array(dp_c2)
        dp_c2a = dp_c2[dp_c2 < dp_sb]
        dp_c2b = dp_c2[dp_c2 > dp_se]

        data.select(ants=ant, pol=pol, scans='track', targets=3)
        dp_c3 = data.dumps
        for i in dp_c3:
            if i not in dp_tt or flags_1ch[i] == True:
                dp_c3 = list(dp_c3)
                dp_c3.remove(i)
                #print 'rm '+str(i)+' from dp_c3'
        dp_c3 = kf.deg_filter(dp_c3, ang_deg, sigma_level, n_iter)
        dp_c3 = np.array(dp_c3)
        dp_c3a = dp_c3[dp_c3 < dp_sb]
        dp_c3b = dp_c3[dp_c3 > dp_se]

        data.select(ants=ant, pol=pol, scans='track', targets=4)
        dp_c4 = data.dumps
        for i in dp_c4:
            if i not in dp_tt or flags_1ch[i] == True:
                dp_c4 = list(dp_c4)
                dp_c4.remove(i)
                #print 'rm '+str(i)+' from dp_c4'
        dp_c4 = kf.deg_filter(dp_c4, ang_deg, sigma_level, n_iter)
        dp_c4 = np.array(dp_c4)
        dp_c4a = dp_c4[dp_c4 < dp_sb]
        dp_c4b = dp_c4[dp_c4 > dp_se]
        # overwrite
        dp_ca = list(dp_c0a)+list(dp_c1a)+list(dp_c2a) + \
            list(dp_c3a)+list(dp_c4a)
        dp_cb = list(dp_c0b)+list(dp_c1b)+list(dp_c2b) + \
            list(dp_c3b)+list(dp_c4b)
        result = dp_ca, dp_cb, dp_c0a, dp_c1a, dp_c2a, dp_c3a, dp_c4a, dp_c0b, dp_c1b, dp_c2b, dp_c3b, dp_c4b
    data.select()  # recover after select!!!
    data.select(ants=ant, pol=pol)
    return result


def cal_dp_u(dp_tt, dp_ss):
    """

    Parameters
    ----------

    Returns
    -------

    """
    dp_u = list(dp_tt)+list(dp_ss)
    dp_u.sort()
    dp_u = np.array(dp_u)
    return dp_u


def label_nd_injection(fname, vis, timestamps, dp_ss, dump_period):

    # Validate filename
    fname = valid_filename(fname)

    ######set a jump limit###############
    f = 10.
    if fname in ['1555793534', '1551055211', '1551037708', '1555775533']:
        f = 10
    if fname == '1556120503':
        f = 2
    if fname == '1556052116':
        f = 15.

    ch_plot0 = 800  # only for edge detection

    lmax = abs(np.nanmax(vis[dp_ss, ch_plot0]))
    lmin = abs(np.nanmin(vis[dp_ss, ch_plot0]))
    lim = (lmax-lmin)/f

    #print lmax,lmin,lim

    mark = []
    nd_1 = []
    nd_1a = []
    nd_1b = []
    for i in range(1, len(timestamps)):
        if (np.abs(vis[i, ch_plot0])-np.abs(vis[i-1, ch_plot0]) > lim  # have a jump
            and vis[i, ch_plot0] > 0 and vis[i-1, ch_plot0] > 0
                and timestamps[i-1]-timestamps[0] - dump_period/2. not in mark):  # not jump the one before

            m = timestamps[i]-timestamps[0]-dump_period/2.
            #print i,m
            mark.append(m)  # for plot
            nd_1.append(i)  # on
            nd_1a.append(i)  # on1
            if i+1 < len(timestamps):
                nd_1.append(i+1)  # on
                nd_1b.append(i+1)  # on2
    return mark, nd_1, nd_1a, nd_1b, lmin, lmax


def gap_list(list):
    gap_list = []
    for i in range(1, len(list)):
        gap_list.append(list[i]-list[i-1])

    gap_list = np.array(gap_list)
    return gap_list


def gap_mode(list):
    gap_list = []
    for i in range(1, len(list)):
        gap_list.append(list[i]-list[i-1])

    gap_list = np.array(gap_list)
    mode = stats.mode(gap_list)[0][0]
    return mode


def cal_nd_wro_0(nd_1a):
    nd_gap_list = gap_list(nd_1a)
    nd_gap_mode = gap_mode(nd_1a)
    nd_wro_0 = np.where(nd_gap_list != nd_gap_mode)[0]
    return nd_gap_mode, nd_wro_0


def cal_t_line(fname, timestamps, nd_set, nd_cycle, dump_period):
    # t_line samples at start of injection, not needed for cal
    t_line = []
    nd_tb = nd_set-timestamps[0]
    nd_te = dump_period*len(timestamps)

    for t in np.arange(nd_tb, nd_te+dump_period, nd_cycle):
        t_line.append(t)

    # Adjustment caused by the diode lead time
    time_idxs = ds.get_diode_info(fname, 'time_range')
    t_line = t_line[time_idxs[0]:time_idxs[1]]
    return t_line


def call_nd_1a_param(fname, ds=None):
    """

    Parameters
    ----------
    fname : str
        Name of file to get diode info for.

    Returns
    -------
    nd_1a_gap : int
        ??

    nd_1a0 : int
        ??
    """
    nd_1a_gap = 10
    nd_1a0 = ds.get_diode_info(fname, 'offset')

    return nd_1a_gap, nd_1a0


def call_nd_1_list(fname, timestamps):
    nd_1a_gap, nd_1a0 = call_nd_1a_param(fname)
    nd_1aa = []
    for i in range(1000):
        a = nd_1a0 + i*nd_1a_gap
        if a >= 0:
            if a >= len(timestamps):
                break
            nd_1aa.append(a)

    nd_1bb = []
    for i in range(1000):
        a = nd_1a0 + 1 + i*nd_1a_gap
        if a >= 0:
            if a >= len(timestamps):
                break
            max = i
            nd_1bb.append(a)
    nd_11 = list(nd_1aa) + list(nd_1bb)
    nd_11.sort()

    nd_00 = []
    for i in range(len(timestamps)):
        if i not in nd_11:
            nd_00.append(i)

    assert(len(nd_00) + len(nd_11) == len(timestamps))

    # 0 = off, 1 = on, 1a = first sample in injection, 1b = second sample in injection
    return nd_1aa, nd_1bb, nd_11, nd_00


def cal_nds_list(dp_ss, nd_1a, nd_1b, nd_1, nd_0):  # dp_ss/dp_s
    nd_s0 = []
    for i in dp_ss:
        if i in nd_0:
            nd_s0.append(i)
    #print np.shape(nd_s0)

    nd_s1 = []
    nd_s1a = []
    nd_s1b = []
    for i in dp_ss:
        if i in nd_1:
            nd_s1.append(i)
        if i in nd_1a:
            nd_s1a.append(i)
        if i in nd_1b:
            nd_s1b.append(i)
    # same as fn above, but only during scan
    return nd_s1a, nd_s1b, nd_s1, nd_s0


def cal_ndt_list(dp_tt, nd_1a, nd_1b, nd_1, nd_0):  # dp_tt/dp_s
    nd_t0 = []
    for i in dp_tt:
        if i in nd_0:
            nd_t0.append(i)
    #print np.shape(nd_t0)

    nd_t1 = []
    nd_t1a = []
    nd_t1b = []
    for i in dp_tt:
        if i in nd_1:
            nd_t1.append(i)
        if i in nd_1a:
            nd_t1a.append(i)
        if i in nd_1b:
            nd_t1b.append(i)
    # same as above, but only during track
    return nd_t1a, nd_t1b, nd_t1, nd_t0


def cal_az_edge_list(az, dp_sb, dp_se):
    dp_az_min = []
    dp_az_max = []
    for i in range(dp_sb, dp_se+1):  # dp_w included
        if az[i] < az[i-1] and az[i] < az[i+1]:
            dp_az_min.append(i)
        if az[i] > az[i-1] and az[i] > az[i+1]:
            dp_az_max.append(i)
    dp_az_edge = list(dp_az_min+dp_az_max)
    dp_az_edge.sort()
    return dp_az_min, dp_az_max, dp_az_edge
