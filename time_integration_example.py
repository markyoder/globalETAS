#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import datetime as dtm
import matplotlib.dates as mpd
import pytz
tzutc = pytz.timezone('UTC')
#
import numpy as np
import os
import sys
#
import matplotlib.pyplot as plt
#
import globalETAS as gep


if __name__=='__main__':
    #
    kwds={}
    pargs=[]
    for arg in sys.argv[1:]:
        arg.replace(',', '')
        # module name is the first argument.
        if '=' in arg:
            kwds.update(dict([arg.split('=')]))
        else:
            pargs+=[arg]
        #
    #
    region = pargs[0]
    
#==============================================================================
# ~~~~~~~~~~ ETAS Parameters and Calc ~~~~~~~~~~
#==============================================================================
    regions = ['nepal', 'chile', 'sichuan', 'tohoku', 'newzealand', 'sumatra', 'iquique', 'swnz', 'hokkaido']
    #regind = 0
    #region = regions[regind]
    
    # This toggles doing time integration (returns aftershock numbers in each cell) 
    #  vs no integration (returns aftershock rate-density at t_now)    
    bTimeInt = 1
    
    if region == 'nepal':
        # End date of catalog (mainshock of study)
        #6:11:25, 7.8
        t_now = dtm.datetime(2015, 4, 25, 6, 13, 00, tzinfo=tzutc)
        #Location (centered on lat0, lon0)
        lat0, lon0  = (28.2305, 84.7314)
    if region == 'chile':
        #06:34:11, 8.8
        t_now = dtm.datetime(2010, 2, 27, 6, 40, 00, tzinfo=tzutc)
        lat0, lon0  = (-35.9089, -72.7327)
    if region == 'sichuan':
        #06:28:01, 7.9
        t_now = dtm.datetime(2008, 5, 12, 6, 35, 00, tzinfo=tzutc)
        lat0, lon0  = (31.586, 104.032)
    if region == 'tohoku':
        #2011-3-11 05:46:24, 9.1
        t_now = dtm.datetime(2011, 3, 11, 5, 48, 00, tzinfo=tzutc)
        lat0, lon0  = (38.297, 142.373)
    if region == 'newzealand':
        #2016-11-13 11:02:56 (UTC), 7.8
        t_now = dtm.datetime(2016, 11, 13, 11, 5, 00, tzinfo=tzutc)
        lat0, lon0  = (-42.737, 173.054)
    if region == 'sumatra':
        #2004-12-26 00:58:53 (UTC), 9.1
        t_now = dtm.datetime(2004, 12, 26, 1, 00, 00, tzinfo=tzutc)
        lat0, lon0  = (3.295, 95.982)
    if region == 'iquique':
        #2014-04-01 23:46:47 (UTC), 8.2
        t_now = dtm.datetime(2014, 4, 1, 23, 48, 00, tzinfo=tzutc)
        lat0, lon0  = (-19.610, -70.769)
    if region == 'hokkaido':
        #2003-09-25 19:50:06 (UTC), 8.3
        t_now = dtm.datetime(2003, 9, 25, 19, 51, 00, tzinfo=tzutc)
        lat0, lon0 = (41.815, 143.910)
    if region == 'ecuador':
        # 2016-04-16 23:58:36 (UTC), 7.8
        t_now = dtm.datetime(2016, 4, 16, 23, 59, 00, tzinfo=tzutc)
        lat0, lon0 = (0.382, -79.922)
    if region == 'illapel':
        # 2015-09-16 22:54:32 (UTC), 8.3
        t_now = dtm.datetime(2015, 9, 16, 22, 56, 00, tzinfo=tzutc)
        lat0, lon0 = (-31.573, -71.674)
    if region == 'gujarat':
        #2001-01-26 03:16:40 (UTC), 7.7
        t_now = dtm.datetime(2001, 1, 26, 3, 18, 00, tzinfo=tzutc)
        lat0, lon0 = (23.419, 70.232)
    if region == 'awaran':
        #2013-09-24 11:29:47 (UTC), 7.7
        t_now = dtm.datetime(2013, 9, 24, 11, 31, 00, tzinfo=tzutc)
        lat0, lon0 = (26.951, 65.501)
    
    # Starting date of catalog (set to one year previous)
    t0 = t_now-dtm.timedelta(365)
    #
    # Check for time integration option; if so, provide an upper time integration limit one month ahead of t_now
    t_int_lim = None
    if bTimeInt:
        # Date out to which we want to integrate the time component of the Omori distribution
        t_int_lim = t_now + dtm.timedelta(30)
    #
    #
    radius_lat=5.
    radius_lon=5.
    lats = [lat0-radius_lat, lat0+radius_lat]
    lons = [lon0-radius_lon, lon0+radius_lon]
    
    eq_prams = {'do_recarray': True, 'D_fract': 1.5,
                   't_0':t0,
                   't_now':t_now,
                   't_int_limit': t_int_lim,
                   'cat_len': None,
                   'lats': lats, 'p': 1.1, 'b1': 1.0, 'mc': 2.5, 'q': 1.5,
                   'lons': lons, 'dmstar': 1.0, 'b2': 1.5, 'd_tau': 2.28,
                   'incat': None, 'etas_fit_factor': 2.0, 'd_lambda': 1.76, 'etas_range_padding':5., 
                   'transform_ratio_max':2.0, 'ab_ratio_expon':1}
    

    this_catalog = None
    
    # This class initializes by fetching the catalog if one is not provided, and then performs all the ETAS calculations.
    #   This is the meat of the ETAS script.  
    #TODO: Fix errors with more than 1 proc.  Only use 1 till then
    
    etas = gep.ETAS_mpp(n_cpu=1, catalog=this_catalog, **eq_prams)
    
    
    #==============================================================================
    # ~~~~~~~~~ ETAS Contour Plot with overlaid large EQ locations ~~~~~~~~~~~
    #==============================================================================
    plt.close('all')
    
    mainshock = etas.catalog[-1]
    for j,eq in enumerate(reversed(etas.catalog)):
        #print('*** ', pytz.utc.localize(eq['event_date'].astype(dtm.datetime)))
        if pytz.utc.localize(eq['event_date'].astype(dtm.datetime))<etas.t_now-dtm.timedelta(days=180): break
        if eq['mag']>mainshock['mag']:
            mainshock = eq
    print("mainshock found")
    
    
    fg=plt.figure(0, figsize=(12,10))
    ax=plt.gca()
    etas.make_etas_contour_map(n_contours=25, fignum=0, map_resolution='f', alpha=.3, ax=ax)
    
    print('ms: coords = ({:0.3f}, {:0.3f}), mag = {}'.format(mainshock['lon'], mainshock['lat'], mainshock['mag']))
    x,y = etas.cm(mainshock['lon'], mainshock['lat'])
    
    #
    # let's get everything m>6 in the last 6 months?
    m6s = [rw for rw in etas.catalog if rw['mag'] >= 5.5 
           and pytz.utc.localize(rw['event_date'].astype(dtm.datetime))>t_now-dtm.timedelta(days=120)]
    #
    # plot mainshock:
    dt = mainshock['event_date'].astype(dtm.datetime)
    dt_str = '{}-{}-{} {}:{}:{}'.format(dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second)
    etas.cm.plot([x], [y], latlon=False, marker='*', color='r', ms=4, zorder=11,
                       label='m={}, {}'.format(mainshock['mag'], dt_str))
    ax.set_title('ETAS, {}\n\n'.format(etas.t_now), size=16)
    for j,m6 in enumerate(m6s):
        dt = m6['event_date'].astype(dtm.datetime)
        dt_str = '{}-{}-{} {}:{}:{}'.format(dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second)
        etas.cm.scatter(m6['lon'], m6['lat'], s=3*(m6['mag']+12.), c='none', marker='o',
                        zorder=11, label='m={}, {}'.format(m6['mag'], dt_str))
        #
    plt.gca().legend()
    
    
    
    #==============================================================================
    # ~~~~~~~~~ File Output ~~~~~~~~~~~~
    #==============================================================================
    
    # Paths for file output
    if bTimeInt:
        f_path = 'etas_outputs/{}_tInt_etas'.format(region)
        f_root = 'etas_ab{}_tInt_{}'.format(str(eq_prams['transform_ratio_max']).replace('.', '-'), region)
    else:
        f_path = 'etas_outputs/{}_rateden_etas_{}'.format(region, etas.t_now)
        f_root = 'etas_ab{}_rateden_{}'.format(str(eq_prams['transform_ratio_max']).replace('.', '-'), region)
    
    if not os.path.isdir(f_path): os.makedirs(f_path)
    
    
    # TODO: we want the datetime part of the filename to come from the etas object itself, for purposes of
    # integrity. BUT, we want this script to be a bit more portable, so we should replace all the etas
    # references/object name to just 'etas'
    #
    etas.export_kml(os.path.join(f_path, '{}.kml'.format(f_root)))
    etas.export_xyz(os.path.join(f_path, '{}.xyz'.format(f_root)))
    fg.savefig(os.path.join(f_path, '{}.png'.format(f_root)))













