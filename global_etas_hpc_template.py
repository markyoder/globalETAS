import matplotlib as mpl
import pylab as plt
mpl.use('Agg')

import datetime as dtm
import matplotlib.dates as mpd
import pytz
tzutc = pytz.timezone('UTC')

#import operator
import math
import random
import numpy
import scipy
import scipy.optimize as spo
from scipy import interpolate
import itertools
import sys
#import scipy.optimize as spo
import os
#import operator
#from PIL import Image as ipp
import multiprocessing as mpp
#
from mpl_toolkits.mplot3d import Axes3D
import json
import pickle
#
import geopy
import geopy.distance
#from geopy.distance import vincenty
#from geopy.distance import great_circle
#
#import shapely.geometry as sgp
# not sure we need this all the time...
#os.environ['PROJ_LIB'] = '{}/anaconda3/share/proj'.format(os.getenv('HOME'))
#os.environ['PROJ_LIB'] = '{}/anaconda3/share/proj'.format(os.getenv('HOME'))
#
from mpl_toolkits.basemap import Basemap as Basemap
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from geographiclib.geodesic import Geodesic as ggp
#

#import ANSStools as atp
from yodiipy import ANSStools as atp
#
import contours2kml
import globalETAS as gep

#import global_etas_auto as ggep

#from eq_params import *
#
#from nepal_figs import *
#import optimizers
#
######
#
#
if __name__ == '__main__':
    # get some input parameters...
    import argparse
    #import multiprocessing as mpp
    #
    # see argparse reff: https://docs.python.org/3.3/library/argparse.html
    #
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('layer_indices', metavar='N', type=int, nargs='+',
                        help='indices of layers to process')
    # set up some parameters, etc.
    #
    # event was some time on the 24th or maybe late the 23rd. this, plus defaults, should find the event:
    #to_dt = dtm.datetime(2016,8,25, tzinfo=pytz.timezone('UTC'))
    to_dt = dtm.datetime.now(pytz.timezone('UTC'))
    #to_dt = dtm.datetime(2019,7,5,0,0,0, tzinfo=pytz.timezone('UTC'))
    #
    Lr_factor = 10.
    # define these from the t_now in the actual etas object, in the event that we load it from pickle,
    #  rather than calc it here.
    #f_path = '/home/myoder/Dropbox/Research/etas/italy_2016_10/etas_{}'.format(to_dt)
    #f_root = 'etas_2016'
    #

    t0 = dtm.datetime.now(pytz.timezone('UTC'))
    t_ms = t0
    #
    # sacramento:
    lat0 = 35.705
    lon0 = -117.506
    #
    #ll_sacramento = (lon0, lat0)

    #m0 = 7.8

    d_lat=2.
    d_lon=2.
    #
    lats = [lat0-d_lat, lat0+d_lat]
    lons = [lon0-d_lon, lon0+d_lon]
    #to_dt = t0-dtm.timedelta(hours=2)
    #to_dt = dtm.datetime.now(pytz.utc)
    #
    #etas = ggep.auto_etas(to_dt=to_dt, Lr_factor=Lr_factor, dt_0=5)
    #italy_prams = {'do_recarray': True, 'D_fract': 1.5,
    #                't_0':dtm.datetime(1990, 1, 1, 0, 0, tzinfo=pytz.timezone('UTC')),
    #                't_now':to_dt,
    #                'lats': [42.,43.5], 'p': 1.1, 'b1': 1.0, 'mc': 2.5, 'q': 1.5,
    #                'lons': [12.,15.], 'dmstar': 1.0, 'b2': 1.5, 'd_tau': 2.28,
    #                'incat': None, 'fit_factor': 2.0, 'd_lambda': 1.76}
    eq_prams = {'do_recarray': True, 'D_fract': 1.5,
                   't_0':dtm.datetime(1990, 1, 1, 0, 0, tzinfo=pytz.timezone('UTC')),
                   't_now':to_dt, 't_future':None ,
                   'lats': lats, 'p_cat': 1.1, 'b1': 1.0, 'mc': 2.5, 'q_cat': 1.5,
                   'p_etas':1.1, 'q_etas':1.5,
                   'lons': lons, 'dmstar': 1.0, 'b2': 1.5, 'd_tau': 2.28,
                   'incat': None, 'fit_factor': 2.0, 'd_lambda': 1.76, 'etas_range_padding':1.5,
                'etas_range_factor':30.0, 'ab_ratio_expon':.25 }
    #eq_prams.update({'mc':3.0, 'd_lat':.25, 'd_lon':.25})
    #
    mycat = None
    #
    #mycat = atp.catfromANSS(lon=lons, lat=lats, minMag=2.5,
    mycat = atp.cat_from_comcat(lon=lons, lat=lats, minMag=2.5,
                                dates0=[dtm.datetime(2005,1,1, tzinfo=pytz.timezone('UTC')),
                                        dtm.datetime.now(pytz.timezone('UTC'))],
                                Nmax=None, fout=None, rec_array=True)
    #                        dates0=[dtm.datetime(2005,1,1, tzinfo=tzutc), None], Nmax=None, fout=None, rec_array=True)
    #
    # NOTE: default behavior is to grab as many CPUs as possible:
    n_cpu = 4
    mycat = gep.make_ETAS_catalog_mpp(incat=mycat, n_cpu=n_cpu)
    #
    # we can adjust paranmeters in the dictionary like this:
    eq_prams['t_now'] = dtm.datetime.now(pytz.timezone('UTC'))
    # eq_prams['lats'] = [lat0 - 1., lat0 + 1.]
    # eq_prams['lons'] = [lon0 - 1., lon0 + 1.]
    #
    # eq_prams['lats'] = [lat0 - 2., lat0 + 2.]
    # eq_prams['lons'] = [lon0 - 2., lon0 + 2.]
    #
    # be careful with this parameter on the tool servers, and other shared resources.
    n_cpu=4
    #n_cpu = 2*mpp.cpu_count()
    #n_cpu=5
    etas = gep.ETAS_mpp(n_cpu=n_cpu, catalog=mycat, **eq_prams)
    #
    # get some data to name this thing:
    event_name = 'Ridgecrest_July_2019'
    #f_path = '/home/myoder/Dropbox/Research/etas/{}/etas_{}'.format(event_name, etas.t_now)
    #f_path = '{}/data_export/Research/etas/{}/etas_{}'.format(os.getenv('HOME'), event_name, etas.t_now)
    #
    #f_path = '{}/Dropbox/Research/etas/{}/etas_{}'.format(os.getenv('HOME'), event_name, etas.t_now)
    f_path = '{}/Mazama_outputs/etas/{}/etas_{}'.format(os.getenv('HOME'), event_name, etas.t_now)
    f_root = 'etas_{}_2019_07'.format(event_name)
    #
    fg=plt.figure(0, figsize=(12,10))
    ax=plt.gca()
    # lats_map= , lons_map=
    etas.make_etas_contour_map(n_contours=25, fignum=0, map_resolution='f', alpha=.3, ax=ax)
    #
    #mainshock = sorted(etas.catalog, key=lambda rw: rw['mag'])[-1]
    #print('mainshock: ', mainshock)
    # get mainshock. it's an m>6 event in the last week or so... this is subjective.
    # if we just look for the biggest event, we get the L'Aquila event, so we'll need to be more creative...
    # or just specify it.

    mainshock = etas.catalog[-1]
    for j,eq in enumerate(reversed(etas.catalog)):
        #print('*** ', pytz.utc.localize(eq['event_date'].astype(dtm.datetime)))
        if pytz.utc.localize(eq['event_date'].astype(dtm.datetime))<etas.t_now-dtm.timedelta(days=180): break
        if eq['mag']>mainshock['mag']:
            mainshock = eq
    #
    #
    #
    print('ms: ', mainshock, mainshock['lon'], mainshock['lat'])
    x,y = etas.cm(mainshock['lon'], mainshock['lat'])
    #
    #print('mm: ', max(etas.catalog['mag']))
    #
    # let's get everything m>6 in the last 6 months?
    m6s = [rw for rw in etas.catalog if rw['mag'] >= 6.
           and pytz.utc.localize(rw['event_date'].astype(dtm.datetime))>to_dt-dtm.timedelta(days=180)]
    #
    # plot mainshock:
    dt = mainshock['event_date'].astype(dtm.datetime)
    dt=t0
    dt_str = '{}-{}-{} {}:{}:{}'.format(dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second)
    #etas.cm.plot([x], [y], latlon=False, marker='*', color='r', ms=16, zorder=11,
    #                   label='m={}, {}'.format(mainshock['mag'], dt_str))
    #etas.cm.plot([lon0], [lat0], latlon=False, marker='*', color='r', ms=16, zorder=11,
    #                   label='m={}, {}'.format(m0, dt_str))

    ax.set_title('ETAS: {}, {}\n\n'.format(event_name, etas.t_now), size=16)
    for j,m6 in enumerate(m6s):
        clr = colors_[j%len(colors_)]
        #
        dt = m6['event_date'].astype(dtm.datetime)
        dt_str = '{}-{}-{} {}:{}:{}'.format(dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second)
        etas.cm.scatter(m6['lon'], m6['lat'], s=2*(m6['mag']+2.), edgecolors=clr,
                        c='none', marker='o', zorder=11, label='m={}, {}'.format(m6['mag'], dt_str))

    #x,y = etas.cm(*ll_sacramento)
    #etas.cm.scatter([x],[y], marker='o', s=18, edgecolors='r', c='r',
    #                    label='Sacramento')
    t_cat = mpd.date2num(etas.t_now-dtm.timedelta(days=15))
    print('tt: ', t_cat, etas.catalog['event_date'][0], type(etas.catalog['event_date'][0]))
    k=0

    # for j,rw in enumerate(etas.catalog):
    #     if mpd.date2num(rw['event_date'].astype(dtm.datetime))<t_cat: continue
    #     k+=1
    #     clr = colors_[k%len(colors_)]
    #     #
    #     dt = rw['event_date'].astype(dtm.datetime)
    #     dt_str = '{}-{}-{} {}:{}:{}'.format(dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second)
    #     #etas.cm.scatter(rw['lon'],rw['lat'], s=3*(rw['mag']+12.), edgecolors=clr,
    #     #                      c='none', marker='o', zorder=11, label='m={}, {}'.format(rw['mag'], dt_str))
    #     etas.cm.plot(rw['lon'],rw['lat'], ms=2.*(rw['mag']+2.), color=clr,
    #                           marker='o', zorder=11, label='m={}, {}'.format(rw['mag'], dt_str), latlon=True)
    plt.gca().legend()
    #
    etas.export_kml(os.path.join(f_path, '{}_{}.kml'.format(f_root, str(etas.t_now).replace(' ', '_'))))
    etas.export_xyz(os.path.join(f_path, '{}_{}.xyz'.format(f_root, str(etas.t_now).replace(' ', '_'))))
    fg.savefig(os.path.join(f_path, '{}_{}.png'.format(f_root, str(etas.t_now).replace(' ', '_'))))
    fg2.savefig(os.path.join(f_path, '{}_{}_with_equakes.png'.format(f_root, str(etas.t_now).replace(' ', '_'))))


    #
    # so this worked, once upon a time, but breaks maybe when the script does not run cleanly all the way through?
    with open (os.path.join(f_path, '{}_{}.pkl'.format(f_root, str(etas.t_now).replace(' ', '_'))), 'wb') as fpkl:
        pickle.dump(etas, fpkl)
