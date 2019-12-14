import datetime as dtm
import matplotlib.dates as mpd
import pytz
tzutc = pytz.timezone('UTC')
import h5py
#
#import operator
import math
import random
import numpy
import scipy
#import scipy.optimize as spo
from scipy import interpolate
import scipy.constants
import itertools
import sys
#import scipy.optimize as spo
import os
#import operator
#from PIL import Image as ipp
import multiprocessing as mpp
#
import matplotlib
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
#import json
#import pickle
#

import geopy.distance
#from geopy.distance import vincenty
#from geopy.distance import great_circle
#
#import shapely.geometry as sgp
os.environ['PROJ_LIB'] = '{}/anaconda3/share/proj'.format(os.getenv('HOME'))
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

from eq_params import *
#
from nepal_figs import *
import optimizers
#
import random
import geopy

from globalETAS import Earthquake
deg2km = 111.11
#
# TODO:
# maybe use numpy.digitize(x,bins) --> indices of bins into which x belong.??
# so set up the bins, then create indices with .digitize(), assign values
# with indices. I think this should be at least comparably as fast as requiring
# and assigning contiguous blocks. One trick will be to be to index on the dinmensionsm,
# so we don't have to evaluate the full 2/3D space (we can evaluate the 2/3 1D dimension
# . arrays).
#
class NETAS_block(Earthquake):
    # Note: it might make more sense to include an Earthquke(), or list of Earthquake()s rather than subclassing.
    #. For batched MPP, for example, it would be possible to pre-construct a complete array of sites for multiple
    #  inputs (Earthquake()s). On the other hand, that introduces some opportunity for memory mis-management, 
    #. mistakes, or requiring indexing (aka, you have earthquakes that are not close together, maybe non-
    #. contiguous lattice sites)...
    def __init__(self, times=None, times_to=None, spatial_intensity_threshold=1e-2, time_intensity_threshold=None, 
                 lon_phase=0., lat_phase=0., d_lon=.1, d_lat=.1,
                 *args, **kwargs):
        '''
        # @times: an iterable of times
        # @spatial_intensity_threshold: compute range by invertng spatial-omori; solve for 1/this.
        # @time_intensity_threshold: same, but for temporal distributionz. again, x_min ~ x0/this
        #
        # removing lon_0, lat_0. We'll map the min addresses (lon, lat) to indices in a parent/calling
        # . function or class. we will substitute these with phase offsets, which we'll usually set to zero.
        # .  note, these should lon_phase < d_lon, lat_phase < d_lat.
        ## @lon_0, lat_0 : zero-point for lon, lat bins, respectively.
        # @lon_phase, @lat_phase: phase to add to lon, lat bins, respectively.
        ## # parent array binning:
        #
        # TODO: what is the best parent binning convention? allow negative bins indices? Maybe we just
        #. leave that up to the parent object. I think that, by itself, negative index labels, centered
        #. on LL=[0,0] is most intuitive.
        # @d_lon, d_lat : lon, lat bin sizes, respectively
        #. bin indices are then j = int((lon - lon_0)/d_lon), k = int( (lat - lat_0)/d_lat)
        # 
        # parent __init__(): __init__(self, dict_or_recarray, transform_type='equal_area', transform_ratio_max=5.,
            ab_ratio_expon=.5, t_1=None, t_2=None)
        #
        '''
        #
        # initialize the Earthquake() part:
        super(NETAS_block, self).__init__(*args, **kwargs)
        #
        # handle intensity threshold inputs:
        # NOTE: None does not make sense with this object. if we want to do an unrestricted, fully connected
        # . nETAS, use a different gridding object (just pass the full lattice dims to all events).
        #if not spatial_intensity_threshold is None:
        spatial_intensity_threshold = min(spatial_intensity_threshold, 1./spatial_intensity_threshold)
        if not time_intensity_threshold is None:
            time_intensity_threshold = min(time_intensity_threshold, 1./time_intensity_threshold)
        #
        self.__dict__.update({key:val for key,val in locals().items() if not key in ('self', '__class__')})
        #
        # now, compute the spatial range and optionally the temporal range, or more specifically
        #   the lon, lat, time dimensions/indices. compute these as pairs: [[k_label, value], ...]
        #.  k_label: the index "label", or the index in a larger (global) lattice.
        #.  an we use their ordering for direct indexing.
        delta_lat = self.etas_range()/deg2km
        delta_lon = delta_lat*numpy.cos(self.lat*scipy.constants.degree)
        #
        # TODO: move this to a module level function?
        lat_min = int((self.lat - delta_lat + lat_phase)/d_lat)*d_lat
        lon_min = int((self.lon - delta_lon + lon_phase)/d_lon)*d_lon
        #
        #self.__dict__.update({'delta_lat':delta_lat, 'delta_lon':delta_lon, 'lat_min':lat_min, 'lon_min':lon_min})
        self.__dict__.update({key:val for key,val in locals().items() if not key in ('self', '__class__')})
        #
        # TODO: write these as functions to save memory? They will typically be called once
        #. and then return a 2D array of intensities.
        #lats = numpy.arange(lat_min, lat_min + 2*delta_lat, d_lat)
        #lats_index_labels = ((lats-lat_0)/d_lat).astype(int)
        #
        #lons = numpy.arange(lon_min, lon_min + 2*delta_lon, d_lon)
        #lons_index_labels = ((lons-lon_0)/d_lon).astype(int)
        #
    #
    @property
    def lats(self):
        # TODO: are we de-modularizing too much here? I think maybe we just compute lats and lons
        #. at this level. we leave the binning to the calling function. 
        #lats = numpy.arange(self.lat_min, self.lat_min + 2*self.delta_lat, self.d_lat)
        #lats_index_labels = ((self.lats-self.lat_0)/self.d_lat).astype(int)
        #
        #return numpy.core.records.fromarrays([lats_index_labels, lats], dtype=[('index', '>i8'), ('lat', '>f8')])
        return numpy.arange(self.lat_min, self.lat_min + 2*self.delta_lat, self.d_lat)
    #                           
    @property
    def lons(self):
        #lons = numpy.arange(self.lon_min, self.lon_min + 2*self.delta_lon, self.d_lon)
        #lons_index_labels = ((self.lons - self.lon_0)/self.d_lon).astype(int)
        #
        #return numpy.core.records.fromarrays([lons_index_labels, lons], dtype=[('index', '>i8'), ('lat', '>f8')])
        return numpy.arange(self.lon_min, self.lon_min + 2*self.delta_lon, self.d_lon)
    #
    def etas_range(self, spatial_intensity_threshold=None):
        spatial_intensity_threshold = spatial_intensity_threshold or self.spatial_intensity_threshold
        if spatial_intensity_threshold is None:
            return None
        #
        # TODO:
        # FIXME: 
        # what are the units for this formulation? i thought km, but...j
        # . this appears to actually be right, just a bit surprising how big r_0 is...
        # anyway, this may be wrong. or not useful, either untis or formula... or maybe just taht r_0 is 
        # . really big -- a lot bigger than rupture length (like 10x for an M7). We probably need to reevaluate
        #   our solutions for r_0. This could be an interesting ML exercise -- to ue the r_0 formulation that
        # . optimized predictability (or something).
        return self.r_0*(spatial_intensity_threshold**(-1./self.q) - 1.)
    #
    def etas_temporal_range(self, time_intensity_threshold=None):
        time_intensity_threshold = time_intensity_threshold or self.time_intensity_threshold
        if time_intensity_threshold is None:
            return None
        #
        return self.t_0*(time_intensity_threshold**(-1./self.p) - 1.)
    #
    def lons_lats_mgrid(self):
        return true_lon_lat(*numpy.meshgrid(self.lons, self.lats) )
    #
    def nETAS_block(self):
        #
        lons, lats = true_lon_lat(*numpy.meshgrid(self.lons, self.lats ) )
        #
        #print('*** DEBUG: shapes:: {}, {}'.format(lons.shape, lats.shape) )
        #print('*** DEBUG: shapes_true:: {}, {}'.format(true_lon(lons).shape, true_lat(lats).shape) )
        #
        return self.local_intensities(ts=self.times, ts_to=self.times_to, lons=lons, lats=lats )
    #
#
class NETAS_collector():
    # the whole enchillada:
    # 1) we'll take a catalog as an input for now.
    # 2) Either from inputs or catalog, determine spatial and temporal extents.
    # 3) Define binning translator, j,k,l = f(lon), g(lat), h(t)
    # 4) create or reference an HDF5, NetCDF, or other container object to aggregate the data
    # 5) Loop through catalog:
    # .   1) Compute a NETAS_block().netas_block() for each event
    # .   2) Bin-Align returned data; aggregate into HDF5-like object
    def __init__(self, catalog=None, h5_file=None,
                 lons=None, lats=None,
                 lon_min=-180., lon_max=180., lat_min=-90., lat_max=90., d_lon=.1, d_lat=.1,
                 times=None, times_to=None
                 ):
        # ... and skip the lon/lat phase; those values are in the min values of the sequences.
        # . lon_bin_phase=0., lat_bin_phase=0.
        # TODO: sort out some default behaviors, like inferring time, lat, lon extents from catalog,
        # . fetching a catalog if lats, lons, times specified, etc...
        # . for now, think modular, and let it break if need be.
        #
        if not os.path.isfile(h5_file): 
            # NOTE: this is not a comprehensive way to build lon, lat domains in a periodic
            # . domain (on a sphere). We probably just need to be specific -- ie, pass the 
            # . lons, lats sequendcs directly.
            if lats is None:
                lats = numpy.arange(lat_min, lat_max+d_lat, d_lat)
            if lons is None:
                lons = numpy.arange(lon_min, lon_max+d_lon, d_lon)
            #
            #
            no_return = create_netas_hdf5(fname_out=h5_file, times=times, times_to=times_to,
                                          lats=lats,
                                          lons=lons)
            #
        #
        # if an h5 file was passed, we'd pull lons, lats from its dimension collection.
        with h5py.File(h5_file, 'r') as h5in:
            lats = numpy.array(h5in['lats'])
            lons = numpy.array(h5in['lons'])
            #
        # NOTE: Do they need to be in order? Is there value in permitting out of order?
        # for now, skip sorting. later, we'll build in some failsafes.
        #lats.sort()
        #lons.sort()
        #
        #print('** DEBUG: ', type(lons) )
        #lon_min, lon_max = min(lons), max(lons)
        #lat_min, lat_max = min(lats), max(lats)
        lon_min, lon_max = lons[0::len(lons)-1]
        lat_min, lat_max = lats[0::len(lats)-1]
        #
        #d_lon = numpy.mean(numpy.diff(lons))
        # TODO: consisder a diagnostic where we evaluate all or a bunch of these
        #. and look for consistency.
        d_lon = lons[1] - lons[0]
        d_lat = lats[1] - lats[0]
            #
        #
        self.__dict__.update({key:val for key,val in locals().items() if not key in ('self', '__class__')})
        #
    #
    def compute_nETAS(self):
        for eq in self.catalog:
            block = NETAS_block(times=ts, **eq.input_parameters, spatial_intensity_threshold=.1)
            #block_netas = block.nETAS_block()
            #
            k_lons, x_lons = self.lon_to_bin(block.lons)
            k_lats, x_lats = self.lat_to_bin(block.lats)
            #
            with h5py.File(self.h5_file, 'a') as h5ap:
                h5ap[:k_lons:k_lats] += block.nETAS_block()
            
    #
    def lon_to_bin(self, lon):
        # return [[lon_bin_index, lon_bin_val], ...]
        #
        # int((self.lon - delta_lon + lon_phase)/d_lon)*d_lon
        #lon_binned = lon
        # from NETAS_block:
        # TODO: figure out how to more easily and consistently synch these up. we need
        # . uniform bin calculations across three interacting classes (hdf5 object, collector,
        # . NETAS_block). Probably, the thing to do is to use NETAS_block() as the primary
        #  (actually use its lat/lon binning function).
        # lat_min = int((self.lat - delta_lat + lat_phase)/d_lat)*d_lat
        # lon_min = int((self.lon - delta_lon + lon_phase)/d_lon)*d_lon
        #pass
        #
        #lon_bin_index = int((lon - self.lon_min)//self.d_lon)
        #
        #lon_bin_value = lon_bin_index*d_lon + self.lon_min
        #
        #return numpy.array([lon_bin_index, lon_bin_index*d_lon + self.lon_min ]).T
        return x_to_bin(lon, x_phase=0, d_x=self.d_lon, x_0=self.lon_min)
    #
    def lat_to_bin(self, lat):
        #lat_bin_index = int((lat - self.lat_min)//self.d_lat)
        ##lon_bin_value = lon_bin_index*d_lon + self.lon_min
        #
        #return numpy.array([lat_bin_index, lat_bin_index*d_lat + self.lat_min ]).T
        return x_to_bin(lat, x_phase=0, d_x=self.d_lat, x_0=self.lat_min)
#
def x_to_bin(x, x_phase=0., d_x=0.1, x_0=0.):
    # TODO: we still need a smart way to synch the lat/lon bins in all of these classes.
    #
    # return index and value, [[lat_index, lat_val], ...]
    # . note re: index, we might be computing the 
    # lat_min = int((self.lat - delta_lat + lat_phase)/d_lat)*d_lat
    #
    #xs = float(int((x+x_phase)/d_x)*d_x)
    #
    # xs*d_x gives us the raw (global) indexing. we might want a smaller
    # . subset. Nominally there are couple of ways to do this -- subtract x_0
    # . before or after we do the binning. Compute speed is an issue, but
    # . this will be performed on the 1D dimensions, so it should be a minor
    # . issue. note the int() operator is transcendental, so might not comply
    # . with standard commutative and associative rules.
    #
    #ks = int((x+x_phase)/d_x) - int((x_0+x_phase)/d_x)
    #
    #print('*** ndim: ', numpy.ndim(x), x)
    if numpy.ndim(x)==0:
        return numpy.array([(int((x+x_phase)/d_x) - int((x_0+x_phase)/d_x)), 
                             float(int((x+x_phase)/d_x)*d_x) ]).T
    else:
        #
        return numpy.core.records.fromarrays([(((x+x_phase)/d_x).astype('>i8') - int((x_0+x_phase)/d_x)).astype('>i8'), 
                             (((x+x_phase)/d_x).astype('>i8')*d_x).astype('>f8') ],
                                             dtype=[('k', '>i8'), ('x', '>f8')])
    # .astype([('k','>i8'),
    #                                      
#
# can we subclass an HDF5... i guess file class?
def create_netas_hdf5(fname_out='my_hdf5.h5', times=None, times_to=None, lats=None, lons=None):
    # NOTE: lats, lons, and time arrays might be better called "dimensions".
    # . they are the unique values.
    # TODO: handle incomplete or improperly formatted times_to inputs.
    if times_to is None:
        times_to = numpy.nan
    times_to = numpy.atleast_1d(times_to)
    #
    with h5py.File(fname_out, 'w') as fout:
        netas = fout.create_dataset("nETAS", ( len(times), len(lats), len(lons) ), dtype='d')
        #lats  = fout.create_dataset('lats', (len(lats),), dtype='>f8') 
        #
        # do we want to embrace the NetCDF 'dimension' convention?
        # dimensions = fout.create_group('dimensions')
        # lats = dimensions.create_dataset('lats', (len(lats),), data=lats) 
        lts     = fout.create_dataset('lats', (len(lats),), data=lats) 
        lns     = fout.create_dataset('lons', (len(lons),), data=lons)
        tmes    = fout.create_dataset('times', (len(times),), data=times)
        tmes_to = fout.create_dataset('times_to', (len(times_to),), data=times_to)
        #
    #
    return None  
#
def create_netas_hdf5_1D(fname_out='my_hdf5.h5', times=None, times_to=None, lats=None, lons=None):
    # NOTE: lats, lons, and time arrays might be better called "dimensions".
    # . they are the unique values.
    # TODO: handle incomplete or improperly formatted times_to inputs.
    if times_to is None:
        times_to = numpy.nan
    times_to = numpy.atleast_1d(times_to)
    #
    with h5py.File(fname_out, 'w') as fout:
        netas = fout.create_dataset("nETAS", ( len(times)*len(lats)*len(lons), ), dtype='d')
        netas.attrs.create('nETAS_shape', numpy.array([ len(times)*len(lats)*len(lons) ]), shape=None, dtype=None)
        #lats  = fout.create_dataset('lats', (len(lats),), dtype='>f8') 
        #
        # do we want to embrace the NetCDF 'dimension' convention?
        # dimensions = fout.create_group('dimensions')
        # lats = dimensions.create_dataset('lats', (len(lats),), data=lats) 
        lts     = fout.create_dataset('lats', (len(lats),), data=lats) 
        lns     = fout.create_dataset('lons', (len(lons),), data=lons)
        tmes    = fout.create_dataset('times', (len(times),), data=times)
        tmes_to = fout.create_dataset('times_to', (len(times_to),), data=times_to)
        #
    #
    return None
#
def true_lon(x):
    # lon can be corrected independently...
    return (x+180.)%360. - 180.
#
def true_lon_lat(X,Y):
    # lat introduces lon corrections...
    #
    # first latitude:
    Y_prime = (Y+90.)%360.
    m = Y_prime//180.
    #print('** DEBUG: ', m)
    #
    #print('*** ', m*180.)
    #print('*** ', (-1.)**m)
    #print('*** ', Y_prime%180.)
    #
    #Y_prime = m*180. + ((-1.)**m)*(Y_prime%180. )- 90.
    #
    #X_prime = ((X+180.) +m*180.)%360 - 180.
    #
    #return X_prime, Y_prime
    return ((X+180.) +m*180.)%360 - 180., m*180. + ((-1.)**m)*(Y_prime%180. )- 90.
#