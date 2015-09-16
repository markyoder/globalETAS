'''
#globalETAS.py:
#author: mark yoder
#email: mryoder@ucdavis.edu, mark.yoder@gmail.com
#institution: dept. of physics, UC Davis
#
# summary:
# these scripts and modules constitute a second generation implementation of the 'mean-field' or
# 'rupture-length' variant of ETAS published by Yoder et al. (2014) in which initial ETAS (aftershock)
# productivity from known earthquake scaling relations, earthquake finite extents, and independent
# observational constraints.
#
# in this version, we adapt the model for large scale, global deployment. we:
# 1) start with the basic ETAS model from Yoder et al. (2014)
# 2) we, again, divide the earth into a large set of discrete locations (basically, a lattice), and calculate
#    the ETAS contribution (aftershock rate) of each earthquake to each lattice site.
# 3) we mitigate the computational requirements of the model by indexing the map production from the earthquake catalog
#    and limiting the spatil range of each earthquake according to its magnitude. for example, we calculate the
#    ETAS contribution of each earthquake to sites up to 5*L_r from the epicenter (and eventually, we'd like to be
#    smarter for large earthquakes with long ruptures, etc.).
#    -- the ETAS contribution drops off pretty quickly with distance, so this should do a pretty good job of reducing
#       computational requirements by skipping the calculation of very small contributions of very small earthquakes
#       a great distance from a given lattice site.
# 4) eventually add temporal indexing as well?
# 5) we're also interested in convolving ETAS with finite source models (aka, map ETAS to faults or some other
#    irregular geometry). we may do something with that here.
'''
#
import datetime as dtm
import matplotlib.dates as mpd
import pytz
tzutc = pytz.timezone('UTC')

import operator
import math
import random
import numpy
import scipy
#import scipy.optimize as spo
#import os
#from PIL import Image as ipp
import multiprocessing as mpp
#
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.mpl as mpl
#
#import shapely.geometry as sgp
#
from mpl_toolkits.basemap import Basemap as Basemap
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from geographiclib.geodesic import Geodesic as ggp
#
import ANSStools as atp
#
days2secs = 60.*60.*24.
year2secs = 60.*60.*24.*365.
deg2km=111.12
deg2rad = 2.0*math.pi/360.
alpha_0 = 2.28
#
calipoly = [[-129.0, 42.75], [-121.5, 29.5], [-113.5, 29.5], [-114.0, 35.6], [-119.3, 39.4], [-119.5, 42.75], [-129.0, 42.75]]
tz_utc = pytz.timezone('UTC')
#
#
class globalETAS_model(object):
	# guts of ETAS. this will contain a catalog, lattice sites, etc. graphical tools, plotting, etc. will likely
	# be kept elsewhere.
	def __init__(self, in_cat=None, lats=[], lons=[], mc=2.5, t_0=dtm.datetime(1990,1,1, tzinfo=tz_utc), t_now=dtm.datetime.now(tzutc)):
		#
		#  basically: if we are given a catalog, use it. try to extract mc, etc. data from catalog if it's not
		# externally provided. otherwise, get a catalog from ANSS or OH based oon the parameters provided.
		#
		lattice_sites = {}	# list of lattice site objects, and let's index it by... probably (i_x/lon, j_y/lat)
							# we might alternativel store this as a simple list [] or a base index/list (that
							# is re-sorting tolerant)) {i:[row]}, and then write indices:
							# index_lat_lon = {(lat, lon):lattice_sites_index}, {(x,y):lattice_sites_index}, etc.
							#
		#
		earthquake_catalog = {}	# and again, let's use an indexed structure, in the event that we are
								# using indices and (inadvertently? sort it?). this may be unnecessary.
								# that said, let's take the next step and return dict. type earthquake entries.

def make_ETAS_catalog(incat=None, lats=[32., 38.], lons=[-117., -114.], mc=2.5, date_range=['1990-1-1', None], grid_size=.1, D_fract=1.5, d_lambda=1.76):
	'''
	# fetch (or receive) an earthquake catalog. for each event, calculate ETAS relevant constants including rupture_length,
	# spatial orientation, elliptical a,b, etc. This catalog will then be provided to other processes to calculate ETAS rates
	'''
	#
	# handle dates:
	if date_range[1]==None: date_range[1] = dtm.datetime.now(pytz.timezone('UTC'))
	#
	for j,dt in enumerate(date_range):
		if isinstance(dt, str):
			date_range[j] = mpd.num2date(mpd.datestr2num(dt))
		#
		# other date-handling....
	#
	start_date = date_range[0]
	end_date = date_range[1]
	#
	if incat==None or (hasattr(incat, '__len__') and len(incat)==0):
		incat = atp.catfromANSS(lon=lons, lat=lats, minMag=mc, dates0=[start_date, end_date], Nmax=None, fout=None, rec_array=True)
	#
	# now, calculate etas properties....
	cols, formats = [list(x) for x in zip(*incat.dtype.descr)]
	
	#
	return incat

class Earthquake(object):
	# an Earthquake object for global ETAS. in parallel operations, treat this more like a bag of member functions than a data container.
	# pass an earthquake catalog list to a process; use Earthquake() to handle each earthquake event row.
	# include "local" ETAS calculation in this object; aka, each earthquake determines its ETAS range based on rupture length and other factors...
	# maybe.
	pass
