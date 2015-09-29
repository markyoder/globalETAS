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
import rtree
from rtree import index
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

def make_ETAS_catalog(incat=None, lats=[32., 38.], lons=[-117., -114.], mc=2.5, date_range=['1990-1-1', None], grid_size=.1, D_fract=1.5, d_lambda=1.76, fit_factor=1.5):
	'''
	# fetch (or receive) an earthquake catalog. for each event, calculate ETAS relevant constants including rupture_length,
	# spatial orientation, elliptical a,b, etc. This catalog will then be provided to other processes to calculate ETAS rates
	'''
	#
	dmstar = 1.0
	dm_tau = 0.0		# constant sort of like dm_star but for temporal component. we usually treat this as 0, but don't really know...
	d_tau = 2.28
	b1=1.0
	b2=1.5
	b=b1
	q=1.5
	p=1.1
	#
	# some constants:
	l_15_factor = (2./3.)*math.log10(1.5)	# pre-calculate some values...
	km2lat = 1./111.1
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
	# first, make a spatial index of the catalog:
	anss_index = index.Index()
	[anss_index.insert(j, (rw['lon'], rw['lat'], rw['lon'], rw['lat'])) for j,rw in enumerate(incat)]
	
	# now, calculate etas properties....
	cols, formats = [list(x) for x in zip(*incat.dtype.descr)]
	#cols += ['L_r', 'r_0', 'chi', 'dt_r', 't_0', 'tau', 'strike_theta', 'strike_epsilon']			# add these parameters so we can switch them out (if we want)
	cols += ['L_r', 'r_0', 'chi', 'dt_r', 't_0', 'tau', 'e0', 'e1', 'v00', 'v01', 'v10', 'v11']
	formats += ['f8' for x in xrange(len(cols)-len(formats))]
	output_catalog = []
	#
	for k_cat, rw in enumerate(incat):
		#
		# tentatively, let's just make a catalog of this and calculate the rupture length, etc.
		# this way, we can specify (different) rupture length, etc. later.
		mag = rw['mag']
		L_r = 10.0**(.5*mag - d_lambda)
		dt_r = 10.0**(.5*mag - d_tau)
		lN_om = b*(mag - dmstar - mc)
		N_om = 10.0**lN_om
		#
		# some short-hand:
		D = D_fract
		#
		# self-similar formulation from yoder et al. 2014/15:
		# start with the (log of the) number of earthquakes/aftershocks in the rupture area:
		lN_chi = (2.0/(2.0 + D))*math.log10(1.0 + D/2.) + D*(mag - dmstar - mc)/(2.+D)
		N_chi = 10.**lN_chi
		#
		# mean linear density and spatial Omori parameters:
		linear_density = 2.0*N_chi/L_r
		r_0 = N_om*(q-1.)/N_chi
		chi = (r_0**(1.-q))/(N_om*(q-1.))
		#
		# temporal Omori parameters:
		rate_max = l_15_factor + d_tau - mag/6. - (dm_tau + mc)/3.
		t_0 = N_om*(p-1.)/rate_max		# note, however, that we can really use just about any value for t_0, so long as we are consistent with tau.
		tau = (t_0**(1.-p))/(N_om*(p-1.))
		#
		# now, guess the earthquake's orientation based on local seismicity. use earthquakes within some distance based on magnitude.
		# use a PCA type method or some sort of line-fit. there should be some good PCA libraries, otherwise it's easy to code.
		area_lat_min = max(-90., rw['lat'] - fit_factor * L_r*km2lat)
		area_lat_max = min(90., rw['lat'] + fit_factor * L_r*km2lat)
		# note: i think we have to catch cases where we cross the dateline (aka, into phi<-180, phi>180) and use two index ranges.
		area_lon_min = rw['lon'] - fit_factor * L_r*km2lat/math.cos(rw['lat']*deg2rad)
		area_lon_max = rw['lon'] + fit_factor * L_r*km2lat/math.cos(rw['lat']*deg2rad)
		#
		#print "rw: ", rw
		#print "window (%f) %s" % (mag, str((area_lat_min, area_lon_min, area_lat_max, area_lon_max)))
		included_indices = list(anss_index.intersection((area_lon_min, area_lat_min, area_lon_max, area_lat_max)))
		if area_lon_min<-180.:
			new_lon_min = area_lon_min%(180.)
			included_indices += list(anss_index.intersection((new_lon_min, area_lat_min, 180., area_lat_max)))
		#
		if area_lon_max>180.:
			new_lon_max = area_lon_max%(-180.)
			included_indices += list(anss_index.intersection((-180., area_lat_min, new_lon_max, area_lat_max)))
		#
		# now, get PCA (or some other method of fitting) for these events. we'll want to convert lat/lon to euclidean (or similar) positions.
		# at some point, we'll want to allow multiple options for fitting and transforms. for now, start by getting PCA eigen-vals/vectors.
		# note: if center_lat/lon=None (default), we subtract mean value to get PCA. if we specify center_lat/lon, we subtract those values.
		if len(included_indices)>10:
			eig_vals, eig_vecs = get_pca_cat(cat=[[incat['lon'][j], incat['lat'][j]] for j in included_indices], center_lat=rw['lat'], center_lon=rw['lon'], xy_transform=True)
		else:
			eig_vals = [1.0, 1.0]
			eig_vecs = [[1.0, 0.], [0., 1.0]]
		#
		# so now what? the linear algebra approach is fast and simple; trig. approach is more compact (6 variables, 2 variables), and
		# we want to use a structured array.
		#
		#output_catalog += [x for x in rw] + [L_r, r_0, chi, dt_r, t_0, tau, strike_theta, strike_epsilon]
		#
		output_catalog += [ [rw[0].tolist()] + list(rw)[1:] + [L_r, r_0, chi, dt_r, t_0, tau] + [x for x in eig_vals] + [x for x in eig_vecs[0]] + [x for x in eig_vecs[1]]]
		# (and strike_theta (or at least something like it -- setting aside conventions) = atan(eig_vecs[0][1]/eig_vecs[0][1])
		#
		#new_row = 
		
		#lt0ssim = 7.0*mag/6.0 - (2.0*mc/3.0) - dtau - dmstar - (2./3.)*math.log10(1.5) + (dmtau/3.0) + math.log10(p-1.0)
		#lr0ssim = mag*((6.0+D)/(4.0+2*D)) - (2.0/(2.0+D))*(dm + mc - math.log10((2.0+D)/2.0)) + math.log10(q-1.0) - 1.76 - math.log10(2.0)
		#r0ssim = 10**lr0ssim
		#t0ssim = 10.**lt0ssim
		#taussim = (t0ssim**(1.0-p))/(Nom*(p-1.0))
	
	#numpy.rec.array((rlist if len(rlist)>0 else [[]]), dtype=[('event_date', 'M8[us]'), ('lat','f8'), ('lon','f8'), ('mag','f8'), ('depth','f8'), ('event_date_float', 'f8')])
	#
	output_catalog = numpy.rec.array((output_catalog if len(output_catalog)>0 else [[]]), dtype = [('event_date', 'M8[us]')] + [(col, type(output_catalog[0][j+1]).__name__) for j,col in enumerate(cols[1:])])
	
	return output_catalog

def get_pca_cat(cat=[], center_lat=None, center_lon=None, xy_transform=True):
	# get pca for the input catalog. if center_lat/lon=None, then subtrac these values as the 'mean' (standard PCA). otherwis,
	# subtract the actual mean.
	# if xy_transform, then transform from lat/lon to x,y
	
	xy_factor=1.0
	if xy_transform:
		xy_factor = deg2km
	#
	if center_lat==None: center_lat = numpy.mean([x[1] for x in cat])*xy_factor
	if center_lon==None:
		center_lon = numpy.mean([x[0]*math.cos(x[1]*deg2rad) for x in cat])*xy_factor
	#
	# subtract off center:
	cat_prime = [[rw[0]*xy_factor*(math.cos(rw[1]*deg2rad) if xy_transform else 1.0)-center_lon, rw[1]*xy_factor-center_lat] for rw in cat]
	#
	# now, get eig_vals, eig_vectors:
	n_dof=len(cat_prime)-1
	cov = numpy.dot(numpy.array(zip(*cat_prime)),numpy.array(cat_prime))/n_dof
	#
	# for now, we can't assume hermitian or symmetric matrices, so use eig(), not eigh()
	#eig_vals, eig_vecs = numpy.linalg.eig(cov)
	return  numpy.linalg.eig(cov)

class Earthquake(object):
	# an Earthquake object for global ETAS. in parallel operations, treat this more like a bag of member functions than a data container.
	# pass an earthquake catalog list to a process; use Earthquake() to handle each earthquake event row.
	# include "local" ETAS calculation in this object; aka, each earthquake determines its ETAS range based on rupture length and other factors...
	# maybe.
	pass
