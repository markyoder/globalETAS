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
import bindex
#
days2secs = 60.*60.*24.
year2secs = 60.*60.*24.*365.
deg2km=111.31
deg2rad = 2.0*math.pi/360.
alpha_0 = 2.28
#
calipoly = [[-129.0, 42.75], [-121.5, 29.5], [-113.5, 29.5], [-114.0, 35.6], [-119.3, 39.4], [-119.5, 42.75], [-129.0, 42.75]]
tz_utc = pytz.timezone('UTC')

# a lookup dictionary for ditance types:
dist_calc_types_dict={}
dist_calc_types_dict.update({key:'spherical' for key in ['sph', 'spherical', 'sphere']})
dist_calc_types_dict.update({key:'cartesian' for key in ['cartesian', 'cart', 'xy', 'rect']})
dist_calc_types_dict.update({key:'geo' for key in ['geo', 'geodesic', 'geodesiclib', 'glib', 'g_lib']})
#
#
class globalETAS_model(object):
	# guts of ETAS. this will contain a catalog, lattice sites, etc. graphical tools, plotting, etc. will likely
	# be kept elsewhere.
	def __init__(self, catalog=None, lats=[32., 38.], lons=[-117., -114.], mc=2.5, d_x=10., d_y=10., bin_x0=0., bin_y0=0., etas_range_factor=5.0, t_0=dtm.datetime(1990,1,1, tzinfo=tz_utc), t_now=dtm.datetime.now(tzutc), transform_type='equal_area', calc_etas=True):
		'''
		#
		#  basically: if we are given a catalog, use it. try to extract mc, etc. data from catalog if it's not
		# externally provided. otherwise, get a catalog from ANSS or OH based oon the parameters provided.
		# note: we want to do this properly in an x-y projection, not lat-lon. our earthquakes will come back to us with lat/lon coordinates
		# so we'll need to translate. we'll also want lat/lon as well as x/y positions for the lattice sites, so we can get even spacing and accurate distances
		'''
		# dx=1., dy=1., x0=0., y0=0., leaf_type=float)
		self.lats=lats
		self.lons=lons
		self.mc=mc
		self.d_x=d_x
		self.d_y=d_y
		self.bin_x0=bin_x0
		self.bin_y0=bin_y0
		self.etas_range_factor = etas_range_factor
		self.t_0=t_0
		self.t_now=t_now
		self.transform_type='equal_area'
		#		
		self.lattice_sites = bindex.Bindex2D(dx=d_x, dy=d_y, x0=0., y0=0.)	# list of lattice site objects, and let's index it by... probably (i_x/lon, j_y/lat)
							# we might alternativel store this as a simple list [] or a base index/list (that
							# is re-sorting tolerant)) {i:[row]}, and then write indices:
							# index_lat_lon = {(lat, lon):lattice_sites_index}, {(x,y):lattice_sites_index}, etc.
							#
		#
		#earthquake_catalog = {}	# and again, let's use an indexed structure, in the event that we are
								# using indices and (inadvertently? sort it?). this may be unnecessary.
								# that said, let's take the next step and return dict. type earthquake entries.
		#
		if catalog==None:
			catalog = make_ETAS_catalog(incat=None, lats=lats, lons=lons, mc=mc, date_range=[t_0, t_now])	# and note there are other variables to consider...
		self.catalog = catalog
		#
		if calc_etas:
			# do the ETAS calculation. for each earthquake, figure out what bins it contributes to and calc. etas for each bin.
			#
			print "calc etas..."
			for quake in catalog:
				x,y = lon_lat_2xy(lon=quake['lon'], lat=quake['lat'])
				delta_x = quake['L_r']*etas_range_factor
				x_min, x_max = x-delta_x, x+delta_x
				y_min, y_max = y-delta_x, y+delta_x
				
				
				#
				# - choose an elliptical transform: equal-area, rotational, etc.
				# - calculate initial rate-density
				# - distribute over relevant sites:
				#    - for each site, D' = D_spherical * R'/R_xy	of course this will be approximately equal to R'.
				#
				for x_site,y_site in [[j,k] for j in numpy.arange(x_min, x_max, self.d_x) for k in numpy.arange(y_min, y_max, self.d_y)]:
					# so we make a square (or maybe a circle later on) and calc. etas at each site. use Bindex() to correct for any misalignment.
					#
					x_bin = self.lattice_sites.get_xbin_center(x_site)	# returns a dict like {'index':j, 'center':x}
					y_bin = self.lattice_sites.get_ybin_center(y_site)
					#
					bin_lonlat = xy2_lon_lat(x_bin['center'], y_bin['center'])
					#
					# now, get distance from earthquake and this bin, calc etas and update this site.
					# the basic model for distances will be like: 1) measure distance using proper spherical transform, 2) get the ellipitical transform from PCA,
					# renormalize
					#
					# currently available distances are {'cartesian', 'spherical', 'geodedic'} (and there is some correction for synonmys). for most accurate
					# distance calculations, use the geodetic option; spherical is a faster approximation; cartesian distances will be used for the elliptical transform.
					#
					eq_obj = Earthquake(quake, transform=self.transform_type)
					distances = dist_to(lon_lat_from=[quake['lon'], quake['lat']], lon_lat_to=bin_lonlat, dist_types=['cart', 'geo'], Rearth = 6378.1)
					#
					self.lattice_sites.add_to_bin(x=x_bin['center'], y=y_bin['center'], z=1.0/distances['geo'])
					
			#
			print "didn't really calculate anything, but spun the catalgo at least..."		
	#
	'''
	# i think these are all handled in the Bindex() object.
	def get_bin_id_x(self, x):
		return self.get_bin_id(x,dx=self.d_x, x0=self.bin_x0)
	#
	def get_bin_id_y(self, x):
		return self.get_bin_id(y,dx=self.d_y, x0=self.bin_y0)
	#
	def get_bin_id(self, x, dx=1., x0=None, bin_0=0):
		# return bin_id/index along a single axis.
		x0 = (x0 or self.bin_x0)
		#
		return int(bin_0) + int(round((x-x0)/dx))

	def bin2x(bin_num, dx=1., x0=0., bin_0=0):
		# return the position x of a bin along one axis.
		return (bin_num - bin_0)*dx + x0x
	'''

def make_ETAS_catalog(incat=None, lats=[32., 38.], lons=[-117., -114.], mc=2.5, date_range=['1990-1-1', None], D_fract=1.5, d_lambda=1.76, d_tau = 2.28, fit_factor=1.5, do_recarray=True):
	'''
	# fetch (or receive) an earthquake catalog. for each event, calculate ETAS relevant constants including rupture_length,
	# spatial orientation, elliptical a,b, etc. This catalog will then be provided to other processes to calculate ETAS rates
	'''
	#
	dmstar = 1.0
	dm_tau = 0.0		# constant sort of like dm_star but for temporal component. we usually treat this as 0, but don't really know...
	#d_tau = 2.28
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
	#
	# first, set up the numpy.recarray column formats:
	cols, formats = [list(x) for x in zip(*incat.dtype.descr)]
	#cols += ['L_r', 'r_0', 'chi', 'dt_r', 't_0', 'tau', 'strike_theta', 'strike_epsilon']			# add these parameters so we can switch them out (if we want)
	#cols += ['L_r', 'r_0', 'chi', 'dt_r', 't_0', 'tau', 'e0', 'e1', 'v00', 'v01', 'v10', 'v11']
	cols += ['L_r', 'r_0', 'chi', 'dt_r', 't_0', 'tau']
	formats += ['f8' for x in xrange(len(cols)-len(formats))]
	#
	my_dtype = [(cols[j], formats[j]) for j in xrange(len(cols))]
	#cols += ['e_vals', 'e_vecs']
	#formats += ['object', 'object']
	my_dtype += [('e_vals', '>f8', 2), ('e_vecs', '>f8', (2,2))]
	
	#print "dtype:",  my_dtype
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
		included_indices = list(anss_index.intersection((area_lon_min, area_lat_min, area_lon_max, area_lat_max)))
		# ... and if we wrap around the other side of the world...
		if area_lon_min<-180.:
			new_lon_min = area_lon_min%(180.)
			included_indices += list(anss_index.intersection((new_lon_min, area_lat_min, 180., area_lat_max)))
		if area_lon_max>180.:
			new_lon_max = area_lon_max%(-180.)
			included_indices += list(anss_index.intersection((-180., area_lat_min, new_lon_max, area_lat_max)))
		#
		# now, get PCA (or some other method of fitting) for these events. we'll want to convert lat/lon to euclidean (or similar) positions.
		# at some point, we'll want to allow multiple options for fitting and transforms. for now, start by getting PCA eigen-vals/vectors.
		# note: if center_lat/lon=None (default), we subtract mean value to get PCA. if we specify center_lat/lon, we subtract those values.
		# note: if there aren't any data, don't bother fitting (though we could just leave this to the PCA code)
		if len(included_indices)>10:
			eig_vals, eig_vecs = get_pca(cat=[[incat['lon'][j], incat['lat'][j]] for j in included_indices], center_lat=rw['lat'], center_lon=rw['lon'], xy_transform=True)
		else:
			eig_vals = [1.0, 1.0]
			eig_vecs = [[1.0, 0.], [0., 1.0]]
		#
		#
		output_catalog += [ list(rw) + [L_r, r_0, chi, dt_r, t_0, tau] + [eig_vals] + [eig_vecs] ]
		
		# (and strike_theta (or at least something like it -- setting aside conventions) = atan(eig_vecs[0][1]/eig_vecs[0][1])
		#
			
	if do_recarray: output_catalog = numpy.core.records.fromarrays(zip(*output_catalog), dtype=my_dtype)
	
	return output_catalog

def get_pca(cat=[], center_lat=None, center_lon=None, xy_transform=True):
	# get pca for the input catalog. if center_lat/lon=None, then subtrac these values as the 'mean' (standard PCA). otherwis,
	# subtract the actual mean.
	# if xy_transform, then transform from lat/lon to x,y
	
	xy_factor=1.0
	if xy_transform:
		xy_factor = deg2km
	#
	if center_lat==None: center_lat = numpy.mean([x[1] for x in cat])*xy_factor
	if center_lon==None:
		center_lon = numpy.mean([x[0]*(math.cos(x[1]*deg2rad) if xy_transform else 1.0) for x in cat])*xy_factor
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

def lon_lat_2xy(lon=0., lat=0.):
	#
	# super simple lat/lon -> xy mapper. this is primarily for purposes of indexing, so don't allow any free parameters (like x0, y0, etc.).
	#
	return [(lon)*math.cos(lat*deg2rad)*deg2km, (lat)*deg2km]


def xy2_lon_lat(x=0., y=0.):
	#
	# super simple xy -> lat/lon mapper
	lat = y/deg2km
	#
	if abs(lat)==90.:
		lon=0.
	else:
		lon = x/(deg2km*math.cos(lat*deg2rad))
	#
	return [lon, lat]

def dist_to(lon_lat_from=[0., 0.], lon_lat_to=[0.,0.], dist_types=['spherical'], Rearth = 6378.1):
	# general purpose distance getter.
	if dist_types in (None, [], 'all'): dist_types = dist_calc_types_dict.keys() 
	types_dict = dist_calc_types_dict
	#
	# spherical, cartesian, geo
	return_distances = {}
	dist_types = [types_dict.get(key, key) for key in dist_types]
	#
	if 'spherical' in dist_types:
		return_distances['spherical'] = spherical_dist(lon_lat_from, lon_lat_to)
	if 'cartesian' in dist_types:
		d_lon    = lon_lat_from[0]-lon_lat_to[0]
		mean_lat = .5*(lon_lat_from[1]+lon_lat_to[1])
		d_lat    = lon_lat_from[1]-lon_lat_to[1]
		return_distances['cartesian'] = ( ((d_lon)*math.cos(mean_lat*deg2rad)*deg2km)**2. + (d_lat*deg2km)**2.)**.5
	if 'geo' in dist_types:
		# .Inverse(lat1, lon1, lat2, lon2)
		g1=ggp.WGS84.Inverse(lon_lat_from[1], lon_lat_from[0], lon_lat_to[1],lon_lat_to[0])
		return_distances['geo'] = g1['s12']/1000.0
	#
	if len(return_distances)==1: return_distances = return_distances.items()[0]
	#
	return return_distances
	

def spherical_dist(lon_lat_from=[0., 0.], lon_lat_to=[0.,0.], Rearth = 6378.1):
	# Geometric spherical distance formula...
	# displacement from inloc...
	# inloc is a vector [lon, lat]
	# return a vector [dLon, dLat] or [r, theta]
	# return distances in km.
	#
	# also, we need to get the proper spherical angular displacement (from the parallel)
	#
	#Rearth = 6378.1	# km
	deg2rad=2.0*math.pi/360.
	#
	# note: i usually use phi-> longitude, lambda -> latitude, but at some point i copied a source where this is
	# reversed. oops. so just switch them here.
	# phi: latitude
	# lon: longitude
	#
	#phif  = inloc[0]*deg2rad
	#lambf = inloc[1]*deg2rad
	#phis  = self.loc[0]*deg2rad
	#lambs = self.loc[1]*deg2rad
	
	phif  = lon_lat_to[1]*deg2rad
	lambf = lon_lat_to[0]*deg2rad
	phis  = lon_lat_from[1]*deg2rad
	lambs = lon_lat_from[0]*deg2rad
	#
	#print 'phif: ', phif
	#print 'lambf: ', lambf
	#
	dphi = (phif - phis)
	dlambda = (lambf - lambs)
	#this one is supposed to be bulletproof:
	sighat3 = math.atan( math.sqrt((math.cos(phif)*math.sin(dlambda))**2.0 + (math.cos(phis)*math.sin(phif) - math.sin(phis)*math.cos(phif)*math.cos(dlambda))**2.0 ) / (math.sin(phis)*math.sin(phif) + math.cos(phis)*math.cos(phif)*math.cos(dlambda))  )
	R3 = Rearth * sighat3
	#
	return R3

def griddata_plot_xyz(xyz, n_x=None, n_y=None):
	#
	# eventually, we're going to need to put this stuff together to produce contours, etc. here's a test script that uses griddata() to, well... grid up the data.
	# also see scipy.interpolate.griddata() -- which i'm guessing is called by matplotlib.mlab.griddata(), but who knows. the calling signature for the scipy library is
	# a little bit different and allows method options ('nearest', 'linear', 'cubic') as opposed to ('nn', 'linear') for mlab version.
	#
	min_x, max_x = min(xyz['x']), max(xyz['x'])
	min_y, max_y = min(xyz['y']), max(xyz['y'])
	#
	padding = .05
	n_x = (n_x or int(math.sqrt(len(xyz))/10.))
	n_y = (n_y or int(math.sqrt(len(xyz))/10.))
	#
	xi = numpy.linspace(min_x-abs(min_x)*padding, max_x + abs(max_x)*padding, n_x)
	yi = numpy.linspace(min_y-abs(min_x)*padding, max_y + abs(max_x)*padding, n_y)
	#
	zi = matplotlib.mlab.griddata(xyz['x'], xyz['y'], numpy.log10(xyz['z']), xi, yi, interp='nn')
	#
	plt.figure(0)
	plt.clf()
	cs = plt.contourf(xi,yi,zi,15,cmap=plt.cm.rainbow,vmax=abs(zi).max(), vmin=-abs(zi).min())
	plt.colorbar()

class Earthquake(object):
	# an Earthquake object for global ETAS. in parallel operations, treat this more like a bag of member functions than a data container.
	# pass an earthquake catalog list to a process; use Earthquake() to handle each earthquake event row.
	# include "local" ETAS calculation in this object; aka, each earthquake determines its ETAS range based on rupture length and other factors...
	# maybe.
	#def __init__(self, event_date=None, lat=0., lon=0., mag=None, eig_vals=[1., 1.], eig_vecs=[[1.,0.], [0., 1.]], transform_type='equal_area'):
	def __init__(self, dict_or_recarray, transform_type='equal_area'):
		if isinstance(dict_or_recarray, dict): self.__dict__.update(dict_or_recarray)
		if hasattr(dict_or_recarray, 'dtype'):
			# recarray.
			self.__dict__.update({key:dict_or_recarray[key] for key in dict_or_recarray.dtype.fields.keys()})
			
		
	#
	def spherical_dist(to_lon_lat=[]):
		return shperical_dist(lon_lat_from=[self.lon, self.lat], lon_lat_to=to_lon_lat)
		
