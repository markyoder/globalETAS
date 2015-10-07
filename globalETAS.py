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
import itertools
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
#
import ANSStools as atp
import bindex
#
import rtree
from rtree import index

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
#dist_calc_types_dict.update({key:'cartesian2' for key in ['cartesian2']})
dist_calc_types_dict.update({key:'geo' for key in ['geo', 'geodesic', 'geodesiclib', 'glib', 'g_lib']})
# and this one to return the x,y components for cartesian measurements:
dist_calc_types_dict.update({key:'dx_dy' for key in ['dxdy', 'dx_dy']})
#
#
# Classes:
class globalETAS_model(object):
	# guts of ETAS. this will contain a catalog, lattice sites, etc. graphical tools, plotting, etc. will likely
	# be kept elsewhere.
	# questions: use equal distance or equal angel bin spacing? we can nominally do a spherical->cartesian transformation like: x = lon*cos(lat)*deg2km
	#  but it behaves sort of funny when we try to calculate disanctes between points... which we don't necessaryily do; we map (x,y) -> (lon, lat) and then
	#  do proper geodetic transformations, but 1) it becomes easy to make mistakes, and 2) the only gain is that we have fewer points in the higher latitudes.
	#  as long as 1) we have sufficient spatial resolution in the lower latitudes and 2) the problem remains computational tractable (equal spacing can reduce
	#  the grid size by maybe halfish??), there's no real harm in just using (a much simpler) lon,lat lattice with equal angular spacing.
	#
	#def __init__(self, catalog=None, lats=[32., 38.], lons=[-117., -114.], mc=2.5, d_x=10., d_y=10., bin_x0=0., bin_y0=0., etas_range_factor=5.0, t_0=dtm.datetime(1990,1,1, tzinfo=tz_utc), t_now=dtm.datetime.now(tzutc), transform_type='equal_area', transform_ratio_max=5., calc_etas=True):
	def __init__(self, catalog=None, lats=[32., 38.], lons=[-117., -114.], mc=2.5, d_lon=.1, d_lat=.1, bin_lon0=0., bin_lat0=0., etas_range_factor=5.0, etas_fit_factor=1.0, t_0=dtm.datetime(1990,1,1, tzinfo=tz_utc), t_now=dtm.datetime.now(tzutc), transform_type='equal_area', transform_ratio_max=5., calc_etas=True):
		'''
		#
		#  basically: if we are given a catalog, use it. try to extract mc, etc. data from catalog if it's not
		# externally provided. otherwise, get a catalog from ANSS or OH based oon the parameters provided.
		# note: we want to do this properly in an x-y projection, not lat-lon. our earthquakes will come back to us with lat/lon coordinates
		# so we'll need to translate. we'll also want lat/lon as well as x/y positions for the lattice sites, so we can get even spacing and accurate distances
		'''
		# dx=1., dy=1., x0=0., y0=0., leaf_type=float)
		# load all the input parameters into class variables (note i believe locals() is all currently defined variables in the function scope,
		# so this syntax must be executed at the very begining of the function to work properly):
		self.__dict__.update(locals())
		#
		# and handle some specific cases...
		if isinstance(t_now, float):
			t_forecast = t_now
		elif isinstance(t_now, numpy.datetime64):
			t_forecast = mpd.date2num(t_now.tolist())
		else:
			t_forecast = mpd.date2num(t_now)
		'''
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
		self.transform_type=transform_type
		self.transform_ratio_max=transform_ratio_max
		'''
		#		
		self.lattice_sites = bindex.Bindex2D(dx=d_lon, dy=d_lat, x0=0., y0=0.)	# list of lattice site objects, and let's index it by... probably (i_x/lon, j_y/lat)
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
			catalog = make_ETAS_catalog(incat=None, lats=lats, lons=lons, mc=mc, date_range=[t_0, t_now], fit_factor=etas_fit_factor)	# and note there are other variables to consider...
		self.catalog = catalog
		#
		if calc_etas:
			# do the ETAS calculation. for each earthquake, figure out what bins it contributes to and calc. etas for each bin.
			#
			print "calc etas..."
			for quake in catalog:
				eq = Earthquake(quake, transform_type=self.transform_type, transform_ratio_max=self.transform_ratio_max)
				# calculate the bins within range of this earthquake:
				# nominally, the proper thing to do now (or at lest one proper approach) is to compute a proper geodesic north and east to get the lat/lon bounds
				# near the poles, however, this will create a problem, that is a bit difficult to track, where the geodesic reaches over the pole and to lat<90
				# on the other side. if we just use a recta-linear approximation, we can use [max(-90, lat-d_lat), min(90, lat+d_lat)]
				#x,y = lon_lat_2xy(lon=quake['lon'], lat=quake['lat'])
				#
				# range of influence:
				delta_lat = eq.L_r*etas_range_factor/deg2km
				if abs(eq.lat)==90.:
					delta_lon=180.
				else:
					delta_lon = delta_lat/math.cos(eq.lat*deg2rad)
				#
				lon_min, lon_max = eq.lon - delta_lon, eq.lon + delta_lon
				lat_min, lat_max = eq.lat - delta_lat, eq.lat + delta_lat
				#
				# - choose an elliptical transform: equal-area, rotational, etc.
				# - calculate initial rate-density
				# - distribute over relevant sites:
				#    - for each site, D' = D_spherical * R'/R_xy	of course this will be approximately equal to R'.
				#
				#for lon_site,lat_site in [[j,k] for j in numpy.arange(lon_min, lon_max, d_lon) for k in numpy.arange(lat_min, lat_max, d_lat)]:
				for lon_site, lat_site in itertools.product(numpy.arange(lon_min, lon_max, d_lon), numpy.arange(lat_min, lat_max, d_lat)):
					# so we make a square (or maybe a circle later on) and calc. etas at each site. use Bindex() to correct for any misalignment.
					#
					lon_bin = self.lattice_sites.get_xbin_center(lon_site)	# returns a dict like {'index':j, 'center':x}
					lat_bin = self.lattice_sites.get_ybin_center(lat_site)
					#
					#print "lon-lat bins: ", lon_bin, lat_bin
					#
					#bin_lonlat = xy2_lon_lat(x_bin['center'], y_bin['center'])
					#bin_lonlat=[lon_bin, lat_bin]
					#
					# now, calculate local ETAS intensity from this eq at this site...
					# (this should be defined in the Earthquake() object:
					# note: at this time, we have only the self-similar type local intensity. eventually, we need to split this up. nominally,
					# this can occur at the Earthquake level, so -- if we so desire, different earthquakes can have differend spatial distributions.
					#
					local_intensity = eq.local_intensity(t=t_forecast, lon=lon_bin['center'], lat=lat_bin['center']) # anything else?
					#
					#self.lattice_sites.add_to_bin(x=x_bin['center'], y=y_bin['center'], z=1.0/distances['geo'])
					#
					#self.lattice_sites.add_to_bin(x=lon_bin['center'], y=lat_bin['center'], z=local_intensity)	# note: if we give x,y instead of the bin index, bindex
					self.lattice_sites.add_to_bin(bin_x=lon_bin['index'], bin_y=lat_bin['index'], z=local_intensity)
					
					
			#
			print("finished calculating ETAS")	
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
#
class Earthquake(object):
	# an Earthquake object for global ETAS. in parallel operations, treat this more like a bag of member functions than a data container.
	# pass an earthquake catalog list to a process; use Earthquake() to handle each earthquake event row.
	# include "local" ETAS calculation in this object; aka, each earthquake determines its ETAS range based on rupture length and other factors...
	# maybe.
	#def __init__(self, event_date=None, lat=0., lon=0., mag=None, eig_vals=[1., 1.], eig_vecs=[[1.,0.], [0., 1.]], transform_type='equal_area'):
	def __init__(self, dict_or_recarray, transform_type='equal_area', transform_ratio_max=5., ab_rato_expon=.5):
		#
		# first, load all indexed elements from the input dict or recarray row as class members:
		if isinstance(dict_or_recarray, dict): self.__dict__.update(dict_or_recarray)
		if hasattr(dict_or_recarray, 'dtype'):
			# recarray.
			self.__dict__.update({key:dict_or_recarray[key] for key in dict_or_recarray.dtype.fields.keys()})
		#
		# elliptical transform bits:
		self.transform_type = transform_type
		self.transform_ratio_max = transform_ratio_max
		self.ab_rato_expon = ab_rato_expon
		#
		# check for float datetime...
		if not hasattr(self, 'event_date_float'): self.event_date_float = mpd.date2num(self.event_date.tolist())
		#self.event_date_float_secs = self.event_date_float*days2secs
		#
		#self.ab_ratio = min(transform_ratio_max, max(self.e_vals)/min(self.e_vals))
		self.ab_ratio_raw = max(self.e_vals)/min(self.e_vals)
		self.set_transform()
		#
		#####
		# now, make some preliminary calculations, namely the peak spatial density and maybe some omori constants? anything we'll calculate
		# again and again...
		#
		self.spatial_intensity_factor=1.0		# corrects for local aftershock density in rotational type transforms.
	#	
	@property
	def lat_lon(self):
		return [self.lat, self.lon]
	@property
	def lon_lat(self):
		return [self.lon, self.lat]		
	#
	def set_transform(self, e_vals=None, e_vecs=None, transform_type=None, transform_ratio_max=None, ab_rato_expon=None):
		'''
		# note: it might be better to not allow custom-runs (aka not allow parameterized calls to this function; always use native member values).
		#
		# define the elliptical transformation for calculating intensities.
		# transform_type: elliptical transformation type = { 'equal_area', 'rotation'}
		# transform_ratio_max: maximum allowable eigenvalue (and transformed axes length) ratio.
		'''
		e_vals = (e_vals or self.e_vals)
		e_vecs = (e_vecs or self.e_vecs)
		transform_ratio_max = float(transform_ratio_max or self.transform_ratio_max)
		transform_type = (transform_type or self.transform_type)
		ab_rato_expon = (ab_rato_expon or self.ab_rato_expon)
		#
		ab_ratio = min(transform_ratio_max, (max(e_vals[0]/e_vals[1], e_vals[1]/e_vals[0]))**ab_rato_expon)	# note: this **.5 on the e0/e1 value is quasi-arbitrary.
		#																								# ok, so what specifically is e0/e1 supposed to be? stdev, var?
		#																								# in any case, it seems to be pretty huge almost all the time, so
		#																								# let's put a fractional power on it...
		self.ab_ratio=ab_ratio
		#
		if transform_type=='equal_area':
			# so that pi*r^2 = pi*ab = pi*b*(ab_ratio)*b
			self.e_vals_n = [math.sqrt(ab_ratio), 1./math.sqrt(ab_ratio)]
			self.spatial_intensity_factor = 1.0
			#
		#	
		elif transform_type=='rotation':
			# set the *small* eigen-value --> 1; scale the large one.
			#self.e_vals_n = [min(transform_ratio_max, x/min(e_vals)) for x in e_vals]
			self.e_vals_n  = [ab_ratio, 1.]
			self.spatial_intensity_factor = min(self.e_vals_n)/max(self.e_vals_n)		# area of "core" (rupture area) is reduced, so intensity increases.
			#
		else:
			return self.get_transform(e_vals=e_vals, e_vecs=e_vecs, transform_type='equal_area')
		#
		E = self.e_vals_n
		self.T = numpy.dot([[E[0], 0.], [0., E[1]]], e_vecs)
		self.T_inverse = numpy.dot([[E[0], 0.], [0., E[1]]], e_vecs.transpose())
	#
	def elliptical_transform_inverse(self, lon=0., lat=0.):
		return self.elliptical_transform(lon=lon, lat=lat, T=self.T_inverse, invert=False)
	#
	def elliptical_transform_inverse_prime(self, lon=0., lat=0.):
		return self.elliptical_transform(lon=lon, lat=lat, T=self.T_inverse, invert=True)
	#
		#
	def elliptical_transform_prime(self, lon=0., lat=0.):
		'''
		#wrapper around elliptical_transform()
		'''
		return self.elliptical_transform(lon=lon, lat=lat, T=None, invert=True)
	#
	def elliptical_transform(self, lon=0., lat=0., T=None, invert=False):
		'''
		# return elliptical transform of position (lon, lat). submit the raw x,y values; we'll subtract this Earthquake's position here...
		# ... and let's just return all the information: {'R':radial_dist (spherical or geo), 'R_prime':transformed distance, 'lon':lon,'lat':, 
		# 'lon_prime':transf_lon, 'lat_prime':tranf_lat}
		#
		# "invert" invert the transformation, specifically False: x' = T dot x
		#                                                  True:  x' = x dot T
		#  basically, to we rotate the frame or the data... i think. in practice, when you do a 1/R field (so the distance R at some point is transformed to R', use False.
		#   if you want to draw a transformed ellipse (aka, relocate points from a circle to their corresponding locations on an ellipse), ust True.
		#
		#use module function:
		#def dist_to(lon_lat_from=[0., 0.], lon_lat_to=[0.,0.], dist_types=['spherical'], Rearth = 6378.1):
		'''
		#
		if T==None: T=self.T
		#
		# first, get the actual distance (and related data) to the point:
		dists = dist_to(lon_lat_from=[self.lon, self.lat], lon_lat_to=[lon, lat], dist_types=['geo', 'xy', 'dx_dy'])
		R = dists['geo']	# the ['s12'] item from Geodesic.Inverse() method, and convert to km...
		#
		dx,dy = dists['dx'], dists['dy']	# cartesian approximations from the dist_to() method; approximate cartesian distance coordinates from e_quake center in general FoR.
		#
		if invert:
			dx_prime, dy_prime = numpy.dot([dx,dy],T)
		else:
			dx_prime, dy_prime = numpy.dot(T, [dx,dy])	# rotated and dilated into the elliptical FoR...
		#R_prime = R*((dx_prime*dx_prime + dy_prime*dy_prime)/(dx*dx+dy*dy))**.5	# use this to mitigate artifacts of the spherical transform: R_prime = R_geo*(R_prime_xy/R_xy)
		
		R_prime = R*numpy.linalg.norm([dx_prime, dy_prime])/numpy.linalg.norm([dx,dy])
		
		#
		dists.update({'R':R, 'R_prime':R_prime, 'dx':dx, 'dy':dy, 'dx_prime':dx_prime, 'dy_prime':dy_prime})
		#return {'R':R, 'R_prime':R_prime, 'dx':dx, 'dy':dy, 'dx_prime':dx_prime, 'dy_prime':dy_prime}
		return dists
	#
	def spherical_dist(self, to_lon_lat=[]):
		return shperical_dist(lon_lat_from=[self.lon, self.lat], lon_lat_to=to_lon_lat)
	#
	def local_intensity(self, t=None, lon=None, lat=None, p=None, q=None):
		'''
		# calculate local ETAS density-rate (aka, earthquakes per (km^2 sec)
		# take time in days.
		'''
		#print "inputs: ", t, lon, lat, p, q
		p = (p or self.p)
		q = (q or self.q)
		# calculate ETAS intensity.
		# note: allow p,q input values to override defaults. so nominally, we might want to use one value of p (or q) to solve the initial ETAS rate density,
		# but then make a map using soemthing else. for example, we might use p=0 (or at least p<p_0) to effectivly modify (eliminate) the time dependence.
		#
		delta_t = (t - self.event_date_float)*days2secs
		if delta_t<0.:
			return 0.
		#
		et = self.elliptical_transform(lon=lon, lat=lat)
		#
		orate = 1.0/(self.tau * (self.t_0 + delta_t)**self.p)
		#
		# for now, just code up the self-similar 1/(r0+r) formulation. later, we'll split this off to allow different distributions.
		# radial density (dN/dr), and as i recall this is normalized so that int(dN/dr)_0^inf --> 1.0
		radial_density = (q-1.0)*(self.r_0**(q-1.0))*(self.r_0 + et['R_prime'])**(-q)
		#
		# ... and this is distributed along an elliptical contour. we could approximate with a circle, but what we really want is to distribute along the ellipse
		# that is orthogonal to our R-R' transformation. fortunately, they have the same radius. calculating the radius of an ellipse is hard. we an exact solution
		# (that uses a "special" function) and an approximation (that should be faster).
		# note major,semi-major axes are R*e_1, R*e_2 (doesn't matter which is which)
		
		#circumf = ellipse_circumference_exact(a=self.e_vals_n[0]*et['R'], b=self.e_vals_n[1]**et['R'])
		circumf = ellipse_circumference_approx1(a=self.e_vals_n[0]*et['R'], b=self.e_vals_n[1]*et['R'])
		#
		spatialdensity = radial_density/circumf
		#
		
		return spatialdensity
#
class Shape(object):
	# helper base class for shapes...
	def __init__(self):
		# anything to do for general case?
		pass
	#
	def plot_poly(self, fignum=0, do_clf=True, poly_len=100, ax=None, **kwargs):
		#
		# some default plotting bits:
		defaults = {'color':'b', 'marker':'.', 'ls':'-', 'lw':1.5}
		for key,val in defaults.iteritems():
			if not kwargs.has_key(key): kwargs[key]=val
		#
		plt.figure(fignum)
		if do_clf: plt.clf()
		if ax==None: ax=plt.gca()
		#
		#if self.poly==None or len(self.poly==0):
		#	self.poly(poly_len)
		#
		ax.plot(*zip(*self.poly()), **kwargs)
#
class Circle(Shape):
	# could be derived from Ellipse,and we could derive from Shape() (which does not exist), but we're just going to code them...
	def __init__(self, R=1.0, x=0., y=0.):
		#
		self.R=R
		self.x=x
		self.y=y
		#
	#
	@property
	def area(self):
		return math.pi*self.R**2.

	@property
	def circumference(self):
		
		return 2.0*math.pi*self.R
	#
	def poly(self, n_points=100):
		d_theta = 2.0*math.pi/n_points
		poly = [[self.x+self.R*math.cos(theta), self.y + self.R*math.sin(theta)] for theta in numpy.arange(0., 2.0*math.pi+d_theta, d_theta)]
		# there's probably a smarter way to do this...
		#print "theta: %f" % (self.theta)
		#
		return poly
	#
#
class Ellipse(Shape):
	def __init__(self, a=1.0, b=.5, ab_ratio=None, theta=0.):
		if a==None and b==None: a,b = 1.0, .5
		if not (a==None and b==None): ab_ratio=a/float(b)
		#
		if a==None: a=b*ab_ratio
		if b==None: b=a/ab_ratio
		#
		self.a = a
		self.b = b
		#
		if theta>1.1*math.pi*2.0: theta*=deg2rad
		self.theta=theta
		#
		self.ab_ratio = ab_ratio
		self.h = ((a-b)/(a+b))**2
		#self.polygon=None
	#
	@property
	def area(self):
		return math.pi*self.a*self.b

	@property
	def circumference_exact(self):
		
		return math.pi*(self.a+self.b)*scipy.special.hyp2f1(-.5, -.5, 1.0, self.h)
	
	@property
	def circumference_approx1(self):
		# there are two good approximations from Ramanujan (see wikipedia); this is one of them...
		#
		return math.pi*(self.a+self.b)*(1. + 3.*self.h/(10 + math.sqrt(4. - 3.*self.h)))
	#
	def poly(self, n_points=100):
		d_theta = 2.0*math.pi/n_points
		poly = [[self.a*math.cos(theta), self.b*math.sin(theta)] for theta in numpy.arange(0., 2.0*math.pi+d_theta, d_theta)]
		# there's probably a smarter way to do this...
		print "theta: %f" % (self.theta)
		if self.theta!=0.:
			#poly = zip(*numpy.dot( [[math.cos(self.theta), -math.sin(self.theta)],[math.sin(self.theta), math.cos(self.theta)]], zip(*poly)))
			poly = numpy.dot(poly, zip(*[[math.cos(self.theta), -math.sin(self.theta)],[math.sin(self.theta), math.cos(self.theta)]]))
		#
		return poly
	#
	def plot_poly(self, fignum=0, do_clf=True, poly_len=100):
		plt.figure(fignum)
		if do_clf: plt.clf()
		#
		#if self.poly==None or len(self.poly==0):
		#	self.poly(poly_len)
		#
		plt.plot(*zip(*self.poly()), color='b', marker='.', ls='-', lw='1.5')
	#	
#

# Working scripts
# make a class out of this; derive from recarray; keep universal valuess like p,q,dmstar, etc. as member variables.
def make_ETAS_catalog(incat=None, lats=[32., 38.], lons=[-117., -114.], mc=2.5, date_range=['1990-1-1', None], D_fract=1.5, d_lambda=1.76, d_tau = 2.28, fit_factor=1.5, p=1.1, q=1.5, dmstar=1.0, b1=1.0,b2=1.5, do_recarray=True):
	'''
	# fetch (or receive) an earthquake catalog. for each event, calculate ETAS relevant constants including rupture_length,
	# spatial orientation, elliptical a,b, etc. This catalog will then be provided to other processes to calculate ETAS rates
	'''
	#
	# save a copy of the input parameters to ammend to the output array as metadata. we can do this in lieu of creating a derived class -- sort of a (sloppy?) shortcut...
	if incat==None:
		# use all the local inputs:
		meta_parameters = locals().copy()
	else:
		meta_parameters = {}
		if hasattr('meta_parameters', incat): incat.meta_parameters.copy()
		meta_parameters.update({'D_fract':D_fract, 'd_lambda':d_lambda, 'd_tau':d_tau, 'fit_factor':fit_factor, 'p':p, 'q':q, 'dmstar':dmstar, 'b1':b1, 'b2':b2})
		
	#
	#dmstar = 1.0
	dm_tau = 0.0		# constant sort of like dm_star but for temporal component. we usually treat this as 0, but don't really know...
	#d_tau = 2.28
	#b1=1.0
	#b2=1.5
	b=b1
	#q=1.5
	#p=1.1
	#
	# some constants:
	l_15_factor = (2./3.)*math.log10(1.5)	# pre-calculate some values...
	km2lat = 1./deg2km		# deg3km defined globally
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
	cols += ['L_r', 'r_0', 'chi', 'dt_r', 't_0', 'tau', 'dmstar', 'p', 'q']
	formats += ['f8' for x in xrange(len(cols)-len(formats))]
	#
	my_dtype = [(cols[j], formats[j]) for j in xrange(len(cols))]
	#cols += ['e_vals', 'e_vecs']
	#formats += ['object', 'object']
	my_dtype += [('e_vals', '>f8', 2), ('e_vecs', '>f8', (2,2)), ('N_eig_cat','>f8')]
	
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
		if len(included_indices)>10 or True:
			eig_vals, eig_vecs = get_pca(cat=[[incat['lon'][j], incat['lat'][j]] for j in included_indices], center_lat=rw['lat'], center_lon=rw['lon'], xy_transform=True)
		else:
			eig_vals = [1.0, 1.0]
			eig_vecs = [[1.0, 0.], [0., 1.0]]
		#
		#
		output_catalog += [ list(rw) + [L_r, r_0, chi, dt_r, t_0, tau, dmstar, p, q] + [eig_vals, eig_vecs, len(included_indices)] ]
		
		# (and strike_theta (or at least something like it -- setting aside conventions) = atan(eig_vecs[0][1]/eig_vecs[0][1])
		#
			
	if do_recarray:
		output_catalog = numpy.core.records.fromarrays(zip(*output_catalog), dtype=my_dtype)
		output_catalog.meta_parameters = meta_parameters
	
	return output_catalog
#
# helper scripts
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
	#print 'cov: ', cov
	try:
		return  numpy.linalg.eig(cov)
	except:
		return numpy.array([1.,1.]), numpy.identity(2)

def lon_lat_2xy(lon=0., lat=0.):
	#
	# super simple lat/lon -> xy mapper. this is primarily for purposes of indexing, so don't allow any free parameters (like x0, y0, etc.).
	# note: this formulation is good for basically indexing lat and lon bins, but we can't accurately calculate cartesian distances directly
	# from this mapping. but i think we can still use the mapping to get equal-spaced lattice sites.
	#
	if hasattr(lon, '__len__'): lon,lat = lon
	#
	return [(lon)*math.cos(lat*deg2rad)*deg2km, (lat)*deg2km]


def xy2_lon_lat(x=0., y=0.):
	#
	# super simple xy -> lat/lon mapper
	#
	if hasattr(x,'__len__'): x,y=x	# accidentally passed an xy=[x,y] list...
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
	#
	'''
	if 'cartesian2' in dist_types or 'dx_dy' in dist_types:
		# note: we can't calculate xy distances directly from the lat_lot2xy mapping, so be careful...
		xy_from = lon_lat_2xy(*lon_lat_from)
		xy_to   = lon_lat_2xy(*lon_lat_to)
		#
		dx = xy_to[0]-xy_from[0]
		dy = xy_to[1]-xy_from[1]
		#
		return_distances['cartesian2']=(dx*dx+dy*dy)**.5
		#
	'''
	#
	if 'cartesian' in dist_types or 'dx_dy' in dist_types:
		d_lon    = lon_lat_to[0]-lon_lat_from[0]
		mean_lat = .5*(lon_lat_from[1]+lon_lat_to[1])
		d_lat    = lon_lat_to[1]-lon_lat_from[1]
		#
		# we might use a geodesic.Direct method to get these distances (and sort of a dumb way to get the parity):
		#dx = (-1. if d_lon<0. else 1.)*ggp.WGS84.Inverse(lon_lat_from[1], lon_lat_from[0], lon_lat_from[1], lon_lat_from[0] + d_lon)['s12']/1000
		#dy = (-1. if d_lat<0. else 1.)*ggp.WGS84.Inverse(lon_lat_from[1], lon_lat_from[0], lon_lat_from[1] + d_lat, lon_lat_from[0])['s12']/1000.
		dx = ((d_lon)*math.cos(mean_lat*deg2rad)*deg2km)
		dy = (d_lat*deg2km)
		if 'cartesian' in dist_types:
			return_distances['cartesian'] = ( dx**2. + dy**2.)**.5
		if 'dx_dy' in dist_types:
			return_distances.update({'dx':dx, 'dy':dy})
		#
	if 'geo' in dist_types:
		# .Inverse(lat1, lon1, lat2, lon2)
		g1=ggp.WGS84.Inverse(lon_lat_from[1], lon_lat_from[0], lon_lat_to[1],lon_lat_to[0])
		#return_distances['geo'] = g1['s12']/1000.0
		return_distances.update({'geo':g1['s12']/1000., 'azi1':g1['azi1'], 'azi2':g1['azi2']})
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

def ellipse_circumference_exact(a=None, b=None):
	h = ((a-b)/(a+b))**2.
	return math.pi*(a+b)*scipy.special.hyp2f1(-.5, -.5, 1.0, h)
	
def ellipse_circumference_approx1(a=None, b=None):
	# there are two good approximations from Ramanujan (see wikipedia); this is one of them...
	#
	h = ((a-b)/(a+b))**2.
	return math.pi*(a + b)*(1. + 3.*h/(10 + math.sqrt(4. - 3.*h)))

def ab_ratio_distribution(e_vals, **kwargs):
	x = [max(x[0]/x[1], x[1]/x[0]) for x in e_vals]
	#
	plt.figure(0)
	plt.clf()
	plt.hist(x, bins=min(250, len(x)/10), histtype='step', **kwargs)

def griddata_plot_xyz(xyz, n_x=None, n_y=None):
	#
	# eventually, we're going to need to put this stuff together to produce contours, etc. here's a test script that uses griddata() to, well... grid up the data.
	# also see scipy.interpolate.griddata() -- which i'm guessing is called by matplotlib.mlab.griddata(), but who knows. the calling signature for the scipy library is
	# a little bit different and allows method options ('nearest', 'linear', 'cubic') as opposed to ('nn', 'linear') for mlab version.
	#
	min_x, max_x = min(xyz['x']), max(xyz['x'])
	min_y, max_y = min(xyz['y']), max(xyz['y'])
	#
	dx = min([x-xyz['x'][j] for j,x in enumerate(xyz['x'][1:]) if x-xyz['x'][j]!=0.])
	dy = min([y-xyz['y'][j] for j,y in enumerate(xyz['y'][1:]) if x-xyz['y'][j]!=0.])
	#
	padding = .05
	n_x = (n_x or int((max_x-min_x)/dx))
	n_y = (n_y or int((max_y-min_y)/dy))
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
#
# testing and development scripts:
#
def test_earthquake_transforms(fignums=(0,1)):
	x1=test_earthquake_transform(fignum=fignums[0], transform_type='equal_area')
	x2=test_earthquake_transform(fignum=fignums[1], transform_type='rotation')
	
	
def test_earthquake_transform(fignum=0,transform_type='equal_area', lats=[33.8, 37.8], lons=[-122.4, -118.4], fit_factor=1.5, etas_mc=4.0, catalog=None):
	# let's do a test script in the Parkfield area...
	# this script will find the parkfield earthquake in the selected catalog, fit local data (via PCA) and apply a transform.
	# we then plot the earthquakes (.), parkfield(*), a curcle, and the elliptical transform (R -> R') of that circle.
	#
	# parkfield={'dt':dtm.datetime(2004,9,28,17,15,24, tzinfo=pytz.timezone('UTC')), 'lat':35.815, 'lon':-120.374, 'mag':5.96}
	if catalog==None:
		catalog = make_ETAS_catalog(incat=None, lats=lats, lons=lons, mc=2.5, date_range=['1990-1-1', None], D_fract=1.5, d_lambda=1.76, d_tau = 2.28, fit_factor=fit_factor, do_recarray=True)
	C = catalog
	#
	# find parkfield (we might also look at san-sim and coalinga):
	for j,rw in enumerate(C):
		if rw['event_date']>dtm.datetime(2004,9,28,17) and rw['event_date']<dtm.datetime(2004,9,28,18,15,24) and rw['mag']>5.9 and rw['mag']<6.1:
			# this should be enough...
			j_pf = j
			break
		#
	#
	E_pf = Earthquake(C[j_pf], transform_type=transform_type)
	#
	# now, make a circle around the earthuqake and do some elliptical transforms...
	circ = Circle(R=E_pf.L_r*2.0/deg2km, x=E_pf.lon, y=E_pf.lat)
	xy_prime = []
	for rw in circ.poly():
		et = E_pf.elliptical_transform_prime(*rw)		# aka, use x' = x dot T, not T dot x
		xy_prime += [[E_pf.lon + et['dx_prime']/(math.cos(E_pf.lat*deg2rad)*deg2km), E_pf.lat+et['dy_prime']/deg2km]]
	
	plt.figure(fignum)
	plt.clf()
	plt.plot(C['lon'], C['lat'], color='b', ls='', marker='.')
	plt.plot(E_pf.lon, E_pf.lat, color='r', ls='', marker='*', ms=18)
	#
	circ.plot_poly(fignum=fignum, do_clf=False, marker='', color='m')
	plt.plot(*zip(*xy_prime), ls='--', color='r', marker='.')
	d_lon1, d_lat1 = E_pf.L_r*E_pf.T[0][0]/(deg2km*math.cos(E_pf.lat*deg2rad)), E_pf.L_r*E_pf.T[0][1]/deg2km
	d_lon2, d_lat2 = E_pf.L_r*E_pf.T[1][0]/(deg2km*math.cos(E_pf.lat*deg2rad)), E_pf.L_r*E_pf.T[1][1]/deg2km
	
	plt.plot([E_pf.lon, E_pf.lon+d_lon1], [E_pf.lat, E_pf.lat+d_lat1], '^-', color='c')
	plt.plot([E_pf.lon, E_pf.lon+d_lon2], [E_pf.lat, E_pf.lat+d_lat2], 'o-', color='m')
	#
	for j,rw in enumerate(circ.poly()):
		# this  makes a cool figure. the transformation is not exactly as simple as it looks; we're not really just rotating and stretching our points; ww're also flipping them.
		# (mirroring across the minor axis) but this won't matter for R -> R' transformations.
		plt.plot(*zip(rw,xy_prime[j]), color='r', ls='--', marker='')
		#plt.plot([rw[0], xy_prime[j][0]], [rw[1], xy_prime[j][1]], color='r', ls='--', marker='')
	
	#
	# now, let's try a (different?) 1/r field:
	# note: we forgot to push some changes, so we'll have conflicts here. maybe we keep both examples...
	plt.figure(2)
	plt.clf()
	d_lat=.05
	d_lon=.05
	field_vals = []
	#
	field_vals = [numpy.zeros(len(numpy.arange(lons[0], lons[1], d_lon))) for rw in numpy.arange(lats[0], lats[1], d_lat)]
	#
	# we could do this in one big itertools.product():
	#for eq, j_lat, k_lon in itertools.product(C, enumerate(numpy.arange(lats[0], lats[1], d_lat)), enumerate(numpy.arange(lons[0], lons[1], d_lon)))
	# however, since we're using the Earthquake class, aka eq=Earthquake(eq), this will probably create a bunch of overhead, so we'll just nest this loop.
	# nominally this might be a case where avoiding the instantiation of a class facilitates a more efficient method (itertools.product() rather than nested loop).
	#
	for eq in C:
		#eq = E_pf
		if eq['mag']<etas_mc:continue
		#if eq!=C[j_pf]: continue
		eq=Earthquake(eq)
		for j_lat, k_lon in itertools.product(enumerate(numpy.arange(lats[0], lats[1], d_lat)), enumerate(numpy.arange(lons[0], lons[1], d_lon))):
		#for j,lat in enumerate(numpy.arange(lats[0], lats[1], d_lat)):
			#field_vals += [[]]	
			#for k, lon in enumerate(numpy.arange(lons[0], lons[1], d_lon)):
			j, lat = j_lat
			k, lon = k_lon
			if True:
				# (maintaining -- for the time being, our indentation structure...)
								
				#r = numpy.linalg.norm(numpy.dot(eq.T, [lon-eq.lon, lat-eq.lat]))
				
				r = eq.elliptical_transform(lon, lat)['R_prime']
				
			
				#field_vals[-1]+= [1./r]
				field_vals[j][k] += 1./(r + 10.**(.5*eq.mag-2.2))**1.5

	plt.contourf(numpy.arange(lons[0], lons[1], d_lon), numpy.arange(lats[0], lats[1], d_lat), numpy.log10(field_vals), 25)
	plt.plot(*zip(*circ.poly()), color='k', ls='-')
	#
	#
	# now, let's do a quick diagnostic plot of earthquak fit angle...
	thetas = [math.acos(V[0][0]/numpy.linalg.norm(V[0])) for V in C['e_vecs']]
	plt.figure(3)
	plt.clf()
	z=plt.hist(thetas, min(150,len(thetas)/25))
	#
	return E_pf, C[j_pf]


if __name__=='__main__':
	# do main stuff...
	pass
else:
	plt.ion()
	


	
		
