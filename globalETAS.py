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
#
# TODO notes:
# 1: revisit Bindex model; i think it will be faster, simpler, and lower memory footprint than rtree, and i think the main obstacles to making it work properly
#    have actually been resolved.
# 2: revise the boundary selection framework for rtree. basically, there is a problem constraining the boundaries when lat/lon parameters are provided to ETAS.
# 3: implement (or confirm implementation of) the brute-force version of ETAS, where we loop-loop over all events,locations. for small ETAS, we probably
#    want to do this anyway and it will be faster because we forego calculating indices.
# 4: revise the spatial limits: instead of just using n*L_r for spatial extetns, actually calculate the distance for z -> x*z0, aka the distance for 90% reduction.
#    this will be similar to the L_r approach, bit will account for proper scaling and should improve the output appearance.
# 5: modify the MPP handler(s) to inherit from Pool() and actually work like a Pool(), rather than explicitly handling Process() objects. the problem is that when
#    we divide up the job, we inevitably give one processor a light job and one processor a heavy job, so we spend a lot of time waiting for one process to finish.
#    yes, there is a bit of overhead in running with n_cpu > n_processors, but it is minor in comparison. also, by dividing the job into lots of small little
#    jobs, we can probably also reduce the real-time memory footprint, so we can run on smaller systems.
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
import scipy.optimize as spo
import itertools
import sys
#import scipy.optimize as spo
import sys
import os
#from PIL import Image as ipp
import multiprocessing as mpp
#
import matplotlib as mpl
import matplotlib.pyplot as plt
import functools
#
# let's see if we can compile some of thes functions.
import numba
#
#import shapely.geometry as sgp
#
from mpl_toolkits.basemap import Basemap as Basemap
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from geographiclib.geodesic import Geodesic as ggp
#
#
#import ANSStools as atp
from yodiipy import ANSStools as atp
import bindex
#import contours2kml.contours2kml as contours2kml
from yodiipy import contours2kml as c2kml
#
import rtree
from rtree import index
import geopy
from geopy.distance import great_circle
#
# python3 vs python2 issues:
# a bit of version 3-2 compatibility:
if sys.version_info.major>=3:
	xrange=range
from eq_params import *
#
# TODO: add a p_map (or p, p0 distinction) variable to distinguish the p value for calculating ETAS parameters (in catalog) and p to calculate
# ETAS itself. this will facilitate time independent ETAS, etc., which we need to proprly evaluate geo-spatial ROC.
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
dist_calc_types_dict.update({key:'euclidian' for key in ['euc', 'euclidian', 'euclid']})
#dist_calc_types_dict.update({key:'cartesian2' for key in ['cartesian2']})
dist_calc_types_dict.update({key:'geo' for key in ['geo', 'geodesic', 'geodesiclib', 'glib', 'g_lib']})
# and this one to return the x,y components for cartesian measurements:
dist_calc_types_dict.update({key:'dx_dy' for key in ['dxdy', 'dx_dy']})
#
# Developer notes and status:
# as of 9 oct 2015: Global_ETAS_model and the various subclasses (which use differend indexing) run (apprently) successfully, but the figures they turn out look a bit
# odd, so i'm working on diagnostic scripts (maybe i just need to run a bigger model... or maybe just one more directly comparable to some standard runs like nor/so CA.
# see also etas_testing.py. it looks like the Earthquake() object products the correct local_intensity(), as compared to the former BASScast/ETASmf model.
# so maybe the rates are not aggregating correctly, or maybe the indexing is screwing up something.
# the Bindex() indexing approach appears to be having some float-integer indexing problems (basically, error in floating point accuracy causing some bins to go the wrong
# way. probably, we need a more explicit approach that maps out the bins by counting (aka, instead of 35.5 --> 35.0/.1 --> j=350, we just need to line up the values we've got
# and count off the bins. note, skipping Bindex and just using the rtree approach... or use Bindex, but pre-calc the domain. maybe also look into using
# numpy.digitize(x,bins) or there are some intesting examples of using hist() to make bins.
#
# Classes:

class Global_ETAS_model(object):
	# guts of ETAS. this will contain a catalog, lattice sites, etc. graphical tools, plotting, etc. will likely
	# be kept elsewhere.
	# questions: use equal distance or equal angel bin spacing? we can nominally do a spherical->cartesian transformation like: x = lon*cos(lat)*deg2km
	#  but it behaves sort of funny when we try to calculate disanctes between points... which we don't necessaryily do; we map (x,y) -> (lon, lat) and then
	#  do proper geodetic transformations, but 1) it becomes easy to make mistakes, and 2) the only gain is that we have fewer points in the higher latitudes.
	#  as long as 1) we have sufficient spatial resolution in the lower latitudes and 2) the problem remains computational tractable (equal spacing can reduce
	#  the grid size by maybe halfish??), there's no real harm in just using (a much simpler) lon,lat lattice with equal angular spacing.
	#
	#def __init__(self, catalog=None, lats=[32., 38.], lons=[-117., -114.], mc=2.5, d_x=10., d_y=10., bin_x0=0., bin_y0=0., etas_range_factor=5.0, t_0=dtm.datetime(1990,1,1, tzinfo=tz_utc), t_now=dtm.datetime.now(tzutc), transform_type='equal_area', transform_ratio_max=5., calc_etas=True):
	def __init__(self, catalog=None, lats=None, lons=None, mc=2.5, mc_etas=None, d_lon=.1, d_lat=.1, bin_lon0=0., bin_lat0=0., etas_range_factor=25.0, 
                          etas_range_padding=1.5, etas_fit_factor=1.0, t_0=dtm.datetime(1990,1,1, tzinfo=tz_utc), t_now=dtm.datetime.now(tzutc), 
                          t_int_limit=None, transform_type='equal_area', transform_ratio_max=2.5, cat_len=5.*365., calc_etas=True, n_contours=15, 
                          etas_cat_range=None, etas_xyz_range=None, p_cat=1.1, q_cat=1.5, ab_ratio_expon=1.0, p_etas=None, D_fract=1.5, **kwargs):
		'''
		#
		#  basically: if we are given a catalog, use it. try to extract mc, etc. data from catalog if it's not
		# externally provided. otherwise, get a catalog from ANSS or OH based oon the parameters provided.
		# note: we want to do this properly in an x-y projection, not lat-lon. our earthquakes will come back to us with lat/lon coordinates
		# so we'll need to translate. we'll also want lat/lon as well as x/y positions for the lattice sites, so we can get even spacing and accurate distances
		#
		# etas_cat_range: range of catalog to do ETAS. mainly, this is to facilitate MPP processing; each process will have a full catalog of earthquakes and
		# (nominally) a full array of etas lattice sites (though we might be able to improve upon this requirement later -- probably the best appraoch is to
		# accept this limitation and write a make_etas() that uses an mpp.Array() object (simplified script).
		#
		# p_,q_cat: p,q values passed to the catalog getter. these will affect not only the p,q values in the catalog, but the
		# other calculated ETAS parameters.
		#
		# p_etas: use this to control p for etas calculations, separate for p calculated in catalog (aka, p to solve initial rates). use p_etas to make time-independent maps.
		'''
		print("begin globalETAS.__init()__")
		# dx=1., dy=1., x0=0., y0=0., leaf_type=float)
		# load all the input parameters into class variables (note i believe locals() is all currently defined variables in the function scope,
		# so this behaves differently if we execute now or down-stream. if we really only want to record the input values, then
		# execute now. if we want to record corrected inputs, then execute after corrections; if we execute at the end, every declared
		# variable becomes a class member.
		#self.__dict__.update(locals())
		#
		# we might just want the last N days, as a consistent standard. note we might, later on, make this a bit more sophisticated
		# by processing the full t_0 -> t_now catalog, but only doing ETAS for the most recent cat_len days. BUT, to do this, we have
		# to enforce in all the do_ETAS() functions
		t_now = (t_now or dtm.datetime.now(pytz.timezone('UTC')))
		if cat_len is not None:
			t_0=t_now - dtm.timedelta(days=cat_len)
			print("Overriding t_0 (etas catalog start date/time) for ETAS calculations. using catalog start, t_0 = t_now - catlen (%f) = %s" % (cat_len, t_0))
		#
		if lats is None and catalog is None: lats = [-89.9, 89.9]
		if lons is None and catalog is None: lons = [-180., 180.]
		#
		# for now, assume the catalog is string-indexed -- aka, recarray, PANDAS,etc.
		if lats is None and not (catalog is None or len(catalog) is 0): lats = [min(catalog['lat']), max(catalog['lat'])]
		if lons is None and not (catalog is None or len(catalog) is 0): lons = [min(catalog['lon']), max(catalog['lon'])]
		if mc   is None and not (catalog is None or len(catalog) is 0): mc = min(catalog['mag'])
		#
		# and handle some specific cases...
		if isinstance(t_now, float):
			self.t_forecast = t_now
		elif isinstance(t_now, numpy.datetime64):
			self.t_forecast = mpd.date2num(t_now.tolist())
		else:
			self.t_forecast = mpd.date2num(t_now)
		#
		if isinstance(t_int_limit, float):
			self.t_int_lim = t_int_limit
		elif isinstance(t_int_limit, numpy.datetime64):
			self.t_int_lim = mpd.date2num(t_int_limit.tolist())
		elif isinstance(t_int_limit, dtm.datetime):
			self.t_int_lim = mpd.date2num(t_int_limit)
		else:
			self.t_int_lim = None
		#
		mc_etas = (mc_etas or mc)	# mc_eats: minimum mag. for etas calculations -- aka, mc for the catalog, but only do etas for m>mc_eatas.
		#
		self.mc_etas = mc_etas
		# inputs massaged; now update class dictionary.
		self.__dict__.update(locals())
		#
		self.d_lat = d_lat
		self.d_lon = d_lon
		#		
		self.latses = numpy.arange(lats[0], lats[1], d_lat)		# note: if we want lats[1], lons[1] inclusive, we need to add +d_lat, +d_lon
		self.lonses = numpy.arange(lons[0], lons[1], d_lon)		# to the range().
		self.n_lat = len(self.latses)
		self.n_lon = len(self.lonses)
		#
		# calculate xyz_range (for mpp applications where we split up the geo-spatial array to processes):
		if etas_xyz_range is None: etas_xyz_range = [0,self.n_lat*self.n_lon]
		etas_xyz_range[0] = (etas_xyz_range[0] or 0)
		etas_xyz_range[1] = (etas_xyz_range[1] or self.n_lat*self.n_lon)
		#
		self.ETAS_array = numpy.array([])
		# [etas_xyz_range[0]:etas_xyz_range[1]]
		self.ETAS_array = numpy.array([[lon, lat, 0.] for j, (lat,lon) in enumerate(itertools.product(self.latses, self.lonses)) if (j>= etas_xyz_range[0] and j<etas_xyz_range[1])])
		self.ETAS_array = numpy.core.records.fromarrays(zip(*self.ETAS_array), dtype = [('x', '>f8'), ('y', '>f8'), ('z', '>f8')])
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
		# this should probably be moved into the specific etas type parts...		
		#self.lattice_sites = bindex.Bindex2D(dx=d_lon, dy=d_lat, x0=bin_lon0, y0=bin_lat0)	# list of lattice site objects, and let's index it by... probably (i_x/lon, j_y/lat)
							# we might alternativel store this as a simple list [] or a base index/list (that
							# is re-sorting tolerant)) {i:[row]}, and then write indices:
							# index_lat_lon = {(lat, lon):lattice_sites_index}, {(x,y):lattice_sites_index}, etc.
							#
		#
		#earthquake_catalog = {}	# and again, let's use an indexed structure, in the event that we are
								# using indices and (inadvertently? sort it?). this may be unnecessary.
								# that said, let's take the next step and return dict. type earthquake entries.
		#
		if catalog is None:
			print("fetch and process catalog for dates: {}".format([t_0, t_now]))
			#catalog = make_ETAS_catalog(incat=None, lats=lats, lons=lons, mc=mc, date_range=[t_0, t_now], fit_factor=etas_fit_factor)	# and note there are other variables to consider...
			catalog = make_ETAS_catalog_mpp(incat=None, lats=lats, lons=lons, mc=mc, date_range=[t_0, t_now], fit_factor=etas_fit_factor, p=p_cat, q=q_cat, D_fract=D_fract)	# and note there are other variables to consider...
			print("catalog fetched and processed.")
		self.catalog = catalog
		#
		# set etas_cat_range as necessary:
		if etas_cat_range is None: etas_cat_range = [0,len(catalog)]
		etas_cat_range[0] = (etas_cat_range[0] or 0)
		etas_cat_range[1] = (etas_cat_range[1] or len(catalog))
		self.etas_cat_range = etas_cat_range
		print("ETAS over etas_cat_range/xyz_range: ", (self.etas_cat_range, self.etas_xyz_range))
		#
		# ... and here is really where we do derived classes for different ETAS forcast objects. so for ETAS_bindex(ETAS), we do:
		if not hasattr(self, 'make_etas'):
			# instantiating directly from this parent class script, or for some other reason we don't get a make_etas() method defined.
			self.make_etas = self.make_etas_bindex
		# ... and for ETAS_brute(ETAS), we'd do:
		# self.make_etas = self.make_etas_all
		#
		if calc_etas: 
			print("make_etas():")
			self.make_etas()
			print("ETAS complete.")
	#
	@property
	def lattice_sites(self):
		X = self.ETAS_array['z']
		X.shape=(len(self.latses), len(self.lonses))
		return X
	#
	def draw_map(self, fignum=0, fig_size=(6.,6.), map_resolution='i', map_projection='cyl', d_lon_range=None, d_lat_range=None, lats_map=None, lons_map=None, ax=None, do_states=True, do_rivers=True, lake_color='blue', lat_label_indices=[1,1,0,0], lon_label_indices=[0,0,1,1]):
		'''
		# TODO: we end up matching up a bunch of procedural calls, which is a big pain. we should write an ETAS_Map() class
		# which includes the contour,etc. figures... but we can keep the variables, like lon_label_indices, etc.
		# in one place...
		#
		# plot contours over a map.
		'''
		lons_map = (lons_map or self.lons)
		lats_map = (lats_map or self.lats)
		#
		# first, get contours:
		#etas_contours = self.calc_etas_contours(n_contours=n_contours, fignum=fignum, contour_fig_file=contour_fig_file, contour_kml_file=contour_kml_file, kml_contours_bottom=kml_contours_bottom, kml_contours_top=kml_contours_top, alpha_kml=alpha_kml, refresh_etas=refresh_etas)
		#
		# now, clear away the figure and set up the basemap...
		#
		d_lon_range = (d_lon_range or 1.)
		d_lat_range = (d_lat_range or 1.)
		#
		if ax==None:
			plt.figure(fignum, fig_size)
			plt.clf()
			ax=plt.gca()
		#
		#lons, lats = self.lons, self.lats
		#cntr = [numpy.mean(lons), numpy.mean(lats)]
		cntr = [numpy.mean(lons_map), numpy.mean(lats_map)]
		#cm = Basemap(llcrnrlon=self.lons[0], llcrnrlat=self.lats[0], urcrnrlon=self.lons[1], urcrnrlat=self.lats[1], resolution=map_resolution, projection=map_projection, lon_0=cntr[0], lat_0=cntr[1])
		cm = Basemap(llcrnrlon=lons_map[0], llcrnrlat=lats_map[0], urcrnrlon=lons_map[1], urcrnrlat=lats_map[1], resolution=map_resolution, projection=map_projection, lon_0=cntr[0], lat_0=cntr[1], ax=ax)
		#
		#cm.drawlsmask(land_color='0.8', ocean_color='b', resolution=map_resolution)
		cm.drawcoastlines(color='gray', zorder=1)
		cm.drawcountries(color='black', zorder=1)
		if do_states: cm.drawstates(color='black', zorder=1)
		if do_rivers: cm.drawrivers(color='blue', zorder=1)
		cm.fillcontinents(color='beige', lake_color=lake_color, zorder=0)
		# drawlsmask(land_color='0.8', ocean_color='w', lsmask=None, lsmask_lons=None, lsmask_lats=None, lakes=True, resolution='l', grid=5, **kwargs)
		#cm.drawlsmask(land_color='0.8', ocean_color='c', lsmask=None, lsmask_lons=None, lsmask_lats=None, lakes=True, resolution=self.mapres, grid=5)
		# lat_label_indices
		#cm.drawmeridians(numpy.arange(int(lons_map[0]/d_lon_range)*d_lon_range, lons_map[1], d_lon_range), color='k', labels=[0,0,1,1])
		#cm.drawparallels(numpy.arange(int(lats_map[0]/d_lat_range)*d_lat_range, lats_map[1], d_lat_range), color='k', labels=[1, 1, 0, 0])
		cm.drawmeridians(numpy.arange(int(lons_map[0]/d_lon_range)*d_lon_range, lons_map[1], d_lon_range), color='k', labels=lon_label_indices)
		cm.drawparallels(numpy.arange(int(lats_map[0]/d_lat_range)*d_lat_range, lats_map[1], d_lat_range), color='k', labels=lat_label_indices)

		#
		return cm
	def make_etas_contour_map(self, n_contours=None, fignum=0, fig_size=(6.,6.), contour_fig_file=None, contour_kml_file=None, kml_contours_bottom=0., kml_contours_top=1.0, alpha=.5, alpha_kml=.5, refresh_etas=False, map_resolution='i', map_projection='cyl', map_cmap='jet', lat_interval=None, lon_interval=None, lats_map=None, lons_map=None, ax=None, do_colorbar=True, do_states=True, do_rivers=True, lake_color='blue' ):
		#
		n_contours = (n_contours or self.n_contours)
		if ax is None:
			fg=plt.figure(fignum)
			ax=fg.gca()
		#
		# mm.draw_map(d_lat_range=10., d_lon_range=20., fignum=0)
		#cm = self.draw_map(fignum=fignum, fig_size=fig_size, map_resolution=map_resolution, map_projection=map_projection)
		cm = self.draw_map(fignum=fignum, fig_size=fig_size, map_resolution=map_resolution, map_projection=map_projection, d_lon_range=lon_interval, d_lat_range=lat_interval, lons_map=lons_map, lats_map=lats_map, ax=ax, do_states=do_states, do_rivers=do_rivers, lake_color=lake_color)
		#
		fg=plt.gcf()
		#
		X,Y = cm(numpy.array(self.lonses), numpy.array(self.latses))
		#print("xylen: ", len(X), len(Y))
		#
		Z = numpy.log10(self.lattice_sites)
		#etas_contours = ax.contourf(X,Y, numpy.log10(self.lattice_sites), n_contours, zorder=8, alpha=alpha, cmap=map_cmap)
		etas_contours = ax.contourf(X,Y, Z, n_contours, zorder=8, alpha=alpha, cmap=map_cmap)
		# ax.colorbar() ??
		if do_colorbar:
			#plt.colorbar(ax)
			plt.colorbar(etas_contours, cax=None, ax=ax, cmap=map_cmap)
			#mpl.colorbar.ColorbarBase(ax=ax, cmap=map_cmap, values=sorted(Z.ravel()), orientation="vertical")
		#
		self.cm=cm
		self.etas_contours = etas_contours
		#
		return cm
		#
	#
	def make_etas_boxy_map(self, n_contours=None, fignum=0, fig_size=(6.,6.), contour_fig_file=None, contour_kml_file=None, kml_contours_bottom=0., kml_contours_top=1.0, alpha=.6, alpha_kml=.5, refresh_etas=False, map_resolution='i', map_projection='cyl', map_cmap='jet'):
		#
		cm = self.draw_map(fignum=fignum, fig_size=fig_size, map_resolution=map_resolution, map_projection=map_projection)
		c_map = plt.get_cmap(map_cmap)
		zs = numpy.log10(self.ETAS_array['z'])
		cNorm = mpl.colors.Normalize(vmin=min(zs), vmax=max(zs))
		scalarMap = mpl.cm.ScalarMappable(norm=cNorm, cmap=c_map)
		#
		for rw in self.ETAS_array:
			plt.fill_between(x=[rw['x']-self.d_lon/2., rw['x']+self.d_lon/2.], y1=[rw['y']-.5*self.d_lat, rw['y']-.5*self.d_lat], y2=[rw['y']+.5*self.d_lat, rw['y']+.5*self.d_lat], color=scalarMap.to_rgba(numpy.log10(rw['z'])))
		#plt.colorbar()		# not sure how to make this work for non-contour plot...
		#
		return cm
		#
	def plot_mainshock_and_aftershocks(self, m0=6.0, n_contours=25, mainshock=None, fignum=0, ax=None):
		#
		map_etas = self.make_etas_contour_map(n_contours=n_contours, fignum=fignum, ax=ax)
		if mainshock is None:
			mainshock = self.catalog[0]
			for rw in self.catalog:
				if rw['mag']>mainshock['mag']: mainshock=rw
		ms=mainshock
		#
		ax = (ax or plt.gca())
		#
		for eq in self.catalog:
			if eq['mag']<m0 or eq['event_date']<ms['event_date']: continue
			if eq==ms:
				x,y = map_etas(eq['lon'], eq['lat'])
				ax.plot([x], [y], 'k*', zorder=7, ms=20, alpha=.8)
				ax.plot([x], [y], 'r*', zorder=8, ms=18, label='mainshock', alpha=.8)
			if eq['event_date']>eq['event_date']:
				x,y = map_etas(eq['lon'], eq['lat'])
				ax.plot([x], [y], 'o', zorder=7, ms=20, alpha=.8)	
		#
		#return plt.gca()
		return ax
	#
	def calc_etas_contours(self, n_contours=None, fignum=0, contour_fig_file=None, contour_kml_file=None, kml_contours_bottom=0., kml_contours_top=1.0, alpha_kml=.5, refresh_etas=False):
		# wrapper for one-stop-shopping ETAS calculations.
		# (and these calc_contours schemes need to be cleaned up a bit. there is a bit of redundancy and disorganization)
		#
		n_contours = (n_contours or self.n_contours)
		#
		if refresh_etas or ('ETAS_array' not in self.__dict__.keys()):
			self.make_etas()
		#
		plt.figure(fignum)
		plt.clf()
		#
		self.etas_contours = plt.contourf(self.lonses, self.latses, numpy.log10(self.lattice_sites), n_contours)
		plt.colorbar()
		#
		if contour_fig_file!=None:
			p_name, f_name = os.path.split(contour_fig_file)
			if not os.path.isdir(p_name): os.makedirs(p_name)
			plt.savefig(contour_fig_file)
			#
		if contour_kml_file!=None:
			# make kml and write to file like:
			# kml_str = kml_from_contours(cset=contours, colorbarname='napa_colorbar.png', open_file=True, close_file=True, contour_labels=None, top=top, bottom=bottom, fname_out=fname_out, alpha_kml=alpha_kml)
			# ... and this is maybe not the best coding framework, but this function will automatically output to a file if we
			# give it fname_out!=None. BUT, since it's sort of a crappy way to write code, and we'll likely 'fix' it later, 
			# let's just return a string and write our own file. we'll even dump the KML to the class dict. maybe... or maybe we should
			# use the built in file-writer, since i think we need/want to also include the colorbar...
			#
			self.contours_kml_str = c2kml.kml_from_contours(cset=self.etas_contours, colorbarname=None, open_file=True, close_file=True, contour_labels=None, top=kml_contours_top, bottom=kml_contours_bottom, alpha_kml=alpha_kml, fname_out=contour_kml_file)
			p_name, f_name = os.path.split(contour_kml_file)
			if not os.path.isdir(p_name): os.makedirs(p_name)
			with open(contour_kml_file, 'w') as f_kml:
				f_kml.write(self.contours_kml_str)
			#pass
			#
		return self.etas_contours
	#
	def export_kml(self, fout='etas_contours.kml', kml_contours_bottom=0., kml_contours_top=1.0, alpha_kml=.5):
		#
		self.contours_kml_str = c2kml.kml_from_contours(cset=self.etas_contours, colorbarname=None, open_file=True, close_file=True, contour_labels=None, top=kml_contours_top, bottom=kml_contours_bottom, alpha_kml=alpha_kml, fname_out=None)
		p_name, f_name = os.path.split(fout)
		if not os.path.isdir(p_name): os.makedirs(p_name)
		with open(fout,'w') as f:
			f.write(self.contours_kml_str)
		#
	#
	def kml_from_contours(self, contours=None, contour_kml_file=None, kml_contours_bottom=0., kml_contours_top=1.0, alpha_kml=.5, refresh_etas=False):
		if contours is None: contours = self.etas_contours
		if (contours is None or refresh_etas): contours = plt.contourf(self.lonses, self.latses, numpy.log10(self.lattice_sites), n_contours)
		#
		self.contours_kml_str = c2kml.kml_from_contours(cset=self.etas_contours, colorbarname=None, open_file=True, close_file=True, contour_labels=None, top=kml_contours_top, bottom=kml_contours_bottom, alpha_kml=alpha_kml, fname_out=contour_kml_file)
		p_name, f_name = os.path.split(contour_kml_file)
		if not os.path.isdir(p_name): os.makedirs(p_name)
		with open(contour_kml_file, 'w') as f_kml:
			f_kml.write(self.contours_kml_str)
		return contours_kml_str
	#
	def export_xyz(self, fout='etas_xyz.xyz'):
		#
		# export ETAS_array as xyx format.
		#
		p_name, f_name = os.path.split(fout)
		#
		if not os.path.isdir(p_name):
			os.makedirs(os.path.split(p_name))
		#
		with open(fout, 'w') as f:
			f.write('#globalETAS xyz (lon, lat, z_etas) export\n')
			f.write('#eventually, add metadata\n')
			f.write('#!x\t\y\tz\n')
			[f.write('\t'.join([str(x) for x in rw]) + '\n') for rw in self.ETAS_array.tolist()]
		#
		
		
	#
	def make_etas_rtree(self):
		# use the same basic framework as etas_all (aka, instantiate a lattice), but map an rtree index to the lattice, then use a delta_lat, delta_lon
		# approach (like etas_bindex) bit loop only over the rtree.intersection() of the delta_lat/lon window.
		#
		print("begin make_etas_rtree()")
		#
		latses=self.latses
		lonses=self.lonses
		n_lat=self.n_lat
		n_lon=self.n_lon
		#
		# note: it's important to do the (latses, lonses) in the same sequence so we have consistency in the indices.
		#     this would make a reasonable case to use a conventional for-loop instead of multiple list comprehensions.
		#print("initialize ETAS array and indices.")
		#if not hasattr(self, 'ETAS_array') or len(self.ETAS_array)==0:
		#	self.ETAS_array = [[lon, lat, 0.] for lat,lon in itertools.product(latses, lonses)]
		#	#self.ETAS_array = numpy.array([[lon, lat, 0.] for lat,lon in itertools.product(self.latses, self.lonses)])
		#	self.ETAS_array = numpy.core.records.fromarrays(zip(*self.ETAS_array), dtype = [('x', '>f8'), ('y', '>f8'), ('z', '>f8')])
		#
		# make a dict of lattice sites and their x,y coordinates.
		#lattice_dict = {i:{'lat':lat, 'lon':lon, 'j_lon':int(i%n_lon), 'j_lat':int(i/n_lon)} for i, (lat,lon) in enumerate(itertools.product(latses, lonses))}
		# Wilson: added the area calculation here
		lattice_dict = {i:{'lat':lat, 'lon':lon, 'j_lon':int(i%n_lon), 'j_lat':int(i/n_lon),
						  'area': (ggp.WGS84.Inverse(lat-self.d_lat/2., lon, lat+self.d_lat/2., lon)['s12']
								  * ggp.WGS84.Inverse(lat, lon-self.d_lon/2., lat, lon+self.d_lon/2.)['s12'])/1e6 #ETAS densities in km^-2
						} for i, (lon,lat,z) in enumerate(self.ETAS_array)}
		self.lattice_dict=lattice_dict
		print("len(local_lattice_dict): ", len(self.lattice_dict))
		#
		# make an rtree index:
		lattice_index = index.Index()
		# build index. we should build directly from ETAS_array, particualrly since we'll be (maybe) building an mpp version that splits up ETAS_array between proceses.
		#[lattice_index.insert(j, (lon, lat, lon, lat)) for j, (lat,lon) in enumerate(itertools.product(latses,lonses))]	# like [[lat0,lon0],[lat0,lon1], [lat0,lon2]...]
		[lattice_index.insert(j, (lon, lat, lon, lat)) for j, (lon, lat, z) in enumerate(self.ETAS_array)]	# like [[lat0,lon0],[lat0,lon1], [lat0,lon2]...]
		#
		print("Indices initiated. begin ETAS :: ", self.etas_cat_range)
		#
		#for quake in self.catalog:
		for quake in self.catalog[self.etas_cat_range[0]:self.etas_cat_range[1]]:
			if quake['mag']<self.mc_etas: continue
			#
			if quake['event_date_float']>self.t_forecast: continue
			#
			eq = Earthquake(quake, transform_type=self.transform_type, transform_ratio_max=self.transform_ratio_max,			ab_ratio_expon=self.ab_ratio_expon)
			#
			# ab_ratio_expon is, indeed, making it successfully to Earthquake(). 
			#print('**debug: eq abexpon: {}'.format(eq.ab_ratio_expon))
			#
			# get lat/lon range:
			# TODO: introduce a new lat/lon range model; use the spatial-omori distribution; calc. to r = r(x=.9x0)
			delta_lat = self.etas_range_padding + eq.L_r*self.etas_range_factor/deg2km
			if abs(eq.lat)==90.:
				delta_lon=180.
			else:
				delta_lon = self.etas_range_padding + delta_lat/math.cos(eq.lat*deg2rad)
			#
			# and let's also assume we want to limit our ETAS map to the input lat/lon:
			# this formulation can get confused if lons, lats, and a catalog are provided separately. look for a smarter way...
			#lon_min, lon_max = max(eq.lon - delta_lon, self.lons[0]), min(eq.lon + delta_lon, self.lons[1])
			#lat_min, lat_max = max(eq.lat - delta_lat, self.lats[0]), min(eq.lat + delta_lat, self.lats[1])
			# for now, just let the space be big.
			lon_min, lon_max = eq.lon - delta_lon, eq.lon + delta_lon
			lat_min, lat_max = eq.lat - delta_lat, eq.lat + delta_lat
			
			#
			#print('rtree indexing: ', quake, lon_min, lat_min, lon_max, lat_max, self.lons, self.lats, delta_lon, delta_lat)
			#
			site_indices = list(lattice_index.intersection((lon_min, lat_min, lon_max, lat_max)))
			# ... and if we wrap around the other side of the world...
			# there's probably a smarter way to do this...
			if lon_min<-180.:
				#new_lon_min = lon_min%(180.)
				#new_lon_min = 180.+(lon_min+180.)
				new_lon_min = 360. + lon_min
				site_indices += list(lattice_index.intersection((new_lon_min, lat_min, 180., lat_max)))
			if lon_max>180.:
				#new_lon_max = lon_max%(-180.)
				#new_lon_max = -180. + lon_max-180.
				new_lon_max = -360. + lon_max
				site_indices += list(lattice_index.intersection((-180., lat_min, new_lon_max, lat_max)))
			#
			#print("LLRange: ", lon_min, lat_min, lon_max, lat_max, len(list(site_indices)))
			#
			# Wilson: set up array for this quake's ETAS for normalization purposes
			this_quake_ETAS = numpy.zeros(len(self.ETAS_array['z']), dtype=float)
			#
			for site_index in site_indices:
				X = lattice_dict[site_index]
				#
				# do we need to do this, or does the Earthquake self-transform?
				et = eq.elliptical_transform(lon=X['lon'], lat=X['lat'])
				#
				# this looks right:
				#local_intensity = 1.0/((75. + et['R_prime'])**1.5)
				#
				local_intensity = eq.local_intensity(t=self.t_forecast, lon=X['lon'], lat=X['lat'], p=self.p_etas, t_int_lim=self.t_int_lim)
				if numpy.isnan(local_intensity):
					#print("NAN encountered: ", site_index, self.t_forecast, X['lon'], X['lat'], eq.lon, eq.lat, eq.__dict__)
					continue
				#
				# Wilson: Multiply by the area of the cell if we're outputting a number instead of a rate-density
				if(self.t_int_lim):
					local_intensity *= X['area']
				
				this_quake_ETAS[site_index] = local_intensity
				#
				#self.lattice_sites.add_to_bin(bin_x=lon_bin['index'], bin_y=lat_bin['index'], z=local_intensity)
				#self.lattice_sites[j_lon][j_lat] += local_intensity
				#
				# TODO: generalize this; write a function like: self.add_to_ETAS(site_index, local_intensity).
				#       for current, SPP or non-memory shared versions, this will be the same as below -- directly add to the array.
				#       however, if we 1) write a c++ extension version or 2) use a shared memory array (or both), we can override the class
				#       definition of self.add_to_ETAS() to be something like, self.shared_ETAS_array[3*site_index+2]+=local_intensity.
				#       that said, it might be possible to accomplish this in the __init__ by declaring the ETAS_array as a list from
				#       the shared array by reff. (??)
			#
			if(self.t_int_lim):
				# Here I'm doing an integration over time from 't_now' to 't_int_lim' (relative to event time).
				# If we're doing the integration, then local_intensity spat out rate, so we need to normalize it 
				# Normalize this_quake_ETAS to N_om.  Values of dmstar and b are the same as those used to calc N_om in the catalog creation,
				# but those are just hardcoded too.
				delta_t = (self.t_forecast - quake['event_date_float'])*days2secs
				t_final = (self.t_int_lim  - quake['event_date_float'])*days2secs
				tau_normed = quake['t_0']**(1-quake['p'])/(quake['p']-1)  # This normalizes the below integration (0-inf), for use after spatial dist noramalized to N_om
				integ_result = scipy.integrate.quad(omori_time, delta_t, t_final, args=(tau_normed, quake['t_0'], quake['p']))
				time_int = integ_result[0]
				
				dmstar = 1.0
				b = 1.0
				lN_om = b*(quake['mag'] - dmstar - self.mc_etas)
				N_om = 10.0**lN_om
				this_quake_ETAS *= (N_om*time_int/numpy.sum(this_quake_ETAS))
			#			
			self.ETAS_array['z'] += this_quake_ETAS
				#
		#
		print('finished calculateing ETAS (rtree). wrap up in recarray and return.')
		# this conversion to the 2d array should probably be moved to a function or @property in the main class scope.
		
		#self.lattice_sites = numpy.array([rw[2] for rw in self.ETAS_array])
		#
		#self.lattice_sites.shape=(len(latses), len(lonses))
		#
		self.ETAS_array = numpy.core.records.fromarrays(zip(*self.ETAS_array), dtype = [('x', '>f8'), ('y', '>f8'), ('z', '>f8')])
	#		
	def make_etas_all(self):
		# TODO: this does not appear to work correctly, or at least not in MPP mode. might be because i'm trying to hijack the 
		#       rtree MPP model; the best approach might be to more directly hijack the rtree model and just circumvent the indexing
		#       or burn a little bit of memory and impose a default index (of all elements).
		# loop-loop over the whole lattice space...
		# first, make an empty rates-lattice:
		#
		#latses = numpy.arange(self.lats[0], self.lats[1]+self.d_lat, self.d_lat)
		#lonses = numpy.arange(self.lons[0], self.lons[1]+self.d_lon, self.d_lon)
		latses = self.latses
		lonses = self.lonses
		n_lat  = self.n_lat
		n_lon  = self.n_lon
		#
		#for quake in self.catalog:
		print("Begin ETAS all/brute :: ", self.etas_cat_range)
		#
		#for quake in self.catalog:
		for quake in self.catalog[self.etas_cat_range[0]:self.etas_cat_range[1]]:
		#for quake in self.catalog[etas_cat_range[0]:etas_cat_range[1]]:
			if quake['mag']<self.mc_etas: continue
			if quake['event_date_float']>self.t_forecast: continue
			#
			eq = Earthquake(quake, transform_type=self.transform_type, transform_ratio_max=self.transform_ratio_max, 		ab_ratio_expon=self.ab_ratio_expon)
			#for j, lat_lon in enumerate(itertools.product(latses, lonses)):
			for j, (lat,lon) in enumerate(itertools.product(latses, lonses)):
				# TODO: for some reason, this is running over the site index (j=len(catalog) or something).
				if j>=len(self.catalog):continue
				#
				#self.ETAS_array[j][2] += eq.local_intensity(t=self.t_forecast, lon=lat_lon[1], lat=lat_lon[0], p=self.p_etas)
				self.ETAS_array[j][2] += eq.local_intensity(t=self.t_forecast, lon=lon, lat=lat, p=self.p_etas, t_int_lim=self.t_int_lim)
			#
			#for lat_tpl,lon_tpl in itertools.product(enumerate(latses), enumerate(lonses)):
			#	local_intensity = eq.local_intensity(t=self.t_forecast, lon=lon_tpl[1], lat=lat_tpl[1])
			#	#
			#	#self.lattice_sites.add_to_bin(bin_x=lon_bin['index'], bin_y=lat_bin['index'], z=local_intensity)
			#	self.lattice_sites[lon_tpl[0]][lat_tpl[0]] += local_intensity
			#	#
			#
		#
		#self.lattice_sites = numpy.array([rw[2] for rw in self.ETAS_array])
		#self.lattice_sites.shape=(len(latses), len(lonses))
		#
		#self.ETAS_array = []
		#
		self.ETAS_array = numpy.core.records.fromarrays(zip(*self.ETAS_array), dtype = [('x', '>f8'), ('y', '>f8'), ('z', '>f8')])
		
	def make_etas_bindex(self):
		# rewrite this version; minimize the use of functional indexing; aka, get initial position from functional index; then just count.
		#
		# do the ETAS calculation. for each earthquake, figure out what bins it contributes to and calc. etas for each bin.
		# (new in p3?) can't seem to set this attribute... maybe in __init__()? generally, i think we're going to handle this a bit differently anyway...
		self.lattice_sites = None 
		bb=bindex.Bindex2D(dx=self.d_lon, dy=self.d_lat, x0=self.bin_lon0, y0=self.bin_lat0)	# list of lattice site objects, and let's index it by...
		self.lattce_sites=bb																			# probably (i_x/lon, j_y/lat)
		# copy class level variables:
		#
		print( "calc bindex etas...")
		#for quake in self.catalog:
		for quake in self.catalog[etas_cat_range[0]:etas_cat_range[1]]:
			if quake['mag']<self.mc_etas: continue
			#
			eq = Earthquake(quake, transform_type=self.transform_type, transform_ratio_max=self.transform_ratio_max, 			ab_ratio_expon=self.ab_ratio_expon)
			# calculate the bins within range of this earthquake:
			# nominally, the proper thing to do now (or at lest one proper approach) is to compute a proper geodesic north and east to get the lat/lon bounds
				# near the poles, however, this will create a problem, that is a bit difficult to track, where the geodesic reaches over the pole and to lat<90
			# on the other side. if we just use a recta-linear approximation, we can use [max(-90, lat-d_lat), min(90, lat+d_lat)]
			#x,y = lon_lat_2xy(lon=quake['lon'], lat=quake['lat'])
			#
			# range of influence:
			delta_lat = self.etas_range_padding + eq.L_r*self.etas_range_factor/deg2km
			if abs(eq.lat) == 90.:
				delta_lon  = 180.
			else:
				delta_lon = delta_lat/math.cos(eq.lat*deg2rad)
			#
			# and let's also assume we want to limit our ETAS map to the input lat/lon:
			lon_min, lon_max = max(eq.lon - delta_lon, self.lons[0]), min(eq.lon + delta_lon, self.lons[1])
			lat_min, lat_max = max(eq.lat - delta_lat, self.lats[0]), min(eq.lat + delta_lat, self.lats[1])
			#
			#print ("lon, lat range: (%f, %f), (%f, %f):: m=%f, L_r=%f, dx=%f/%f" % (lon_min, lon_max, lat_min, lat_max, eq.mag, eq.L_r, eq.L_r*self.etas_range_factor, eq.L_r*self.etas_range_factor/deg2km))
			#
			# - choose an elliptical transform: equal-area, rotational, etc.
			# - calculate initial rate-density
			# - distribute over relevant sites:
			#    - for each site, D' = D_spherical * R'/R_xy	of course this will be approximately equal to R'.
			#
			#for lon_site,lat_site in [[j,k] for j in numpy.arange(lon_min, lon_max, d_lon) for k in numpy.arange(lat_min, lat_max, d_lat)]:
			for lon_site, lat_site in itertools.product(numpy.arange(lon_min+self.bin_lon0, lon_max+self.bin_lon0, self.d_lon), numpy.arange(lat_min+self.bin_lat0, lat_max+self.bin_lat0, self.d_lat)):
				# so we make a square (or maybe a circle later on) and calc. etas at each site. use Bindex() to correct for any misalignment.
				#
				lon_bin = self.lattice_sites.get_xbin_center(lon_site)	# returns a dict like {'index':j, 'center':x}
				lat_bin = self.lattice_sites.get_ybin_center(lat_site)
				#
				#print("lon-lat bins: ", lon_bin, lat_bin)
				#
				#bin_lonlat = xy2_lon_lat(x_bin['center'], y_bin['center'])
				#bin_lonlat=[lon_bin, lat_bin]
				#
				# now, calculate local ETAS intensity from this eq at this site...
				# (this should be defined in the Earthquake() object:
				# note: at this time, we have only the self-similar type local intensity. eventually, we need to split this up. nominally,
				# this can occur at the Earthquake level, so -- if we so desire, different earthquakes can have differend spatial distributions.
				#
				local_intensity = eq.local_intensity(t=self.t_forecast, lon=lon_bin['center'], lat=lat_bin['center'], p=self.p_etas, t_int_lim=self.t_int_lim) # anything else?
				#
				#self.lattice_sites.add_to_bin(x=x_bin['center'], y=y_bin['center'], z=1.0/distances['geo'])
				#
				#self.lattice_sites.add_to_bin(x=lon_bin['center'], y=lat_bin['center'], z=local_intensity)	# note: if we give x,y instead of the bin index, bindex
				self.lattice_sites.add_to_bin(bin_x=lon_bin['index'], bin_y=lat_bin['index'], z=local_intensity)
				
					
			#
		self.ETAS_array = self.lattice_sites.to_array()
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
		return int(round(bin_0) + int(round(round((x-x0)/dx))

	def bin2x(bin_num, dx=1., x0=0., bin_0=0):
		# return the position x of a bin along one axis.
		return (bin_num - bin_0)*dx + x0x
	'''
#
class ETAS_bindex(Global_ETAS_model):
	# default case...
	def __init__(self, *args, **kwargs):
		self.make_etas=self.make_etas_bindex
		super(ETAS_bindex,self).__init__(*args, **kwargs)
#
class ETAS_brute(Global_ETAS_model):
	def __init__(self, *args, **kwargs):
		self.make_etas = self.make_etas_all
		super(ETAS_brute, self).__init__(*args, **kwargs)
		#
	#
#
class ETAS_rtree(Global_ETAS_model):
	def __init__(self, *args, **kwargs):
		self.make_etas = self.make_etas_rtree
		super(ETAS_rtree, self).__init__(*args, **kwargs)
		#
	#
#
class ETAS_rtree_mpp(ETAS_rtree, mpp.Process):
	''' 
	# an mpp version of ETAS_rtree (to be run as an mpp process, aka not an mpp handler, but an mpp worker**).
	'''	
	def __init__(self, pipe_r, *args, **kwargs):
		self.make_etas = self.make_etas_rtree
		kwargs['calc_etas']=False
		self.pipe_r = pipe_r
		#print('initting worker, xyz_range = ', locals().get('etas_xyz_range', 'none'), kwargs.get('etas_xyz_range', 'nnonee'))
		super(ETAS_rtree_mpp, self).__init__(*args, **kwargs)
		mpp.Process.__init__(self)
		#
	#
	def run(self):
		# the mpp run bit.
		# i think we need to disable the 'execute make etas' bit in __init__ (super); we need to kick off ETAS in run()
		self.make_etas()
		#
		# now, pipe back the results
		print('etas complete (from mpp_rtree run() loop); now pipe back(%s)' % self.etas_cat_range)
		self.pipe_r.send(self.ETAS_array)		# might need to pipe back in a simpler form. eventually, we want to use a shared memory mpp.Array()...
		self.pipe_r.close()
#
class ETAS_mpp_handler(Global_ETAS_model):
	# a sub/container class to manage and collect results from a bunch of mpp_etas instances.
	# Note: this version works, but can be super memory-intensive for large catalogs. its cousin class, ETAS_mpp_handler_xyz
	# is recommended for almost all operations (to the extent that this version should probalby be depricated and removed).
	#
	def __init__(self, n_cpu=None, *args, **kwargs):
		#self.make_etas = self.make_etas_rtree
		#
		self.etas_kwargs = {key:val for key,val in kwargs.items() if key not in ['n_proocesses']}		# or any other kwargs we want to skip. note: might need to handle
																										# the ETAS_range parameter carefully...
		#
		#if n_cpu is None: n_cpu = max(1, mpp.cpu_count()-1)
		n_cpu = (n_cpu or mpp.cpu_count())
		self.n_cpu = n_cpu
		#
		# go ahead and run the base class __init__. this basically handles lats, lons, forecast_time, and some other bits and then kicks off make_etas(), where we'll
		#kwargs['make_etas']=False
		self.make_etas = self.make_etas_mpp
		#
		# now, when we run __init__, it will build the catalog (unless one is provided). we'll build up a set of etas_rtree objects. we need to sort out how
		# to generalize. nominally, we set up a sort of parallel-inherited class structure for mpp_in_general and mpp_specific inheritance.
		#
		# define the ETAS worker class:
		self.ETAS_worker = ETAS_rtree_mpp
		#
		#
		super(ETAS_mpp_handler, self).__init__(*args, **kwargs)
	#
	def make_etas_mpp(self):
		#
		#
		cat_len = len(self.catalog)
		#proc_len = cat_len/self.n_cpu
		proc_len = (self.etas_cat_range[1]-self.etas_cat_range[0])/self.n_cpu
		# first, gather a bunch of rtree processes.
		#
		etas_workers = []
		etas_pipes   = []
		prams = self.etas_kwargs.copy()
		prams['catalog']=self.catalog		# .copy() ??
		for j in range(self.n_cpu):
			pipe_r, pipe_s = mpp.Pipe()
			prams['etas_cat_range'] = [int(self.etas_cat_range[0]+j*proc_len), int((j+1)*proc_len)]		# (sort out these indices to be sure we always get the whole catalog...)
			prams['pipe_r'] = pipe_r
			etas_pipes += [pipe_s]
			#
			etas_workers += [self.ETAS_worker(**prams)]
		#
		# now, go through the list again and start each intance (this can probably be done when they're instantiated):
		for j,p in enumerate(etas_workers):
			p.start()
			#etas_workers[j].start()
			#
		#
		# and join() them? do we need to join if we're doing send()-recv()? (see vq code for examples):
		# in fact, i think we don't need this, and be careful with the join() syntax (see notes in vq_analyzer code).

		# now, processes are finished; they return ETAS_array[] like [[x,y,z]]... ETAS_array should be initiated. we can aggregate them by index, but they should all
		# have the same indexing (they should all be a complete set of lats, lons. so we should check for that here... and i think actually we should have an empty copy
		# from __init__(), but for now, just add them row by row (maybe sort first). in later versions (maybe this one??) all processes will write directly to an mpp.Array()
		# shared memory object.
		#
		print('now gather sub-arrays...')
		for j,(etas,pp) in enumerate(zip(etas_workers, etas_pipes)):
			ary = pp.recv()
			for k,rw in enumerate(ary): self.ETAS_array['z'][k]+=rw['z']
			pp.close()
		#
		# still not sure if we need this...
		#for j,p in enumerate(etas_workers):
		#	p.join()
		#
		del etas_workers			
#
#
class ETAS_mpp_handler_xyz(Global_ETAS_model):
	# a sub/container class to manage and collect results from a bunch of mpp_etas instances.
	# this is a semi-memory friendly version in which the forecast lattice is divided amongst the processes,
	# as opposed to ETAS_mpp_handler() in which the catalog is split up, which uses a lot of memory (and probably lots
	# of time piping results) for large maps.
	#
	# TODO: Figure out how to derive this from mpp.Pool() and, rather than manually managing processes, run as a Pool().
	# the main problem is that we divide the job into jobs based on geography (aka, slice up the map and give each processor
	# a section of map), but this is not an accurate representation of the actual compute requirements. seismically active
	# regions end up doing way more flops than quiescent regions. a good, simple approach, then, is to divide into, say
	# 2*n_cpu() jobs, which we process n_cpu() at a time using a Pool(). this way, non-intensive jobs will be discarded and replaced
	# quickly, and the compute intensive jobs will be smaller, as a product of smaller geometry.
	#
	def __init__(self, n_cpu=None, worker_class=ETAS_rtree_mpp, *args, **kwargs):
		#self.make_etas = self.make_etas_rtree
		#
		self.etas_kwargs = {key:val for key,val in kwargs.items() if key not in ['n_proocesses']}		# or any other kwargs we want to skip. note: might need to handle
																										# the ETAS_range parameter carefully...
		#
		if n_cpu is None: n_cpu = max(1, mpp.cpu_count()-1)
		n_cpu = (n_cpu or mpp.cpu_count())
		self.n_cpu = int(numpy.ceil(n_cpu))
		#
		# go ahead and run the base class __init__. this basically handles lats, lons, forecast_time, and some other bits and then kicks off make_etas(), where we'll
		#kwargs['make_etas']=False
		self.make_etas = self.make_etas_mpp
		#
		# now, when we run __init__, it will build the catalog (unless one is provided). we'll build up a set of etas_rtree objects. we need to sort out how
		# to generalize. nominally, we set up a sort of parallel-inherited class structure for mpp_in_general and mpp_specific inheritance.
		#
		# define the ETAS worker class:
		#self.ETAS_worker = ETAS_rtree_mpp
		self.ETAS_worker = worker_class
		#
		#
		super(ETAS_mpp_handler_xyz, self).__init__(*args, **kwargs)
	#
	def make_etas_mpp(self):
		#
		cat_len = len(self.catalog)
		##proc_len = cat_len/self.n_cpu
		#proc_len = (self.etas_cat_range[1]-self.etas_cat_range[0])/self.n_cpu
		#
		xyz_len = self.n_lat*self.n_lon/self.n_cpu
		#
		# first, gather a bunch of rtree processes.
		#
		etas_workers = []
		etas_pipes   = []
		prams = self.etas_kwargs.copy()
		prams['etas_cat_range']=None
		prams['catalog']=self.catalog		# .copy() ??
		for j in range(self.n_cpu):
			pipe_r, pipe_s = mpp.Pipe()
			prams['etas_xyz_range'] = [int(j*xyz_len), int((j+1)*xyz_len)]		# (sort out these indices to be sure we always get the whole catalog...)
			print('etas_mpp worker xyz_range: ', prams['etas_xyz_range'])
			prams['pipe_r'] = pipe_r
			etas_pipes += [pipe_s]
			#
			etas_workers += [self.ETAS_worker(**prams)]
		#
		# now, go through the list again and start each intance (this can probably be done when they're instantiated):
		for j,p in enumerate(etas_workers):
			p.start()
			#
		#
		# and join() them? do we need to join if we're doing send()-recv()? (see vq code for examples):
		# in fact, i think we don't need this, and be careful with the join() syntax (see notes in vq_analyzer code).

		# now, processes are finished; they return ETAS_array[] like [[x,y,z]]... ETAS_array should be initiated. we can aggregate them by index, but they should all
		# have the same indexing (they should all be a complete set of lats, lons. so we should check for that here... and i think actually we should have an empty copy
		# from __init__(), but for now, just add them row by row (maybe sort first). in later versions (maybe this one??) all processes will write directly to an mpp.Array()
		# shared memory object.
		#
		self.ETAS_array = []
		#self.arys = []
		#
		print('now gather sub-arrays...')
		'''
		for j,(etas,pp) in enumerate(zip(etas_workers, etas_pipes)):
			#ary = pp.recv()
			#self.arys += [ary]
			#for k,rw in enumerate(ary): self.ETAS_array['z'][k]+=rw['z']
			#self.ETAS_array += list(ary)
			self.ETAS_array += list(pp.recv())
			pp.close()
			# can we delete each worker here? might be necessary to reduce max memory footprint. no. instead, use an mpp.Queue() (i think) or a while-loop
			# and delete off the top or bottom of the list after each pipe.close()
		'''
		for j in range(len(etas_workers)):
			self.ETAS_array += list(etas_pipes[0].recv())
			etas_pipes[0].close()
			del etas_pipes[0]
			del etas_workers[0]
		#
		self.ETAS_array.sort(key=lambda rw: (rw[1],rw[0]))
		self.ETAS_array = numpy.core.records.fromarrays(zip(*self.ETAS_array), dtype = [('x', '>f8'), ('y', '>f8'), ('z', '>f8')])
		
		# still not sure if we need this...
		#for j,p in enumerate(etas_workers):
		#	p.join()
		#
		del etas_workers			
#
class ETAS_mpp(ETAS_mpp_handler_xyz):
	# container for default mpp handler.
	def __init__(self, *args, **kwargs):
		super(ETAS_mpp,self).__init__(*args, **kwargs)
		# in other words, just be the parent class.
#
class Earthquake(object):
	# an Earthquake object for global ETAS. in parallel operations, treat this more like a bag of member functions than a data container.
	# pass an earthquake catalog list to a process; use Earthquake() to handle each earthquake event row.
	# include "local" ETAS calculation in this object; aka, each earthquake determines its ETAS range based on rupture length and other factors...
	# maybe.
	#def __init__(self, event_date=None, lat=0., lon=0., mag=None, eig_vals=[1., 1.], eig_vecs=[[1.,0.], [0., 1.]], transform_type='equal_area'):
	def __init__(self, dict_or_recarray, transform_type='equal_area', transform_ratio_max=2.0, ab_ratio_expon=1):
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
		self.ab_ratio_expon = ab_ratio_expon
		#
		# yoder: be sure eigen-tuples are properly sorted (eventually anyway).
		# TODO: test this eigen-value/vector sorting bit. we might be able to strip out a clumsy "if" logic down the way.
		# be sure we're not changing something we 'fixed' by doing this the right way.
		##e_vals, e_vecs = list(zip(*sorted([[lamb, evec] for lamb, evec in zip(e_vals, e_vecs.T)], key=lambda rw:rw[0])))
		
		e_vals, e_vecs = list(zip(*reversed(sorted(list(zip(self.e_vals, self.e_vecs.T)), key=lambda rw:rw[0]))))
		self.e_vecs = numpy.array(e_vecs).T		# put the eigen-vectors back into column matrix form.
		self.e_vals = numpy.array(e_vals)
		#
		# check for float datetime...
		if not hasattr(self, 'event_date_float'): self.event_date_float = mpd.date2num(self.event_date.tolist())
		#self.event_date_float_secs = self.event_date_float*days2secs
		#
		# yoder 2016-8-7:
		# so actually, the ab_ratio should be sqrt(\lambda_1/\lambda_2).
		# in PCA, the eigen-values are basically the magnitude of variance (aka, stdev**2) in that eigen-direction, so the approximate length
		# of the spatial vectors is the sqrt() of eigen-value. note then that for the equal-area transform, we take another sqrt() so that
		# a_ea = r*(\lambda_1 / \lambda_2)**.25
		# b_ea = r*(\lambda_2 / \lambda_1)**.25
		#
		#self.ab_ratio = min(transform_ratio_max, max(self.e_vals)/min(self.e_vals))
		#self.ab_ratio_raw = max(self.e_vals)/min(self.e_vals)
		#
		#self.ab_ratio_raw = math.sqrt(abs(max(self.e_vals)/min(self.e_vals)))
		#self.ab_ratio_raw = (math.sqrt(abs(max(self.e_vals)/min(self.e_vals))) if not min(self.e_vals)==0. else transform_ratio_max)
		#
		# yoder: but we don't need to set ab_ratio (in any form) here. let's do all of this in set_transform()
		#self.ab_ratio_raw = ((abs(max(self.e_vals)/min(self.e_vals)))**.5 if not min(self.e_vals)==0. else transform_ratio_max)
		self.set_transform()
		#
		#####
		# now, make some preliminary calculations, namely the peak spatial density and maybe some omori constants? anything we'll calculate
		# again and again...
		#
		# this gets done in set_transform()
		#self.spatial_intensity_factor=1.0		# corrects for local aftershock density in rotational type transforms.
	#	
	@property
	def lat_lon(self):
		return [self.lat, self.lon]
	@property
	def lon_lat(self):
		return [self.lon, self.lat]		
	#
	@property
	def theta_rad(self):
		return self.eq_theta(deg_rad='rad')
	@property
	def theta_deg(self):
		return self.eq_theta(deg_rad='deg')
	@property
	def theta(self):
		return self.eq_theta(deg_rad='deg')
	#
	def eq_theta(self, deg_rad='deg'):
		dr_factor = deg2rad
		if deg_rad=='rad': dr_factor=1.0
		#
		return math.atan(self.e_vecs[0][0]/self.e_vecs[0][1])/deg2rad
	#
	def dist_to(self, lon=None, lat=None, dist_types=['spherical'], Rearth = 6378.1, *args, **kwargs):
		if hasattr(lon, '__len__'):
			lon = lon[0]
			lat = lat[1]
		#
		return dist_to(lon_lat_from=self.lon_lat, lon_lat_to=[lon, lat], dist_types=dist_types, Rearth=Rearth, *args, **kwargs)
	#
	def set_transform(self, e_vals=None, e_vecs=None, transform_type=None, transform_ratio_max=None, ab_ratio_expon=None):
		'''
		# note: it might be better to not allow custom-runs (aka not allow parameterized calls to this function; always use native member values).
		#
		# define the elliptical transformation for calculating intensities.
		# transform_type: elliptical transformation type = { 'equal_area', 'rotation'}
		# transform_ratio_max: maximum allowable eigenvalue (and transformed axes length) ratio.
		#
		#... and let's clean up a bit, removing some of the inverse, etc. variants.
		#
		# TODO: this sub-script seems to conflict a bit with __init__, where we also set ab_ratio... right? for now, we've adjusted the code to be consistent.
		# in the very, very near fugure, we need to clean it up to be not redundant.
		#
		'''
		#
		# get transform variables. all raw, proper values, so the eigenvalues are the variances, on par
		# with a/b ~ sqrt(lambda_1/lambda_2)
		e_vals = (e_vals or self.e_vals)
		e_vecs = (e_vecs or self.e_vecs)
		transform_ratio_max = transform_ratio_max or self.transform_ratio_max #float(transform_ratio_max or self.transform_ratio_max)
		transform_type = (transform_type or self.transform_type)
		ab_ratio_expon = (ab_ratio_expon or self.ab_ratio_expon)
		ab_ratio_expon = (ab_ratio_expon or .5)
		#
		# now, sort eigenvectors by eigenvalue:
		# actually, let's do this at the __init__() level.
		# note: eigen vectors are columns of e_vecs array.
		# TODO: test this eigen-value/vector sorting bit. we might be able to strip out a clumsy "if" logic down the way.
		# be sure we're not changing something we 'fixed' by doing this the right way.
		##e_vals, e_vecs = list(zip(*sorted([[lamb, evec] for lamb, evec in zip(e_vals, e_vecs.T)], key=lambda rw:rw[0])))
		#e_vals, e_vecs = list(zip(*reversed(sorted([list(zip(e_vals, e_vecs.T)), key=lambda rw:rw[0]))))
		#e_vecs = e_vecs.T		# put the eigen-vectors back into column matrix form.
		#print('**debug, evals, evecs: ', e_vals, e_vecs)
		#
		# notes on ab_ratio: in the strictest sense, ab_ratio expon. should be 0.5, in the sense that the 'singular values' of the decomposition are equal to
		# the sqrt(eigen_values) of the covariance (which makes sense; the basis lengths are approximately the standard deviation in some direction;
		# the covariance eigenvalues are variance). 
		# let's catch the cases where an eigenvalue=0...
		#
		# are we getting negative eigenvalues (i think we shouldn't, but f we do, just keep the magnitude)
		# and strip off the imaginary component.
		abs_evals = numpy.abs(numpy.real(e_vals))
		# then, biggest lambda / smallest lambda, unless the small e-val is 0.
		# ... and i think this is carried over from anther way of doing this, long long ago. the better approach is to just get
		# the ratio of the eigen-values, limit it to too big/too small, and apply in a proper linear transformation... later
		if min(abs_evals)==0.:
			ab_ratio = transform_ratio_max**ab_ratio_expon
		elif transform_ratio_max == None:
			ab_ratio = (max(abs_evals)/min(abs_evals))**ab_ratio_expon
		else:
			ab_ratio = transform_ratio_max**ab_ratio_expon #min(transform_ratio_max**ab_ratio_expon, (max(abs_evals)/min(abs_evals))**ab_ratio_expon)
		#
		self.ab_ratio=ab_ratio
		#print('**debug: self.ab_ratio: ', self.ab_ratio)
		#
		if transform_type=='equal_area':
			# so that pi*r^2 = pi*ab = pi*b*(ab_ratio)*b
			# 
			# "normalized eigenvalues" is not the best name, but it's what we have right now... but now, these are the lengths
			# of the ellipitical axes.
			#self.e_vals_n = list(reversed(sorted([abs(ab_ratio), 1./abs((ab_ratio))])))
			#
			# this works, but it's sloppy; clearly, i'm missing someting (in my current foggy condition) about the right way to sort
			# and dot these vectors. the right thing to do is, i think, to just do a proper linear transformation,
			# M = diag( (lambda1/lambda2)**b, (lambda2/lambda1)**b ), but also handle the large/small value exceptions. maybe use
			#  a min/max inside: norm_eval_0 = max(1/x0, min(x0, lambda1/lambda2)), norm_eval_1 = max(1/x0, min(x0, lambda1/lambda2)),
			# or maybe sort the eigen-vectors by eigen_values (but be careful we're still getting the correct transformation)
			# when we consider that we have to catch the 0 valued eigenvalue as well as extreme values, this approach is maybe
			# not so bad...
			#
			# this works, except we need to (properly) incorporate the min/max ab_ratio filtering.
			'''
			if min(abs_evals)==0:
				self.e_vals_n = [transform_ratio_max, 1./transform_ratio_max]
			else:
				#self.e_vals_n = [(abs_evals[1]/abs_evals[0])**ab_ratio, (abs_evals[0]/abs_evals[1])**ab_ratio]
				self.e_vals_n = [max(1./transform_ratio_max, min(transform_ratio_max, (abs_evals[1]/abs_evals[0])))**ab_ratio, max(1./transform_ratio_max, min(transform_ratio_max, (abs_evals[0]/abs_evals[1])))**ab_ratio]
			'''
			#
			self.e_vals_n = [1./abs(ab_ratio), abs((ab_ratio))]
			#if abs_evals[0]<abs_evals[1]:
			#	#self.e_vals_n = [abs(ab_ratio), 1./abs((ab_ratio))]
			#	self.e_vals_n = [ab_ratio, 1./ab_ratio]
			#else:
			#	#self.e_vals_n = [1./abs(ab_ratio), abs((ab_ratio))]
			#	self.e_vals_n = [1./ab_ratio, ab_ratio]
			#
			self.spatial_intensity_factor = 1.0
			#
		#	
		elif transform_type=='rotation':
			# TODO: this rotation projection is not properly tested just yet... and it should be.
			#
			# set the *small* eigen-value --> 1; scale the large one.
			# remember (i think) that this geometry is inverted. we leave the short eigen-vector=1, so when we calc r', it is unchanged in this direction.
			# in the orothogonal direction, r'>>r, and we have etas halos elongated along the "short" (rupture) axis).
			#
			#self.e_vals_n = [min(transform_ratio_max, x/min(e_vals)) for x in e_vals]
			self.e_vals_n  = [1., ab_ratio]
			#if abs_evals[0]<abs_evals[1]:
			#	self.e_vals_n  = [ab_ratio, 1.]
			#else:
			#	self.e_vals_n  = [1., ab_ratio]
			#
			self.spatial_intensity_factor = min(self.e_vals_n)/max(self.e_vals_n)
			#
		else:
			return self.set_transform(e_vals=e_vals, e_vecs=e_vecs, transform_type='equal_area')
		#
	#
	def elliptical_transform(self, lon=0., lat=0.):
		#
		# ... no idea if this works, but it should be about right. perhaps a bit redundant in the LT's ??
		# simple elliptical transform:
		# get v = [lon,lat]-[self.lon, self.lat]
		# rotate v with transpose: v_prime = numpy.dot(self.e_vecs.transpose(),v)
		#    -- x_prime = v dot x_prime
		#    -- y_prime = v dot y_prime
		#
		# yoder:
		# TODO: it may be necessary to handle cross-international date line modulus operations here. self.lat/lon should be fine,
		# but we might get a less than letal lat/lon in some cases (which should be handled upstream, but why not here as well?).
		# let's just have a sloppy go at it for now:
		# this can be done with a modulus operator (maybe abs(lon)%180 ?), but we still have to stitch together the parity case.
		# ... and now we need to see how that affects dist_to(), and in fact this might be the better place to handle these issues.
		#if lon>180:
		#	lon -= 180
		#if lon<-180:
		#	lon += 180
		#
		# TODO: are we using this? let's clean it up and probably not use the 'geo' distance... in fact, not calculate it
		# (i think we have actually removed this, for the reason that the geo distance calc. is expensive.).
		dists = dist_to(lon_lat_from=[self.lon, self.lat], lon_lat_to=[lon, lat], dist_types=['geo', 'xy', 'dx_dy'])
		R = dists['geo']	# the ['s12'] item from Geodesic.Inverse() method, (and converted to km in dist_to() )...
		dx,dy = dists['dx'], dists['dy']	# cartesian approximations from the dist_to() method; approximate cartesian distance coordinates from e_quake center in
		#
		#v = [dx, dy]
		#v_prime = numpy.dot(self.e_vecs.transpose(),v)
		#
		# this should probably be reorganized to look like a prober singular value decomposition (SVD).
		# 1) rotate to pca frame (v' = numpy.dot(e_vecs.T, v)
		# 2) elongate: v'' = numpy.dot(numpy.diag(reversed(self.e_vals_n)), v' ) ## check reversed, etc.
		# 3) now, rotate back to original frame: v''' = numpy.dot(e_vecs, v'')
		# TODO: is this right, or should we be dotting to e_vecs.not_transpose()?
		#dx_prime = self.e_vals_n[1]*numpy.dot(v,self.e_vecs.transpose()[0])
		#dy_prime = self.e_vals_n[0]*numpy.dot(v,self.e_vecs.transpose()[1])		# ... but isn't this the rotation transform?
		dx_prime = self.e_vals_n[0]*numpy.dot([dx, dy],self.e_vecs.T[0])
		dy_prime = self.e_vals_n[1]*numpy.dot([dx, dy],self.e_vecs.T[1])		# ... but isn't this the rotation transform?
		#
		R_prime = R * numpy.linalg.norm([dx_prime,dy_prime])/numpy.linalg.norm([dx,dy])
		
		dists.update({'R':R, 'R_prime':R_prime, 'dx':dx, 'dy':dy, 'dx_prime':dx_prime, 'dy_prime':dy_prime})
		#
		return dists
		#
	#
	def elliptical_transform_depricated(self, lon=0., lat=0., T=None):
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
		# T is the transformation matrix from the PCA.
		if T is None: T=self.T
		E_hat=self.E_hat		# like (e_vals_n x I)
		
		#T = (T or self.T)	# note: this syntax sometimes fails for numpy.array() objects; gives a "multiple truth" error, in other words, " if each element is None, then..."
		
		#T = numpy.dot([[self.e_vals_n[0], 0.],[0., self.e_vals_n[1]]], self.e_vecs.transpose())
		#
		# diagnostics:
		#T = self.e_vecs.transpose()
		#
		# first, get the actual distance (and related data) to the point:
		
		dists = dist_to(lon_lat_from=[self.lon, self.lat], lon_lat_to=[lon, lat], dist_types=['geo', 'xy', 'dx_dy'])
		R = dists['geo']	# the ['s12'] item from Geodesic.Inverse() method, (and converted to km in dist_to() )...
		dx,dy = dists['dx'], dists['dy']	# cartesian approximations from the dist_to() method; approximate cartesian distance coordinates from e_quake center in general FoR.
		
		#dx = (self.lon-lon)*111.3
		#dy = (self.lat-lat)*111.3
		#R=(dx*dx + dy*dy)**.5
		
		#print('dists: ', dists)
		#
		# some diagnostics:
		#dists = dist_to(lon_lat_from=[self.lon, self.lat], lon_lat_to=[lon, lat], dist_types=['geo', 'euc', 'xy', 'dx_dy'])
		#R = dists['euclidian']
		#dx,dy = dists['dx_e'], dists['dy_e']		
		#
		dx_prime, dy_prime = numpy.dot(T, [dx,dy])	# rotated and dilated into the elliptical FoR...
		#dx_prime, dy_prime = numpy.dot(self.e_vecs.inverse(), [dx,dy])
		#
		#R_prime = R*((dx_prime*dx_prime + dy_prime*dy_prime)/(dx*dx+dy*dy))**.5	# use this to mitigate artifacts of the spherical transform: R_prime = R_geo*(R_prime_xy/R_xy)
		#
		R_prime = R*numpy.linalg.norm([dx_prime, dy_prime])/numpy.linalg.norm([dx,dy])
		#R_prime = numpy.linalg.norm([dx_prime, dy_prime])
		#
		dists.update({'R':R, 'R_prime':R_prime, 'dx':dx, 'dy':dy, 'dx_prime':dx_prime, 'dy_prime':dy_prime})
		#return {'R':R, 'R_prime':R_prime, 'dx':dx, 'dy':dy, 'dx_prime':dx_prime, 'dy_prime':dy_prime}
		#
		return dists
	#
	def spherical_dist(self, to_lon_lat=[]):
		return shperical_dist(lon_lat_from=[self.lon, self.lat], lon_lat_to=to_lon_lat)
	#
	def local_intensity(self, t=None, lon=None, lat=None, p=None, q=None, t0_prime=None, t_int_lim=None):
		'''
		# t0_prime: use this to re-scale the temporal component. we'll transform the initial rate (and effectively minimum delta_t)
		# to avoid near-field artifacts (where a recent earthquake dominates the ETAS).
		# calculate local ETAS density-rate (aka, earthquakes per (km^2 sec)
		# take time in days.
		#
		# TODO: we can speed up large ETAS by pre-calculating (or more specifically, calculating only once) the temporal rate.
		#       this could be done pre-mpp, but since the remaining processes would then have to wait for this to finish, 
		#       not much gain would be realized; we're probably better off just computing the rates on each process (keep it simple...)
		'''
		#print("inputs: ", t, lon, lat, p, q)
		t = (t or mpd.date2num(dtm.datetime.now(pytz.timezone('UTC'))))
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
		# get elliptial transformation for this earthquake:
		#et = self.elliptical_transform_prime(lon=lon, lat=lat)		# use "inverse" transformation (see code), so x' = numpy.dot(x,T)
		#et = self.elliptical_transform(lon=lon, lat=lat, T=self.T_inverse)
		et = self.elliptical_transform(lon=lon, lat=lat)
		#
		# in some cases, let's adjust t0 so maybe we dont' get near-field artifacts...
		'''
		if t0_prime!=None:
			t_0_prime = 60.*10.	# ten minutes...
			tau_prime = (self.tau*self.t_0**self.p)/(t_0_prime)	#this looks wrong; note we have tau in there twice...
			#
			#orate = 1.0/(tau_prime * (t_0_prime + delta_t)**p)
		else:
			#orate = 1.0/(self.tau * (self.t_0 + delta_t)**p)
			tau_prime = self.tau
			t_0_prime = self.t_0
		'''
		#
		# note: everything from here down can (i think) be compiled with numba.jit, if we are so inclined.
		#
		#orate = 1./(tau_prime * (t_0_prime + delta_t)**p)
        
		# Wilson: If we're doing a time integration, we do that in the make_etas_rtree method
		if not t_int_lim:
			orate = 1./(self.tau * (self.t_0 + delta_t)**p)
		else:
			orate = 1.0
			
		#
		# for now, just code up the self-similar 1/(r0+r) formulation. later, we'll split this off to allow different distributions.
		# radial density (dN/dr). the normalization constant comes from integrating N'(r) --> r=inf = N_om (valid, of course, ony for q>1).
		#
		radial_density = (q-1.0)*(self.r_0**(q-1.0))*((self.r_0 + et['R_prime'])**(-q))
		#
		# ... and this is distributed along an elliptical contour. we could approximate with a circle, but what we really want is to distribute along the ellipse
		# that is orthogonal to our R-R' transformation. fortunately, they have the same radius. calculating the radius of an ellipse is hard. we an exact solution
		# (that uses a "special" function) and an approximation (that should be faster).
		# note major,semi-major axes are R*e_1, R*e_2 (doesn't matter which is which)
		#
		# radial normalization:
		# notes:
		# - using an R_prime geometry (circle or ellipse) produces nice, elliptical looking contours. normalizing on the physical distances/geometries R,
		# produces diamond like shapes. so which one is the artifact?
		# - normalize with r -> r0 + r_prime, (aka, the distance as "seen" by the earthquake). if we normalize with geometric distances, we get singularities around the earthquakes. 
		#circumf = ellipse_circumference_approx1(a=self.e_vals_n[0]*et['R'], b=self.e_vals_n[1]*et['R'])
		#circumf = ellipse_circumference_approx1(a=self.e_vals_n[0]*et['R_prime'], b=self.e_vals_n[1]*et['R_prime'])
		#
		# ummm... this is area? i think this arose from some desperate trouble-shooting to find the singular behavior
		#  (that is resolved by including the r0 term in the effective radius).
		#circumf = ellipse_circumference_approx1(a=self.e_vals_n[0]*(self.r_0 + et['R_prime']), b=self.e_vals_n[1]*(et['R_prime'] + self.r_0))
		#circumf = math.pi*et['R']**2.
		#circumf = math.pi*(self.r_0#etas = gep.ETAS_mpp(n_cpu=2*mpp.cpu_count(), catalog=this_catalog, **eq_prams) + et['R_prime'])**2.
		
		circumf = 2.*math.pi*(self.r_0 + et['R_prime'])
		#
		# @spatial_intensity_factor: corrects for local aftershock density in rotational type transforms.
		#  aka, in rotations, space is effectively compressed (by trigonometry), so intensity is boosted. for equal-area
		#  transforms, spatil_intensity_factor=1.0
		spatialdensity = self.spatial_intensity_factor*radial_density/circumf
		#		
		return spatialdensity*orate
		#return orate*self.spatial_intensity_factor*radial_density/circumf		
	#
	# some diagnostics:
	def plot_linear_density(self, r_max=None, r_max_factor=5., fignum=0, n_points=1000):
		'''
		# a diagnostic plot of thie earthquake's linear density distribution. r_max is, as it sounds, the max value of r to plot. alternatively, use r_max_factor, 
		# for r_max = r_max_factor * l_r
		'''
		#
		r_max=(r_max or r_max_factor*self.L_r)
		print("r_max: ", r_max)
		#
		plt.figure(fignum)
		plt.clf()
		X=numpy.arange(0., r_max, r_max/float(n_points))
		ax=plt.gca()
		ax.set_xscale('log')
		ax.set_yscale('log')
		ax.plot(X, [1.0/(self.chi*(self.r_0 + r)**self.q) for r in X], '.-', lw=2.)
		#	
	#
	def sketch_evecs(self, fignum=0):
		plt.figure(fignum)
		plt.clf()
		#
		# axes:
		x_max = max(self.e_vals_n)
		plt.plot([x*x_max for x in [-1., 1.]], [0.,0.], 'k-', lw=2)
		plt.plot( [0.,0.], [x*x_max for x in [-1., 1.]], 'k-', lw=2)
		plt.plot([-1., 1.], [0.,0.], 'k-', lw=2)
		plt.plot( [0.,0.], [-1., 1.], 'k-', lw=2)
		#
		plt.plot([0., self.e_vecs[0][0]], [0., self.e_vecs[0][1]], 'mo-')
		plt.plot([0., self.e_vecs[1][0]], [0., self.e_vecs[1][1]], 'co-')
		#
		T = numpy.dot(self.E_hat, self.e_vecs)
		plt.plot([0., T[0][0]], [0., T[0][1]], 'bo-')
		plt.plot([0., T[1][0]], [0., T[1][1]], 'go-')
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
		for key,val in defaults.items():
			#if not kwargs.has_key(key): kwargs[key]=val
			#if not key in kwargs.keys(): kwargs[key]=val
			kwargs[key] = kwargs.get(key, val)
		#
		plt.figure(fignum)
		if do_clf: plt.clf()
		if ax is None: ax=plt.gca()
		#
		#if self.poly is None or len(self.poly==0):
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
		if a is None and b is None: a,b = 1.0, .5
		if not (a is None and b is None): ab_ratio=a/float(b)
		#
		if a is None: a=b*ab_ratio
		if b is None: b=a/ab_ratio
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
		print( "theta: %f" % (self.theta))
		if self.theta!=0.:
			#poly = zip(*numpy.dot( [[math.cos(self.theta), -math.sin(self.theta)],[math.sin(self.theta), math.cos(self.theta)]], zip(*poly)))
			poly = numpy.dot(poly, [[math.cos(self.theta), math.sin(self.theta)],[-math.sin(self.theta), math.cos(self.theta)]])
		#
		return poly
	#
	def plot_poly(self, fignum=0, do_clf=True, poly_len=100):
		plt.figure(fignum)
		if do_clf: plt.clf()
		#
		#if self.poly is None or len(self.poly==0):
		#	self.poly(poly_len)
		#
		plt.plot(*zip(*self.poly()), color='b', marker='.', ls='-', lw='1.5')
	#	
#
# Working scripts
#
def make_ETAS_catalog_mpp(incat=None, lats=[32., 38.], lons=[-117., -114.], mc=2.5, date_range=['1990-1-1', None], D_fract=1.5, d_lambda=1.76, d_tau = 2.28, fit_factor=1.5, p=1.1, q=1.5, dmstar=1.0, b1=1.0,b2=1.5, do_recarray=False, n_cpus=None):
	# multiprocessing wrapper for make_ETAS_catalog(). note, this would be more efficient and faster if we can figure out how to use shared memory
	# (don't have to make multiple copies of catalog; don't have to pickle back entire catalog). maybe we need to learn python mpi, rather
	# than multiprocessing?
	# ... and this totally does not work (yet?)
	# ... dunno; throwing some sort of "attr must be string" error on the .get().
	#
	etas_prams = locals().copy()
	etas_prams.__delitem__('n_cpus')
	etas_prams['incat']=None
	#etas_prams.__delitem__('self')
	print('etas_prams: ', etas_prams)
	# handle dates:
	if date_range[1] is None: date_range[1] = dtm.datetime.now(pytz.timezone('UTC'))
	#
	for j,dt in enumerate(date_range):
		if isinstance(dt, str):
			date_range[j] = mpd.num2date(mpd.datestr2num(dt))
	#
	if incat is None or (hasattr(incat, '__len__') and len(incat)==0):
		# try to catch network errors:
		n_tries_max = 10
		t_sleep = 5
		n_tries=0
		while n_tries<=n_tries_max:
			try:
				incat = atp.catfromANSS(lon=lons, lat=lats, minMag=mc, dates0=date_range, Nmax=None, fout=None, rec_array=True)
				n_tries = n_tries_max + 1
			except:
				print("network failed, or something. trying again (%d)" % n_tries)
				n_tries+=1
	#
	etas_prams['incat']=incat
	if n_cpus is None:
		#n_cpus = max(1, mpp.cpu_count()-1)
		n_cpus = mpp.cpu_count()
		#
	#
	P=mpp.Pool(n_cpus)
	cat_len=len(incat)
	# yoder (2016-12-4):
	#n = cat_len/n_cpus
	n = int(numpy.ceil(cat_len/n_cpus))
	#
	pool_handlers = []
	pool_results = []
	for k in range(n_cpus):
		# yoder (2016-12-4: int() and ceil() teh catalog_range terms. otherwise, this breaks when we try to process
		# a length 1 catalog.
		#etas_prams['catalog_range']=[k*n, min((k+1)*n, cat_len)]
		catalog_range = [int(k*n), int(numpy.ceil(min((k+1)*n, cat_len)))]
		etas_prams['catalog_range']=catalog_range
		#
		#print("*debug** parameters: ", etas_prams)
		#print("*debug** catalog_range: ", etas_prams['catalog_range'])
		#
		# TODO: convert thsi call to use "kwds" instead of "args"
		# i think, however, we could skip all the lats, lons, etc. and just pass the catalog, catatalog_range, and the omori parameters.
		#my_etas_prams = {'incat':incat, 'lats':lats, 'lons':lons, 'mc':mc, 'date_range':['1990-1-1', None], 'D_fract':1.5, 'd_lambda':1.76, 'd_tau': 2.28, 'fit_factor':1.5, 'p':1.1, 'q':1.5, 'dmstar':1.0, 'b1':1.0, 'b2':1.5, 'do_recarray':True, 'catalog_range':catalog_range}
		my_etas_prams = {'incat':incat, 'lats':lats, 'lons':lons, 'mc':mc, 'date_range':date_range, 'D_fract':D_fract, 'd_lambda':d_lambda, 'd_tau':d_tau, 'fit_factor':fit_factor, 'p':p, 'q':q, 'dmstar':dmstar, 'b1':b1, 'b2':b2, 'do_recarray':True, 'catalog_range':catalog_range}
		#
		# new syntax (older syntax being depricated? undepricated? keep an eye on this in the future...)
		#pool_handlers += [P.apply_async(make_ETAS_catalog(**my_etas_prams))]
		pool_handlers += [P.apply_async(make_ETAS_catalog, kwds=my_etas_prams)]
		#
	P.close()
	#
	for R in pool_handlers:
		#print("*R: ", R)
		R_return = R.get()
		if hasattr(R_return, '__len__'):
			pool_results+=[R_return]
			#pool_results += [numpy.core.records.fromarrays(zip(*R_return[1:]), dtype = R_return[0])]
			
		#else:
		#	pool_results+=[None]
		#pool_results+=[R.get()[0]]
	
	print("results fetched.")
	pool_results = functools.reduce(numpy.append, [pr for pr in pool_results if not pr is None])
	
	#pool_results = functools.reduce(numpy.append, pool_results)
	pool_results.sort(order='event_date_float')
	
	#output_catalog.sort(order='event_date_float')
	#return numpy.core.records.fromarrays(zip(*output_catalog), dtype = dtype)
	return pool_results
	#
	#
# make a class out of this; derive from recarray; keep universal valuess like p,q,dmstar, etc. as member variables.
def make_ETAS_catalog(incat=None, lats=[32., 38.], lons=[-117., -114.], mc=2.5, date_range=['1990-1-1', None], D_fract=1.5, d_lambda=1.76, d_tau = 2.28, fit_factor=1.5, p=1.1, q=1.5, dmstar=1.0, b1=1.0,b2=1.5, do_recarray=True, catalog_range=None, N_min_ellip_fit=5):
	'''
	# fetch (or receive) an earthquake catalog. for each event, calculate ETAS relevant constants including rupture_length,
	# spatial orientation, elliptical a,b, etc. This catalog will then be provided to other processes to calculate ETAS rates
	# @catalog_range: use this for mpp operations. to parallelize this script, create n_cpu instances of make_ETAS_catalog; each instance
	# contains the whole catalog (or work out a shared memory approach), but each instance only fits/processes [catalog_range] of the catalog.
	# then, aggregate the catalogs after they finish processing.
	'''
	#
	# save a copy of the input parameters to ammend to the output array as metadata. we can do this in lieu of creating a derived class -- sort of a (sloppy?) shortcut...
	#
	if incat is None:
		# use all the local inputs:
		meta_parameters = locals().copy()
	else:
		meta_parameters = {}
		if hasattr(incat, 'meta_parameters'): incat.meta_parameters.copy()
		meta_parameters.update({'D_fract':D_fract, 'd_lambda':d_lambda, 'd_tau':d_tau, 'fit_factor':fit_factor, 'p':p, 'q':q, 'dmstar':dmstar, 'b1':b1, 'b2':b2})
		#print('**debug: source catalog provided... {}'.format(len(incat)))
	#
	#dmstar = 1.0
	dm_tau = 0.0		# constant sort of like dm_star but for temporal component. we usually treat this as 0, but don't really know...
	#d_tau = 2.28
	b=b1
	#
	# some constants:
	l_15_factor = (2./3.)*math.log10(1.5)	# pre-calculate some values...
	km2lat = 1./deg2km		# deg3km defined globally
	#
	# handle dates:
	if date_range[1] is None: date_range[1] = dtm.datetime.now(pytz.timezone('UTC'))
	#
	for j,dt in enumerate(date_range):
		if isinstance(dt, str):
			date_range[j] = mpd.num2date(mpd.datestr2num(dt))
		#
		# other date-handling....
	#
	start_date = date_range[0]
	end_date   = date_range[1]
	#
	if incat is None or (hasattr(incat, '__len__') and len(incat)==0):
		n_tries = 0
		try_wait = 5
		max_tries=10
		while n_tries<max_tries:
			try:
				incat = atp.catfromANSS(lon=lons, lat=lats, minMag=mc, dates0=[start_date, end_date], Nmax=None, fout=None, rec_array=True)
				n_tries = max_tries+1
			except:
				# sometimes the network blips. {wait try_wait} seconds and try again.
				print('Network is down, or something. waiting {} sec to try again.)'.format(try_wait))
				n_tries+=1
				time.sleep(try_wait)
		# catalog_range:
	# catalog_range parameter, for MPP applications (usually None).
	if catalog_range is None or len(catalog_range)==0:
		catalog_range = [0,len(incat)]
	#
	#print('**debug: catalog_processing range: {}/{}'.format(catalog_range, len(incat)))
	#
	# first, make a spatial index of the catalog:
	anss_index = index.Index()
	[anss_index.insert(j, (rw['lon'], rw['lat'], rw['lon'], rw['lat'])) for j,rw in enumerate(incat)]
	#
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
	#
	output_catalog = []
	#
	#for k_cat, rw in enumerate(incat):
	for k_cat, rw in enumerate(incat[catalog_range[0]:catalog_range[1]]):
		#
		# right now, we don't use k_cat for anything, but note it's the index of the sub-catalog that we're calculating.
		# we actually pull the entire catalog, so we can do orientation fitting.
		#
		# tentatively, let's just make a catalog of this and calculate the rupture length, etc.
		# this way, we can specify (different) rupture length, etc. later.
		#
		#lNas = (2.0*(2.0+D))*math.log10(1.+D/2.) + (D/(2.+D))*(self.mag - self.dm - self.mc)
		#lNprimemax = lNas - lrupture + math.log10(2.0)	# log(2) from: 0 <r < L/2 -- assume mainshock is centered.   (aka, N/(L/2))
		#Nas=10.0**lNas
		#Nprimemax = 10.0**lNprimemax
		#r0 = Nom*(q-1.0)/Nprimemax
		#chi = (r0**(1.0-q))/(Nom*(q-1.0))
		#
		# note: look (one day...) for a compiled external (not in this class) function that makes these calcuations.
		# 1) it will be portable, 2) by compiling, say, with numba.jit, we might be able to get some speed. we might even
		# be able to port it to a GPU module...
		mag = rw['mag']
		L_r = 10.0**(.5*mag - d_lambda)
		dt_r = 10.0**(.5*mag - d_tau)		# approximate duration of rupture
		lN_om = b*(mag - dmstar - mc)
		N_om = 10.0**lN_om
		#
		# some short-hand:
		D = D_fract
		#
		# self-similar formulation from yoder et al. 2014/15:
		# start with the (log of the) number of earthquakes/aftershocks in the rupture area:
		# (in yoder et al. 2015, this is "N_as", or "number of aftershocks (inside rupture area)".
		#
		lN_chi = (2.0/(2.0 + D))*numpy.log10(1.0 + D/2.) + D*(mag - dmstar - mc)/(2.+D)
		N_chi  = 10.**lN_chi
		#
		# mean linear density and spatial Omori parameters:
		linear_density = 2.0*N_chi/L_r		# (linear density over rupture area, or equivalently the maximum spatial (linear) density of aftershocks.
		#
		## this version is in use with BASScast:
		# ... and these appear to result in different initial rate-densities, though the maps are qualitatively similar beyond that.
		# not sure at this point which version was actually published...
		# ... looks like this is not the publised version, and that it's not the best approach to compting these parameters.
		# go with the apporoach from Yoder et al. (2015).
		#l_r0 = mag*((6.0+D)/(4.0+2*D)) - (2.0/(2.0+D))*(dmstar+ + mc - math.log10((2.0+D)/2.0)) + math.log10(q-1.0) - d_lambda - math.log10(2.0)
		# r_0 = 10.**l_r0
		#
		# from yoder et al. 2015 and sort of from BASScast:
		r_0 = N_om*(q-1.0)/linear_density		
		#
		## let's try this formulation; sort out details later.
		#lr_0 = mag*((6.0+D)/(4.0+2*D)) - (2.0/(2.0+D))*(dmstar + mc - math.log10((2.0+D)/2.0)) + math.log10(q-1.0) - d_lambda - math.log10(2.0)
		#r_0 = 10.**lr_0
		#
		chi = (r_0**(1.-q))/(N_om*(q-1.0))
		#radialDens = (q-1.0)*(r0ssim**(q-1.0))*(r0ssim + rprime)**(-q)
		#
		# temporal Omori parameters:
		# (something is not correct here; getting negative rate_max (i think))
		# ... but isn't this supposed to be the exponent (log)?
		rate_max = 10.**(l_15_factor + d_tau - mag/6. - (dm_tau + mc)/3.)   
		t_0 = N_om*(p-1.)/rate_max		# note, however, that we can really use just about any value for t_0, so long as we are consistent with tau.
		# something is wrong with this tau calc; we're getting t_0<0. needs fixin...
		tau = (t_0**(1.-p))/(N_om*(p-1.))
		if numpy.isnan(tau): print("nan tau: ", t_0, p, N_om, p)
		#
		# now, guess the earthquake's orientation based on local seismicity. use earthquakes within some distance based on magnitude.
		# use a PCA type method or some sort of line-fit. there should be some good PCA libraries, otherwise it's easy to code.
		# print('**debug: fit_factor = ', fit_factor)
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
		#		
		if len(included_indices)>N_min_ellip_fit:
			# since we're just getting the eigen- vals, vecs, it might be faster to just do the covariance calculation in place here.
			# note that in general, whether or not it's faster to re-calc cov(X) in place or save a copy will depend on the size of the array,
			# but since we're always doing dim(cov)=(2,2), it should be faster to make a copy.
			# cov_pca = numpy.cov([[incat['lon'][j], incat['lat'][j]] for j in included_indices])
			# eig_vals, eig_vecs = numpy.eigh(numpy.dot(cov_pca.transpose(), cov_pca ))
			#
			#print('***debug: gettng pca: ', incat['lon'], incat['lat'])
			# TODO: not sure this is happening correctly in some cases. maybe we just code the pca from scratch? otherwise, check get_pca()
			# this cov, then eig_vals/vecs works, so there's a problem with get_pca()... also, using numpy.cov() will be faster
			# than our all-python method, so let's keep it for a while.
			# note: numpy.eig() returns eigen_vectors as *columns* in the eigen_vector matrix. so the j'th e-vectors is eigen_vecs.T[j]
			cov = numpy.cov(incat['lon'], incat['lat'])
			eig_vals, eig_vecs = numpy.linalg.eig(cov)
			#eig_vals, eig_vecs = get_pca(cat=[[incat['lon'][j], incat['lat'][j]] for j in included_indices], center_lat=rw['lat'], center_lon=rw['lon'], xy_transform=True)
		else:
			eig_vals = [1.0, 1.0]
			eig_vecs = [[1.0, 0.], [0., 1.0]]
			# or: 
			#eig_vecs = numpy.identity(2)
		#
		#
		output_catalog += [ list(rw) + [L_r, r_0, chi, dt_r, t_0, tau, dmstar, p, q] + [eig_vals, eig_vecs, len(included_indices)] ]
		
		# (and strike_theta (or at least something like it -- setting aside conventions) = atan(eig_vecs[0][1]/eig_vecs[0][1])
		#
	if len(output_catalog)==0:
		# if we give parameters that return no events or if we mpp a catalog with only a few events, we get an
		# empty catalog, and the rec_array conversion will pitch a fit. so let's kick back a None and handle it on the return.
		return None
	#	
	if do_recarray:
		output_catalog = numpy.core.records.fromarrays(zip(*output_catalog), dtype=my_dtype)
		output_catalog.meta_parameters = meta_parameters
	else:
		# 2016_12_10 yoder: subcats don't seem to be piping back correctly. might be that we need to change our call signature
		# from apply_async(function, args, kwargs) to apply_async(function(*args, **kwargs)) for mpp calls... and recarrays nor porting back?
		# so for list-returns, let's return the first row/two rows as the dtype parameters:
		#print('***debug: returning list with dtype appended...')
		output_catalog = [my_dtype] + output_catalog
	
	return output_catalog
#
def elliptical_transform_test(theta = 0., x0=0., y0=0., n_quakes=100, m=6.0, max_ratio=2.0, fignum=0):
	# script to test elliptical transforms...
	# keep this simple. we'll line up n_quakes along a line with angle theta, length L.
	# then, run them through catalog pre-processing then etas and see what we get... halos should line up along the line of events.
	#
	if abs(theta)>math.pi*2.1: theta*=deg2rad
	dx = math.cos(theta)
	dy = math.sin(theta)
	Rx = random.Random()
	Ry = random.Random()
	R_l = random.Random()
	#
	L = (10.**(.5*m-1.76))/deg2km
	#
	# test catalog should have a fromat like:
	# test_cat_like = incat = atp.catfromANSS(lon=lons, lat=lats, minMag=mc, dates0=date_range, Nmax=None, fout=None, rec_array=True)
	#
	# which returns rows like:
	# datetime.datetime(2005, 1, 1, 22, 20, 16, 950000), 30.471, 141.204, 4.1, 40.0, 731947.9307517362)
	#
	# then process the catalog like:
	# processed_catalog = make_ETAS_catalog_mpp(incat=None, lats=[32., 38.], lons=[-117., -114.], mc=2.5, date_range=['1990-1-1', None], D_fract=1.5, d_lambda=1.76, d_tau = 2.28, fit_factor=1.5, p=1.1, q=1.5, dmstar=1.0, b1=1.0,b2=1.5, do_recarray=True, n_cpus=None)
	#
	x1 = x0 - L*math.cos(theta)
	y1 = y0 - L*math.sin(theta)
	#
	t_now = dtm.datetime.now(pytz.timezone('UTC'))
	catalog = []
	#
	for j in range(n_quakes):
		#x,y = (x1+
		#
		#x = x1 + 2.*L*dx*Rx.random()
		#y = y1 + 2.*L*dy*Ry.random()
		L_prime = L*2.0*R_l.random()
		x = x1 + L_prime*math.cos(theta)
		y = y1 + L_prime*math.sin(theta)
		#pass
		#
		dt = t_now - dtm.timedelta(days=10)
		#
		catalog += [[dt, y, x, m-2.0, 42., mpd.date2num(dt)]]
		#
	#
	dt = dtm.datetime.now(pytz.timezone('UTC'))
	# grab a dummy catalog --  a stupid way to get a dtype array.
	ct = atp.catfromANSS(dates0=[dtm.datetime.now(pytz.timezone('UTC'))-dtm.timedelta(days=10), dtm.datetime.now(pytz.timezone('UTC'))], lat=[31., 33.], lon=[-117., -115.])
	catalog += [[dt, y0, x0, m, 21., mpd.date2num(dt)]]
	catalog = numpy.core.records.fromarrays(zip(*catalog), dtype=ct.dtype)
	#
	print("cat_len: ", len(catalog))
	#
	ct = make_ETAS_catalog(incat=catalog, lats=[min(catalog['lat']), max(catalog['lat'])], lons=[min(catalog['lon']), max(catalog['lon'])], mc=2.5, date_range=['1990-1-1', None], D_fract=1.5, d_lambda=1.76, d_tau = 2.28, fit_factor=1.5, p=1.1, q=1.5, dmstar=1.0, b1=1.0,b2=1.5, do_recarray=True, catalog_range=None)
	print("len new catalog: ", len(ct))
	#
	#ct = make_ETAS_catalog(incat=catalog, lats=None, lons=None, mc=2.5, date_range=['1990-1-1', None], D_fract=1.5, d_lambda=1.76, d_tau = 2.28, fit_factor=1.5, p=1.1, q=1.5, dmstar=1.0, b1=1.0,b2=1.5, do_recarray=True, catalog_range=None)
	#
	#etas = ETAS_rtree(catalog=ct, t_0 = t_now-dtm.timedelta(days=2), t_now = t_now+dtm.timedelta(days=1), cat_len=None, lats=[min(catalog['lat']), max(catalog['lat'])], lons=[min(catalog['lon']), max(catalog['lon'])], transform_ratio_max=15.)
	etas = ETAS_rtree(catalog=ct, t_0 = t_now-dtm.timedelta(days=12), t_now = t_now+dtm.timedelta(days=15), cat_len=None, lats=[y0-2.0*L, y0+2.0*L], lons=[x0-2.0*L, x0+2.0*L], transform_ratio_max=max_ratio)
	#
	#
	plt.figure(fignum)
	plt.clf()
	ax1=plt.subplot('121')
	ax2=plt.subplot('122')
	etas.make_etas_contour_map(fignum=fignum, ax=ax1)
	#plt.plot(*zip(*[[rw[2], rw[1]] for rw in catalog]), marker='o', ls='')
	#plt.plot(*zip(*[[rw[2], rw[1]] for rw in catalog[-1:]]), marker='*', color='r', ms=18., ls='')
	x,y = etas.cm(catalog['lon'], catalog['lat'])
	plt.plot(x,y, marker='o', ls='')
	x,y = etas.cm(catalog['lon'][-1], catalog['lat'][-1])
	plt.plot([x],[y], marker='*', ms=18., ls='')
	#
	#plt.figure(1)
	#plt.clf()
	ms = Earthquake(ct[-1])
	print('*********************\n*******************')
	print('evecs: ', ms.e_vals, ms.e_vecs)
	# TODO: do this properly as a matrix; remember that eig() returns eigen_vectors as columns.
	# but remember also, that it does not really matter because 1) vector direction is not important, 2) if we do the rotation wrong,
	# we fix it with the a/b ratio (vs b/a ratio). really, we need the direction and relative ratio of eigen-values.
	# ... but this is sort of a silly way to express this linear algebra:
	#x_prime = numpy.dot(ms.e_vecs.transpose(), [1.0,0.])
	#y_prime = numpy.dot(ms.e_vecs.transpose(), [0., 1.])
	# this is better.
	x_prime, y_prime = numpy.dot([[1.0,0.], [0., 1.]], ms.e_vecs)
	#
	#  so we should be able to do this, but check to see we have not messed up the parity somewhere.
	#x_prime, y_prime = numpy.dot(numpy.identity(2), ms.e_vecs)
	#
	print('xprime: ', x_prime, y_prime)
	#print('xprime2: ', x_prime_2, y_prime_2)
	ax2.plot(*zip([0.,0.],x_prime), ls='-', lw=2., marker='s')		# ... so we want to rotate with the transpose matrix, then adjust basis vector lengths...
	plt.plot(*zip([0.,0.],y_prime), ls='-', lw=2., marker='s')		# ... so we want to rotate with the transpose matrix, then adjust basis vector lengths...
	#
	#E = [[min(ms.e_vals[1], 2.0),0.], [0., min(ms.e_vals[0], 2.0)]]
	e1 = min(2.0, max(ms.e_vals[1], .5))
	e2 = min(2.0, max(ms.e_vals[0], .5))
	#
	E = [[e2,e1], [e2, e1]]
	#x_pp = numpy.dot(E,x_prime)
	#y_pp = numpy.dot(E,y_prime)
	x_pp = e1*x_prime
	y_pp = e2*y_prime
	# right... so we don't want to dot E dot x_prime, we want to just to scalar multiply the prime vectors. maybe if we reverse the order so it's like R*E*v ???,
	# but for now, let's just do the rotation and multiply the new basis vectors (prime coordinates) by the corresponding eigen-vectors.
	#e_vals_n = numpy.linalg.norm(self.e_vals) # ... and enforce a maximum aspect ratio of between 2 and 5 or so...
	#
	ax2.plot(*zip([0.,0.],x_pp), ls='--', lw=2., marker='D')		# ... so we want to rotate with the transpose matrix, then adjust basis vector lengths...
	ax2.plot(*zip([0.,0.],y_pp), ls='--', lw=2., marker='D')		# ... so we want to rotate with the transpose matrix, then adjust basis vector lengths...
	#
	print('x_prime dot y_prime: ', numpy.dot(x_prime, y_prime), numpy.dot(y_prime, x_prime), x_pp, y_pp, numpy.dot(x_pp, y_pp))
	
	ax2.set_xlim(-1.1, 1.1)
	ax2.set_ylim(-1.1, 1.1)
	#
	x0 = .4
	ax1.set_xlim(-x0, x0)
	ax1.set_ylim(-x0, x0)
	
	
#
# helper scripts
# can we numba.jit compile this? it looks like... no. we can't probably because of the compiled numpy.linalg bits.
# we might try useing scipy linear algebra; rumor has it that scipy is better compiled.
#@numba.jit
def get_pca(cat=[], center_lat=None, center_lon=None, xy_transform=True):
	# get pca (principal component analysis) for the input catalog. 
	# ... but recognizing that this notaiton is a little sloppy, since the specifics of the pca() transform are taylored to this
	# specific geospatial application.
	#
	# if center_lat/lon != None, then subtract these values as
	# the 'mean' (standard PCA). otherwise, subtract the actual mean.
	# if xy_transform, then transform from lat/lon to x,y
	# we really need to find a good PCA library or write an extension for this. maybe a good D project?
	# ... or maybe since most of the work gets done in numpy, it's ok speed-wise, but it would still make a nice D project.
	#
	if len(cat)<2:
		return numpy.array([1.,1.]), numpy.identity(2)
	#
	xy_factor=1.0
	if xy_transform:
		xy_factor = deg2km
	#
	if center_lat is None: center_lat = numpy.mean([x[1] for x in cat])*xy_factor
	if center_lon is None:
		center_lon = numpy.mean([x[0]*(math.cos(x[1]*deg2rad) if xy_transform else 1.0) for x in cat])*xy_factor
	#
	# subtract off center:
	cat_prime = [[rw[0]*xy_factor*(math.cos(rw[1]*deg2rad) if xy_transform else 1.0)-center_lon, rw[1]*xy_factor-center_lat] for rw in cat]
	#
	# now, get eig_vals, eig_vectors:
	# note: we could use numpy.cov() for the covariance matrix, except that it assumes data-CM-centering; here, we can arbitrarily center,
	# for example around the mainshock, not the CM of the data points.
	#n_dof=max(1., len(cat_prime)-1)
	cov = numpy.dot(numpy.array(list(zip(*cat_prime))),numpy.array(cat_prime))/max(1., len(cat_prime)-1)
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
		#return_distances['spherical'] = spherical_dist(lon_lat_from, lon_lat_to)
		return_distances['spherical'] = great_circle(reversed(lon_lat_from), reversed(lon_lat_to)).km
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
	# add an Euclidian calc. for diagnostics:
	if 'euclidian' in dist_types:
		dx_e = (lon_lat_to[0]-lon_lat_from[0])*deg2km
		dy_e = (lon_lat_to[1]-lon_lat_from[1])*deg2km
		#if 'cartesian' in dist_types:
		#
		return_distances['euclidian'] = ( dx_e**2. + dy_e**2.)**.5
		#if 'dx_dy' in dist_types:
		#
		return_distances.update({'dx_e':dx_e, 'dy_e':dy_e})
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
	if len(return_distances)==1: return_distances = list(return_distances.items())[0]
	#
	return return_distances
	
@numba.jit
def spherical_dist(lon_lat_from=[0., 0.], lon_lat_to=[0.,0.], Rearth = 6378.1):
	# ... also, we might start using geopy.distances.great_circle() -- or something like that. it's geopy anyway...
	# Geometric spherical distance formula...
	# displacement from inloc...
	# inloc is a vector [lon, lat]
	# return a vector [dLon, dLat] or [r, theta]
	# return distances in km.
	#
	# TODO: this may be responsible for wonky elliptical transforms; at one point i know i'd reversd the phi,lambda variables
	# ... so just double-check it. it would be nice to use this instead of the 'geo' model; also look at geopy for a spherical
	# distance formula... for example:
	# >>> from geopy.distance import great_circle
	# >>> newport_ri = (41.49008, -71.312796)
	# >>> cleveland_oh = (41.499498, -81.695391)
	# >>> print(great_circle(newport_ri, cleveland_oh).miles)
	#
	# also, we need to get the proper spherical angular displacement (from the parallel)
	#
	#Rearth = 6378.1	# km
	#deg2rad=2.0*math.pi/360.
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
	#print ('phif: ', phif)
	#print('lambf: ', lambf)
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
	#x = [max(x[0]/x[1], x[1]/x[0]) for x in e_vals]
	x = [math.sqrt(max(x[0]/x[1]), math.sqrt(x[1]/x[0])) for x in e_vals]
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
	if not hasattr(xyz, 'dtype'):
		# we got a list (or something), not a recarray.
		xyz = numpy.core.records.fromarrays(zip(*xyz)[0:3], dtype=[('x','>f8'), ('y', '>f8'), ('z', '>f8')])
	#
	min_x, max_x = min(xyz['x']), max(xyz['x'])
	min_y, max_y = min(xyz['y']), max(xyz['y'])
	#
	padding = .05
	if n_x is None:
		dx = min([abs(x-xyz['x'][j]) for j,x in enumerate(xyz['x'][1:]) if x-xyz['x'][j]!=0.])
		n_x = (n_x or int(round((max_x-min_x)/dx)))
	if n_y is None:
		dy = min([abs(y-xyz['y'][j]) for j,y in enumerate(xyz['y'][1:]) if x-xyz['y'][j]!=0.])
		n_y = (n_y or int(round((max_y-min_y)/dy)))
	#
	xi = numpy.linspace(min_x-abs(min_x)*padding, max_x + abs(max_x)*padding, n_x)
	yi = numpy.linspace(min_y-abs(min_x)*padding, max_y + abs(max_x)*padding, n_y)
	#
	zi = mpl.mlab.griddata(xyz['x'], xyz['y'], numpy.log10(xyz['z']), xi, yi, interp='nn')
	#
	plt.figure(0)
	plt.clf()
	cs = plt.contourf(xi,yi,zi,15,cmap=plt.cm.rainbow,vmax=abs(zi).max(), vmin=-abs(zi).min())
	plt.colorbar()
#
def griddata_brute_plot(xyz, xrange=None, yrange=None, dx=None, dy=None):
	# do a brute force grid-n-contour of an xyz input.
	# ... in the future, use plt.imshow()
	#
	#
	if not hasattr(xyz, 'dtype'):
		# we got a list (or something), not a recarray.
		xyz = numpy.core.records.fromarrays(zip(*xyz)[0:3], dtype=[('x','>f8'), ('y', '>f8'), ('z', '>f8')])
	#
	min_x, max_x = min(xyz['x']), max(xyz['x'])
	min_y, max_y = min(xyz['y']), max(xyz['y'])
	#
	#
	padding = .05
	if dx is None:
		dx = min([abs(x-xyz['x'][j]) for j,x in enumerate(xyz['x'][1:]) if x-xyz['x'][j]!=0.])
		#n_x = (n_x or int(round((max_x-min_x)/dx))
	if dy is None:
		dy = min([abs(y-xyz['y'][j]) for j,y in enumerate(xyz['y'][1:]) if x-xyz['y'][j]!=0.])
		#n_y = (n_y or int(round((max_y-min_y)/dy))
	#
	index_x = lambda x: int(round((x-min_x)/dx))
	index_y = lambda y: int(round((y-min_y)/dy))
	# now, make a big (empty) 2D array for the data.
	xy = numpy.array([[None for j in numpy.arange(min_x, max_x+dx, dx)] for k in numpy.arange(min_y, max_y+dy, dy)])
	#xy=xy.transpose()
	for rw in xyz:
		#print(round("indexing: ", rw['x'], index_x(rw['x']), rw['y'], index_y(rw['y']))
		xy[index_y(rw['y'])][index_x(rw['x'])]=math.log10(rw['z'])
	#
	plt.figure(0)
	plt.clf()
	plt.contourf(xy, 15)
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
	if catalog is None:
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
	plt.figure(fignum+2)
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

def etas_diagnostic_1(lons=[-118., -114.], lats=[31., 38.], mc=5.0, date_range=[dtm.datetime(2000,1,1, tzinfo=pytz.timezone('UTC')), dtm.datetime.now(pytz.timezone('UTC'))], etas_fit_factor=1.5, gridsize=.05):
	# make_ETAS_catalog(incat=None, lats=lats, lons=lons, mc=mc, date_range=[t_0, t_now], fit_factor=etas_fit_factor)
	# make a test catalog or two. make a real catalog; just keep a couple big earthquakes.
	#
	# status: looks like the r_0 parameter is way too small.
	#
	cat0 = make_ETAS_catalog(incat=None, lats=lats, lons=lons, mc=mc, date_range=date_range, fit_factor=etas_fit_factor)
	l_mags = sorted(cat0['mag'])
	#
	m8 = l_mags[int(.8*len(l_mags))]
	m_max = l_mags[-1]
	#
	rw=cat0[-1]
	#
	rw = [rw for rw in cat0 if rw['mag']==m_max][0]
	eq = rw.copy()
	eq_obj = Earthquake(eq)
	#
	Lr = 10.**(.5*rw['mag'] - 1.76)
	d_lon = 5.*Lr*math.cos(deg2rad*rw['lat'])/deg2km
	d_lat = 5.*Lr/deg2km
	etas = ETAS_rtree(catalog=[rw], lons=[rw['lon']-d_lon, rw['lon']+d_lon], lats=[rw['lat']-d_lat, rw['lat']+d_lat], d_lat=gridsize, d_lon=gridsize)
	etas.calc_etas_contours(contour_fig_file='temp/temp.png', contour_kml_file='temp/temp.kml', fignum=0)
	plt.figure(0)
	plt.title('etas_1')
	#
	# note: check aggregation
	# outcom: aggregation seems to be working (test with two identical earthquakes; we get exactly x2=2*x1 for the first 10 (and other) elements.
	etas2 = ETAS_rtree(catalog=[rw,rw], lons=[rw['lon']-d_lon, rw['lon']+d_lon], lats=[rw['lat']-d_lat, rw['lat']+d_lat], d_lat=gridsize, d_lon=gridsize)
	etas.calc_etas_contours(contour_fig_file='temp/temp.png', contour_kml_file='temp/temp.kml', fignum=5)
	plt.figure(5)
	plt.title('etas_5')
	#
	# so now, 1) get a linear density plot
	# and 2) create a catalog with multiple entries of the same event; see that they summ together correctly.
	# g1=ggp.WGS84.Inverse(lon_lat_from[1], lon_lat_from[0], lon_lat_to[1],lon_lat_to[0])
	# return_distances.update({'geo':g1['s12']/1000., 'azi1':g1['azi1'], 'azi2':g1['azi2']})
	#
	print("eq: ", eq)
	print('r_0:', eq.r_0, eq.chi)
	#dists = [[ggp.WGS84.Inverse(rw['lat'], rw['lon'], X['y'], X['x'])['s12']/1000., X['z']] for X in etas.ETAS_array]
	
	#dists = [[111.2*math.sqrt((eq['lat']-X['y'])**2. + ((eq['lon']-X['x'])*math.cos(rw['lat']*deg2rad))**2.), X['z']] for X in etas.ETAS_array] 
	dists = [[ggp.WGS84.Inverse(rw['lat'], rw['lon'], X['y'], X['x'])['s12']/1000., X['z']] for X in etas.ETAS_array]
	#
	plt.figure(1)
	plt.clf()
	ax=plt.gca()
	ax.set_yscale('log')
	ax.set_xscale('log')
	plt.plot(*zip(*dists), marker='.', ls='', lw=1.5)
	#plt.plot(*[[x[0], math.log10(x[1])] for x in dists], marker='.', ls='-', lw=1.5)
	plt.plot([Lr/2., Lr/2.], [min([rw[1] for rw in dists]), max([rw[1] for rw in dists])], 'o-', lw=2.)
	plt.xlabel('distance $r$')
	plt.ylabel('ETAS rate-density $z$')
	#
	xbin=5
	#X_bins = numpy.zeros(int(max(zip*(dists[0]))-min(zip*(dists)[0])/xbin), 5.)
	#N = int((max([rw[0] for rw in dists])-min([rw[0] for rw in dists]))/xbin)+1
	N = int((max([rw[0] for rw in dists]))/xbin) + 1
	z_bins = numpy.zeros(N)
	for dist, z in dists:
		z_bins[int(dist/xbin)]+=z*gridsize*gridsize*math.cos(deg2rad*numpy.mean(lats))*deg2km*deg2km
	#
	plt.figure(2)
	plt.clf()
	plt.clf()
	ax=plt.gca()
	ax.set_yscale('log')
	ax.set_xscale('log')
	ax.plot([xbin*j for j,b in enumerate(z_bins)], z_bins, '.-')
	#
	plt.figure(3)
	plt.clf()
	ax=plt.gca()
	ax.set_yscale('log')
	ax.set_xscale('log')
	f_omori = lambda r,r_0,chi,q: 1.0/(chi*(r_0 + r)**q)
	X = numpy.arange(0., eq.L_r*5., .01)
	ax.plot(X, [f_omori(x, eq.r_0, eq.chi, eq.q) for x in X], '.-')
	#
	'''
	plt.figure(4)
	plt.clf()
	rmap = []
	
	Y = list(numpy.arange(lats[0], lats[1], gridsize))
	X = list(numpy.arange(lons[0], lons[1], gridsize))
	
	n_lats = len(Y)
	n_lons = len(X)
	
	for lat, lon in itertools.product(Y,X):
		et=eq_obj.elliptical_transform(lat=lat, lon=lon)
		rmap += [-math.log10(et['R_prime']+75.)]
		#
	#
	rmap = numpy.array(rmap)
	print("rmap: ", rmap.size, n_lats, n_lons)
	#rmap.shape=(n_lats, n_lons)
	rmap.shape=(n_lats, n_lons)
	
	plt.contourf(X,Y,rmap, 25)
	'''
	#
	return etas, etas2
#
def etas_diagnostic_r0(lons=[-118., -114.], lats=[31., 38.], mc=5.0, date_range=[dtm.datetime(2000,1,1, tzinfo=pytz.timezone('UTC')), dtm.datetime.now(pytz.timezone('UTC'))], etas_fit_factor=1.5, gridsize=.05, D_fract=1.5, b=1.0):
	'''
	# another diagnostic:
	# start by checking out the r_0(m) scaling:
	
	'''
	cat0 = make_ETAS_catalog(incat=None, lats=lats, lons=lons, mc=mc, date_range=date_range, fit_factor=etas_fit_factor, D_fract=D_fract, b1=b)
	#l_mags = sorted(cat0['mag'])
	f_lin = lambda x,a,b: a + b*x
	#
	lstsq, cov = spo.curve_fit(f_lin, xdata=numpy.array(cat0['mag']), ydata=numpy.log10(cat0['r_0']), p0=numpy.array([0.,0.]))
	#print("lstsq: ", lstsq)
	#print('a=%f, b=%f' % tuple(lstsq.tolist()))
	#
	plt.figure(0)
	plt.clf()
	ax = plt.gca()
	ax.set_yscale('log')
	ax.set_xscale('linear')
	#
	ax.plot(cat0['mag'], cat0['r_0'], '.')
	#
	X = [min(cat0['mag']), max(cat0['mag'])]
	ax.plot(X, [10**(f_lin(X[0], *lstsq)), 10**(f_lin(X[1], *lstsq))], '-', label='a=%f, b=%f' % tuple(lstsq))
	plt.legend(loc=0, numpoints=1)
	ax.set_xlabel('magnitude $m$')
	ax.set_ylabel('$r_0$')
	#

def ellipse_test(fignume=0, N=1000):
	plt.figure(0)
	plt.clf()
	
	E = Ellipse(a=10., b=5., theta=45.)
	xy = E.poly()
	plt.plot(*zip(*xy), marker='.', ls='-', color='b')
	gs=.1
	xy_prime = [[float(gs*int(x/gs)), float(gs*int(y/gs))] for x,y in xy]
	plt.plot(*zip(*xy_prime), marker='x', ls='-', color='g') 
	#
	return E
#
def get_eq_properties(m=None, mc=None, b=1.0, d_lambda=1.76, d_tau=2.3, D=1.5, dm=1.0, p=1.1, q=1.5, dm_tau=0):
	#mag = rw['mag']
	r_vals = locals().copy()
	#
	# bridging some new and old syntax:
	mag=m
	dmstar=dm
	l_15_factor = (2./3.)*numpy.log10(1.5)
	#
	L_r = 10.0**(.5*mag - d_lambda)
	dt_r = 10.0**(.5*mag - d_tau)		# approximate duration of rupture
	lN_om = b*(mag - dmstar - mc)
	N_om = 10.0**lN_om
	#
	# self-similar formulation from yoder et al. 2014/15:
	# start with the (log of the) number of earthquakes/aftershocks in the rupture area:
	# (in yoder et al. 2015, this is "N_as", or "number of aftershocks (inside rupture area)".
	#
	lN_chi = (2.0/(2.0 + D))*numpy.log10(1.0 + D/2.) + D*(mag - dmstar - mc)/(2.+D)
	N_chi  = 10.**lN_chi
	#
	# mean linear density and spatial Omori parameters:
	linear_density = 2.0*N_chi/L_r		# (linear density over rupture area, or equivalently the maximum spatial (linear) density of aftershocks.
	#
	## this version is in use with BASScast:
	# ... and these appear to result in different initial rate-densities, though the maps are qualitatively similar beyond that.
	# not sure at this point which version was actually published.
	#l_r0 = mag*((6.0+D)/(4.0+2*D)) - (2.0/(2.0+D))*(dmstar+ + mc - math.log10((2.0+D)/2.0)) + math.log10(q-1.0) - d_lambda - math.log10(2.0)
	# r_0 = 10.**l_r0
	#
	# from yoder et al. 2015 and sort of from BASScast:
	r_0 = N_om*(q-1.0)/linear_density		
	#
	## let's try this formulation; sort out details later.
	#lr_0 = mag*((6.0+D)/(4.0+2*D)) - (2.0/(2.0+D))*(dmstar + mc - math.log10((2.0+D)/2.0)) + math.log10(q-1.0) - d_lambda - math.log10(2.0)
	#r_0 = 10.**lr_0
	#
	chi = (r_0**(1.-q))/(N_om*(q-1.0))
	#radialDens = (q-1.0)*(r0ssim**(q-1.0))*(r0ssim + rprime)**(-q)
	#
	# temporal Omori parameters:
	# (something is not correct here; getting negative rate_max (i think))
	# ... but isn't this supposed to be the exponent (log)?
	rate_max = 10.**(l_15_factor + d_tau - mag/6. - (dm_tau + mc)/3.)
	t_0 = N_om*(p-1.)/rate_max		# note, however, that we can really use just about any value for t_0, so long as we are consistent with tau.
	# something is wrong with this tau calc; we're getting t_0<0. needs fixin...
	tau = (t_0**(1.-p))/(N_om*(p-1.))
	#
	r_vals.update({'tau':tau, 'chi':chi, 't_0':t_0, 'r_0':r_0})
	#
	# initial etas rate density:
	# (1/N_om)*omori_rate*omori_density/circumf
	z0=(10**(dm+mc-m))*1./(tau*(t_0**p)*chi*(r_0**(q+1))*2*numpy.pi)
	r_vals.update({'z0':z0})
	
	return r_vals
#
def omori_time(t_var, tau, t_0, p):
	return 1./(tau * (t_0 + t_var)**p)
#	
if __name__=='__main__':
	# do main stuff...
	pass
else:
	plt.ion()
	


	
		
