#
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
import itertools
import sys
#import scipy.optimize as spo
import os
#from PIL import Image as ipp
import multiprocessing as mpp
#
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.mpl as mpl
#import functools
#
#import shapely.geometry as sgp
#
from mpl_toolkits.basemap import Basemap as Basemap
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from geographiclib.geodesic import Geodesic as ggp
#
#
import ANSStools as atp
#import bindex
import contours2kml
import globalETAS as gep
from eq_params import *
import random
#
#import rtree
#from rtree import index
#

sischuan_prams = {'to_dt':dtm.datetime(2008,6,12, tzinfo=pytz.timezone('UTC')), 'mainshock_dt':dtm.datetime(2008,5,13, tzinfo=pytz.timezone('UTC')), 'lat_center':31.021, 'lon_center':103.367, 'Lr_map_factor':4.0, 'mc':4.0, 'mc_0':None, 'dm_cat':2.0, 'gridsize':.1, 'fnameroot':'etas_auto_sichuan', 'catlen':10.0*365., 'd_lambda':1.76, 'doplot':True}
#sischuan_prams['to_dt'] = dtm.datetime(2014,4,20, tzinfo=pytz.timezone('UTC'))
sischuan_prams['to_dt'] = dtm.datetime.now(pytz.timezone('UTC'))

def nepal_etas_roc():
	# def __init__(self, catalog=None, lats=[32., 36.], lons=[-117., -114.], mc=2.5, mc_etas=None, d_lon=.1, d_lat=.1, bin_lon0=0., bin_lat0=0., etas_range_factor=10.0, etas_range_padding=.25, etas_fit_factor=1.0, t_0=dtm.datetime(1990,1,1, tzinfo=tz_utc), t_now=dtm.datetime.now(tzutc), transform_type='equal_area', transform_ratio_max=5., cat_len=2.*365., calc_etas=True, n_contours=15,**kwargs)
	#
	nepal_etas_fc = get_nepal_etas_fc()
	nepal_etas_test = get_nepal_etas_test()
	#
	# get mainshock:
	ms = nepal_etas_fc.catalog[0]
	for rw in nepal_etas_fc.catalog:
		if rw['mag']>ms['mag']: ms=rw
	#
	z1 = plot_mainshock_and_aftershocks(etas=nepal_etas_fc, m0=6.0, fignum=0)
	#
	z2 = plot_mainshock_and_aftershocks(etas=nepal_etas_test, m0=6.0, fignum=1)
	#
	# now, normalize z1,z2. we can normalize a number of different ways, namely 1) normalize so that sum(z)=1, 2) normailze to max(z).
	# nominally, if we want to compare a big catalog to a small catalog and we don't want to mess around with time dependence, we just normalize to sum(z)=1.
	#
	return nepal_etas_fc, nepal_etas_test

def get_nepal_etas_fc():
	np_prams = {key:nepal_ETAS_prams[key] for key in ['lats', 'lons', 'mc']}
	np_prams.update({'d_lat':0.1, 'd_lon':0.1, 'etas_range_factor':10.0, 'etas_range_padding':.25, 'etas_fit_factor':1.5, 't_0':dtm.datetime(1990,1,1, tzinfo=tz_utc), 't_now':dtm.datetime(2015,5,7,tzinfo=tzutc), 'transform_type':'equal_area', 'transform_ratio_max':2., 'cat_len':5.*365., 'calc_etas':True, 'n_contours':15})
	#nepal_etas_fc = gep.ETAS_rtree(**np_prams)
	#
	#return gep.ETAS_rtree(**np_prams)
	return gep.ETAS_mpp_handler_xyz(**np_prams)

def get_nepal_etas_test(**pram_updates):
	# pram_updates: any earthquake parameters (aka, np_prams) we might want to specify, like "q"...
	#
	# 
	#
	np_prams = {key:nepal_ETAS_prams[key] for key in ['lats', 'lons', 'mc']}
	np_prams.update({'d_lat':0.1, 'd_lon':0.1, 'etas_range_factor':10.0, 'etas_range_padding':.25, 'etas_fit_factor':1.5, 't_0':dtm.datetime(1990,1,1, tzinfo=tz_utc), 't_now':dtm.datetime(2015,5,7,tzinfo=tzutc), 'transform_type':'equal_area', 'transform_ratio_max':2., 'cat_len':5.*365., 'calc_etas':False, 'n_contours':15})
	#
	np_prams_test = np_prams
	np_prams_test.update({'t_now':dtm.datetime(2015,5,21,tzinfo=tzutc), 't_0':dtm.datetime(2015,5,8,tzinfo=tzutc)})
	np_prams_test.update(pram_updates)
	#
	#nepal_etas_test = gep.ETAS_rtree(**np_prams_test)
	#
	# there's a more proper way to do this, probably to create an inherited class of globalETAS and/or Earthquake for a "stationary" output (p=0 for calculating ETAS, but remember NOT for
	# the initial rate calculation (which is done in the catalog construction bits) ). for now, create the ETAS object; then spin through the catalog and set p=0 for all the earthquakes.
	# remember also that etas.catalog is a recarray, not an array of Earthquake objects.
	#
	#etas = gep.ETAS_rtree(**np_prams_test)
	etas = gep.ETAS_mpp_handler_xyz(**np_prams_test)
	for j,eq in enumerate(etas.catalog):
		etas.catalog['p'][j] = 0.0
		# ... and sort of a sloppy way to do this as well...
		for key,val in pram_updates.items(): etas.catalog[key]=val
	#
	etas.make_etas()
	return etas

class Toy_etas(object):
	def __init__(self, etas_in, mainshock={'mag':7.3, 'lon':84.698, 'lat':28.175}):
		# nepal_epi_lon = 84.698
		# nepal_epi_lat = 28.175
		self.__dict__.update(etas_in.__dict__)
		self.__dict__.update(locals())
		self.lattice_sites = etas_in.lattice_sites
		#
		# now, replace all the ETAS z values with 1/r to epicenter... let's do 1/R+L_r/2, so we don't get singularities (it won't matter since the comaprison will
		# be rank ordered).
		#
	
	def normalize(self):
		norm_factor = numpy.sum(self.ETAS_array['z'])
		self.ETAS_array['z']/=norm_factor
		#
	#
class  Toy_etas_invr(Toy_etas):
	def __init__(self, *args, **kwargs):
		super(Toy_etas_invr,self).__init__(*args, **kwargs)
		self.L_r = 10.**(.5*self.mainshock['mag']-1.76)
		#
		for j,rw in enumerate(self.ETAS_array):
			g1=ggp.WGS84.Inverse(self.mainshock['lat'], self.mainshock['lon'], rw['y'], rw['x'])
			#r_prime = (g1['s12']/1000.) + .5*L_r
			self.ETAS_array['z'][j] = 1./((g1['s12']/1000.) + .5*self.L_r)
		#
		self.normalize()
#
class Toy_etas_random(Toy_etas):
	def __init__(self, *args, **kwargs):
		super(Toy_etas_random, self).__init__(None,*args, **kwargs)
		R=random.Random()
		for j,rw in enumerate(self.ETAS_array):
			self.ETAS_array['z'][j] = R.random()
		#
		self.normalize()
#class Toy_etas_fromxyz(Toy_etas):
class Toy_etas_fromxyz(object):
	# maybe inherit from Toy_etas later...
	# (still needs testing and development)
	def __init__(self, fname='global_map_20151129.xyz', *args, **kwargs):
		#super(Toy_etas_fromxyz, self).__init__(*args, **kwargs)
		#self.ETAS_array=[]
		with open(fname,'r') as f:
			self.ETAS_array = [[float(x) for x in rw] for rw in f if rw[0]!='#']
		#
		self.ETAS_array = numpy.core.records.fromarrays(zip(*self.ETAS_array), names=('x', 'y', 'z'), formats = ('double', 'double', 'double'))
		#
		self.lons = [min(self.ETAS_array['x']), max(self.ETAS_array['x'])]
		self.lats = [min(self.ETAS_array['y']), max(self.ETAS_array['y'])]
		#
		#
	#
		
	#
'''	
def toy_etas(ETAS_array=None, lats=None, lons=None, l_lon=.1, d_lat=.1, epicenter=None, mc=None):
	# make a toy object that contains all the necessary bits to look like an etas object.
	#
	# ... but actually, maybe a better way to do this is to just clone (most of) another etas objec, like:
	# obj.update(etas_template.__dict__)
	if ETAS_array==None:
		ETAS_array = [[x,y] for x,y in itertools.product(numpy.arange(lats[0], lats[1]+d_lat, d_lat), numpy.arange(lons[0], lons[1]+d_lon, d_lon))]
	#
	epicen_lat = 28.147
	epicen_lon = 84.708
	#
	for j,rw in enumerate(ETAS_array):
		if len(rw)<3:
			ETAS_array[j]+=[0]
			g1=ggp.WGS84.Inverse(epicen_lat, epicen_lon, rw[1], rw[0])
			ETAS_array[j][2] = 1.0/(g1['s12']/1000.)
			#
	my_etas = object()
	my_etas.lattice_xy = ETAS_array
	my_etas.lons=[min([rw[0] for rw in ETAS_array]), max([rw[0] for rw in ETAS_array])]
	my_etas.lats=[min([rw[1] for rw in ETAS_array]), max([rw[1] for rw in ETAS_array])]
	#
'''	

def roc_normalses(etas_fc, test_catalog=None, to_dt=None, cat_len=120., mc_rocs=[4.0, 5.0, 6.0, 7.0], fignum=1, do_clf=True, roc_ls='-'):
	#
	# make a set of ROCs.
	plt.figure(fignum)
	if do_clf: plt.clf()
	ax=plt.gca()
	FHs=[]
	#
	for mc in mc_rocs:
		# ... we should probalby modify roc_normal() so we can pass a catalog (for speed optimization), but we'll probably only run this a few times.
		print('roc for %f', mc)
		FH = roc_normal(etas_fc, test_catalog=None, to_dt=None, cat_len=120., mc_roc=mc, fignum=0)
		ax.plot(*zip(*FH), marker='', ls=roc_ls, lw=2.5, alpha=.8, label='$m_c=%.2f$' % mc)
		FHs += [[mc,FH]]
		#
	#
	ax.plot(range(2), range(2), ls='--', marker='', lw=2.75, alpha=.7, zorder=1)
	plt.figure(fignum)
	ax.legend(loc=0, numpoints=1)
	ax.set_ylim([-.1,1.15])
	ax.set_xlim([-.1,1.15])
	ax.set_title('ROC Analysis', size=18)
	ax.set_xlabel('False Alarm Rate $F$', size=18)
	ax.set_ylabel('Hit Rate $H$', size=18)
	plt.draw()
	#
	return FHs

#
class ROC_base(object):
	# a base class object for MPP ROC stuff. we'll inherit this for both the worker and handerl mpp classes.
	#
	def __init__(self, fc_xyz, test_catalog=None, from_dt=None, to_dt=None, dx=None, dy=None, cat_len=120., mc=5.0):
		#
		if isinstance(fc_xyz, str):
			# if we're given a filename...
			with open(fc_xyz, 'r') as froc:
				fc_xyz= [[float(x) for x in rw.split()] for rw in froc if rw[0] not in('#', ' ', '\t', '\n')]
		#
		if to_dt == None:
			to_dt = dtm.datetime.now(pytz.timezone('UTC'))
		if from_dt == None or from_dt>to_dt:
			from_dt = to_dt - dtm.timedelta(days=cat_len)
	
			#
		#
		#return fc_xyz
		if not hasattr(fc_xyz, 'dtype'):
			fc_xyz = numpy.core.records.fromarrays(zip(*fc_xyz), names=('x','y','z'), formats=['>f8', '>f8', '>f8'])
		#
		lats = [min(fc_xyz['y']), max(fc_xyz['y'])]
		lons = [min(fc_xyz['x']), max(fc_xyz['x'])]
		#
		#mc   = mc_roc
		X_set = list(set(fc_xyz['x']))
		Y_set = list(set(fc_xyz['y'])) 
		d_lon = (dx or abs(X_set[1] - X_set[0]))
		d_lat = (dy or abs(Y_set[1] - Y_set[0]))
		nx = len(X_set)
		ny = len(Y_set)
		#
		# (for this application, we can also just get nasty and to a loop-loop with geodetic distancing).
		# (maybe ought to set this as a regular class function?)
		get_site = lambda x,y: int(round((x-lons[0]+.5*d_lon)/d_lon)) + int(round((y-lats[0]+.5*d_lat)/d_lat))*nx
		#
		
		print("get cataog: ", lons, lats, mc, from_dt, to_dt)
		if test_catalog==None: test_catalog = atp.catfromANSS(lon=lons, lat=lats, minMag=mc, dates0=[from_dt, to_dt])
		print("catlen: ", len(test_catalog))
		#
		# note: this will reduce the number of z values a bit, but what we really need to do is bin them into N bins... or maybe
		# we do explicitly want to remove each site, one at a time to get a proper ROC...
		#Zs = list(set(sorted(list(fc_xyz['z'].copy()))))
		Zs = sorted(list(fc_xyz['z'].copy()))
		ROCs = [[0,0,0,0]]
		# nominally, we'd make a copy of the whole catalog, but for memory conservation, just index the lsit.
		#eq_ks = [[get_site(eq['lon'], eq['lat']), eq] for eq in test_catalog]
		eq_site_indices = [get_site(eq['lon'], eq['lat']) for eq in test_catalog]
		#
		# basic bits and pieces established. now, handler class will allocate sub-sets of Zs to workers; workers will aggregate roc.
		#
	#
	def calc_ROCs(self):
		#for j_z, z0 in enumerate(Zs['z']):
		ROCs = [[0,0,0,0]]
		for j_z, z0 in enumerate(Zs):
			# z0 is the threshold z for predicted=True/False
			#
			for j_eq, eq in enumerate(test_catalog):
				# TODO: and parse it up for MPP.
				#k = get_site(eq['lon'], eq['lat'])
				k = eq_site_indices[j_eq]
				#print('site: ', k)
				z_val = fc_xyz['z'][k]
				#
				if z_val>=z0:
					# predicted!
					ROCs[-1][0]+=1
					#
					# ... and subtract from falsies; in the end, we'll assume all sites>z were false alarms:
					ROCs[-1][1]-=1
				else:
					# missed it
					ROCs[-1][2]+=1
					ROCs[-1][3]-=1		# an earthquake occurred in this site. it did not correctly predict non-occurrence.
				#
			n_gt = float(len([z for z in fc_xyz['z'] if z>=z0]))
			n_lt = float(len([z for z in fc_xyz['z'] if z<z0]))
			#
			ROCs[-1][1]+=n_gt
			ROCs[-1][3]+=n_lt
			#
			ROCs += [[0,0,0,0]]
		#
		self.ROCs=ROCs
		#
		# return ROCs
	#
	def calc_HF(self):
		self.calc_ROCs()
		#
		#Hs2=[]
		#Fs2=[]
		# note: this migh make a nice practice problem for mpp.Array() ....
		HF = [[float(rw[0])/float(len(test_catalog)), roc[1]/float(len(fc_xyz))] for rw in ROCs[:-1]]
		#
		self.HF = HF
		#for roc in ROCs[:-1]:
		#	roc=[float(x) for x in roc]
		#	#
		#	Hs += [roc[0]/float(len(test_catalog))]
		#	Fs += [roc[1]/float(len(fc_xyz))]
	#
	def plot_HF(self, fignum=0):
		plt.figure(fignum)
		if do_clf: plt.clf()
		#plt.plot(Fs,Hs, '-', label='ROC_approx.', lw=2., alpha=.8)
		plt.plot(*zip(*self.HF), ls='-', label='ROC_approx.', lw=2., alpha=.8)
		#plt.plot(Fs2, Hs2, '-', label='ROC', lw=2., alpha=.8)
		plt.plot(range(2), range(2), 'r--', lw=2.5, alpha=.6)


#
class ROC_mpp_handler(ROC_base):
	def __init__(self, fc_xyz, test_catalog=None, from_dt=None, to_dt=None, dx=None, dy=None, cat_len=120., mc=5.0, n_procs=None):
		super(ROC_mpp_handler,self).__init__(**{key:val for key,val in locals().items() if key!='self'})
		#
		if n_procs==None: n_procs = min(1, mpp.cpu_count()-1)
		#
		
	#
class ROC_worker(ROC_base):
	# worker class for ROC calculations.
	def __init__(self, fc_xyz, test_catalog=None, from_dt=None, to_dt=None, dx=None, dy=None, cat_len=120., mc=5.0):
		super(ROC_mpp_handler,self).__init__(**{key:val for key,val in locals().items() if key!='self'})
		pass	


	
#
def roc_normal_from_xyz(fc_xyz, test_catalog=None, from_dt=None, to_dt=None, dx=None, dy=None, cat_len=120., mc=5.0, fignum=0, do_clf=True):
	#
	# roc from an xyz forecast input. eventually, convolve this with roc_normal() which takes etas_fc, an etas type, object as an input.
	# dx, dy are grid-sizes in the x,y direction. if none, we'll figure them out.
	#
	if isinstance(fc_xyz, str):
		# if we're given a filename...
		with open(fc_xyz, 'r') as froc:
			fc_xyz= [[float(x) for x in rw.split()] for rw in froc if rw[0] not in('#', ' ', '\t', '\n')]
	#
	if to_dt == None:
		to_dt = dtm.datetime.now(pytz.timezone('UTC'))
	if from_dt == None or from_dt>to_dt:
		from_dt = to_dt - dtm.timedelta(days=cat_len)
	
		#
	#
	#return fc_xyz
	if not hasattr(fc_xyz, 'dtype'):
		fc_xyz = numpy.core.records.fromarrays(zip(*fc_xyz), names=('x','y','z'), formats=['>f8', '>f8', '>f8'])
	#
	lats = [min(fc_xyz['y']), max(fc_xyz['y'])]
	lons = [min(fc_xyz['x']), max(fc_xyz['x'])]
	#
	#mc   = mc_roc
	X_set = list(set(fc_xyz['x']))
	Y_set = list(set(fc_xyz['y'])) 
	dx = (dx or abs(X_set[1] - X_set[0]))
	dy = (dy or abs(Y_set[1] - Y_set[0]))
	#
		
	print("get cataog: ", lons, lats, mc, from_dt, to_dt)
	if test_catalog==None: test_catalog = atp.catfromANSS(lon=lons, lat=lats, minMag=mc, dates0=[from_dt, to_dt])
	print("catlen: ", len(test_catalog))
	#
	# note: this will reduce the number of z values a bit, but what we really need to do is bin them into N bins... or maybe
	# we do explicitly want to remove each site, one at a time to get a proper ROC...
	#Zs = list(set(sorted(list(fc_xyz['z'].copy()))))
	Zs = sorted(list(fc_xyz['z'].copy()))
	#Zs.sort(order='z')
	#
	d_lat = dy
	d_lon = dx
	nx = len(X_set)
	ny = len(Y_set)
	#
	# (for this application, we can also just get nasty and to a loop-loop with geodetic distancing).
	get_site = lambda x,y: int(round((x-lons[0]+.5*d_lon)/d_lon)) + int(round((y-lats[0]+.5*d_lat)/d_lat))*nx
	#
	'''
	#test:
	print('testing get_site:')
	Rx = random.Random()
	Ry = random.Random()
	for k in range(10):
		x = lons[0] + Rx.random()*(lons[1]-lons[0])
		y = lats[0] + Ry.random()*(lats[1]-lats[0])
		#
		j_lattice = get_site(x,y)
		print("for %f, %f: " % (x,y), etas_fc.ETAS_array[j_lattice])
	#
	'''
	#
	'''
	# hits (observed yes, forecast yes):
	roc_A = [0]
	# falsies (over-predict) (observed no, fc. yes)
	roc_B = [0]
	# misses (observed yes, fc no):
	roc_C = [0]
	# didn't happen (observed no, fc no):
	roc_D = [0]
	'''
	#
	ROCs = [[0,0,0,0]]
	# nominally, we'd make a copy of the whole catalog, but for memory conservation, just index the lsit.
	#eq_ks = [[get_site(eq['lon'], eq['lat']), eq] for eq in test_catalog]
	eq_site_indices = [get_site(eq['lon'], eq['lat']) for eq in test_catalog]
	#print("eZs: ", etas_fc.ETAS_array['z'][0:10], len(etas_fc.ETAS_array['z']))
	#
	# (and here, we really need to figure out how to MPP this).
	#for j_z, z0 in enumerate(Zs['z']):
	for j_z, z0 in enumerate(Zs):
		# z0 is the threshold z for predicted=True/False
		#
		for j_eq, eq in enumerate(test_catalog):
			# TODO: and parse it up for MPP.
			#k = get_site(eq['lon'], eq['lat'])
			k = eq_site_indices[j_eq]
			#print('site: ', k)
			z_val = fc_xyz['z'][k]
			#
			if z_val>=z0:
				# predicted!
				ROCs[-1][0]+=1
				#
				# ... and subtract from falsies; in the end, we'll assume all sites>z were false alarms:
				ROCs[-1][1]-=1
			else:
				# missed it
				ROCs[-1][2]+=1
				ROCs[-1][3]-=1		# an earthquake occurred in this site. it did not correctly predict non-occurrence.
			#
		n_gt = float(len([z for z in fc_xyz['z'] if z>=z0]))
		n_lt = float(len([z for z in fc_xyz['z'] if z<z0]))
		#
		ROCs[-1][1]+=n_gt
		ROCs[-1][3]+=n_lt
		#
		ROCs += [[0,0,0,0]]
	#
	
	Hs=[]
	Fs=[]
	
	Hs2=[]
	Fs2=[]
	# note: this migh make a nice practice problem for mpp.Array() ....
	for roc in ROCs[:-1]:
		#try:
		if True:
			roc=[float(x) for x in roc]
			#h=float(len(etas_fc.ETAS_array))
			#f=roc[1]/float(len(etas_fc.ETAS_array))
			#
			Hs += [roc[0]/float(len(test_catalog))]
			Fs += [roc[1]/float(len(fc_xyz))]
			#
			## diagnostic formulation:
			#Hs2 += [roc[0]/(roc[0]+roc[2])]
			#Fs2 += [roc[1]/(roc[1]+roc[3])]
			
		#except:
		#	print('ROC error, probably div/0: ', roc, len(test_catalog), len(etas_fc.ETAS_array), roc[0]/float(len(test_catalog)), roc[1]/float(float(len(etas_fc.ETAS_array))) )
		#	
	#
	if fignum!=None:
		plt.figure(fignum)
		if do_clf: plt.clf()
		plt.plot(Fs,Hs, '-', label='ROC_approx.', lw=2., alpha=.8)
		#plt.plot(Fs2, Hs2, '-', label='ROC', lw=2., alpha=.8)
		plt.plot(range(2), range(2), 'r--', lw=2.5, alpha=.6)
	#
	return list(zip(Fs,Hs))
#	
#
def roc_normal(etas_fc, test_catalog=None, from_dt=None, to_dt=None, cat_len=120., mc_roc=5.0, fignum=0, do_clf=True):
	#
	if from_dt==None:
		from_dt=max([dt.tolist() for dt in etas_fc.catalog['event_date']])
	if to_dt==None:
		to_dt = from_dt + dtm.timedelta(days=cat_len)
		#
	#
	lats = etas_fc.lats
	lons = etas_fc.lons
	mc   = etas_fc.mc
	print("get cataog: ", lons, lats, mc_roc, from_dt, to_dt)
	if test_catalog==None: test_catalog = atp.catfromANSS(lon=lons, lat=lats, minMag=mc_roc, dates0=[from_dt, to_dt])
	print("catlen: ", len(test_catalog))
	#
	Zs = etas_fc.ETAS_array.copy()
	Zs.sort(order='z')
	#
	d_lat = etas_fc.d_lat
	d_lon = etas_fc.d_lon
	nx,ny = etas_fc.lattice_sites.shape	# should be (i think)= int((max(lon)-min(lon))/d_lon)
	#
	#lat0 = min(etas_fc.ETAS_array['y'])
	#lon0 = min(etas_fc.ETAS_array['x'])
	#
	# (for this application, we can also just get nasty and to a loop-loop with geodetic distancing).
	get_site = lambda x,y: int(round((x-lons[0]+.5*d_lon)/d_lon)) + int(round((y-lats[0]+.5*d_lat)/d_lat))*nx
	#get_site = lambda x,y: int(round(x-lons[0])/d_lon) + int((y-lats[0])/d_lat)*nx
	#
	'''
	#test:
	print('testing get_site:')
	Rx = random.Random()
	Ry = random.Random()
	for k in range(10):
		x = lons[0] + Rx.random()*(lons[1]-lons[0])
		y = lats[0] + Ry.random()*(lats[1]-lats[0])
		#
		j_lattice = get_site(x,y)
		print("for %f, %f: " % (x,y), etas_fc.ETAS_array[j_lattice])
	#
	'''
	#
	'''
	# hits (observed yes, forecast yes):
	roc_A = [0]
	# falsies (over-predict) (observed no, fc. yes)
	roc_B = [0]
	# misses (observed yes, fc no):
	roc_C = [0]
	# didn't happen (observed no, fc no):
	roc_D = [0]
	'''
	#
	ROCs = [[0,0,0,0]]
	#print("eZs: ", etas_fc.ETAS_array['z'][0:10], len(etas_fc.ETAS_array['z']))
	#
	for j_z, z0 in enumerate(Zs['z']):
		# z0 is the threshold z for predicted=True/False
		#
		for eq in test_catalog:
			k = get_site(eq['lon'], eq['lat'])
			#print('site: ', k)
			z_val = etas_fc.ETAS_array['z'][k]
			#
			if z_val>=z0:
				# predicted!
				ROCs[-1][0]+=1
				#
				# ... and subtract from falsies; in the end, we'll assume all sites>z were false alarms:
				ROCs[-1][1]-=1
			else:
				# missed it
				ROCs[-1][2]+=1
				ROCs[-1][3]-=1		# an earthquake occurred in this site. it did not correctly predict non-occurrence.
			#
		n_gt = float(len([z for z in etas_fc.ETAS_array['z'] if z>=z0]))
		n_lt = float(len([z for z in etas_fc.ETAS_array['z'] if z<z0]))
		#
		ROCs[-1][1]+=n_gt
		ROCs[-1][3]+=n_lt
		#
		ROCs += [[0,0,0,0]]
	#
	
	Hs=[]
	Fs=[]
	
	Hs2=[]
	Fs2=[]
	# note: this migh make a nice practice problem for mpp.Array() ....
	for roc in ROCs[:-1]:
		#try:
		if True:
			roc=[float(x) for x in roc]
			Hs2 += [roc[0]/(roc[0]+roc[2])]
			#h=float(len(etas_fc.ETAS_array))
			#f=roc[1]/float(len(etas_fc.ETAS_array))
			#
			Hs += [roc[0]/float(len(test_catalog))]
			Fs2 += [roc[1]/(roc[1]+roc[3])]
			Fs += [roc[1]/float(len(etas_fc.ETAS_array))]
			
		#except:
		#	print('ROC error, probably div/0: ', roc, len(test_catalog), len(etas_fc.ETAS_array), roc[0]/float(len(test_catalog)), roc[1]/float(float(len(etas_fc.ETAS_array))) )
		#
	#
	# now, make a heavy-sizde forecast (aka, "there will be N earthquakes", assume within some alpha*L_r.
	
	#
	if fignum!=None:
		plt.figure(fignum)
		if do_clf: plt.clf()
		plt.plot(Fs,Hs, '-', label='ROC_approx.', lw=2., alpha=.8)
		plt.plot(Fs2, Hs2, '-', label='ROC', lw=2., alpha=.8)
		plt.plot(range(2), range(2), 'r--', lw=2.5, alpha=.6)
	#
	return list(zip(Fs,Hs))

def nepal_roc_normal_script(fignum=0):
	# full, one stop shopping script for nepal ROC analysis.
	#
	# first, get nepal ETAS objects:
	etas_fc, etas_test = nepal_etas_roc()
	#
	ROC_n = roc_normalses(etas_fc, test_catalog=None, to_dt=None, cat_len=120., mc_rocs=[4.5, 5.0, 6.0, 7.0], fignum=fignum, do_clf=True)
	#
	# now, make a toy catalog:
	etas_toy = Toy_etas_invr(etas_in=etas_fc, mainshock={'mag':7.3, 'lon':84.698, 'lat':28.175})
	#
	ROC_t = roc_normalses(etas_toy, test_catalog=None, to_dt=None, cat_len=120., mc_rocs=[4.5, 5.0, 6.0, 7.0], fignum=fignum, do_clf=False, roc_ls='--')
	#
	# now, some random catalogs:
	for j in range(25):
		this_etas = Toy_etas_random(etas_in=etas_fc)
		FH = roc_normal(this_etas, fignum=None)
		plt.plot(*zip(*FH), marker='.', ls='', alpha=.6)
		
	
	#
def etas_roc_geospatial_fcset(q_fc_min=1.1, q_fc_max=2.5, q_test_min=1.1, q_test_max=2.5, do_log=True, dq_fc=.1, dq_test=.1, fignum=0):
	#
	#if etas_fc==None: 
	etas_fc=get_nepal_etas_fc()
	#if etas_test==None: 
	etas_test = get_nepal_etas_test() #... and we don't really need to calc eatas here, so later on maybe clean this up.
	
	FH=[]
	#
	for q_fc in numpy.arange(q_fc_min, q_fc_max, dq_fc):
		for j,rw in enumerate(etas_fc.catalog): etas_fc.catalog['q'][j] = q_fc
		etas_fc.make_etas()
		#
		fh = etas_roc_geospatial_set(etas_fc=etas_fc, etas_test=etas_test, q_test_min=q_test_min, q_test_max=q_test_max, do_log=do_log, dq=dq_test, fignum=None)
		for rw in fh: FH += [[q_fc] + rw]
	#
	plt.figure(fignum)
	plt.clf()
	plt.plot([rw[2] for rw in FH], [rw[3] for rw in FH], 'o')
	#
	with open('roc_geospatial.csv', 'w') as fout:
		fout.write('#roc output.\n#q_fc\tq_test\tF\tH\n')
		for rw in FH:
			fout.write('%s\n' % '\t'.join([str(x) for x in rw]))
			#
		#
	#
	return FH
#
def etas_roc_geospatial_set(etas_fc=None, etas_test=None, do_log=True, q_test_min=1.1, q_test_max=2.0, dq=.1, fignum=0):
	# compare ETAS for a bunch of different q. right now this is just the "test" q. we'll probably want to vary the forecast as well, but of course taht will be expensive.
	if etas_fc==None: etas_fc=get_nepal_etas_fc()
	if etas_test==None: etas_test = get_nepal_etas_test() #... and we don't really need to calc eatas here, so later on maybe clean this up.
	FH=[]
	#
	for q in numpy.arange(q_test_min, q_test_max, dq):
		#etas_test = get_nepal_etas_test(q=q)		# and we should sort it out to keep a copy of the catalog, or just re-calc etas with new q...
		for j,rw in enumerate(etas_test.catalog): etas_test.catalog['q'][j]=q
		etas_test.make_etas()
		
		FH += [[q] + list(analyze_etas_roc_geospatial(etas_fc=etas_fc, etas_test=etas_test, do_log=do_log))]
	#
	if fignum != None:
		plt.figure(fignum)
		plt.plot([rw[0] for rw in FH], [rw[2]-rw[1] for rw in FH], 'bo-', lw=2.5)
		plt.xlabel('Test catalog scaling exponent $q_{test}$', size=18)
		plt.ylabel('Skill, $H-F$', size=18)
		#
		plt.figure(fignum+1)
		plt.plot([rw[1] for rw in FH], [rw[2] for rw in FH], 'bo-', lw=2.5)
		plt.plot(range(2), range(2), 'r--', lw=3., alpha=.8)
		plt.plot([0., 0.], [0.,1.], 'k-', lw=2)
		plt.plot([0.,1.], [0.,0.], 'k-', lw=2)
		plt.xlabel('False Alarm Rate $F$', size=18)
		plt.ylabel('Hit Rate $H$', size=18)
		#
	return FH
#		
#
def roc_plots_from_gsroc(FH, fignum=0):
	if len(FH[0])==4: cols = {key:val for key,val in zip(['q_fc', 'q_t', 'F', 'H'], range(4))}
	if len(FH[0])==3: cols = {key:val for key,val in zip(['q_t', 'F', 'H'], range(3))}
	#
	skl = [rw[-1]-rw[-2] for rw in FH]
	#
	fg=plt.figure(fignum)
	plt.clf()
	lft=.05
	btm=.05
	dx=.4
	dy=.4
	ax_ll = fg.add_axes([lft, btm, dx, dy])
	ax_lr = fg.add_axes([lft+dx + .05, btm, dx, dy])
	ax_ul = fg.add_axes([lft, btm+dy + .05, dx, dy])
	ax_ur = fg.add_axes([lft+dx + .05, btm+dy+.05,dx,dy])
	#
	ax_ll.plot([x[cols['F']] for x in FH], [x[cols['H']] for x in FH], 'o', lw=2.)
	ax_ll.plot([FH[0][cols['F']]], [FH[0][cols['H']]], 'rd', lw=2.)
	ax_ll.plot([FH[-1][cols['F']]], [FH[-1][cols['H']]], 'cs', lw=2.)
	ax_ll.plot(range(2),range(2), 'r--', lw=3, alpha=.8)
	ax_ll.set_xlabel('False Alarm Rate $F$')
	ax_ll.set_ylabel('Hit Rate $H$')
	ax_ll.set_title('ROC')
	#
	ax_lr.plot(skl, 'o-')
	ax_lr.plot([0], [skl[0]], 'rd')
	ax_lr.plot([len(skl)-1], [skl[-1]], 'cs')
	ax_lr.set_ylabel('Skill Score, $H-F$')
	ax_lr.set_xlabel('$q_{fc}, q_{test}$')
	ax_lr.set_title('Skill Score')
	#
	if len(FH[0])==4:
		# contour:
		X = sorted(list(set([x[cols['q_fc']] for x in FH])))
		Y = sorted(list(set([x[cols['q_t']] for x in FH])))
		#
		zs = numpy.array([x[cols['H']]-x[cols['F']] for x in FH])
		zs.shape=(len(X), len(Y))
		#
		#ax_ul.contourf(X,Y,zs, 15, alpha=.75)
		ax_ul.contourf(Y,X,zs,15,alpha=.75)
		#
		ax_ul.set_ylabel('$q_{fc}$')
		ax_ul.set_xlabel('$q_{test}$')
	#
	# best skill:
	best_skill = max(skl)
	for rw in FH:
		if rw[-1]-rw[-2]==best_skill: print("best skill: ", rw)
	#besk_skill = [rw for rw in skl if skl[-1]==best_skill][0]
	#print('Best Skill: ', best_skill)
	#	
#
def analyze_etas_roc_geospatial(etas_fc=None, etas_test=None, do_log=True, diagnostic=False):
#def analyze_etas_roc_geospatial(etas_fc=None, etas_test=None, do_log=True):
	#
	if etas_fc   == None: etas_fc   = get_nepal_etas_fc()
	if etas_test == None: etas_test = get_nepal_etas_test()
	#
	etas_fc.make_etas_contour_map(fignum=0)
	etas_test.make_etas_contour_map(fignum=1)
	#
	lon_vals = sorted(list(set(etas_fc.ETAS_array['x'])))
	lat_vals = sorted(list(set(etas_fc.ETAS_array['y'])))
	#
	# we need normalization here...
	z_fc_norm = etas_fc.ETAS_array['z'].copy()
	z_test_norm = etas_test.ETAS_array['z'].copy()
	#
	if do_log:
		z_fc_norm   = numpy.log10(z_fc_norm)
		z_test_norm = numpy.log10(z_test_norm)
	#
	z_fc_norm -= min(z_fc_norm)
	z_test_norm -= min(z_test_norm)
	#
	norm_fc   = sum(z_fc_norm)
	norm_test = sum(z_test_norm)
	#
	z_fc_norm /= norm_fc
	z_test_norm /= norm_test
	#
	z1 = z_fc_norm
	z2 = z_test_norm
	#
	#
	# [z1, z2, diff, h, m, f(predicted, didn't happen)
	#diffs = [[z1, z2, z1-z2, max(z1, z2), -min(z1-z2,0.), max(z1-z2,0.)] for z1,z2 in zip(z_fc_norm, z_test_norm)] 
	# hits: accurately predicted; min(z1,z2)
	# misses: prediction deficite, or excess events: min(z2-z1,0.)
	# falsie: excess prediction: min(z1-z2,0.)
	# then rates: H = hits/sum(z2), F =falsies/sum(z1)
	#diffs = [[z1, z2, z1-z2, min(z1, z2), max(z2-z1,0.), max(z1-z2, 0.)] for z1,z2 in zip(z_fc_norm, z_test_norm)]
	#
	# so we can test this properly, we'll want to move diffs offline to a function call (eventually)...
	diffs = [[z1, z2, z1-z2, min(z1, z2), max(z2-z1,0.), max(z1-z2, 0.)] for z1,z2 in zip(z1, z2)]
	diffs_lbls = ['z_fc', 'z_test', 'z1-z2', 'hits: min(z1,z2)','misses:min(z2-z1,0)', 'falsie: min(z1-z2,0)']
	diffs_lbl_basic = ['z_fc', 'z_test', 'z1-z2', 'hits','misses', 'falsie']
	#

	# to plot contours, we'll want to use the shape from: etas.lattice_sites.shape
	#
	sh1 = etas_fc.lattice_sites.shape
	sh2 = etas_test.lattice_sites.shape
	#
	print('shapes: ', sh1, sh2)
	#
	zs_diff, h, m, f = list(zip(*diffs))[2:]
	#
	# and ROC bits:
	H = sum(h)/sum(z2)
	F = sum(f)/sum(z1)
	#
	#for z in [zs_diff, h, m, f]:
	for j,z in enumerate(list(zip(*diffs))):
		plt.figure(j+2)
		plt.clf()
		#
		zz=numpy.array(z)
		zz.shape=sh1
		#plt.contourf(list(set(etas_fc.ETAS_array['x'])), list(set(etas_fc.ETAS_array['y'])), zz, 25)
		#plt.contourf(numpy.log10(zz), 25)
		plt.contourf(lon_vals, lat_vals, zz, 25)
		plt.title(diffs_lbls[j])
		plt.colorbar()
	#plt.figure(j+3)
	#plt.clf()
	#gzintas.shape=sh1
	#fgz = plt.contourf(lon_vals, lat_vals, numpy.log10(gzintas), 15)
	#plt.title('z_fc/z_cat')
	#plt.colorbar()
	#
	if diagnostic:
		return [diffs_lbls] + diffs
	else:
		return F,H
	#return F,H
	#
#
def nepal_linear_roc():
	# production figure script (almost... just script ranges).
	diffs = analyze_etas_roc_geospatial(etas_fc=None, etas_test=None, do_log=True, diagnostic=True)
	AA=eap.roc_gs_linear_figs(diffs)
#
def roc_gs_linear_figs(diffs, fignum=0):
	# test the roc_gs bit. basically, take two xyz arrays, do the gs_roc thing;
	# plot out the various arrays like time-series. show the various H,F, etc. in time series.
	#pass
	cols = diffs[0]
	diffs = diffs[1:]
	print('cols: ', cols)
	#
	diffs = numpy.core.records.fromarrays(zip(*diffs),  names=cols, formats=['float' for c in cols])
	#
	f = plt.figure(fignum)
	plt.clf()
	ax0 = f.add_axes([.1,.1,.4,.4])
	ax1 = f.add_axes([.1,.55,.4,.4], sharex=ax0)
	ax2 = f.add_axes([.55,.1,.4,.4], sharex=ax0, sharey=ax1)
	ax3 = f.add_axes([.55,.55, .4,.4], sharex=ax0, sharey=ax1)
	#
	plt.figure(fignum+1)
	plt.clf()
	ax_main=plt.gca()
	
	#
	X=numpy.arange(len(diffs))
	#
	for ax in (ax0, ax_main):
		ax.plot(X, diffs['z_fc'], '-', lw=3., label='forecast')
		ax.plot(X, diffs['z_test'], '-', lw=3., label='test')
		ax.fill_between(X, y1=numpy.zeros(len(diffs['z_fc'])), y2=diffs['hits'], color='c', alpha=.25)
		ax.fill_between(X, y1=diffs['z_test'], y2=diffs['z_fc'], where=[True if zfc<ztest else False for zfc,ztest in zip(diffs['z_fc'],diffs['z_test'])], label='misses', color='r', alpha=.8)
		ax.fill_between(X, y1=diffs['z_test'], y2=diffs['z_fc'], where=[True if zfc>ztest else False for zfc,ztest in zip(diffs['z_fc'],diffs['z_test'])], label='falsies', color='m', alpha=.25)
		ax.legend(loc=0, numpoints=1)
		ax.set_title('Hits, False Alarms, Misses')
	
	ax1.plot(X, diffs['hits'], 'm-',   lw=2., alpha=.8, label='hits')
	ax1.plot(X, diffs['falsie'], 'r-', lw=2., alpha=.8, label='falsies')
	ax1.plot(X, diffs['misses'], 'b-', lw=2., alpha=.8, label='misses')
	ax1.legend(loc=0, numpoints=1)
	ax1.set_title('Hit, falsies, misses')
	#
	ax2.plot(X, diffs['z_fc'], '-', lw=2., label='forecast')
	ax2.plot(X, diffs['z_test'], '-', lw=2., label='test')
	ax2.fill_between(X, y1=diffs['z_fc'],y2=diffs['z_test'], where=[True if zfc>ztest else False for zfc,ztest in zip(diffs['z_fc'],diffs['z_test'])], label='falsies', color='b')
	ax2.set_title('False Alarms')
	ax2.legend(loc=0, numpoints=1)
	#
	ax3.plot(X, diffs['z_fc'], '-', lw=2., label='forecast')
	ax3.plot(X, diffs['z_test'], '-', lw=2., label='test')
	ax3.fill_between(X, y1=diffs['z_fc'],y2=diffs['z_test'], where=[True if zfc<ztest else False for zfc,ztest in zip(diffs['z_fc'],diffs['z_test'])], label='misses', color='r')
	ax3.set_title('Misses')
	ax3.legend(loc=0, numpoints=1)
	
#
def plot_mainshock_and_aftershocks(etas, m0=6.0, mainshock=None, fignum=0):
	
	map_etas = etas.make_etas_contour_map(n_contours=25, fignum=fignum)
	if mainshock==None:
		mainshock = etas.catalog[0]
		for rw in etas.catalog:
			if rw['mag']>mainshock['mag']: mainshock=rw
	ms=mainshock
	#
	for eq in etas.catalog:
		if eq['mag']<m0 or eq['event_date']<ms['event_date']: continue
		if eq==ms:
			x,y = map_etas(eq['lon'], eq['lat'])
			plt.plot([x], [y], 'k*', zorder=7, ms=20, alpha=.8)
			plt.plot([x], [y], 'r*', zorder=8, ms=18, label='mainshock', alpha=.8)
		if eq['event_date']>eq['event_date']:
			x,y = map_etas(eq['lon'], eq['lat'])
			plt.plot([x], [y], 'o', zorder=7, ms=20, alpha=.8)	
	
