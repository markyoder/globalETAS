import pylab as plt
import numpy
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
import random
import multiprocessing as mpp
import sys
import os
import datetime as dtm
import pytz
from geographiclib.geodesic import Geodesic as ggp

import globalETAS
import etas_analyzer
import roc_generic
import contours2kml
import ANSStools as atp

from eq_params import *

colors_ =  mpl.rcParams['axes.color_cycle']

nepal_mainshock = {'mag':7.8, 'lon':84.708, 'lat':28.147}

class Map_drawer(object):
	# container class to draw maps and stuff like that.
	def __init__(self, xyz='data/global_xyz_20151129.xyz'):
		#self.__dict__.update(locals())
		#self.fignum=fignum
		if isinstance(xyz, str):
			with open(xyz, 'r') as f:
				xyz = [[float(x) for x in rw.split()] for rw in f if rw[0]!='#']
			#
		#
		self.XYZ = xyz
		#
		self.lonses=sorted(list(set([rw[0] for rw in self.XYZ])))
		self.latses=sorted(list(set([rw[1] for rw in self.XYZ])))
		self.lons = (min(self.lonses), max(self.lonses))
		self.lats = (min(self.latses), max(self.latses))
		#
		
			
	def draw_map(self, fignum=0, fig_size=(6.,6.), map_resolution='i', map_projection='cyl', d_lon_range=1., d_lat_range=1.):
		'''
		# plot contours over a map.
		'''
		#
		# first, get contours:
		#etas_contours = self.calc_etas_contours(n_contours=n_contours, fignum=fignum, contour_fig_file=contour_fig_file, contour_kml_file=contour_kml_file, kml_contours_bottom=kml_contours_bottom, kml_contours_top=kml_contours_top, alpha_kml=alpha_kml, refresh_etas=refresh_etas)
		#
		# now, clear away the figure and set up the basemap...
		#		
		plt.figure(fignum, fig_size)
		plt.clf()
		#
		lons, lats = self.lons, self.lats
		cntr = [numpy.mean(lons), numpy.mean(lats)]
		cm = Basemap(llcrnrlon=self.lons[0], llcrnrlat=self.lats[0], urcrnrlon=self.lons[1], urcrnrlat=self.lats[1], resolution=map_resolution, projection=map_projection, lon_0=cntr[0], lat_0=cntr[1])
		#
		#cm.drawlsmask(land_color='0.8', ocean_color='b', resolution=map_resolution)
		cm.drawcoastlines(color='gray', zorder=1)
		cm.drawcountries(color='black', zorder=1)
		cm.drawstates(color='black', zorder=1)
		#cm.drawrivers(color='blue', zorder=1)
		cm.fillcontinents(color='beige', lake_color='blue', zorder=0)
		# drawlsmask(land_color='0.8', ocean_color='w', lsmask=None, lsmask_lons=None, lsmask_lats=None, lakes=True, resolution='l', grid=5, **kwargs)
		#cm.drawlsmask(land_color='0.8', ocean_color='c', lsmask=None, lsmask_lons=None, lsmask_lats=None, lakes=True, resolution=self.mapres, grid=5)
		#
		#
		cm.drawmeridians(numpy.arange(int(lons[0]), int(lons[1]), d_lon_range), color='k', labels=[0,0,1,1])
		cm.drawparallels(numpy.arange(int(lats[0]), int(lats[1]), d_lat_range), color='k', labels=[1, 1, 0, 0])
		#
		return cm
	def make_etas_contour_map(self, n_contours=None, fignum=0, fig_size=(6.,6.), contour_fig_file=None, contour_kml_file=None, kml_contours_bottom=0., kml_contours_top=1.0, alpha=.6, alpha_kml=.5, refresh_etas=False, map_resolution='i', map_projection='cyl', map_cmap='spectral'):
		n_contours = (n_contours or self.n_contours)
		#
		cm = self.draw_map(fignum=fignum, fig_size=fig_size, map_resolution=map_resolution, map_projection=map_projection)
		#X,Y = cm(numpy.array(self.lonses), numpy.array(self.latses))
		X,Y = sorted(list(set([rw[0] for rw in self.XYZ]))), sorted(list(set([rw[1] for rw in self.XYZ])))
		Zs = numpy.array(numpy.log10([rw[2] for rw in self.XYZ]))
		#
		Zs.shape = (len(Y), len(X))
		#
		etas_contours = plt.contourf(X,Y, Zs, n_contours, zorder=8, alpha=alpha)
		plt.colorbar()
		#
		self.cm=cm
		#
		return cm
		#

def draw_global_etas_contours(xyz='data/global_xyz_20151129.xyz', fignum=0, n_conts=15):
	mm = Map_drawer(xyz=xyz)
	mm.draw_map(d_lat_range=10., d_lon_range=20., fignum=0)
	#return mm
	#
	lns, lts = mm.lonses, mm.latses
	Zs = numpy.log10([rw[2] for rw in mm.XYZ])
	Zs.shape=(len(lts), len(lns))
	#
	plt.figure(fignum)
	# plt.cm.coolwarm
	plt.contourf(lns, lts, Zs, n_conts, alpha=.65, zorder=11, cmap=plt.cm.hot)
	plt.colorbar()

def roc_random(n_events=100, n_fc=10000, n_rocs=100, n_cpus=None, ax=None, n_bins=100, line_color='m', shade_color='m'):
	R=random.Random()
	n_cpus = (n_cpus or mpp.cpu_count())
	#
	if ax==None:
		plt.figure()
		plt.clf()
		ax=plt.gca()
	#
	dx=1./n_bins
	j_bin = lambda x: int(x/dx)
	x_min_max = [[(j+1.)/float(n_bins), 1., 0.] for j in range(n_bins)]
	#
	for j in range(n_rocs):
		Z_events=[R.random() for j in range(n_events)]
		Z_fc=sorted([R.random() for j in range(n_fc)])
		#
		C = roc_generic.ROC_mpp(n_procs=n_cpus, Z_events=Z_events, Z_fc=Z_fc)
		roc = C.calc_roc()
		#
		#ax=plt.plot(C.F, C.H, '-', lw=2.0, alpha=.5)
		#
		for k,(f,h) in enumerate(zip(C.F, C.H)):
			bin = j_bin(f)
			while bin>=len(x_min_max): x_min_max += [[x_min_max[-1][0]+dx, 0.,0.]]	# sometimes we get a little integer overhang
			#print('bin: ', bin, f)
			x_min_max[bin][1]=min(x_min_max[bin][1], h)
			x_min_max[bin][2]=max(x_min_max[bin][2], h)
		#
	#
	X,mn, mx = zip(*x_min_max)
	#return X,mn,mx
	ax.plot(X,mn, color=line_color, ls='-', lw=2., alpha=.7)
	ax.plot(X,mx, color=line_color, ls='-', lw=2., alpha=.7)
	ax.fill_between(X,mn,mx, color=shade_color, alpha=.3)
	#
	return X,mn,mx
	

def nepal_roc_script():
	#
	etas_fc = etas_analyzer.get_nepal_etas_fc()
	#nepal_etas_test = get_nepal_etas_test()
	Xs = sorted(list(set(etas_fc.ETAS_array['x'])))
	Ys = sorted(list(set(etas_fc.ETAS_array['y'])))
	get_site = lambda x,y: int(numpy.floor((x-lons[0])/d_lon)) + int(numpy.floor((y-lats[0])/d_lat))*nx
	
	A=etas_analyzer.roc_normalses(etas_fc, test_catalog=None, to_dt=None, cat_len=120., mc_rocs=[4.0, 5.0, 6.0, 7.0], fignum=1, do_clf=True, roc_ls='-')
	#
	ax=plt.gca()
	# now, get roc for a 1/r map (i think there's a script for that)
	etas_toy = etas_analyzer.Toy_etas_invr(etas_in=etas_fc, mainshock=nepal_mainshock)
	B=etas_analyzer.roc_normalses(etas_toy, test_catalog=None, to_dt=None, cat_len=120., mc_rocs=[4.0, 5.0, 6.0, 7.0], fignum=1, do_clf=False, roc_ls='--') 
	#
	#
	# and draw in roc for random...
	bins, mins, maxes = roc_random(n_events=100, n_fc=10000, n_rocs=100, n_cpus=None, ax=ax, n_bins=100, line_color='m', shade_color='m')
#
def global_roc():
	roc_global = etas_analyzer.roc_normal_from_xyz(fc_xyz='data/global_xyz_20151129.xyz', test_catalog=None, from_dt=None, to_dt=None, dx=None, dy=None, cat_len=120., mc=5.0, fignum=0, do_clf=True)
	return roc_global
#
def global_roc1(fc_xyz='global/global_xyz_20151129.xyz', n_cpu=None, fnum=0, m_cs=[4.0, 4.5, 5.0, 6.0, 7.0]):
	# a global ROC script; we may have some competing candidates for this right now. this seems to work. can we duplicate it with generic_roc?
	#
	# this script produced a really nice ROC, so let's clean it up a bit.
	#A=etas_analyzer.ROC_mpp_handler(n_procs=8, fc_xyz='global/global_xyz_20151129.xyz', from_dt = eap.dtm.datetime(2015, 11, 30, tzinfo=eap.pytz.timezone('UTC')), to_dt=dtm.datetime.now(eap.pytz.timezone('UTC')), mc=5.5)
	#
	n_cpu = (n_cpu or mpp.cpu_count())
	print('cpu_count: ', n_cpu)
	#
	etas_end_date = dtm.datetime(2015,11,30, tzinfo=pytz.timezone('UTC'))		# ran ETAS sometime on 29 Nov. so we'll start our test period after the 30th.
	fc_len=120	# test period is 120 days.
	#
	# something's not right; maybe a memory state problem, so let's try just re-delcaring/creating this object for each mc.
	# this will be expensive, but create the object here and then poach its components.
	roc0=etas_analyzer.ROC_mpp_handler(n_procs=n_cpu, fc_xyz=fc_xyz, from_dt = etas_end_date, to_dt=etas_end_date + dtm.timedelta(days=120), mc=min(m_cs))
	#
	# now, make a 1/r forecast. we'll need to specify a mainshock (or in general, center) lat,lon.
	# (note: for the sake of conserving memory (of which this can use a lot), we should do the regualr roc then the 1/r roc afterwards. we might also use a
	# faster distance algorithm. .Inverse() can take an option to use a simple spherical solution -- i think. we might also parallelize the loop through the array.
	# ... of course, unless we can come up with a more general 1/r (aka, not 1/r from gorkah, this is all nonesense. so what this basically will amount to 
	# (if we do it correctly) is comparing a bunch of time-dependent global (or other) fc to a static fc (with p=0)... so save it for another paper.
	'''
	print('creating 1/r forecast...')
	x0,y0,m0 = nepal_mainshock['lon'], nepal_mainshock['lat'], nepal_mainshock['mag']
	L_r = 10**(.5*m0-1.76)
	xyz=[]
	# ... (eventually) design this to first read the xyz array, then split it up to a special function (here) which returns [[x,y,1/r]]
	if n_cpu>1 and False:
		# ... eventually, make the mpp work...
		#with open(fc_xyz,'r') as f:
		#	xy = [[float(x) for x in rw.split()[0:2]] for rw in f if rw[0]!='#']
		xy = listt(zip(roc.fc_xyz['x'].tolist(), roc.fc_xyz['y'].tolist()))
		P=mpp.Pool(processes=n_cpu)
		P_results = [P.apply_async(inv_dist_to, kwds={'xy':xyz[j:int(numpy.ceil(len(xyz)/n_cpu))], 'x0':x0, 'y0':y0}) for j in range(n_cpu)]
		P.close()
		xyz=[]
		for r in P_results:
			xyz += r.get()
		P.join()
		#
		xyz.sort(key=lambda rw: (rw[0], rw[1]))
		#
	else:
		# do spp
		# a little slower than a .tolist(), or something, but we explicitly get list-o-lists.
		xyz = [[rw[0], rw[1], 1./(globalETAS.spherical_dist(lon_lat_from=[x0,y0], lon_lat_to=[rw[0], rw[1]], Rearth = 6378.1) + .5*L_r)] for rw in roc.fc_xyz]
	#
	print('xyz[0:5]: ', xyz[0:5])
	xyz = numpy.core.records.fromarrays(zip(*xyz), dtype=[('x','<f8'), ('y', '<f8'), ('z','<f8')])
	#
	print('1/r forecast created. now make 1/r roc object.')
	roc_r=etas_analyzer.ROC_mpp_handler(n_procs=n_cpu, fc_xyz=xyz, from_dt = etas_end_date, to_dt=etas_end_date + dtm.timedelta(days=120), mc=min(m_cs))
	'''
	#
	print('... and start doing roc calcs...')
	#
	fg1=plt.figure(fnum)
	plt.clf()
	for j,mc in enumerate(m_cs):
		print('calcing roc for m_c=%f' % mc)
		clr = colors_[j%len(colors_)]
		test_cat = numpy.core.records.fromarrays(zip(*[rw for rw in roc0.test_catalog if rw['mag']>=mc]), dtype=roc0.test_catalog.dtype)
		roc=etas_analyzer.ROC_mpp_handler(n_procs=n_cpu, fc_xyz=roc0.fc_xyz, test_catalog=test_cat, from_dt = etas_end_date, to_dt=etas_end_date + dtm.timedelta(days=120), mc=mc)
		X=roc.calc_ROCs(n_procs=n_cpu, m_c=mc)
		#roc.plot_HF(fignum=fnum, do_clf=False)
		plt.figure(fnum)
		plt.plot(*zip(*roc.FH), ls='-', marker='', lw=3., color=clr, label='m_c=%f' % mc)
		plt.figure()
		plt.plot(*zip(*roc.FH), ls='-', marker='', lw=3., color=clr, label='m_c=%f' % mc)
		plt.figure(fnum)
		
	#
	plt.plot(range(2),range(2), 'r--', alpha=.7, lw=2.5)
	plt.legend(loc=0, numpoints=1)
	plt.title('Global ROC', size=18)
	plt.xlabel('False Alarm Rate $F$', size=18)
	plt.ylabel('Hit Rate $H$', size=18)
	plt.savefig('global_roc1.png')
	#
	return roc

def global_roc3(fc_xyz='global/global_xyz_20151129.xyz', n_cpu=None, fnum=0, m_cs=[4.0, 5.0, 6.0, 6.5], test_catalog=None, fc_start_date=None, fc_end_date=None, cat_len=120, fc_frac=1.0, fout='global_roc3.png'):
	# ok, so all of this is a mess. the roc bits need to be consolidated, cleaned up, and better modularized. we'll do some of that here, at lest by example.
	# @fc_frac: fraction of fc sites to use (aka, skip the first (1-frac)*N (low z) sites.
	# fetch teh fc_xyz, generate a test catalog, then use generic_roc tools (mpp implementations) to do the roc.
	# with the test catalog: 1) get the test catalog for the region, time domain (maybe write  script to do this)
	# then, get the fc_z-values, but keep those as pairs like [[z, mag], ...]. this way, we can quickly make z_events lists.
	#
	n_cpu = (n_cpu or mpp.cpu_count())
	#
	# for now, with this script, we're assuming that we are using this specific file, but we might pre-load it.
	if isinstance(fc_xyz, str):
		with open(fc_xyz, 'r') as froc:
			fc_xyz= [[float(x) for x in rw.split()] for rw in froc if rw[0] not in('#', ' ', '\t', '\n')]
		#
	#
	#return fc_xyz
	fc_start_date = (fc_start_date or dtm.datetime(2015,11,30, tzinfo=pytz.timezone('UTC')))
	fc_end_date   = (fc_end_date or fc_start_date + dtm.timedelta(days=cat_len))
	#
	if not hasattr(fc_xyz, 'dtype'):
		fc_xyz = numpy.core.records.fromarrays(zip(*fc_xyz), names=('x','y','z'), formats=['>f8', '>f8', '>f8'])
	#
	#lons=[-180., 180.]
	#lats=[-89., 89.]
	mc=min(m_cs)
	#
	X_set = sorted(list(set(fc_xyz['x'])))
	Y_set = sorted(list(set(fc_xyz['y'])))
	nx = len(X_set)
	ny = len(Y_set)
	lons = [min(X_set), max(X_set)]
	lats = [min(Y_set), max(Y_set)]
	d_lon = abs(X_set[1] - X_set[0])
	d_lat = abs(Y_set[1] - Y_set[0])
	get_site = lambda x,y: int(round((x-lons[0]+.5*d_lon)/d_lon)) + int(round((y-lats[0]+.5*d_lat)/d_lat))*nx
	#	
	print("get cataog: ", lons, lats, mc, fc_start_date, fc_end_date)
	if test_catalog==None:
		test_catalog = atp.catfromANSS(lon=lons, lat=lats, minMag=min(m_cs), dates0=[fc_start_date, fc_end_date])
	print("catlen: ", len(test_catalog))
	#
	# forecast z-values.
	#Zs = sorted(list(fc_xyz['z'].copy()))
	Zs = sorted(list(fc_xyz['z']))
	#
	# now, get both the eq magnitudes and site_z values, so we can change the mc threshold later.
	eq_site_zs = [[fc_xyz['z'][get_site(eq['lon'], eq['lat'])], eq['mag']] for eq in test_catalog]
	#eq_site_zs.sort(key=lambda x: x[1])
	#
	plt.figure(fnum)
	plt.clf()
	plt.plot(range(2), range(2), ls='--', color='m', lw=3., alpha=.75, zorder=2)
	FHs = {}		# we'll use mc as a key, FH as a val: {mc:[FH]...}
	#				# it won't key well because of floating point error (aka, FHs[5.5] will not be reliable. but it will make a decent container.
	#
	for j,mc in enumerate(m_cs):
		print('doing ROC for mc=%f' % mc)
		#
		# n_procs,Z_events, Z_fc, h_denom=None, f_denom=None, f_start=0., f_stop=None
		roc = roc_generic.ROC_mpp(n_procs=n_cpu, Z_events=[z for z,m in eq_site_zs if m>=mc], Z_fc=Zs, h_denom=None, f_denom=None, f_start=0, f_stop=None)
		a=roc.calc_roc()		# no parameters, and in fact no return value, but it's never a bad idea to leave a place-holder for one.
		#
		clr = colors_[j%len(colors_)]
		plt.plot(roc.F, roc.H, ls='-', color=clr, marker='', lw=2.5, label='$m_c=%.2f$' % mc)
		#FHs[mc]=[[f,h] for f,h in zip(roc.F, roc.H)]
		#
		plt.show()	# just in case...
	bins, mins, maxes = roc_random(n_events=100, n_fc=10000, n_rocs=100, n_cpus=None, ax=plt.gca(), n_bins=100, line_color='m', shade_color='m')
	#plt.fill_between(bins, mins, maxes, color='m', alpha=.3)
	#plt.plot(bins, mins, '-', lw=1.5, alpha=.8)
	#plt.plot(bins, maxes, '-', lw=1.5, alpha=.8)
	plt.legend(loc=0, numpoints=1)
	plt.title('Global ROC', size=18)
	plt.xlabel('False Alarm Rate $F$', size=18)
	plt.ylabel('Hit Rate $H$', size=18)
	plt.savefig(fout)
	#
	#
	return FHs
#
def global_roc_comparison(fc_xyz='global/global_xyz_20151129.xyz', n_cpu=None, fnum=0, mc=6.0, roc_fracs=[1.0, .8, .5], test_catalog=None, fc_start_date=None, fc_end_date=None, fc_frac=1.0):
	#
	# (turns out that i don't think we really need this script. its apparent value appears to have resulted from an error
	# in the global_roc scripts.
	n_cpu = (n_cpu or mpp.cpu_count())
	#
	# for now, with this script, we're assuming that we are using this specific file, but we might pre-load it.
	if isinstance(fc_xyz, str):
		with open(fc_xyz, 'r') as froc:
			fc_xyz= [[float(x) for x in rw.split()] for rw in froc if rw[0] not in('#', ' ', '\t', '\n')]
		#
	#
	fc_start_date = (fc_start_date or dtm.datetime(2015,11,30, tzinfo=pytz.timezone('UTC')))
	fc_end_date   = (fc_end_date or fc_start_date + dtm.timedelta(days=120))
	#
	if not hasattr(fc_xyz, 'dtype'):
		fc_xyz = numpy.core.records.fromarrays(zip(*fc_xyz), names=('x','y','z'), formats=['>f8', '>f8', '>f8'])
	#
	X_set = sorted(list(set(fc_xyz['x'])))
	Y_set = sorted(list(set(fc_xyz['y'])))
	nx = len(X_set)
	ny = len(Y_set)
	lons = [min(X_set), max(X_set)]
	lats = [min(Y_set), max(Y_set)]
	d_lon = abs(X_set[1] - X_set[0])
	d_lat = abs(Y_set[1] - Y_set[0])
	get_site = lambda x,y: int(round((x-lons[0]+.5*d_lon)/d_lon)) + int(round((y-lats[0]+.5*d_lat)/d_lat))*nx
	#	
	print("get cataog: ", lons, lats, mc, fc_start_date, fc_end_date)
	if test_catalog==None:
		test_catalog = atp.catfromANSS(lon=lons, lat=lats, minMag=mc, dates0=[fc_start_date, fc_end_date])
	print("catlen: ", len(test_catalog))
	#
	# forecast z-values.
	#Zs = sorted(list(fc_xyz['z'].copy()))
	Zs = sorted(list(fc_xyz['z']))
	#
	# now, get both the eq magnitudes and site_z values, so we can change the mc threshold later.
	#eq_site_zs = [[Zs[get_site(eq['lon'], eq['lat'])], eq['mag']] for eq in test_catalog if eq['mag']>=mc]
	eq_site_zs = [fc_xyz['z'][get_site(eq['lon'], eq['lat'])] for eq in test_catalog if eq['mag']>=mc]
	#eq_site_zs.sort(key=lambda x: x[1])
	#
	plt.figure(fnum)
	plt.clf()
	plt.plot(range(2), range(2), ls='--', color='r', lw=3., alpha=.75, zorder=2)
	#
	for j,frac in enumerate(roc_fracs):
		#
		print('calcing roc for mc=%f, frac=%f' % (mc, frac))
		j_roc = int(len(fc_xyz)*(1.-frac))
		Zs = sorted(list(fc_xyz['z'][j_roc:].copy()))
		#
		# n_procs,Z_events, Z_fc, h_denom=None, f_denom=None, f_start=0., f_stop=None
		roc = roc_generic.ROC_mpp(n_procs=n_cpu, Z_events=eq_site_zs, Z_fc=Zs, h_denom=None, f_denom=None, f_start=0, f_stop=None)
		a=roc.calc_roc()		# no parameters, and in fact no return value, but it's never a bad idea to leave a place-holder for one.
		#
		plt.plot(roc.F, roc.H, ls='-', label='$m_c=%f, frac=%f' % (mc, frac))
	#
	plt.title('ROC Comparison')
	plt.xlabel('False Alarm Rate $F$')
	plt.ylabel('Hit Rate $H$')
#
def nepal_roc_script(fignum=0, mcs = [4., 5., 6., 7.]):
	# this needs to be rewritten a bit to:
	# 1) use the same color for each magnitude
	# 2) should probably use the roc_generic class; see _rocs3()
	#
	# full, one stop shopping script for nepal ROC analysis.
	#
	# first, get nepal ETAS objects:
	etas_fc, etas_test = etas_analyzer.nepal_etas_roc()
	test_catalog = etas_test.catalog
	#
	x0 = nepal_epi_lon
	y0 = nepal_epi_lat
	L_r = .5*10**(.5*nepal_ETAS_prams['mag'] - 1.76)
	xyz = etas_fc.ETAS_array
	#
	X_set = sorted(list(set(xyz['x'])))
	Y_set = sorted(list(set(xyz['y'])))
	nx = len(X_set)
	ny = len(Y_set)
	get_site = lambda x,y: int(round((x-etas_fc.lons[0]+.5*etas_fc.d_lon)/etas_fc.d_lon)) + int(round((y-etas_fc.lats[0]+.5*etas_fc.d_lat)/etas_fc.d_lat))*nx
	#
	xyz_r = xyz.copy()
	for j,(x,y,z) in enumerate(xyz_r):
		rw['z']=1./dist_to(x,y,x0,y0,.5*L_r)
	#
	Zs = sorted(list(fc_xyz['z'].copy()))
	Zs_r = sorted(list(fc_xyz_r['z'].copy()))
	#
	eq_site_zs = [[fc_xyz['z'][get_site(eq['lon'], eq['lat'])], eq['mag']] for eq in test_catalog]
	#
	plt.figure(fnum)
	plt.clf()
	plt.plot(range(2), range(2), ls='--', color='m', lw=3., alpha=.75, zorder=2)
	FHs = {}		# we'll use mc as a key, FH as a val: {mc:[FH]...}
	for j,mc in enumerate(mcs):
		clr = colors_[j%len(colors_)]
		roc = roc_generic.ROC_mpp(n_procs=n_cpu, Z_events=[z for z,m in eq_site_zs if m>=mc], Z_fc=Zs, h_denom=None, f_denom=None, f_start=0, f_stop=None)
		a=roc.calc_roc()		# no parameters, and in fact no return value, but it's never a bad idea to leave a place-holder for one.
		#
		plt.plot(roc.F, roc.H, ls='-', color=clr, marker='', lw=2.5, label='$m_c=%.2f$' % mc)
		#FHs[mc]=[[f,h] for f,h in zip(roc.F, roc.H)]
		#
		roc = roc_generic.ROC_mpp(n_procs=n_cpu, Z_events=[z for z,m in eq_site_zs if m>=mc], Z_fc=Zs_r, h_denom=None, f_denom=None, f_start=0, f_stop=None)
		roc.calc_roc()
		plt.plot(roc.F, roc.H, ls='--', color=clr, marker='', lw=2.5, label='$m_c=%.2f$' % mc)
		#
		plt.show()	# just in case...
	
	
	#ROC_n = roc_normalses(etas_fc, test_catalog=None, to_dt=None, cat_len=120., mc_rocs=[4.5, 5.0, 6.0, 7.0], fignum=fignum, do_clf=True)
	#
	# now, make a toy catalog:
	#etas_toy = Toy_etas_invr(etas_in=etas_fc, mainshock={'mag':7.3, 'lon':84.698, 'lat':28.175})
	#
	#ROC_t = roc_normalses(etas_toy, test_catalog=None, to_dt=None, cat_len=120., mc_rocs=[4.5, 5.0, 6.0, 7.0], fignum=fignum, do_clf=False, roc_ls='--')
	#
	# now, some random catalogs:
	for j in range(25):
		this_etas = Toy_etas_random(etas_in=etas_fc)
		FH = roc_normal(this_etas, fignum=None)
		plt.plot(*zip(*FH), marker='.', ls='', alpha=.6)

#
def inv_dist_to(xy,x0,y0,r0):
	return [[x,y, 1./(globalETAS.spherical_dist(lon_lat_from=[x0,y0], lon_lat_to=[x, y], Rearth = 6378.1) + r0)] for x,y in xy]
def dist_to(x,y,x0,y0,r0):
	return globalETAS.spherical_dist(lon_lat_from=[x0,y0], lon_lat_to=[x, y], Rearth = 6378.1) + r0

#
def global_roc1_single(fc_xyz='global/global_xyz_20151129.xyz', n_cpu=None, fnum=0, mc=6.0):
	# a global ROC script; we may have some competing candidates for this right now. this seems to work. can we duplicate it with generic_roc?
	#
	# this script produced a really nice ROC, so let's clean it up a bit.
	#A=etas_analyzer.ROC_mpp_handler(n_procs=8, fc_xyz='global/global_xyz_20151129.xyz', from_dt = eap.dtm.datetime(2015, 11, 30, tzinfo=eap.pytz.timezone('UTC')), to_dt=dtm.datetime.now(eap.pytz.timezone('UTC')), mc=5.5)
	#
	n_cpu = (n_cpu or mpp.cpu_count())
	#
	etas_end_date = dtm.datetime(2015,11,30, tzinfo=pytz.timezone('UTC'))		# ran ETAS sometime on 29 Nov. so we'll start our test period after the 30th.
	fc_len=120	# test period is 120 days.
	fc_end_date = etas_end_date + dtm.timedelta(days=fc_len)
	#
	roc=etas_analyzer.ROC_mpp_handler(n_procs=n_cpu, fc_xyz=fc_xyz, from_dt = etas_end_date, to_dt=fc_end_date, mc=mc)
	X=roc.calc_ROCs(n_procs=n_cpu, m_c=6.0)
	roc.plot_HF(fignum=fnum)
	#
	return roc

def global_etas_and_roc(fc_len=120, fout_xyz='figs/global_etas.xyz', fnum=0, m_cs=[4.0, 5.0, 6.0, 6.5]):
	# a soup-to-nuts global ETAS and roc bit. calculate a global ETAS up to fc_len days ago (fc_len+1?); then do ROC on that data set.
	#
	if not os.path.isdir(os.path.split(fout_xyz)[0]): os.makedirs(os.path.split(fout_xyz)[0])
	#
	lats=[-89., 89.]
	lons=[-180., 180.]
	mc=3.0
	d_lon=.1
	d_lat=.1
	etas_range_factor=15.
	etas_range_padding=.5
	etas_fit_factor=1.5
	t_now=dtm.datetime.now(globalETAS.tzutc)-dtm.timedelta(days=fc_len-1)
	cat_len=3650.
	#
	etas = globalETAS.ETAS_mpp(lats=lats, lons=lons, mc=mc, d_lon=d_lon, d_lat=d_lat, etas_range_factor=etas_range_factor, etas_range_padding=etas_range_padding, etas_fit_factor=etas_fit_factor, t_now=t_now, cat_len=cat_len)
	mp = etas.make_etas_contour_map(fignum=fnum, map_resolution='f')
	plt.savefig('%s/etas_contours_%s.png' % (os.path.split(fout_xyz)[0], str(t_now)))
	#
	with open(fout_xyz,'w') as fout:
		fout.write('#global ETAS\n#lats={lats!s}\tlons={lons!s}\tmc={mc!s}\td_lon={dlon!s}\td_lat={dlat!s}\tetas_range_factor={erf!s}\tetas_range_padding={erp!s}\tetas_fit_factor={eff!s}\tt_now={tnow!s}\tcat_len={catlen!s}\n'.format(lats=lats, lons=lons, mc=mc, dlon=d_lon, dlat=d_lat, erf=etas_range_factor, erp=etas_range_padding, eff=etas_fit_factor,tnow=t_now,catlen=cat_len))
		#
		[fout.write('\t'.join([str(x) for x in rw])+'\n') for j,rw in enumerate(etas.ETAS_array)]
		#
	#
	roc_glob = global_roc3(fc_xyz=etas.ETAS_array, n_cpu=None, fnum=fnum+1, m_cs=m_cs, test_catalog=None, fc_start_date=t_now+dtm.timedelta(days=1), fc_end_date=t_now+dtm.timedelta(days=121))
	plt.savefig('%s/etas_global_roc_a__%s.png' % (os.path.split(fout_xyz)[0], str(t_now)))
	#
	return{'etas':etas, 'roc':roc_glob}
	
	
