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
from optimizers import roc_tools

import contours2kml
import yodiipy.ANSStools as atp

from eq_params import *

#colors_ =  mpl.rcParams['axes.color_cycle']
colors_ = ['b','g','r','c','m','y','k']
nepal_mainshock = {'mag':7.8, 'lon':84.708, 'lat':28.147}
#
# TODO: replace all ROCgeneric with optimizers.roc_tools
#  there may be problems with the ROCgeneric ROC calcs... and also roc_tools is simpler and maybe faster (though
#  we need to investigate to see if we can get a gain from mpp).

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
		
			
	def draw_map(self, fignum=0, fig_size=(6.,6.), map_resolution='i', map_projection='cyl', d_lon_range=1., d_lat_range=1., lats=None, lons=None):
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
		#lons, lats = self.lons, self.lats
		lons = (lons or self.lons)
		lats = (lats or self.lats)
		#	
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

#
def toy_gs_roc(fignum=0):
	z1=list(range(10))
	z2=reversed(list(range(10)))
	diffs = etas_analyzer.get_gs_diffs(z1,z2)
	#
	X=list(range(len(diffs)))
	plt.figure(fignum)
	plt.clf()
	#plt.plot(*zip(*stepify(list(zip(X, diffs['z_fc'])))), marker='.', ls='-', label='foreccast')
	#plt.plot(*zip(*stepify(list(zip(X, diffs['z_test'])))), marker='.', ls='-', label='test_cat')
	plt.plot(X, diffs['z_fc'], marker='.', ls='-', label='foreccast')
	plt.plot(X, diffs['z_test'], marker='.', ls='-', label='test_cat')
	#
	plt.plot(X, diffs['hits'], '--', label='hits')
	plt.plot(X, diffs['misses'], '--', label='misses')
	plt.plot(X, diffs['falsie'], '--', label='falsies')
	plt.legend(loc=0, numpoints=1)
	#
	return diffs
#
def stepify(xy):
	# lets move (or copy) stepify function(s) to a generic repo.
	# also, the main purpose of this little function was to make molchan calculations look like ROC. we have better
	# ROC scripts, and maybe we're better off just letting molchan have slopes, so we can visually distinguish it from ROC.
	#
	xy_prime = [[x,y] for x,y in xy] + [[xy[j][0],y] for j,(x,y) in enumerate(xy[1:])]
	
	xy_prime.sort(key=lambda rw:rw[0])
	return xy_prime

#
def nepal_roc_script(fignum=0, mcs = [4., 5., 6., 7.], n_cpu=None):
	return Nepal_ROC_script(**locals())
class nepal__ROC_script(object):
	def __init__(self, fignum=0, mcs = [4., 5., 6., 7.], n_cpu=None):
		# this needs to be rewritten a bit to:
		# 1) use the same color for each magnitude
		# 2) should probably use the roc_generic class; see _rocs3()
		#    **** correction **: REMOVE roc_generic, use optimziers.roc_tools ROC tools instead. see tons of new examples.
		#     ... also, this script has been rewritten in the nepal-revisions notebook, so before too long, replace this
		#     implementation with the new, shiny, code.
		#
		# full, one stop shopping script for nepal ROC analysis.
		#
		# first, get nepal ETAS objects:
		etas_fc, etas_test = etas_analyzer.nepal_etas_roc()
		test_catalog = etas_test.catalog
		#
		x0 = nepal_epi_lon
		y0 = nepal_epi_lat
		mag=7.8
		L_r = .5*10**(.5*mag - 1.76)
		xyz = etas_fc.ETAS_array
		#
		# now, replace all of this "get x,y and ROC" stuff with the optimizers.roc_tools equivalents.
		X_set = sorted(list(set(xyz['x'])))
		Y_set = sorted(list(set(xyz['y'])))
		nx = len(X_set)
		ny = len(Y_set)
		get_site = lambda x,y: int(round((x-etas_fc.lons[0]+.5*etas_fc.d_lon)/etas_fc.d_lon)) + int(round((y-etas_fc.lats[0]+.5*etas_fc.d_lat)/etas_fc.d_lat))*nx
		#
		xyz_r = xyz.copy()
		for j,(x,y,z) in enumerate(xyz_r):
			xyz_r['z'][j]=1./(dist_to(x,y,x0,y0) + .5*L_r)
		#
		Zs = sorted(list(xyz['z'].copy()))
		Zs_r = sorted(list(xyz_r['z'].copy()))
		#
		eq_site_zs = [[xyz['z'][get_site(eq['lon'], eq['lat'])], eq['mag']] for eq in test_catalog]
		#
		plt.figure(fignum)
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
			this_etas = etas_analyzer.Toy_etas_random(etas_in=etas_fc)
			FH = etas_analyzer.roc_normal(this_etas, fignum=None)
			plt.plot(*zip(*FH), marker='.', ls='', alpha=.6)
		#
		self.__dict__.update(locals())
		#

#
def inv_dist_to(xy,x0,y0,r0):
	return [[x,y, 1./(globalETAS.spherical_dist(lon_lat_from=[x0,y0], lon_lat_to=[x, y], Rearth = 6378.1) + r0)] for x,y in xy]
def dist_to(x,y,x0,y0):
	return globalETAS.spherical_dist(lon_lat_from=[x0,y0], lon_lat_to=[x, y], Rearth = 6378.1)
##################################
#################################
# good stuff, that i think we can keep, happening below this line...
#
# yoder: in revisions... let's try to consolidate and/or rewrite the ROC and Molchan codes. they also may need better unit testing, so set up a well
# defined unit test. note: for some implementations, keep it simple; one earthquake, one bin (especially for small lattices). if we make this approximation,
# we can run much faster algorithms than explicit ROC. note that if we run into cases where this is not true, particularly for coarse grain maps, we can --
# somewhat ironically, improve performance by fine-meshing the lattice. we can also devise a more sophisticated fast approach that counts the number of events
# in each bin/at each level.
#
# yoder 2016-07-07: this is a new, shiny, working roc-from-etas_file script. hook it up wiht "calc global ROC, and we're good to go...'
# and before very long, we need to tie in the date fields...
#
# ... and this needs to be morphed into a general etas,roc handler function... or class...
def global_etas_and_roc(fout_xyz='global_etas.xyz', fc_len=120, out_path = 'figs', fignum=0, m_cs=[4.0, 5.0, 6.0, 6.5], n_cpu=None, t_now=None):
	# a soup-to-nuts global ETAS and roc bit. calculate a global ETAS up to fc_len days ago (fc_len+1?); then do ROC on that data set.
	#
	# t_now for original globalETAS paper draft: dtm.datetime(2015,11,30, tzinfo=pytz.timezone('UTC'))
	# note: etas end time should be in the .xyz header file, so we can get it from there as well...
	#
	if not os.path.isdir(out_path): os.makedirs(out_path)
	fout_xyz = os.path.join(out_path, fout_xyz)
	#
	n_cpu = (n_cpu or mpp.cpu_count())
	if not os.path.isdir(os.path.split(fout_xyz)[0]): os.makedirs(os.path.split(fout_xyz)[0])
	#
	lats=[-89., 89.]
	lons=[-180., 180.]
	mc=3.0
	d_lon=.1
	d_lat=.1
	etas_range_factor=15.
	etas_range_padding=1.0
	etas_fit_factor=1.5
	#
	# date to which etas is calculated:
	t_now =  (t_now or dtm.datetime.now(globalETAS.tzutc)-dtm.timedelta(days=fc_len))
	cat_len=3650.
	#
	plt.figure(fignum, figsize=(12,10))
	etas = globalETAS.ETAS_mpp(lats=lats, lons=lons, mc=mc, d_lon=d_lon, d_lat=d_lat, etas_range_factor=etas_range_factor, etas_range_padding=etas_range_padding, etas_fit_factor=etas_fit_factor, t_now=t_now, cat_len=cat_len, n_cpu=n_cpu)
	mp = etas.make_etas_contour_map(fignum=fignum, map_resolution='f', lat_interval=20, lon_interval=20)
	plt.savefig('%s/etas_contours_%s.png' % (os.path.split(fout_xyz)[0], str(t_now)))
	#
	# not sure if this will take an array as an input. if not, it should...
	#draw_global_etas_contours(xyz=etas.ETAS_array, fignum=0, n_conts=15, cmap=plt.cm.jet)
	#
	with open(fout_xyz,'w') as fout:
		fout.write('#global ETAS\n#lats={lats!s}\tlons={lons!s}\tmc={mc!s}\td_lon={dlon!s}\td_lat={dlat!s}\tetas_range_factor={erf!s}\tetas_range_padding={erp!s}\tetas_fit_factor={eff!s}\tt_now={tnow!s}\tcat_len={catlen!s}\n'.format(lats=lats, lons=lons, mc=mc, dlon=d_lon, dlat=d_lat, erf=etas_range_factor, erp=etas_range_padding, eff=etas_fit_factor,tnow=t_now,catlen=cat_len))
		#
		[fout.write('\t'.join([str(x) for x in rw])+'\n') for j,rw in enumerate(etas.ETAS_array)]
		#
	#
	roc_marker=''
	roc_ls = '-'
	roc_lw=2.
	roc_glob = global_roc_from_optimizer(fc_xyz=etas.ETAS_array, fignum=fignum+1, etas_end_date=t_now+dtm.timedelta(days=1), mcs=[4.,5., 6.], fc_len=fc_len, ls=roc_ls, marker=roc_marker, lw=roc_lw)
	#
	#plt.savefig('%s/etas_global_roc_a__%s.png' % (os.path.split(fout_xyz)[0], str(t_now)))
	etas_png_fpath = os.path.join(out_path, 'etas_global_roc_a_{}.png'.format(t_now))
	plt.savefig(etas_png_fpath)
	#
	# ... and maybe we need to save the ROC data as well?
	#
	return{'etas':etas, 'roc':roc_glob}
#
# 2016-04-12 12:52:58.803348+00:00
def global_roc_from_optimizer(fc_xyz='global/global_xyz_20151129.xyz', fignum=0, etas_end_date = None, mcs=6.0, fc_len=120, ls='-', marker='.', lw=2.5, x_scale='linear', y_scale='linear'):
	#yoder, 2016_08_01:
	# note, of course, this is for a specific run of a global ETAS, so get this all stitched together as soon as possible...
	#
	#etas_end_date = dtm.datetime(2015,11,30, tzinfo=pytz.timezone('UTC'))		# ran ETAS sometime on 29 Nov. so we'll start our test period after the 30th.
	
	# this is still a bit of a mess, but this default behavior is consistent with the calc_global_etas behavior, above.
	# dtm.datetime(2015,11,30, tzinfo=pytz.timezone('UTC'))
	etas_end_date = (etas_end_date or dtm.datetime.now(globalETAS.tzutc)-dtm.timedelta(days=fc_len-1))
	
	fc_end_date = etas_end_date + dtm.timedelta(days=fc_len)
	if not hasattr(mcs, '__getitem__'): mcs = [mcs]
	#
	if isinstance(fc_xyz,str):
		with open(fc_xyz,'r') as f:
			fc_xyz = [[float(x) for x in rw.split()] for rw in f if rw[0] not in ('#', ' ', '\t', '\n')]
			fc_xyz = numpy.core.records.fromarrays(zip(*fc_xyz), dtype=[('x','float'), ('y','float'),('z','float')])
		#
	#	
	X_set = sorted(list(set(fc_xyz['x'])))
	Y_set = sorted(list(set(fc_xyz['y'])))
	lons = [min(X_set), max(X_set)]
	lats = [min(Y_set), max(Y_set)]
	#
	d_lon = abs(X_set[1] - X_set[0])
	d_lat = abs(Y_set[1] - Y_set[0])
	nx = len(X_set)
	ny = len(Y_set)
	#
	mc0 = min(mcs)
	print("get cataog: ", lons, lats, mc0, etas_end_date, fc_end_date)
	test_catalog = atp.catfromANSS(lon=lons, lat=lats, minMag=mc0, dates0=[etas_end_date, fc_end_date])
	#test_catalog = atp.catfromANSS(lon=lons, lat=lats, minMag=mc, dates0=[etas_end_date, fc_end_date])
	print("catlen: ", len(test_catalog))
	#
	FHs = []
	for mc in mcs:
		events_xyz = [[rw['lon'], rw['lat'], rw['mag']] for rw in test_catalog if rw['mag']>=mc]
		#
		#FH = roc_tools.calc_roc(Z_fc, Z_ev)
		roc_obj = roc_tools.ROC_xyz_handler(fc_xyz=fc_xyz, events_xyz=events_xyz)
		FHs += [[mc, roc_obj.calc_roc()]]
	#	
	if not fignum is None:
		plt.figure(fignum)
		plt.clf()
		ax=plt.gca()
		ax.set_xscale(x_scale)
		ax.set_yscale(y_scale)
		#		
		for j,(mc, FH) in enumerate(FHs):
			ax.plot(*zip(*FH), marker=marker, ls=ls, lw=lw, label='$m_c={:.2f}$'.format(mc), zorder=5)
		# plot a bunch of points along the H=F line in case we log-transform...
		ax.plot(numpy.linspace(0.,1., 250), numpy.linspace(0.,1., 250), color='r', ls='--', lw=2., label='$H=F$', zorder=4)
		ax.legend(loc=0, numpoints=1)
		#
		# draw random roc:
		roc_random(n_events=1000, n_fc=10000, n_rocs=100, ax=ax, n_bins=100, line_color='m', shade_color='m', zorder=1)
		#
		ax.set_xlabel('False Alarm rate $F$', size=18)
		ax.set_ylabel('Hit Rate $H$', size=18)
		ax.set_title('Golbal ETAS ROC\n{} + {} days'.format(etas_end_date, fc_len))
	#
	ax.set_xlim(0., 1.05)
	ax.set_ylim(0.,1.05)
	if len(FH)==1:
		return FH[0][1]
	else:
		return FHs
####################


def q_q_skill_figs(data='data/roc_geospatial_nepal_q11_24_11_24.csv'):
	# make some figures for geospatial roc
	#
	with open(data,'r') as f:
		#X = [[x,y, h-f] for x,y,f,h in rw.split() for rw in f if rw[0]!='#']
		X = []
		for rw in f:
			if rw[0]=='#': continue
			x,y,f,h=[float(x) for x in rw.split()]
			X += [[x,y,h-f]]
	#
	#
	fg = plt.figure(0)
	plt.clf()
	ax = fg.add_subplot(111, projection='3d')
	ax.set_xlabel('$q_{fc}$', size=18)
	ax.set_ylabel('$q_{test}$', size=18)
	ax.set_zlabel('skill = $H-F$', size=18)
	#
	ax.scatter(*zip(*X), marker='.')
	#
	Xs = list(set([rw[0] for rw in X]))
	Ys = list(set([rw[1] for rw in X]))
	n_x=len(Xs)
	n_y=len(Ys)
	Z=numpy.array([rw[2] for rw in X])
	Z.shape=(n_x, n_y)
	#plt.clf()
	xx = numpy.array([rw[0] for rw in X])
	yy = numpy.array([rw[1] for rw in X])
	xx.shape = Z.shape
	yy.shape = Z.shape
			
	#ax.plot_wireframe(xx,yy,Z)
	ax.plot_surface(xx,yy,Z, cmap='jet')
	ax.plot_trisurf(*zip(*X), cmap='coolwarm', lw=.5)
	cset=ax.contourf(Xs, Ys, list(zip(*Z)), 25, zdir='z', offset=.67, cmap=mpl.cm.coolwarm)
	#
	plt.figure(1)
	plt.clf()
	plt.contourf(Xs,Ys,list(zip(*Z)), 25, cmap=mpl.cm.coolwarm)
	#plt.contourf(Xs,Ys,Z, 25, cmap=mpl.cm.coolwarm)
	plt.xlabel('$q_{fc}$', size=18)
	plt.ylabel('$q_{test}$', size=18)
	plt.title('Continuum ROC Skill, $H-F$')
	plt.colorbar()
	
	return X


def draw_global_etas_contours(xyz='data/global_xyz_20151129.xyz', fignum=0, n_conts=15, cmap=plt.cm.jet):
	plt.figure(fignum)
	plt.clf()
	mm = Map_drawer(xyz=xyz)
	mm.draw_map(d_lat_range=10., d_lon_range=20., fignum=fignum)
	#return mm
	#
	lns, lts = mm.lonses, mm.latses
	Zs = numpy.log10([rw[2] for rw in mm.XYZ])
	Zs.shape=(len(lts), len(lns))
	#
	plt.figure(fignum)
	# plt.cm.coolwarm
	print('cmap: ', cmap)
	plt.contourf(lns, lts, Zs, n_conts, alpha=.65, zorder=11, cmap=cmap)
	plt.colorbar()

def roc_random(n_events=100, n_fc=10000, n_rocs=100, ax=None, n_bins=100, line_color='m', shade_color='m', zorder=1):
	# calculate a bunch of random ROCs (random events, random forecasts) and plot an envelope figure.
	# by passing ax={matplotlib.axes (axis?)} object, we can use this script to plot a "random-envelope" onto an 
	# independently computed ROC figure.
	#
	R=random.Random()
	#
	if ax==None:
		plt.figure()
		plt.clf()
		ax=plt.gca()
	#
	dx=1./n_bins
	j_bin = lambda x: int(x/dx)
	#x_min_max = [[(j+1.)/float(n_bins), 1., 0.] for j in range(n_bins)]
	x_min_max = [[(j)/float(n_bins), 1., 0.] for j in range(n_bins+1)]
	#
	for j in range(n_rocs):
		Z_events=[R.random() for j in range(n_events)]
		Z_fc=sorted([R.random() for j in range(n_fc)])
		#
		roc_FH = roc_tools.calc_roc(Z_fc, Z_events, f_denom=None, h_denom=None)
		#
		#ax=plt.plot(C.F, C.H, '-', lw=2.0, alpha=.5)
		#ax.plot(*zip(*roc_FH), marker='.', ls='')
		#
		# bin up the F,H values (you know, maybe a better way is to use scipy.interpolate...
		#for k,(f,h) in enumerate(zip(C.F, C.H)):
		for k, (f,h) in enumerate(roc_FH):
			bin = j_bin(f)
			while bin>=len(x_min_max): 
				print('extending...')
				x_min_max += [[x_min_max[-1][0]+dx, 0.,0.]]	# sometimes we get a little integer overhang
			#print('bin: ', bin, f)
			x_min_max[bin][1]=min(x_min_max[bin][1], h)
			x_min_max[bin][2]=max(x_min_max[bin][2], h)
		#
	#
	X,mn, mx = zip(*x_min_max)
	#return X,mn,mx
	ax.plot(X,mn, color=line_color, ls='-', lw=2., alpha=.7, zorder=zorder)
	ax.plot(X,mx, color=line_color, ls='-', lw=2., alpha=.7, zorder=zorder)
	ax.fill_between(X,mn,mx, color=shade_color, alpha=.3, zorder=zorder)
	#
	return X,mn,mx
	
	
def roc_fig_geospatial_fast_raw():
	# roc figure for q/q range(s).use pre-compiled data.
	#
	FH_fast = []
	FH_raw  = []
	qfc_max=2.4
	qtest_max=2.4
	#
	with open('data/roc_geospatial_nepal_raw.csv', 'r') as f:
		for rw in f:
			if rw[0]=='#': continue
			qfc,qtest,f,h = [float(x) for x in rw[:-1].split()]
			if qfc>qfc_max or qtest>qtest_max: continue
			FH_raw += [qfc,qtest,f,h,h-f]
			#
		#
	#
	with open('data/roc_geospatial_nepal_q11_24_11_24.csv', 'r') as f:
		for rw in f:
			if rw[0]=='#': continue
			qfc,qtest,f,h = [float(x) for x in rw[:-1].split()]
			if qfc>qfc_max or qtest>qtest_max: continue
			FH_raw += [qfc,qtest,f,h,h-f]
	#
	plt.figure(0)
	plt.clf()
	plt.plot(range(2), range(2), ls='--', lw=3.5, alpha=.8, color='r', marker='', label='Random')
	plt.plot(*zip(*[[rw[2], rw[4]] for rw in FH_raw]), marker='o', color='b', ls='',label='ROC (raw)', zorder=4)
	plt.plot(*zip(*[[rw[2], rw[4]] for rw in FH_fast]), marker='s', color='g', ls='',label='ROC (raw)', zorder=5)
	plt.legend(loc=0, numpoints=1)
	plt.title('Nepal ROC-geospatial Analysis')
	plt.xlabel('False Alarm Rate $F$', size=18)
	plt.ylabel('Hit Rate $H$', size=18)
	#
	FH_raw.sort(key = lambda rw: rw[-1])
	FH_fast.sort(key=lambda rw: rw[-1])
	#
	print('fast:')
	print(FH_fast[-4:])
	print('raw:')
	print(FH_raw[-4:])
	
