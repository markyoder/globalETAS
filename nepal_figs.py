import pylab as plt
import numpy
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
import random
import multiprocessing as mpp
import sys
import datetime as dtm
import pytz
from geographiclib.geodesic import Geodesic as ggp

import globalETAS
import etas_analyzer
import roc_generic
import contours2kml

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

def draw_global_etas_contours(xyz='data/global_xyz_20151129.xyz', fignum=0, n_conts=12):
	mm = Map_drawer(xyz=xyz)
	mm.draw_map(d_lat_range=10., d_lon_range=10., fignum=0)
	#return mm
	#
	lns, lts = mm.lonses, mm.latses
	Zs = numpy.log10([rw[2] for rw in mm.XYZ])
	Zs.shape=(len(lts), len(lns))
	#
	plt.figure(fignum)
	plt.contourf(lns, lts, Zs, n_conts, alpha=.5, zorder=11, cmap=plt.cm.coolwarm)

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
	plt.plot(X,mn, color=line_color, ls='-', lw=2., alpha=.7)
	plt.plot(X,mx, color=line_color, ls='-', lw=2., alpha=.7)
	plt.fill_between(X,mn,mx, color=shade_color, alpha=.3)
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
def global_roc1(fc_xyz='global/global_xyz_20151129.xyz', n_cpu=None, fnum=0, m_cs=[4.0, 5.0, 6.0, 7.0]):
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
	roc=etas_analyzer.ROC_mpp_handler(n_procs=n_cpu, fc_xyz=fc_xyz, from_dt = etas_end_date, to_dt=etas_end_date + dtm.timedelta(days=120), mc=min(m_cs))
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
		'''	
		xyz=[]	
		with open(fc_xyz,'r') as f:
			for rw in f:
				if rw[0]=='#': continue
				#rws = [float(x) for x in rw.split()]
				#
				xyz += [[float(x) for x in rw.split()]]
				
				# this geodessic  method is CRAZY slow. for global catalogs, definitely use the spherical. geodesic maybe if we can get parallel working...
				#g1=ggp.WGS84.Inverse(y0, x0, rw[1], rw[0])
				#xyz[-1][-1] = 1./1./((g1['s12']/1000.) + .5*L_r)
				xyz[-1][-1] = 1./(globalETAS.spherical_dist(lon_lat_from=[x0,y0], lon_lat_to=[xyz[-1][0], xyz[-1][1]], Rearth = 6378.1) + .5*L_r)
			
			#xyz = numpy.core.records.fromarrays(zip(*[(float(x) for x in rw.split()) for rw in f if rw[0]!='#']), dtype=[('x','<f8'), ('y', '<f8'), ('z','<f8')])
		'''
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
	ax1=plt.gca()
	#fg2=plt.figure(fnum+1)
	#plt.clf()
	#ax2=plt.gca()
	#plt.figure(fnum)
	for j,mc in enumerate(m_cs):
		print('calcing roc for m_c=%f' % mc)
		clr = colors_[j%len(colors_)]
		X=roc.calc_ROCs(n_procs=n_cpu, m_c=mc)
		#roc.plot_HF(fignum=fnum, do_clf=False)
		ax1.plot(*zip(*roc.FH), ls='-', marker='', lw=3., color=clr, label='m_c=%f' % mc)
		#ax2.plot(*zip(*roc.FH), ls='-', marker='', lw=3., color=clr, label='m_c=%f' % mc)
		#
		#X2 = roc_r.calc_ROCs(n_procs=n_cpu, m_c=mc)
		#ax1.plot(*zip(*roc_r.FH), ls='--', marker='', lw=3., color=clr, label='m_c=%f' % mc)
		#ax2.plot(*zip(*roc_r.FH), ls='--', marker='', lw=3., color=clr, label='m_c=%f' % mc)
	#
	ax1.plot(range(2),range(2), 'r--', alpha=.7, lw=2.5)
	#ax2.plot(range(2),range(2), 'r--', alpha=.7, lw=2.5)
	ax1.legend(loc=0, numpoints=1)
	#ax2.legend(loc=0, numpoints=1)
	return roc
#
def inv_dist_to(xy,x0,y0):
	return [[x,y, 1./(globalETAS.spherical_dist(lon_lat_from=[x0,y0], lon_lat_to=[x, y], Rearth = 6378.1) + .5*L_r)] for x,y in xy]

#
def global_roc1_single(fc_xyz='global/global_xyz_20151129.xyz', n_cpu=None, fnum=0):
	# a global ROC script; we may have some competing candidates for this right now. this seems to work. can we duplicate it with generic_roc?
	#
	# this script produced a really nice ROC, so let's clean it up a bit.
	#A=etas_analyzer.ROC_mpp_handler(n_procs=8, fc_xyz='global/global_xyz_20151129.xyz', from_dt = eap.dtm.datetime(2015, 11, 30, tzinfo=eap.pytz.timezone('UTC')), to_dt=dtm.datetime.now(eap.pytz.timezone('UTC')), mc=5.5)
	#
	n_cpu = (n_cpu or mpp.cpu_count())
	#
	etas_end_date = dtm.datetime(2015,11,30, tzinfo=pytz.timezone('UTC'))		# ran ETAS sometime on 29 Nov. so we'll start our test period after the 30th.
	fc_len=120	# test period is 120 days.
	#
	roc=etas_analyzer.ROC_mpp_handler(n_procs=n_cpu, fc_xyz=fc_xyz, from_dt = etas_end_date, to_dt=etas_end_date + dtm.timedelta(days=120), mc=5.5)
	X=roc.calc_ROCs(n_procs=n_cpu, m_c=6.0)
	roc.plot_HF(fignum=fnum)
	#
	return roc
	
