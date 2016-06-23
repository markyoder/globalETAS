import globalETAS as gep
import BASScast as bcp

import numpy
import scipy
import scipy.interpolate
import math
import pylab as plt
import datetime as dtm
import pytz
import random
import itertools

import matplotlib.dates as mpd
import matplotlib as mpl

days2secs = 60.*60.*24.
year2secs = 60.*60.*24.*365.
deg2km=111.31
deg2rad = 2.0*math.pi/360.
alpha_0 = 2.28

class equake_test1(object):

	def __init__(self, mc=3.0, L_r_factor=5., r_N=25, r_range=2., n_grid=None):
		#
		# use the parkfield earthquake to do some comparisons between the new Earthquake() object from globalETAS an the older
		# one from BASScast or ETASmf. at the end, run a bunch of "local intensity" calculations for random locations,times
		# about parkfield and our original "target".
		# summary (so far): the difference between the two calculations is pretty small 
		# (mean, median, max, min) error in the exponents from the two methods are ~ .0276, .0281, .1055, 8.24e-6 respectively.
		#
		if n_grid==None: n_grid = 1 + int(r_N**.5)
		#
		self.__dict__.update(locals())
		#
		# get a catalog around parkfield:
		d_lambda=1.76
		parkfield={'dt':dtm.datetime(2004,9,28,17,15,24, tzinfo=pytz.timezone('UTC')), 'lat':35.815, 'lon':-120.374, 'mag':5.96}
		L_r = 10.**(parkfield['mag']/2. - d_lambda)
		d_lat = L_r_factor * L_r/111.1
		d_lon = L_r_factor * L_r*math.cos(deg2rad*parkfield['lat'])/111.1
		#
		lats = [parkfield['lat']-d_lat, parkfield['lat']+d_lat]
		lons = [parkfield['lon']-d_lon, parkfield['lon']+d_lon]
		#
		t_0 = dtm.datetime(2000,1,1,0,0,0, tzinfo=pytz.timezone('UTC'))
		t_now = dtm.datetime.now(pytz.timezone('UTC'))
		etas_fit_factor=1.0
		#
		C=gep.make_ETAS_catalog(incat=None, lats=lats, lons=lons, mc=mc, date_range=[t_0, t_now], fit_factor=etas_fit_factor)
		#
		for c in C: 
			if c['mag']>5.92 and c['mag']<5.99:
				pf=c
		#
		# from globalETAs:
		eq1 = gep.Earthquake(pf, transform_type='equal_area', transform_ratio_max=2.)
		#
		# from older ETAS:
		# note: def __init__(self, mag=5.0, loc=[0.0,0.0], evtime=1.0, mc=2.0, drho=5.25, p=1.05, dmstar=1.2, dm=1.0, alpha = 2.28, eqtheta=0.0, eqeps=1.0, rtype='ssim', p_map=None)
		eq2 = bcp.earthquake(mag=pf['mag'], loc=[pf['lon'], pf['lat']], evtime=pf['event_date'].tolist(), mc=mc, drho=5.25, p=pf['p'], dmstar=pf['dmstar']+.2, dm=pf['dmstar'], alpha=2.28, eqtheta=None, eqeps=None, rtype='ssim', p_map=pf['p'])
		# BASScast earthquakes get their epsilon,theta fit during the BASScast initialization, so
		# we don't get that here automatically. for this test, we can fake it.
		eq2.eqtheta=-45.
		eq2.eqeps=math.sqrt(2.)
		#
		target = {'time':mpd.date2num(dtm.datetime(2005,1,1, tzinfo=pytz.timezone('UTC'))), 'lon':-120.0, 'lat':35.0}
		#
		z1 = eq1.local_intensity(t=target['time'], lon=target['lon'], lat=target['lat'])
		print("target: ", target)
		print("z1: ", z1)
		print(target)
		z2 = eq2.localIntensity(inloc=[target['lon'],target['lat']] , intime=target['time']*days2secs)
		#
		print("Zs: ", z1, z2)
		#
		# and do some random tests:
		R1=random.Random()
		R2=random.Random()
		R3=random.Random()
		lat_factor=r_range
		lon_factor=r_range
		rlat = lambda x: (R1.random()-.5)*lat_factor + x
		rlon = lambda x: (R2.random()-.5)*lat_factor + x
		rlons_lats = [[rlat(parkfield['lon']), rlat(parkfield['lat'])] for x in range(r_N)]
		mean_err=0.
		errs = []
		j=1
		xyzs = []
		#
		for lon,lat in rlons_lats:
			t = target['time'] + (R3.random()-.5)*50.		# days
		
			z1 = math.log10(eq1.local_intensity(t=t, lon=lon, lat=lat))
			z2 = math.log10(eq2.localIntensity(inloc=[lon,lat] , intime=t*days2secs))
			#
			# add to xyz arrays:
			xyzs+= [[lon, lat, z1,z2]]
			#
			#print("(%f,%f,%f) :: z_1: %f, z_2: %f (%f)" % (lon, lat, t, z1, z1, 2.*(z1-z2)/(z1+z2)))
			err = 2.*abs((z1-z2)/(z1+z2))
			mean_err+=err
			errs += [err]
			j+=1
		#
		print("mean_err: ", mean_err/float(j), numpy.median(errs), max(errs), min(errs))
		plt.figure(11)
		plt.clf()
		plt.plot(pf['lon'], pf['lat'], 'r*', ms=18, zorder=4)
		plt.plot(*zip(*rlons_lats), marker='.', ls='', color='b',zorder=3)
		
		#
		#return {'eq11':eq1, 'eq2':eq2, 'errs':errs}
		X,Y = zip(*xyzs)[0:2]
		Z1,Z2 = zip(*xyzs)[2:4]
		#
		#xi = numpy.linspace(min(X), max(X), num=n_grid)
		#yi = numpy.linspace(min(Y), max(Y), num=n_grid)
		#zi_1 = scipy.interpolate.griddata(zip(X,Y),Z1, itertools.product(xi,yi)) # named variables...
		#zi_2 = scipy.interpolate.griddata(zip(X,Y),Z2, itertools.product(xi,yi)) # named variables...
		
		# there are some differences in these figures, but they're actually an improvement,
		# (note the r' vs elliptical perim. fix, etc.) so all and all, it looks pretty good.
		plt.figure(12)
		plt.clf()
		plt.contourf(*(list(grid(X,Y,Z1, n_grid, n_grid)) + [15]))
		plt.colorbar()
		##plt.imshow(zi_1, extent=(min(X), max(X), min(Y), max(Y)), origin='lower')
		
		plt.figure(13)
		plt.clf()
		plt.contourf(*(list(grid(X,Y,Z2, n_grid, n_grid)) + [15]))
		plt.colorbar()
		#
		# ... and let's just dump everything into the local class dictionary:
		self.__dict__.update(locals())
		#
		return None
#
class equake_test2(object):
	# we should probalby be making inherited classes, but let's just slog through this...
	#
	def __init__(self, mc=3.0, L_r_factor=5., r_N=25, r_range=2., gridsize=.1, n_quakes=1):
		#
		# use the parkfield earthquake to do some comparisons between the new Earthquake() object from globalETAS an the older
		# one from BASScast or ETASmf. at the end, run a bunch of "local intensity" calculations for random locations,times
		# about parkfield and our original "target".
		# summary (so far): the difference between the two calculations is pretty small 
		# (mean, median, max, min) error in the exponents from the two methods are ~ .0276, .0281, .1055, 8.24e-6 respectively.
		#
		self.__dict__.update(locals())
		#
		# get a catalog around parkfield:
		d_lambda=1.76
		parkfield={'dt':dtm.datetime(2004,9,28,17,15,24, tzinfo=pytz.timezone('UTC')), 'lat':35.815, 'lon':-120.374, 'mag':5.96}
		L_r = 10.**(parkfield['mag']/2. - d_lambda)
		d_lat = L_r_factor * L_r/111.1
		d_lon = L_r_factor * L_r*math.cos(deg2rad*parkfield['lat'])/111.1
		#
		lats = [parkfield['lat']-d_lat, parkfield['lat']+d_lat]
		lons = [parkfield['lon']-d_lon, parkfield['lon']+d_lon]
		#
		t_0 = dtm.datetime(2000,1,1,0,0,0, tzinfo=pytz.timezone('UTC'))
		t_now = dtm.datetime.now(pytz.timezone('UTC'))
		etas_fit_factor=1.0
		#
		C=gep.make_ETAS_catalog(incat=None, lats=lats, lons=lons, mc=mc, date_range=[t_0, t_now], fit_factor=etas_fit_factor)
		#
		C.sort(order='mag')
		#
		target = {'time':mpd.date2num(dtm.datetime(2005,1,1, tzinfo=pytz.timezone('UTC'))), 'lon':-120.0, 'lat':35.0}
		t=target['time']
		#
		quakes_1=[]
		quakes_2=[]
		# C[::-1] creates a reversed "view" of C...
		xs = numpy.arange(lons[0], lons[1], .1)
		ys = numpy.arange(lats[0], lats[1], .1)
		Z1 = numpy.zeros((len(xs), len(ys)))
		Z2 = numpy.zeros((len(xs), len(ys)))
		for c in C[::-1][0:n_quakes]: 
			#if c['mag']>5.92 and c['mag']<5.99:
			#	pf=c
			#
			# from globalETAs:
			#eq1 = gep.Earthquake(pf, transform_type='equal_area', transform_ratio_max=2.)
			quakes_1 += [gep.Earthquake(c, transform_type='equal_area', transform_ratio_max=2.)]
			#
			# from older ETAS:
			# note: def __init__(self, mag=5.0, loc=[0.0,0.0], evtime=1.0, mc=2.0, drho=5.25, p=1.05, dmstar=1.2, dm=1.0, alpha = 2.28, eqtheta=0.0, eqeps=1.0, rtype='ssim', p_map=None)
			quakes_2 += [bcp.earthquake(mag=c['mag'], loc=[c['lon'], c['lat']], evtime=c['event_date'].tolist(), mc=mc, drho=5.25, p=c['p'], dmstar=c['dmstar']+.2, dm=c['dmstar'], alpha=2.28, eqtheta=None, eqeps=None, rtype='ssim', p_map=c['p'])]
			# BASScast earthquakes get their epsilon,theta fit during the BASScast initialization, so
			# we don't get that here automatically. for this test, we can fake it.
			quakes_2[-1].eqtheta=-45.
			quakes_2[-1].eqeps=math.sqrt(2.)
			#
			for j,lat in enumerate(ys):
				for k,lon in enumerate(xs):
					Z1[k][j]+=[quakes_1[-1].local_intensity(t=t, lon=lon, lat=lat)]
					Z2[k][j]+=[quakes_2[-1].localIntensity(inloc=[lon,lat],intime=t*days2secs)]
				
		#
		Z1 = numpy.log10(Z1)
		Z2 = numpy.log10(Z2)
		#
		plt.figure(11)
		plt.clf()
		plt.contourf(Z1, 15)
		
		plt.figure(12)
		plt.clf()
		plt.contourf(Z2, 15)
		
#
def grid(x, y, z, resX=100, resY=100):
    "Convert 3 column data to matplotlib grid"
    xi = numpy.linspace(min(x), max(x), resX)
    yi = numpy.linspace(min(y), max(y), resY)
    Z = mpl.mlab.griddata(x, y, z, xi, yi)
    X, Y = numpy.meshgrid(xi, yi)
    return X, Y, Z
