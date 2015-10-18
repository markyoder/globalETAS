import globalETAS as gep
import BASScast as bcp

import numpy
import math
import pylab as plt
import datetime as dtm
import pytz
import random

import matplotlib.dates as mpd

days2secs = 60.*60.*24.
year2secs = 60.*60.*24.*365.
deg2km=111.31
deg2rad = 2.0*math.pi/360.
alpha_0 = 2.28

def equake_test1(mc=3.0, L_r_factor=5., r_N=25, r_range=2.):
	#
	# use the parkfield earthquake to do some comparisons between the new Earthquake() object from globalETAS an the older
	# one from BASScast or ETASmf. at the end, run a bunch of "local intensity" calculations for random locations,times
	# about parkfield and our original "target".
	# summary (so far): the difference between the two calculations is pretty small 
	# (mean, median, max, min) error in the exponents from the two methods are ~ .0276, .0281, .1055, 8.24e-6 respectively.
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
	#
	target = {'time':mpd.date2num(dtm.datetime(2005,1,1, tzinfo=pytz.timezone('UTC'))), 'lon':-120.0, 'lat':35.0}
	#
	z1 = eq1.local_intensity(t=target['time'], lon=target['lon'], lat=target['lat'])
	print "target: ", target
	print("z1: ", z1)
	print target
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
	for lon,lat in rlons_lats:
		t = target['time'] + (R3.random()-.5)*50.		# days
		print "t: ", t
		
		z1 = math.log10(eq1.local_intensity(t=t, lon=lon, lat=lat))
		z2 = math.log10(eq2.localIntensity(inloc=[lon,lat] , intime=t*days2secs))
		#
		print("(%f,%f,%f) :: z_1: %f, z_2: %f (%f)" % (lon, lat, t, z1, z1, 2.*(z1-z2)/(z1+z2)))
		err = 2.*abs((z1-z2)/(z1+z2))
		mean_err+=err
		errs += [err]
		j+=1
	#
	print "mean_err: ", mean_err/float(j), numpy.median(errs), max(errs), min(errs)
	plt.figure(11)
	plt.clf()
	plt.plot(pf['lon'], pf['lat'], 'r*', ms=18, zorder=4)
	plt.plot(*zip(*rlons_lats), marker='.', ls='', color='b',zorder=3)
		
	#
	return {'eq11':eq1, 'eq2':eq2, 'errs':errs}
