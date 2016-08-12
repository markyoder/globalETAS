import matplotlib
import numpy
import pylab as plt
import math
import random

from  matplotlib.mlab import PCA
import pca_tools as ptp
import ANSStools as atp
import datetime as dtm
import pytz

deg2rad = math.pi*2.0/360.

def parkfield_pca(L_r_factor=3.0):
	# a test example using the parkfield earthquake. let's throw in some rtree as well.
	# catfromANSS(lon=[135., 150.], lat=[30., 41.5], minMag=4.0, dates0=[dtm.datetime(2005,01,01, tzinfo=tzutc), None], Nmax=None, fout=None, rec_array=True):
	#
	d_lambda   = 1.76
	#
	parkfield={'dt':dtm.datetime(2004,9,28,17,15,24, tzinfo=pytz.timezone('UTC')), 'lat':35.815, 'lon':-120.374, 'mag':5.96}
	L_r = 10.**(parkfield['mag']/2. - d_lambda)
	d_lat = L_r_factor * L_r/111.1
	d_lon = L_r_factor * L_r*math.cos(deg2rad*parkfield['lat'])/111.1
	
	print("d_lat, d_lon: ", d_lat, d_lon)
	#
	x_pf=parkfield['lon']
	y_pf=parkfield['lat']
	parkfield_cat_prams = {'lon':[x_pf-d_lon, x_pf+d_lon], 'lat':[y_pf-d_lat, y_pf+d_lon], 'minMag':1.5, 'dates0':[dtm.datetime(2004,9,28, tzinfo=pytz.timezone('UTC')), dtm.datetime(2010,9,28, tzinfo=pytz.timezone('UTC'))], 'Nmax':None, 'fout':None, 'rec_array':True}
	#
	cat = atp.catfromANSS(**parkfield_cat_prams)
	#
	# i still don't get what this does...
	#my_pca = PCA(numpy.array(zip(cat['lon'], cat['lat'])))
	my_pca = ptp.yoda_pca(list(zip(cat['lon'], cat['lat'])))		# returns (eig_vals, eig_vecs)
	e_vals = my_pca[0]
	e_vecs = numpy.array(my_pca[1])
	#
	#e_vals_n = e_vals/min(e_vals)
	e_vals_n = numpy.array([min(4.0, x/min(e_vals)) for x in e_vals])
	#
	print("e_vecs:", e_vecs[0][0], e_vecs[0][1], e_vecs[1][0], e_vecs[1][1])
	circle_xy = simple_circle(x=x, y=y, r=L_r*L_r_factor/111.1)

	#T = numpy.array([[e_vals_n[j]*x for x in rw] for j,rw in enumerate(e_vecs)])
	#circle_xy_prime = numpy.dot(circle_xy,T.transpose())	# note: this syntax will add [x,y] to all members of circle_xy_prime like [[a+x,b+y], [a+x,b+y],...]
	
	T = numpy.dot([[e_vals_n[0],0.],[0., e_vals_n[1]]], e_vecs)
	circle_xy_prime = numpy.dot(circle_xy, T)
	
	#circle_xy_prime = numpy.dot(circle_xy, [[e_vals_n[0],0.],[0., e_vals_n[1]]])
	#circle_xy_prime = numpy.dot(circle_xy_prime, e_vecs)
	
	#circle_xy = numpy.array(circle_xy)+numpy.array([x,y])
	circle_xy_prime = [[j+x_pf, k+y_pf] for j,k in circle_xy_prime]
	circle_xy = [[j+x_pf, k+y_pf] for j,k in circle_xy]
	#print "T: ", T
	#
	# a rotational transformation:
	#  elliptical distribution; a = r, b = {something < a}. this is distinct from the equal-area transform, in which a>r.
	#  note that for this transform, the initial rate-density is adjusted to the new (surface projection) area.
	#mu_x, mu_y = [numpy.mean(col) for col in zip(*circle_xy)]
	#print "means: ", mu_x, mu_y
	#circle_xy_prime = [[rw[0]-mu_x, rw[1]-mu_y] for rw in circle_xy]
	#
	#circle_xy_prime = [[numpy.dot(rw,e_vecs[0])*e_vals_n[0], numpy.dot(rw, e_vecs[1])*e_vals_n[1]] for rw in circle_xy_prime]
	#circle_xy_prime = numpy.dot([[rw[0] + mu_x, rw[1]+mu_y] for rw in circle_xy_prime], zip(*e_vecs))
	#circle_xy_prime = numpy.dot(circle_xy_prime, e_vecs.transpose())
	#circle_xy_prime = numpy.dot(circle_xy_prime, e_vecs)
	#circle_xy_prime = [[j+mu_x, k+mu_y] for j,k in circle_xy_prime]

	#circle_xy_prime = numpy.dot(circle_xy_prime, zip(*e_vecs))
	#
	plt.figure(0)
	plt.clf()
	plt.plot(cat['lon'], cat['lat'], '.', label='parkfield cat')
	plt.plot([x_pf], [y_pf], 'r*', ms=15)
	plt.legend(loc=0, numpoints=1)
	#
	plot_axes = (L_r/111.2)*numpy.array([e_vecs[0], [0.,0.], e_vecs[1]]) + numpy.array([x_pf, y_pf])
	plt.plot(*zip(*plot_axes), ls='-', marker='o')
	#
	Lry = abs(L_r_factor * L_r/111.1)
	Lrx = abs(Lry*math.cos(y_pf*deg2rad))
	
	#
	Wts = [xx/max(e_vals) for xx in e_vals]
	print("Wts: ", Wts)
	#e_vecs = e_vecs.transpose()
	
	
	#plt.plot([x, x + Lrx*e_vecs[0][0]], [y, y + Lry*e_vecs[0][1]], ls='-', marker='o', color='r')
	#plt.plot([x, x + Lrx*e_vecs[1][0]], [y, y + Lry*e_vecs[1][1]], ls='-', marker='^', color='m')
	
	'''
	#<<<<<<< HEAD
	plt.plot([x, x + Wts[0]*Lry*e_vecs[0][0]], [y, y + Wts[0]*Lry*e_vecs[0][1]], ls='-', marker='o', color='r')
	plt.plot([x, x + Wts[1]*Lry*e_vecs[1][0]], [y, y + Wts[1]*Lry*e_vecs[1][1]], ls='-', marker='^', color='m')
	plt.plot(*list(zip(*circle_xy)), ls='-', marker='', lw=2.)
	plt.plot(*list(zip(*circle_xy_prime)), ls='--', color='r', marker='', lw=1.5, alpha=.7, zorder=11)
	#=======
	'''
	#ax1 = numpy.dot(
	
	plt.plot([x_pf, x_pf + Wts[0]*Lry*e_vecs[0][0]], [y_pf, y_pf + Wts[0]*Lry*e_vecs[0][1]], ls='-', marker='o', color='r')
	plt.plot([x_pf, x_pf + Wts[1]*Lry*e_vecs[1][0]], [y_pf, y_pf + Wts[1]*Lry*e_vecs[1][1]], ls='-', marker='^', color='m')
	
	plt.plot(*zip(*circle_xy), ls='-', marker='', lw=2.)
	plt.plot(*zip(*circle_xy_prime), ls='--', color='r', marker='', lw=1.5, alpha=.7, zorder=11)	
	
	#plt.plot([x, x+Lrx*e_vecs[0][0]], [y, y + Lry*e_vecs[1][0]], ls='-', marker='o', color='r')
	#plt.plot([x, x+Lrx*e_vecs[0][1]], [y, y + Lry*e_vecs[1][1]], ls='-', marker='^', color='m')
	
	plt.figure(1)
	plt.clf()
	#plt.plot([0., Lrx*e_vecs[0][0]], [0., Lry*e_vecs[0][1]], '.-')
	#plt.plot([0., Lrx*e_vecs[1][0]], [0., Lry*e_vecs[1][1]], '.--')
	#plt.plot([0., e_vecs[0][0]], [0., e_vecs[0][1]], '.-')
	#plt.plot([0., e_vecs[1][0]], [0., e_vecs[1][1]], '.--')
	plt.plot(*zip([0.,0.],e_vecs[0]), marker='.', ls='--')
	plt.plot(*zip([0.,0.],e_vecs[1]),marker='^', ls='-')
	
	plt.plot(*zip([0.,0.],T[0]), marker='.', ls='--')
	plt.plot(*zip([0.,0.],T[1]),marker='^', ls='-')
	
	
	#
	return my_pca
	

def lat_lon_to_xy(lat=0., lon=0.):
	return {'x':lon*111.1/math.cos(lat), 'y':lat*111.1}

def xy_to_lat_lon(x,y):
	return {'lat':y/111.1, 'lon':x*math.cos(lat)/111.1}

def simple_circle(x=0., y=0., r=0., n_points=250):
	#
	d_theta = 2.0*math.pi/float(n_points)
	return [[x+r*math.cos(theta), y+r*math.sin(theta)] for theta in numpy.arange(0., math.pi*2.0+d_theta, d_theta)]
	
