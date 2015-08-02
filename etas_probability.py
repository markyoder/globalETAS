
import pytz
import os
import math
import numpy
import scipy
import scipy.special
import itertools


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

plt.ion()

class Earthquake(object):
	#
	def __init__(self, m=None, mc=3.0, t0=1.0, p=1.05, D=1.5, q=1.5, dm=1.0):
		self.m=m
		self.mc=mc
		self.t0=t0
		self.p=p
		self.D=D
		self.q=q
		self.dm=dm
		#
		self.tau = self.get_tau()
	
	def get_tau(self, t0=None, p=None, Nom=None, m=None, mc=None, dm=None):
		#
		if Nom==None: 
			m  = (m or self.m)
			mc = (mc or self.mc)
			dm = (dm or self.dm)
			#
			Nom=self.n_omori(m=m, mc=mc, dm=dm)
		#
		t0 = (t0 or self.t0)
		p =  (p or self.p)
		#
		return (t0**(1.-p))/((p-1.)*Nom)
	#
	def omori_rate(self, t=0., t0=None, tau=None, p=None):
		t0=(t0 or self.t0)
		tau=(tau or self.tau)
		p=(p or self.p)
		#
		return 1.0/(tau*(t0+t)**p)
	#
	def omori_N(self, t2, t1=0., t0=None, tau=None, p=None):
		'''
		# integrated omori rate; expected count of earthquakes between t1,t2.
		'''
		t0=(t0 or self.t0)
		tau=(tau or self.tau)
		p=(p or self.p)
		#
		if p==1.0:
			return (1.0/tau)*numpy.log((t0+t2)/(t0+t1))
		else:
			return -((t0+t2)**(1.-p) - (t0 + t1)**(1.-p))/(tau*(p-1.))
		#
	#
	def n_omori(self, m=None, mc=None, dm=1.0):
		if m==None:  m  = self.m
		if mc==None: mc = self.mc
		if dm==None: dm = self.dm
		#
		return 10.0**(m-dm-mc)
	#
	def omori_prob(self, m_target=None, mc=None, delta_m=0., t0=None, t2=None, t1=None, p=1.05):
		'''
		# calc probability of one or more m>m_0 earthquakes based on Omori statistics. use NHHP process;
		# this ends up being an integration of the rate from t1 to t2 in the NHHP distribution. note that this
		# can be shown to be the right approach using baysean logic: P(a|b) = P(b|a)P(a)/P(b), where P(a) is
		# the probability of (not) an event up until t2; P(b) is prob. of (not) an event until t1, and P(b|a)=1.
		'''
		#
		# now code it up...
		pass
	#
	def observed_omori_prob(self, rate_in, m_target=None, mc=None, delta_m=0., t0=None, delta_t=0.0, p=1.05):
		'''
		# estimate probability for an observed seismicity rate. if t0!=None, use an Omori model to estimate decaying rate
		# use m_target, m_c to adjust the rate to account for magnitudes of interest (this can also be done at the input)
		# in this case, t0 is the current/starting time, and equivalently t0 for Omori.
		#
		# note: this function does not use any class variables, so it can be moved outside the class scope.
		# we can still use the Earthquake() class member functions; do we need to declare an instance to do so?
		''' 
		# def get_tau(self, t0=None, p=None, Nom=None, m=None, mc=None, dm=None)
		if m_target!=None and mc!=None:
			delta_m = m_target-mc
		delta_m = (0.0 or delta_m)
		#
		rate_factor = 10**(-delta_m)
		effective_rate = rate_factor*rate_in
		print "eff_rate: ", effective_rate
		#
		if t0==None or t0<=0.:
			# return a stationary Poisson probability.
			return 1.0 - numpy.exp(-rate_in*delta_t)
		#
		# otherwise, use an Omori-NHHP model.
		# first solve for tau in Omori:
		tau = (t0**(-p))/effective_rate
		N   = self.omori_N(t2=delta_t, tau=tau, t1=0., t0=t0, p=p)
		print "N: ", N
		#
		return 1.0 - numpy.exp(-N)

def etas_prob_test(r0=500, p_min=.01, p_max=1.5, p_step=.01, t0=0., t1=0., t2=1000., dt2=1.):
	#
	# trying to figure out how to calc an event probability from observed etas rates.
	# start with etas; assume p>1, and a Non-homogeneous Poisson Process (NHPP).
	# so, P = 1 - exp(ingegral[r(t)] ), and r(t) is Omoir.
	#
	# in principle, we get R = int(r) = r0*(p-1)*[t0 + t2]**(1-p) - (t0+t1)**(1-p)]
	# and the idea will be to assume t1=0 and t0 -> 0,
	# so
	# R -> r0*(p-1)*t2^(1-p)
	#
	
	f=plt.figure(0)
	plt.clf()
	ax3d = f.add_subplot(111, projection='3d')
	#
	X = numpy.arange(.001, t2, dt2)
	Y = numpy.arange(p_min, p_max, p_step)
	#Z = R_omori(r0=r0, p=Y, t0=t0, t1=0., t2=X)
	#
	XYZ = [[x,y,R_omori(r0=r0,p=y, t0=t0, t1=0.,t2=x)] for x,y in itertools.product(X,Y)]
	
	# let's do R(p,t)
	
	ax3d.plot(zip(*XYZ)[0],zip(*XYZ)[1], zip(*XYZ)[2],  '.')
	ax3d.set_xlabel('time $t$')
	ax3d.set_ylabel('Omori scaling exponent $p$')
	ax3d.set_zlabel('Integrated Omori Rate $R$')
	#
	#
	f2=plt.figure(1)
	plt.clf()
	ax3d2 = f2.add_subplot(111, projection='3d')
	XYZ = [[x,y,poisson_cum_R(r0=r0,p=y, t0=t0, t1=0.,t2=x)] for x,y in itertools.product(X,Y)]
	
	# let's do R(p,t)
	
	ax3d2.plot(zip(*XYZ)[0],zip(*XYZ)[1], zip(*XYZ)[2],  '.')
	ax3d2.set_xlabel('time $t$')
	ax3d2.set_ylabel('Omori scaling exponent $p$')
	ax3d2.set_zlabel('Poisson Probability based on $R$')
#
def basic_exp_prob(t, t1=0.0, r=1.0):
	'''
	# basic exponential (Poisson, k=1) probability.
	# t is the independent variable, t1 is the initial time (value of the ind. variable), r is the rate, so P = exp(-r*t1) - exp(-r*t2)
	'''
	#
	return exp(-r*t1) - exp(-r*t2)

def omori_prob(t2, t1=0., t0=1., tau=1., p=1.0, c=None):
	# note: the NHHP solution for Omori over time interval t -> t+\Delta t is (for p=1, so int(f) = ln(f):
	# F(\Delta t, t) = ( (c + t + \Delta t)/(c + t) ) ** (-c/\tau)
	# and c = t_0, tau' = tau*t_0**p  (aka, initial interval or rate or whatever).
	# if c is given, assume the (1+t/c) format:
	if c!=None:
		t0=c
		tau *= t0**(-p)
	#
	if p==1.0:
		# use logarithmic value:
		F = ((t0 + t2 + (t2-t1))/(t0+t1))**(-1./tau)
	elif p!=1.0:
		F = numpy.exp(((t0 + t2)**(1.-p) - (t0+t1)**(1.-p))/(tau*(p-1.)))
	#
	return 1.0-F

##################
# modifications to rate/spatial density to account for large aftershocks falling outside the rupture area. eventually, this will require some sort of
# renormalization... or maybe just normalization, to not get into an argument about language.
def omori_demos(X=None, x1=1.0, x0=1.0, chi=1.0, q=1.5, y_min=None, **kwargs):
	if X==None: X = numpy.arange(x0/1000., x0*20., x0/1000.)
	#
	kwargs['x_scale'] = kwargs.get('x_scale', 'log')
	kwargs['y_scale'] = kwargs.get('y_scale', 'log')
	#
	plt.figure(0)
	plt.clf()
	ax = plt.gca()
	ax.set_xscale(kwargs['x_scale'])
	ax.set_yscale(kwargs['y_scale'])
	#
	Y_omori = f_omori(X=X, x0=x0, chi=chi, q=q)
	if y_min==None: y_min = min(Y_omori)
	#y_min = 0.
	#
	plt.plot(X, Y_omori, '-', label='omori')
	Y_exp = f_omori_exp(X=X, x0=x0, x1=x1, chi=chi, q=q, y_min=y_min)
	plt.plot(X, Y_exp, '-', label='omori_exp')
	plt.plot(X, f_omori_inv_gamma(X=X, x0=x0, x1=x1, chi=chi, q=q, y_min=y_min), '-', label='omori_gamma')
	#
	plt.plot([x0, x0], [min(Y_exp), 1.1], '--')
	
	plt.legend(loc=0, numpoints=1)

def f_omori(X, x0, chi, q):
	# can be used for reguar space or time omori.
	return (1.0/chi)*(x0+X)**-q
	
def f_omori_exp(X, x1=1.0, x0=1.0, chi=1.0, q=1.5, y_min=0.):
	'''
	# omori-exponential rate distribution
	# x1: exponential factor (exp(-x/x1))
	# x0: omori factor (1/(x0+x))
	# omori function is:
	# f_omori = (1/chi)(x0 + x)**-q
	#
	# this is a simple way to "hollow out" the omiri distribution, but f_omori_inf_gamma() is probably more appropriate (for what reason?), though -- again,
	# since we really don't know anything about r<L_r, at this point there is a lot of operator discretion. note also 	# that when we integrate this function,
	# we normalize with a gamma function.
	'''
	#
	f_in =  lambda r: 1.0 - numpy.exp(-r/x1)
	f_out = lambda r: ((x0 + r)**(-q))/chi
	f = lambda r: f_in(r)*f_out(r) + y_min*(x1-r)*(r<=x1)
	#
	return f(numpy.array([f(x) for x in X]))

def f_omori_inv_gamma(X, x0=1.0, x1=None, chi=1.0, q=1.5, y_min=0.):
	'''
	# this is the "inverse gamma" distribution; see Malamud et al. 2004, "LANDSLIDE INVENTORIES AND THEIR STATISTICAL PROPERTIES"
	# this may be a "proper" or expected distribution for this sort of work. it seems to fit well with landslide data, but f_omori_exp() might fit well too.
	# since we don't necessarily understand the r<L_r domain, and since both -> regular_omori, it shouldn't matter too much what we use.
	# see wikipedia: https://en.wikipedia.org/wiki/Inverse-gamma_distribution
	'''
	if x1==None: x1=chi
	#
	alpha = q-1.0
	#return ((1.0/chi)**alpha)*(1.0/scipy.special.gamma(alpha))*((x0+X)**(-q))*numpy.exp(-x1/X)
	return [((1.0/chi)**alpha)*(1.0/scipy.special.gamma(alpha))*((x0+x)**(-q))*numpy.exp(-x1/x) + y_min*(x1-x)*(x<=x1) for x in X]

def F_omori_exp(x, x1=1.0, x0=1.0, chi=1.0, q=1.5):
	'''
	# cumulative omori-exponential (aka, integrated f_omori_exp)
	#  - number of events in an omori-exponential process, probably used for Non-homogeneous Poisson calculations.
	#  - this might warrant some cleaning up, namely to include an x1, x_final range, so we don't always integrate from zero.
	#
	# probability density of omori-exponential distribution
	# x1: exponential factor (exp(-x/x1))
	# x0: omori factor (1/(x0+x))
	# omori function is:
	# f_omori = (1/chi)(x0 + x)**-q
	'''
	#
	if not hasattr(x, '__len__'): x=numpy.array([x])
	# for readability, break this out in components, then return product... or comment and combine.
	f1 = ((x0 + x)**(-q))/(chi*(q-1.0))
	f2 = (q-1.0)*x1*math.exp(x0/x1)*(((x0+x)**q)/x1)
	f3 = scipy.special.gamma((1.0-q), (x0+x)/x1)
	#
	f_diff = ((x0)**(-q))/(chi*(q-1.0)) * ((q-1.0)*x1*math.exp(x0/x1)*(((x0)**q)/x1)*scipy.special.gamma(1.0-q, numpy.array([x0/x1])) - x0 )
	#
	# ... though i think we'll just get rid of this. since min_y applies to a rate, it's ok if it -> 0.
	min_y_correction = (x<x1)*.5*y_min*(x0**2. - (x0-x)**2.)
	#
	return f1*(f2*f3 - x0 - x) - f_diff + min_y_correction
	#return f1*(f2*f3 - x0 - x)
#
def F_omori_inv_gamma(x, x0=1.0, x1=None, chi=1.0, q=1.5, y_min=0.):
	'''
	#
	'''
	if x1==None: x1 = chi
	#
	alpha = q-1.0
	
	F = scipy.special.gamma(alpha, (1./chi)*x)/scipy.special.gamma(alpha)
	#
	# min_y correction(s):
	min_y_correction = (x<x1)*.5*y_min*(x0**2. - (x0-x)**2.)
	
	#
	return F + min_y_correction
	#
#
def big_mag_distribution(r0=1.0, r1=1.0, chi=1.0, q=1.5, r_min=0., r_max=100., nits=1000, fnum=0, x_scale='log', y_scale='log'):
	'''
	# plots for "Big aftershocks are farther away" distribution 
	'''
	f_in =  lambda r: 1.0 - numpy.exp(-r/r1)
	f_out = lambda r: ((r0 + r)**(-q))/chi
	f = lambda r: f_in(r)*f_out(r)
	#
	X = numpy.arange(r_min, r_max, float(r_max)/float(nits))
	plt.figure(fnum)
	plt.clf()
	ax=plt.gca()
	#
	ax.set_yscale(x_scale)
	ax.set_xscale(y_scale)
	#
	ax.plot(X, f_in(X), '-', label='exponential P')
	ax.plot(X, f_out(X), '-', label='Omori')
	ax.plot(X, f(X), label='product')
	plt.legend(loc=0, numpoints=1)
	#
	

def R_omori(r0=0., p=1.1, t0=0., t1=0., t2=0.):
	# these throwing errors because 0.0 cannot be raised to negative power, so just trap this case.
	#return r0*(p-1.)*((t0+t2)**(1.-p) - (t0+t1)**(1.-p))
	return r0*(p-1.)*( ((t0+t2) if (t0+t2)!=0.0 else 1e-10)**(1.-p) - (((t0+t1) if (t0+t1)!=0.0 else 1e-10)**(1.-p)))

def poisson_cum_R(r0=0., p=1.1, t0=0., t1=0., t2=0.):
	# there's a better way to do this than **locals(), but since we've not set any variables yet, it will work.
	return 1.0 - numpy.exp(R_omori(**locals()))
