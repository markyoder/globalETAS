import matplotlib
import numpy
import pylab as plt
import math
import random
#
from  matplotlib.mlab import PCA
#
def pca1(theta=math.pi/6., a=1.0, b=.5, x0=0., y0=0., N=1000):
	#
	Rx = random.Random()
	Ry = random.Random()
	#
	# make a blob of data; the PCA should find the center, axes, etc.
	#
	#XY = [rotate_ccw([x0 + a*Rx.random(), y0 + b*Ry.random()], theta, x0=x0, y0=y0) for n in xrange(N)]
	XY = [[x0 + a*Rx.random(), y0 + b*Ry.random()] for n in range(N)]
	print(("variances on raw matrix: ", numpy.var(zip(*XY)[0]), numpy.var(zip(*XY)[1])))
	print(("std on raw matrix: ", numpy.std(zip(*XY)[0]), numpy.std(zip(*XY)[1])))
	XY = [rotate_ccw(rw, theta, x0=x0, y0=y0) for rw in XY]
	# first, let's just visualize the rotation matrix:
	#
	pcas = PCA(numpy.array(XY))
	#
	# and do it manually:
	x_mean, y_mean = numpy.mean(zip(*XY)[0]), numpy.mean(zip(*XY)[1])
	my_mu = numpy.array([x_mean, y_mean])
	print(("x_mean, y_mean: ", x_mean, y_mean))		# these check with pcas.mu
	dXY = [[(x-x_mean), (y-y_mean)] for x,y in XY]
	my_cov = numpy.dot(numpy.array(list(zip(*dXY))),numpy.array(dXY))/(float(N)-1.)
	# use eigh() for symmetric or hermitian matrices. it's faster and avoids some round-off errors that can result in complex/imaginary valued elements.
	# ... and sort by eigen-values
	eig_vals, eig_vecs = numpy.linalg.eigh(my_cov)		# ... and there's something weird about eig() vs eigh() they return somewhat different eigen values.
	#eig_vecs = numpy.array(zip(*eig_vecs))	# transpose so each 'row' is an eigen vector.
	eig_vals, eig_vecs = list(zip(*sorted(zip(eig_vals, eig_vecs.transpose()), key=lambda x:x[0])))		# and we can use reverse=True, or just remember that we're reverse order sorted.
	eig_vals_norm = numpy.linalg.norm(eig_vals)
	#axis_weights = [math.sqrt(e/max(eig_vals)) for e in eig_vals]		# i think this is pretty standard, but it seems that normalizing these vectors so that e1*82 + e2**2 = 1
	#																						# also makes sense.
	axis_weights = [e/eig_vals_norm for e in eig_vals]
	print(("normed(eig_vals): %f (normed-check (should be 1.0): %f)" % (eig_vals_norm, numpy.linalg.norm(axis_weights))))
	#	
	print(("eigs: ", eig_vals, eig_vecs))
	frac_x, frac_y = [x/(numpy.linalg.norm(eig_vals)**2.) for x in eig_vals]
	#
	# primed axes?:
	x_prime = numpy.dot(numpy.array(eig_vecs), numpy.array([1., 0.]))*frac_x**.5 #+ numpy.array([x_mean, y_mean])
	y_prime = numpy.dot(numpy.array(eig_vecs), numpy.array([0., 1.]))*frac_y**.5 #+ numpy.array([x_mean, y_mean])
	#
	plt.figure(0)
	plt.clf()
	#
	plt.plot(zip(*XY)[0], zip(*XY)[1], 'b.')
	plt.plot(zip(*pcas.a)[0], zip(*pcas.a)[1], 'r.')
	print(("means again: ", x_mean, y_mean))
	plt.plot(x_mean, y_mean, 'co', ms=9, zorder=4)
	#plt.plot(*pcas.mu, color='r', ms=5, marker='o',zorder=5)
	#
	#plt.plot(numpy.array([0., x_prime[0]]) + x_mean, numpy.array([0., x_prime[1]]) + y_mean, 'o-')
	#plt.plot(numpy.array([0., y_prime[0]]) + x_mean, numpy.array([0., y_prime[1]]) + y_mean, 's-')
	#
	
	#
	ax_fact = -.5*axis_weights[0]
	plt.plot([x_mean, x_mean+eig_vecs[0][0]*ax_fact], [y_mean, y_mean+eig_vecs[1][0]*ax_fact], 'ro-', zorder=5)
	ax_fact = -.5*axis_weights[1]
	plt.plot([x_mean, x_mean+eig_vecs[0][1]*ax_fact], [y_mean, y_mean+eig_vecs[1][1]*ax_fact], 'gs-')
	#
	plt.figure(1)
	plt.clf()
	plt.plot(*list(zip(*pcas.Y)), marker='.', ls='', color='g')

	#return pcas
	return [eig_vals, eig_vecs]

class PCA_transform(object):
	def __init__(self, data_in=None, N=None, theta=None):
		if data_in==None:
			N=(1000 or N)
			theta = (theta or math.pi/6.)
			#
			data_in = make_test_data(theta=theta, N=100, x0=0., y0=0., a=2., b=1.)
		#
		self.data = data_in
		self.calc_pca(data_in=data_in, do_return=False)
	
	def calc_pca(self, data_in=None, do_return=False):
		'''
		# calc pca and assign class scope variables.
		# assume data are in rows, X = [[x0,x1,x2,x3], [x0,x1,x2,x3],...], or for example: X = [[x0,y0,z0], [x1,y1,z1], ...]
		#
		'''
		data = (data_in or self.data)
		#
		# first, get means:
		mus = [numpy.mean(col) for col in zip(*data)]		# mean values for each column. note, this will also constitute a vector 
															# between our transoformed origin and the origina origin.
		data_centered = [[x-mus[j] for j,x in enumerate(rw)] for rw in data]
		cov_normed = numpy.dot(numpy.array(list(zip(*data_centered))),numpy.array(data_centered))/(float(len(data_centered))-1.)
		#
		# get eigen-values/vectors:
		# note: for a symmetric or hermitian matrix, we should use numpy.linalg.eigh(); it's faster, more accurate, and will quash complex value rounding errors.
		#       it might be necessary, down the road, to just spin through these and dump any complex components.
		eig_vals, eig_vecs = numpy.linalg.eig(cov_normed)
		#
		if data_in==None or data_in==self.data:			
			#
			self.max_eig = max(eig_vals)
			self.max_eig_vec_index = eig_vals.tolist().index(self.max_eig)		# this might be a stupid way to do this;
			self.axis_weights = [math.sqrt(e/self.max_eig) for e in eig_vals]
			#
			# assign stuff.
			self.mus = numpy.array(mus)
			self.data_centered = data_centered
			self.eig_vals = eig_vals
			self.eig_vecs = eig_vecs
			#
			#self.axis_weights = axis_weights
			#self.max_eig = max_eig
			#self.max_eig_vec_index = max_eig_vec_index
			self.cov_normed =cov_normed
			#
			self.to_PCA_rotation = numpy.array([[x*self.axis_weights[j] for j,x in enumerate(rw)] for rw in eig_vecs])
		#
		if do_return: return [eig_vals, eig_vecs]
		#
	def to_PCA(self, v):
		#x0 = self.M_to_PCA.dot(v)
		#x1 = self.mus
		#print "adding: ", x0, x1, x0+x1
		return numpy.array(self.to_PCA_rotation.dot(v)) + numpy.array(self.mus)
	#
	@property
	def primary_axis(self):
		return zip(*self.to_PCA_rotation)[self.max_eig_vec_index]
	#
	def plot_data(self, fignum=0, axes=[0,1]):
		'''
		# some diagnostic plotting...
		# for now, just 2d plotting.
		'''
		#
		plt.figure(fignum)
		plt.clf()
		#
		plt.plot(zip(*self.data)[axes[0]], zip(*self.data)[axes[1]], marker='.', color='b', zorder=2, ls='')
		#
		#print "vx: ", self.to_PCA([1.,0.])
		vx = numpy.array([self.mus, self.to_PCA([1.,0.])])
		vy = numpy.array([self.mus, self.to_PCA([0.,1.])])
		
		print(("vx, vy: ", vx, vy))
		#
		plt.plot(*self.mus, marker='*', color='r', ms=15)
		plt.plot(*list(zip(*vx)), color='r', lw=1.5, ls='-', marker='s')
		plt.plot(*list(zip(*vy)), color='g', lw=1.5, ls='-', marker='s')
		#
		#[plt.plot(*(self.mus + v), color='m', ls='--', marker='*') for v in zip(*self.eig_vecs)]		# needs to be corrected to use axes[]
		[plt.plot(*(self.mus + v), color='m', ls='--', ms=10, marker='*') for v in zip(*self.to_PCA_rotation)]
		
		
		
def pca_test2(theta=math.pi/6., N=1000, x0=0., y0=0., fignum=0):
	#
	my_data = make_test_data(theta=theta, N=N, x0=x0, y0=y0)
	my_pca = yoda_pca(my_data)
	#
	print(("my_pca: ", my_pca))
	#return my_pca
	ax_x, ax_y = list(zip(*my_pca[1]))
	axis_weights = [math.sqrt(e/max(my_pca[0])) for e in my_pca[0]]
	ax_x = numpy.array(ax_x)*axis_weights[0]
	ax_y = numpy.array(ax_y)*axis_weights[1]
	#
	plt.figure(fignum)
	plt.clf()
	plt.plot(*list(zip(*my_data)), marker='.', color='b', zorder=2, ls='')
	plt.plot(*list(zip([0,.0], ax_x)), marker='s', ls='-', lw=2, color='r')
	plt.plot(*list(zip([0,.0], ax_y)), marker='s', ls='-', lw=2, color='g')
	#
	xprime = my_pca[1].dot([1.,0.]) + numpy.array([.5,0.])
	yprime = my_pca[1].dot([0.,1.]) + numpy.array([.5,0.])
	#
	plt.plot(*list(zip([.5,0.], xprime)), marker='s', ls='--', lw=2, color='r')
	plt.plot(*list(zip([.5,0.], ypprime)), marker='s', ls='--', lw=2, color='g')
#
#
def yoda_pca(data_in):
	# we'll rename this later. for now, this is just a super simple PCA approach to finding principal axes. we'll leave them in their original order
	# and leave them all in tact. it's basically a fitting algorithm so we can construct a transformation matrix between two frames.
	# so, basically, do normal PCA (mean value, eigen-vals/vecs, etc.
	# assume data are like [[x0, x1, x2, x3], [x0,x1,x2,x3],...]
	#
	mus = [numpy.mean(col) for col in zip(*data_in)]
	#print(mus)
	centered = [[x-mus[j] for j,x in enumerate(rw)] for rw in data_in]
	my_cov = numpy.dot(numpy.array(list(zip(*centered))),numpy.array(centered))/(float(len(data_in))-1.)
	#
	# for now, we can't assume hermitian or symmetric matrices, so use eig(), not eigh()
	eig_vals, eig_vecs = numpy.linalg.eig(my_cov)
	#axis_weights = [math.sqrt(e/max(eig_vals)) for e in eig_vals]
	#
	# note: the transformed axes are the column vectors in eig_vecs; the main axis of interest is the one with the largest eigen vector, but we also want to preserve the
	# geometry of the system.
	#
	
	#
	#return centered
	return [eig_vals, eig_vecs]

def rotate_ccw(v,theta, x0=0., y0=0.):
	# rotates a vector about it's tail. if no tail is given, vector is in [x,y] format, we effectively assume vector at origin.
	# note that this is not a good generalized format; if we input v like v=[[0.,0.], [0.,1]], it returns a v'=[.866,.5] (aka, drops the tail).
	# however, it will be a useful format for these little testing scripts.
	#
	tail = [x0, y0]
	tip = v
	if hasattr(v[0], '__len__'):
		# we have a vector-vector, aka, not from origin. assume we're rotating about the tail
		tail = v[0]
		tip  = [v[-1][j] - v[0][j] for j in range(len(v[0]))]
		#
		print(("tip, tail: ", tip, tail))
	v_prime = numpy.dot([[math.cos(theta), -math.sin(theta)],[math.sin(theta), math.cos(theta)]], tip)
	#
	return [x+tail[j] for j,x in enumerate(v_prime)]

def make_test_data_gen(thetas=[math.pi/6.], N=100, X0=[0., 0.], delta_X=[2.,1.]):
	# eventually, interpret thetas to be a rotation in the [x_j, x_j+1] plane about [x_j+2]. but we'll get to that later. for now, just produce a generalized arary.
	#
	while len(delta_X)<len(X0): delta_X += [1.]
	while len(thetas)<(len(X0)-1): thetas+=0.
	#
	Rs = [random.Random() for x in X0]
	XY = [[x0 + dx*r.random() for x0,dx,r in zip(X0, delta_X, Rs)] for n in range(N)]
	#
	for k,theta in enumerate(thetas):
		for rw in XY:
			xprime, yprime = rotate_ccw(v=[rw[k], rw[k+1]], theta=theta, x0=X0[k], y0=X0[k+1])
			#
		#
	#
	return XY

def make_test_data(theta=math.pi/6., N=100, x0=0., y0=0., a=2., b=1.):
	#
	# for now, limit to 2D.
	Rx,Ry = [random.Random() for x in [0,1]]
	#
	#XY = [rotate_ccw([x0 + a*Rx.random(), y0 + b*Ry.random()], theta, x0=x0, y0=y0) for n in xrange(N)]
	XY = [[x0 + a*Rx.random(), y0 + b*Ry.random()] for n in range(N)]
	#print "variances on raw matrix: ", numpy.var(zip(*XY)[0]), numpy.var(zip(*XY)[1])
	#print "std on raw matrix: ", numpy.std(zip(*XY)[0]), numpy.std(zip(*XY)[1])
	XY = [rotate_ccw(rw, theta, x0=x0, y0=y0) for rw in XY]
	#
	return XY

def lzip(X):
	return list(zip(X))

if __name__=='__main__':
	pass
else:
	plt.ion()	
