'''
# some generic ROC tools, and playing a bit with shared memory mpp as well.
'''
###


#import datetime as dtm
#import matplotlib.dates as mpd
#import pytz
#tzutc = pytz.timezone('UTC')

#import operator
import math
import random
import numpy
import scipy
#import scipy.optimize as spo
import itertools
#import sys
#import scipy.optimize as spo
#import os
#import operator
#from PIL import Image as ipp
import multiprocessing as mpp
#
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.mpl as mpl
#import functools
#
#
#
#import ANSStools as atp
#import bindex
#import contours2kml
#import globalETAS as gep
#from eq_params import *

class Array_test(object):
	def __init__(self, N=1000):
		# just a toy, test script to work with mpp.Array() and other mpp stuff.
		#
		R1 = random.Random()
		R2 = random.Random()
		#
		self.X1 = [R1.random() for x in range(int(N))]
		self.X2 = [R2.random() for x in range(int(N))]
		#
		self.Z_spp = [x1 + x2 for x1,x2 in zip(self.X1, self.X2)]
	#
	def add_arrays(self, X1=None, X2=None):
		if X1==None: X1=self.X1
		if X2==None: X2=self.X2
		#
		return [x1 + x2 for x1,x2 in zip(X1,X2)]

class Array_test_worker(mpp.Process):
	def __init__(self, ary1, ary2, aryout, j_start=0, j_stop=None):
		super(Array_test_worker,self).__init__()
		j_stop = (j_stop or len(ary1))
		self.ary1 = ary1
		self.ary2 = ary2
		self.aryout=aryout
		self.j_start=(j_start or 0)
		self.j_stop=(j_stop or len(ary1))
		#self.pipe=pipe
		#
		pass
	def add_arrays(self):
		print("begin adding for %d, %d" % (self.j_start,self.j_stop))
		#for j in range(self.j_start,self.j_stop):
		#	self.aryout[j] += (self.ary1[j] + self.ary2[j])
		#
		self.aryout[self.j_start:self.j_stop] = [x1+x2 for x1,x2 in zip(self.ary1[self.j_start:self.j_stop], self.ary2[self.j_start:self.j_stop])]
		#
	def run(self):
		self.add_arrays()
		#self.pipe.send(self.aryout)
		#self.pipe.close()
	#
class Array_test_mpp(Array_test):
	def __init__(self, N=1000, n_procs=None):
		super(Array_test_mpp, self).__init__(N=N)
		self.n_procs = (n_procs or min(2,mpp.cpu_count()-1))
		lk = mpp.Lock()
		#
		self.X1=mpp.Array('d', self.X1, lock=lk)
		self.X2=mpp.Array('d', self.X2, lock=lk)
		self.sum_mpp = mpp.Array('d', numpy.zeros(len(self.X1)), lock=lk)
		#
	#
	def add_arrays(self):
		workers = []
		#pipes   = []
		#
		d_len = int(numpy.ceil(len(self.X1))/self.n_procs)
		for n in range(self.n_procs):
			#p1,p2 = mpp.Pipe()
			# this does work:
			workers += [Array_test_worker(self.X1, self.X2, self.sum_mpp, j_start = n*d_len, j_stop=min(len(self.X1), (n+1)*d_len))]
			#
			# this does not work:
			#j_start = n*d_len
			#j_stop=min(len(self.X1), (n+1)*d_len)
			#workers += [Array_test_worker(self.X1[j_start:j_stop], self.X2[j_start:j_stop], self.sum_mpp[j_start:j_stop], j_start = 0, j_stop=None)]
			#pipes += [p1]
			workers[-1].start()
		for j,w in enumerate(workers):
			w.join()
		#
def test_test_arrays(N=1e7):
	# ... and this all checks out. a couple questions:
	# do we need to use the full_array and index notation: pass the whole array and j_start, j_stop,
	# or can we pass do_it(ary=mpp.Array('d', lst[j:k])
	# ... and it looks like this works just fine. still need to test on n>2 cpus and manually check some of the values (just using the
	# test script validation here).
	print('testing shared memory array addition')
	T = Array_test_mpp(N=N)
	print('object made. now do mpp addition...')
	T.add_arrays()
	#
	print('added. now make a new copy and check...')
	#
	XX = [x1 + x2 for x1,x2 in zip(T.X1, T.X2)]
	for j,x in enumerate(XX):
		if XX[j]!=T.sum_mpp[j] or T.Z_spp[j]!=T.sum_mpp[j]: print("difference: ", T.sum_mpp[j], T.Z_spp[j], x)
	#
	print('done.')
	#
	return T	
#
class ROC_generic(object):
	def __init__(self, Z_events, Z_fc, h_denom=None, f_denom=None, f_start=0., f_stop=None):
		#
		print('initializing ROC_generi')
		self.Z_events = Z_events
		self.Z_fc     = Z_fc
		self.h_denom = float(h_denom or len(Z_events))
		self.f_denom = float(f_denom or len(Z_fc))
		self.f_stop  = (f_stop or len(Z_fc))
		self.f_start = (f_start or 0)
		#
		# separating H,F arrays to facilitate easier mpp.Array() access... though we CAN count to 2, and i think
		# mpp.Array() can use a Point() data type, which appears to be a 2D array or array-dict.
		# we'll look into keeping [hits, misses, falsies, correct-non-events] at another time...
		self.H = []
		self.F = []
	#	
	
	# there are lots of ways to calc ROC; lots of shortcuts and approximations for speed optimization. here are a few...
	def roc_simple_sparse(self, h_denom=None, f_denom=None, f_start=0, f_stop=None):
		# simple ROC calculator.
		h_denom = (h_denom or self.h_denom)
		f_denom = (f_denom or self.f_denom)
		f_start = (f_start or self.f_start)
		f_stop  = (f_stop  or self.f_stop)
		Z_events = self.Z_events
		Z_fc = self.Z_fc		# these just make a reference, but it still writes a copy of a variable, so it might be a bit faster to just reff self.x, etc.
		#
		h_denom = (h_denom or len(Z_events))
		f_denom = (f_denom or len(Z_fc))
		#
		# more precuse, but not quite as fast:
		
		# direct calc. approach. there's a loop, so maybe a bit slow.
		'''
		r_XY = []
		# direct way, but with a regular for-loop...
		for j,z_fc in enumerate(Z_fc):
			n_h = sum([(z_ev>=z_fc) for z_ev in Z_events])
			# falsies: approximately number of sites>z0-n_hits. for large, sparse maps, f~sum((z_fc>z0)), but this is a simple correction.
			r_XY += [[(f_start + j - n_h)/f_denom, n_h/h_denom]]
		'''
		#
		# list comprehensions, but 2 loops... not sure which should be faster.
		self.H = [sum([float(z_ev>=z_fc) for z_ev in self.Z_events])/h_denom for j,z_fc in enumerate(self.Z_fc)]
		self.F = [(f_start + j - self.H[j]*h_denom)/f_denom for j,x in enumerate(self.Z_fc)]
		#
		#return r_XY
	#
	def roc_sparse_approx(self):
		# a minor approximation when not sparse:
		self.F, self.H = zip(* [[(self.f_start + j)/self.f_denom, sum([(z_ev>=z_fc) for z_ev in self.Z_events])/self.h_denom] for j,z_fc in enumerate(self.Z_fc)] )
	#
	'''
	def roc_hits_sparse(self, h_denom=None):
		h_denom = (h_denom or self.h_denom)
		# return just the hits. for mpp operations, this might be faster, since the Fs array is just [j/f_denom for j,x in enumerate(z_fc)] (or equivalent).
		# it'll taka longer to pipe it back than to calc it.
		# ... and note, this isn't quite right; it assumes that all z values are unique... which they may not be. that said, it is nominally a reasonable 
		# approximation and it should be much faster than if we actually check z_fc like f = sum([z>=z0 for z in Z_fc]) and easier than if we have to 
		# 'manually' check for unequal values.
		return [sum([float(z_ev>=z_fc) for z_ev in self.Z_events])/h_denom for j,z_fc in enumerate(self.Z_fc)]
	'''
	#
#
class ROC_generic_mpp(ROC_generic):
	def __init__(self, n_procs=None, *args, **kwargs):
		super(ROC_generic_mpp, self).__init__(*args,**kwargs)
		if n_procs==None: n_procs = mpp.cpu_count()
		self.n_procs=n_procs
		#
		lk = mpp.Lock()
		self.Z_events = mpp.Array('d', self.Z_events, lock=lk)
		self.Z_fc     = mpp.Array('d', self.Z_fc, lock=lk)
		self.H        = mpp.Array('d', len(self.Z_fc), lock=lk)
		self.F        = mpp.Array('d', len(self.Z_fc), lock=lk)

		#
	#
	def calc_roc(self):
		workers = []
		d_len = int(numpy.ceil(len(self.Z_fc))/self.n_procs)
		for n in range(self.n_procs):
			workers += [ROC_generic_worker(H=self.H, F=self.F, Z_events=self.Z_events, Z_fc=self.Z_fc, h_denom=self.h_denom, f_denom=self.f_denom, f_start = n*d_len, f_stop=min(len(self.Z_fc), (n+1)*d_len))]
			#
			workers[-1].start()
		for j,w in enumerate(workers):
			w.join()
		#
	def plot_FH(self, fignum=0):
		plt.figure(fignum)
		plt.clf()
		#
		plt.plot(self.F, self.H, '.', lw=2.)
		plt.plot(range(2), range(2), 'r-', lw=2.5)
 
	
	#
class ROC_generic_worker(ROC_generic, mpp.Process):
	def __init__(self, H=None, F=None, pipe_r=None, *args, **kwargs):
		'''
		# worker for generic ROC mpp operations. for shared memory mpp.Array(), pass mpp.Array() objects (addresses) for H,F.
		# for piping approach, pass one end of a Pipe() into pipe_r. we'll check for pipe_r when we finish calculating and
		# send back if it exists.
		'''
		super(ROC_generic_worker, self).__init__(*args, **kwargs)
		mpp.Process.__init__(self)
		self.H = (H or [])
		self.F = (F or [])	# can we set these like self.F = F[f_start:f_stop] ? will this assign by reference?
							# then, we can use f_start to augment f_denom and assign to F[j] instead of F[j+f_start].
							
		self.pipe_r = pipe_r
	def run(self):
		self.roc_simple_sparse()
		if self.pipe_r!=None:
			# ... but actually, this will be a little more complicated. the f_start/f_stop parameters get a bit confused here.
			# 
			self.pipe_r.send(zip(self.F, self.H)
			self.pipe_r.close()
	#
	def roc_simple_sparse(self, h_denom=None, f_denom=None, f_start=0, f_stop=None):
		# simple ROC calculator.
		print('calcingn roc...')
		h_denom = (h_denom or self.h_denom)
		f_denom = (f_denom or self.f_denom)
		f_start = (f_start or self.f_start)
		f_stop  = (f_stop  or self.f_stop)
		print('calcingn roc...', f_start, f_stop)
		#
		Z_events = self.Z_events
		Z_fc = self.Z_fc		# these just make a reference, but it still writes a copy of a variable, so it might be a bit faster to just reff self.x, etc.
		#
		h_denom = (h_denom or len(Z_events))
		f_denom = (f_denom or len(Z_fc))
		#
		# more precuse, but not quite as fast:
		
		# direct calc. approach. there's a loop, so maybe a bit slow.
		
		r_XY = []
		# direct way, but with a regular for-loop...
		n_total = len(self.Z_fc)
		for j,z_fc in enumerate(Z_fc[f_start:f_stop]):
			n_h = sum([(z_ev>=z_fc) for z_ev in Z_events])
			# falsies: approximately number of sites>z0-n_hits. for large, sparse maps, f~sum((z_fc>z0)), but this is a simple correction.
			#r_XY += [[(f_start + j - n_h)/f_denom, n_h/h_denom]]
			#
			n_h = sum([(z_ev>=z_fc) for z_ev in Z_events])
			self.H[j+f_start] = n_h/h_denom
			self.F[j+f_start] = (n_total - f_start - j - n_h)/f_denom
		

		#
def ROC_test1(N1=100, N2=10000, n_procs=None):
	R=random.Random()
	A=ROC_generic_mpp(n_procs=n_procs, Z_events=[R.random() for j in range(N1)], Z_fc=sorted([R.random() for j in range(N2)]))
	roc=A.calc_roc()
	A.plot_FH()
	
	return A

def ROC_test_n_bench(N1=1000, N2=100000, n_tests=10):
	# TODO: run a handfull of test and benchmarks:
	# make a data set (random catalog, random forecast.)
	# for these data, run n_tests of: spp, mpp_shared_Array, mpp_piped and compare times to complete.
	pass
	
#
if __name__=='__main__':
	pass
else:
	plt.ion()
		
#

