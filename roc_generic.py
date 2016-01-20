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
		j_stop = (j_stop or len(ary))
		self.ary1 = ary1
		self.ary2 = ary2
		self.aryout=aryout
		self.j_start=j_start
		self.j_stop=j_stop
		#self.pipe=pipe
		#
		pass
	def add_arrays(self):
		print("begin adding for %d, %d" % (self.j_start,self.j_stop))
		for j in range(self.j_start,self.j_stop):
			self.aryout[j] += (self.ary1[j] + self.ary2[j])
		#
	def run(self):
		self.add_arrays()
		#self.pipe.send(self.aryout)
		#self.pipe.close()
	#
class Array_test_mpp(Array_test):
	def __init__(self, N=1000, n_procs=None):
		super(Array_test_mpp, self).__init__(N=N)
		self.n_procs = (n_procs or max(2,min(1,mpp.cpu_count()-1)))
		lk = mpp.Lock()
		#
		self.X1=mpp.Array('d', self.X1, lock=lk)
		self.X2=mpp.Array('d', self.X2, lock=lk)
		self.sum_mpp = mpp.Array('d', numpy.zeros(len(self.X1)), lock=lk)
		#
	#
	def add_arrays(self):
		workers = []
		pipes   = []
		#
		d_len = int(numpy.ceil(len(self.X1))/self.n_procs)
		for n in range(self.n_procs):
			#p1,p2 = mpp.Pipe()
			# ary1, ary2, aryout, j_start=0, j_stop=None, pipe=None):
			#workers += [Array_test_worker(self.X1, self.X2, self.sum_mpp, j_start = n*d_len, j_stop=min(len(self.X1), (n+1)*d_len))]
			j_start = n*d_len
			j_stop=min(len(self.X1), (n+1)*d_len))
			workers += [Array_test_worker(self.X1[j_start:j_stop], self.X2[j_start:j_stop], self.sum_mpp[j_start:j_stop], j_start = 0, j_stop=None]
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
	
	
		

class ROC_generic(object):
	def __init__(self, Z_events, Z_fc, h_denom=None, f_denom=None, f_start=0.):
		#self.__dict__.update({key:val for key,val in locals().items if key not in ('self', 'class')})
		self.Z_events = Z_events
		self.Z_fc     = Z_fc
		self.h_denom = float(h_denom or len(Z_fc))
		self.f_denom = float(f_denom or len(Z_events))
		#
		self.H = []
		self.F = []
	#	
	
	# there are lots of ways to calc ROC; lots of shortcuts and approximations for speed optimization. here are a few...
	def roc_simple_sparse(self, h_denom=None, f_denom=None, f_start=0):
		# simple ROC calculator.
		h_denom = (h_denom or self.h_denom)
		f_denom = (f_denom or self.f_denom)
		f_start = (f_start or self.f_start)
		Z_events = self.Z_events
		Z_fc = self.Z_fc		# these just make a reference, but it still writes a copy of a variable, so it might be a bit faster to just reff self.x, etc.
		#
		h_denom = (h_denom or len(Z_events))
		f_denom = (f_denom or len(Z_fc))
		#
		# more precuse, but not quite as fast:
		r_XY = []
		for j,z_fc in enumerate(Z_fc):
			n_f = sum([(z_ev>=z_fc) for z_ev in Z_events])
			r_XY += [[(f_start + j - n_f)/f_denom, n_f/h_denom]]
		#
		return r_XY
	#
	def roc_sparse_approx(self):
		# a minor approximation when not sparse:
		return [[(self.f_start + j)/self.f_denom, sum([(z_ev>=z_fc) for z_ev in self.Z_events])/self.h_denom] for j,z_fc in enumerate(self.Z_fc)]
	#
	def roc_hits_sparse(self, h_denom=None):
		h_denom = (h_denom or self.h_denom)
		# return just the hits. for mpp operations, this might be faster, since the Fs array is just [j/f_denom for j,x in enumerate(z_fc)] (or equivalent).
		# it'll taka longer to pipe it back than to calc it.
		# ... and note, this isn't quite right; it assumes that all z values are unique... which they may not be. that said, it is nominally a reasonable 
		# approximation and it should be much faster than if we actually check z_fc like f = sum([z>=z0 for z in Z_fc]) and easier than if we have to 
		# 'manually' check for unequal values.
		return [sum([float(z_ev>=z_fc) for z_ev in self.Z_events])/h_denom for j,z_fc in enumerate(self.Z_fc)]

