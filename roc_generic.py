'''
# some generic ROC tools, and playing a bit with shared memory mpp as well.
'''
###
#
import random
import numpy
import scipy
import itertools
import multiprocessing as mpp
#
import matplotlib.pyplot as plt
import matplotlib.mpl as mpl
import time
#import functools
#
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
	#def __init__(self, Z_events, Z_fc, h_denom=None, f_denom=None, f_start=0., f_stop=None, do_normalize=True):
	def __init__(self, Z_events, Z_fc, h_denom=None, f_denom=None, f_start=0., f_stop=None):
		#
		#print('initializing ROC_generi')
		self.Z_events = sorted(Z_events)
		self.Z_fc     = numpy.array(sorted(list(set(Z_fc))))
		self.h_denom = float(h_denom or len(Z_events))
		self.f_denom = float(f_denom or len(Z_fc))
		self.f_stop  = (f_stop or len(Z_fc))
		self.f_start = (f_start or 0)
		#self.do_normalziation = (do_normalize or True)
		#
		# separating H,F arrays to facilitate easier mpp.Array() access... though we CAN count to 2, and i think
		# mpp.Array() can use a Point() data type, which appears to be a 2D array or array-dict.
		# we'll look into keeping [hits, misses, falsies, correct-non-events] at another time...
		self.H = []
		self.F = []
	#	
	
	# there are lots of ways to calc ROC; lots of shortcuts and approximations for speed optimization. here are a few...
	#def roc_simple_sparse(self, h_denom=None, f_denom=None, f_start=0, f_stop=None, do_normalize=None):
	def roc_simple_sparse(self, h_denom=None, f_denom=None, f_start=0, f_stop=None):
		# simple ROC calculator.
		#@do_normalize: normalise F,H. for mpp calculations, it might be simpler to return the total sums and normalize in the parent class/calling funct.
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
		# this is an attempt to simplify the roc process a bit, but i'm not sure it's the right approach.
		#do_normalize=(do_normalize or self.do_normalize)
		#do_normalize=(do_normalize or True)
		#if do_normalize:
		#	h_denom = 1.0
		#	f_denom = 1.0
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
	def plot_FH(self, fignum=0):
		plt.figure(fignum)
		plt.clf()
		#
		plt.plot(self.F, self.H, '.', lw=2.)
		plt.plot(range(2), range(2), 'r-', lw=2.5)
		ax=plt.gca()
		ax.set_xlabel('False Alarm Rate $F$')
		ax.set_ylabel('Hit Rate $H$')
	#
#
class ROC_generic_mpp_pipes(ROC_generic):
	# an mpp ROC handler that pipes data to/from workers, instead of using shared memory mpp.Array() objects.
	# this will most likely be the default ROC handler; as it turns out, it is much faster than its shared memory counterpart.
	#
	def __init__(self, n_procs=None, *args, **kwargs):
		super(ROC_generic_mpp_pipes, self).__init__(*args,**kwargs)
		if n_procs==None: n_procs = mpp.cpu_count()
		self.n_procs=n_procs
		#
		#self.H = numpy.zeros(len(self.Z_fc))
		#self.F = numpy.zeros(len(self.Z_fc))
		self.H = []
		self.F = []
		
	#
	def calc_roc(self):
		workers = []
		pipes   = []
		d_len = int(numpy.ceil(len(self.Z_fc))/self.n_procs)
		f_j_stop  = lambda n: min(len(self.Z_fc), (n+1)*d_len)
		f_j_start = lambda n: n*d_len 
		#
		for n in range(self.n_procs):
			#workers += [ROC_generic_worker(H=self.H, F=self.F, Z_events=self.Z_events, Z_fc=self.Z_fc, h_denom=self.h_denom, f_denom=self.f_denom, f_start = n*d_len, f_stop=min(len(self.Z_fc), (n+1)*d_len))]
			#j_start = n*d_len
			#j_end   = min(len(self.Z_fc), (n+1)*d_len)
			p1, p2 = mpp.Pipe()
			#
			workers += [ROC_generic_worker(H=None, F=None, Z_events=self.Z_events, Z_fc=self.Z_fc[f_j_start(n):f_j_stop(n)], h_denom=self.h_denom, f_denom=self.f_denom, pipe_r=p2, f_start_index = f_j_start(n), len_fc = len(self.Z_fc))]
			# , f_start = n*d_len, f_stop=min(len(self.Z_fc), (n+1)*d_len)
			pipes += [p1]
			#
			workers[-1].start()
		#
		FH = []
		for j,p in enumerate(pipes):
			#FH += [zip(*p.recv())]
			
			FH = p.recv()
			F,H = zip(*FH)
			#
			F  = numpy.array(F)
			self.F += F.tolist()
			self.H += H
			
			#F = numpy.array(F)
			#F += f_j_stop(j)/self.f_denom
			#FH = zip(F,H)
			p.close()
		for j,w in enumerate(workers):
			w.join()
			
#
class ROC_generic_mpp_shared(ROC_generic):
	def __init__(self, n_procs=None, *args, **kwargs):
		super(ROC_generic_mpp_shared, self).__init__(*args,**kwargs)
		if n_procs==None: n_procs = mpp.cpu_count()
		self.n_procs=n_procs
		#
		lk = mpp.Lock()
		# so i think that since we're only reading Z_events and Z_fc, we can set lock=False (that might not be the right syntax).
		self.Z_events = mpp.Array('d', self.Z_events, lock=False)
		self.Z_fc     = mpp.Array('d', self.Z_fc, lock=False)
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
#
#
class ROC_mpp(ROC_generic_mpp_pipes):
	# a name-wrapper for the "standard" ROC_mpp handler:
	def __init__(self, *args, **kwargs):
		super(ROC_mpp, self).__init__(*args, **kwargs)
		#
	#

#
class ROC_generic_worker(ROC_generic, mpp.Process):
	def __init__(self, H=None, F=None, pipe_r=None, f_start_index = None, len_fc=None, *args, **kwargs):
		'''

		# worker for generic ROC mpp operations. for shared memory mpp.Array(), pass mpp.Array() objects (addresses) for H,F.
		# for piping approach, pass one end of a Pipe() into pipe_r. we'll check for pipe_r when we finish calculating and
		# send back if it exists.
		#
		# in this version of the worker, we can typically pass H,F as None values.
		'''
		super(ROC_generic_worker, self).__init__(*args, **kwargs)
		mpp.Process.__init__(self)
		#self.H = (H or [])		# ... and we might change these defaults to be the length of the Z_fc array.
		#self.F = (F or [])	# can we set these like self.F = F[f_start:f_stop] ? will this assign by reference?
							# then, we can use f_start to augment f_denom and assign to F[j] instead of F[j+f_start].
		#
		# (this "if none" sytax might not work because H,F are arrays).
		self.H = (H or numpy.zeros(len(self.Z_fc)))
		self.F = (F or numpy.zeros(len(self.Z_fc)))
		self.f_start_index = (f_start_index or self.f_start)
		self.len_fc = (len_fc or len(self.Z_fc))
							
		self.pipe_r = pipe_r
	def run(self):
		self.roc_simple_sparse()
		if self.pipe_r!=None:
			# ... but actually, this will be a little more complicated. the f_start/f_stop parameters (indices) get a bit confused here.
			# 
			self.pipe_r.send(zip(self.F, self.H))
			self.pipe_r.close()
	#
	def roc_simple_sparse(self, h_denom=None, f_denom=None, f_start=0, f_stop=None, f_start_index=None, len_fc=None):
		# simple ROC calculator.
		#print('calcingn roc..., f_start_index = ', f_start_index, self.f_start_index)
		# @f_start_index: the F starting index, which might be different than f_start. f_start is the actual starting index
		# for the local array, typically as used by a shared memory (mpp.Array() ) approach. f_start_index indicates the starting
		# index, from the parent list, of the list passed to this object (so for the first half of X, f_start_index=0, for the second
		# half, f_start_index = len(X)/2
		#
		h_denom = (h_denom or self.h_denom)
		f_denom = (f_denom or self.f_denom)
		f_start = (f_start or self.f_start)
		f_start = (f_start or 0)
		f_stop  = (f_stop  or self.f_stop)
		f_stop  = (f_stop  or len(self.Z_fc))
		#
		f_start_index = (f_start_index or self.f_start_index)
		f_start_index = (f_start_index or f_start)
		len_fc = (len_fc or self.len_fc)
		len_fc = (len_fc or len(self.Z_fc))
		#print('calcing roc...', f_start, f_stop, f_start_index, ' **', len_fc)
		#
		Z_events = self.Z_events
		Z_events.sort()
		Z_fc = self.Z_fc		# these just make a reference, but it still writes a copy of a variable, so it might be a bit faster to just reff self.x, etc.
		#
		h_denom = (h_denom or len(Z_events))
		f_denom = (f_denom or len(Z_fc))
		#
		# more precuse, but not quite as fast:
		# direct calc. approach. there's a loop, so maybe a bit slow.
		#
		r_XY = []
		j_events = 0
		len_ev = len(Z_events)
		# direct way, but with a regular for-loop...
		#n_total = len(self.Z_fc)
		for j,z_fc in enumerate(Z_fc[f_start:f_stop]):
			#n_h = sum([(z_ev>=z_fc) for z_ev in Z_events])
			# falsies: approximately number of sites>z0-n_hits. for large, sparse maps, f~sum((z_fc>z0)), but this is a simple correction.
			#r_XY += [[(f_start + j - n_h)/f_denom, n_h/h_denom]]
			#
			#n_h = sum([(z_ev>=z_fc) for z_ev in Z_events])
			
			#n_h = len([z_ev for z_ev in Z_events if z_ev>=z_fc])	# this should be faster. we might also look into dropping th part of z_ev that we know is <z_fc.
			# ... and this should be much faster than both... but needs testing. (... dunno)
			while Z_events[j_events]<z_fc and j_events<(len_ev-1): j_events+=1
			n_h = len_ev-j_events
			#
			self.H[j+f_start] = n_h/h_denom
			self.F[j+f_start] = (len_fc - f_start_index - j - n_h)/f_denom
		#
def ROC_test1(N1=100, N2=10000, n_procs=None):
	R=random.Random()
	#
	# so based on this test, it looks like the piped approach is much, much faster (2-3 times faster) than using mpp.Array(),
	# even having made a holy mess of the final handling of the F array. so we'll want to clean up the return sequences (maybe
	# pass an additinal variable to properly normalize the F array in the child processes.
	#
	Z_events=[R.random() for j in range(N1)]
	Z_fc=sorted([R.random() for j in range(N2)])
	#
	t0=time.time()
	print('starting: ', t0)
	A=ROC_generic_mpp_shared(n_procs=n_procs, Z_events=Z_events, Z_fc=Z_fc)
	roc=A.calc_roc()
	t1=time.time()
	#A.plot_FH(fignum=0)
	#
	
	print('finished 1: ', t1, " :: ", t1-t0)
	
	t2=time.time()
	B = ROC_generic_mpp_pipes(n_procs=n_procs, Z_events=Z_events, Z_fc=Z_fc)
	roc = B.calc_roc()
	t3=time.time()
	B.plot_FH(fignum=1)
	
	
	print("finished 2: ", t2, " :: ", t3, " :: ", t3-t2)
	
	t4 = time.time()
	C = ROC_mpp(n_procs=n_procs+1, Z_events=Z_events, Z_fc=Z_fc)
	roc = C.calc_roc()
	t5 = time.time()
	C.plot_FH(fignum=2)
	
	print("finished 3: ", t4, " :: ", t5, " :: ", t5-t4)
	#
	
	#return A

def roc_bench_1(n_cycles=10, cpu_range=[2,10],N1=100, N2=10000):
	#
	R=random.Random()
	XY = []
	for j,k in itertools.product(range(n_cycles), range(cpu_range[0],cpu_range[1]+1)):
		Z_events=[R.random() for j in range(N1)]
		Z_fc=sorted([R.random() for j in range(N2)])
		t0 = time.time()
		C = ROC_mpp(n_procs=k, Z_events=Z_events, Z_fc=Z_fc)
		roc = C.calc_roc()
		print('max,min(F,H): ', min(C.F), max(C.F), min(C.H), max(C.H))
		t1 = time.time()
		#
		dt = t1-t0
		print("completed (%d,%d), dt = %f" % (j,k,dt))
		XY += [[k,dt]]
	#
	plt.figure(0)
	plt.clf()
	plt.plot(*zip(*XY), marker='o', ls='')
	plt.ylabel('Time interval $\\Delta t$', size=18)
	plt.xlabel('N_threads, $n$', size=18)
	plt.title('n_cpu bench test, actual n_cpu()=%d' % mpp.cpu_count())
	#
	return XY
	


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

