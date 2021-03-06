'''
# DEPRICATION: This module represents a good investigation, development, and benchmarking of some mpp methods, but a
# much simpler, faster ROC can be found in the optimizers/roc_tools.py module.
#
# some generic ROC tools, and playing a bit with shared memory mpp as well.
# (eventually, this needs to be moved into a general-tools repository). also,
# lets see if we can't @numba.jit  compile some of these bits...
#
# also, it might be smarter to write the ROC processes just as procedural functions,
# namely because 1) they don't depend on class-scope variables (input-output),
# 2) we can compile them with @numba.jit, and 3) as per 1,2, they may be more
# portable that way.
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
import matplotlib as mpl
import time
#import functools
try:
	import numba
	have_numba=True
	print('numba imported')
except:
	have_numba=False
	print('numba not available. you should install it.')
#
print('***************** roc_generic.py ***********************')
print('*****************\nDEPRICATION WARNING:\nThis module is being depricated; look at yodiipy.optimizers.roc_tools()')
print('and possibly something like etas_roc_tools.py in the globalETAS folder.')
print('this module contains some working code... and some not working code, so be very very careful,')
print('particularly when running the more optimized codes.')
print('NOTE: the Aray() (shared memory) vs piped mpp tests are probably worth keeping, particularly the Array() code sample.')
#
#
class Array_test(object):
	def __init__(self, N=1000):
		# just a toy, test script to work with mpp.Array() and other mpp stuff.
		# this (shared memory) approach to Python mpp significatnly underperformed standard piping models in earlier tests.
		# however, it might be worth looking at it for the newer, single-pass ROC methods (see optimizers.roc_tools.py).
		# because this approach is single-pass, its spp implementation is faster than piped mpp, but shared memory (Array() ),
		# even in Python, might be faster, especially for a simplified metric, like Molchan vs true ROC...
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
		# (so this function can probably be compiled using numba, but maybe not inside the class object?
		# ... and why not use numpy.array()
		if X1==None: X1=self.X1
		if X2==None: X2=self.X2
		#
		return [x1 + x2 for x1,x2 in zip(X1,X2)]

# here's how to compile with numba:
# ... but how do we make just the compilation part conditional?
if have_numba:
	@numba.jit
	def add_arrays(X1=None, X2=None):
		# (so this function can probably be compiled using numba, but maybe not inside the class object?
		# ... and why not use numpy.array()
		#
		return [x1 + x2 for x1,x2 in zip(X1,X2)]
else:
	def add_arrays(X1=None, X2=None):
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
	# TODO: needs to be further validated...
	# ... it seems to be working but everything needs to be checked...
	#
	# generic, multi-purpose ROC class. this is super simple until multi-processing comes in, and then it gets
	# a little bit tricky to make sure all the denominators are correct. there are a couple of approaches to take.
	# 1) pass the denominator(s) to the workers (aka, the parent thread knows the full list size, so pass that to 
	# the workers; the workers can then do the division on-process; all workers are identical, and we can return
	# (pipe back) smaller data sets, 2) just do sorting and aggregation in-process; re-sort, aggregate, etc.
	# in the parent-process (aka, shift more work to the parent process).
	#
	# ######
	# also, this needs to be benchmarked for speed. tests on another ROC code-base definitely showed that mpp did NOT
	# make it run faster... in fact much slower. basically, ROC is a single pass through the data, and if we pipe
	# it to a process, well that's two passes. so it should be investigated whether or not piping off multiple ROC jobs
	# can be faster, but i'm going to guess that it isn't.
	# for more, see the "optimizers" repository... where this code should eventually end up anyway.
	#
	# at this time, we've chosen to optionally pass denominators to the process. this is a bit more complicated,
	# but it makes the tool more versatile, and if we just pass None objects into the denoms, we get an expected
	# default behavior.
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
	def roc_simple_sparse(self, h_denom=None, f_denom=None, f_start=0, f_stop=None, do_sort=True):
		# ... actually, this needs to be re-done or removed. i don't think we can actually solve for the correct F with the information provided.
		# simple ROC calculator.
		# ... this is not quite right, however, because it does not account for multiple events in a single bin...
		# but it's close for fine lattices.
		#@do_normalize: normalise F,H. for mpp calculations, it might be simpler to return the total sums and normalize in the parent class/calling funct.
		h_denom = (h_denom or self.h_denom)
		f_denom = (f_denom or self.f_denom)
		f_start = (f_start or self.f_start)
		f_stop  = (f_stop  or self.f_stop)
		Z_events = self.Z_events
		Z_fc = self.Z_fc		# these just make a reference, but it still writes a copy of a variable, so it might be a bit faster to just reff self.x, etc.
		#
		print('doing ss, f_start={}'.format(f_start))
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
		#
		
		# direct calc. approach. there's a loop, so maybe a bit slow.
		if do_sort:
			Z_fc.sort()
			Z_events.sort()
		N_f = len(Z_fc)
		#
		FH = []
		#
		FH=[[0., 0.]]
		# fix these params for now:
		do_sort=True
		j_eq0 = 0	# holder
		n_eq=float(len(Z_events))
		#
		k_eq = 0
		for j,z_fc in enumerate(Z_fc):	
			#
			if k_eq<n_eq and z_fc>=Z_events[k_eq]:
				while k_eq<n_eq and z_fc>=Z_events[k_eq]:
					#
					#FH += [[j+1,z_fc, k_eq+1, Z_ev[k_eq]]]
					#
					#FH += [[(j+1 + j_fc0)/f_denom, (k_eq + j_eq0 +1)/h_denom]]
					#print('ev: ', k_eq, Z_ev[k_eq], z_fc)
					k_eq+=1
				#
				FH += [[(j+1 + f_start)/f_denom, (k_eq + j_eq0 )/h_denom]]
		self.F, self.H = numpy.array(list(zip(*FH)))
		#
		# list comprehensions, but 2 loops... not sure which should be faster.
		# (and this can be compiled with numba.jit if we are so inclinded -- which we should be).
		#self.H = [sum([float(z_ev>=z_fc) for z_ev in self.Z_events])/h_denom for j,z_fc in enumerate(self.Z_fc)]
		#self.F = [(f_start + j - self.H[j]*h_denom)/f_denom for j,x in enumerate(self.Z_fc)]
		#
		return FH
	#
	def roc_sparse_approx(self):
		# a minor approximation when not sparse:
		# but note we repeatedly sum over the events array, so this gets really slow for large catalogs.
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
			print('fstart: ', f_j_start(n))
			#workers += [ROC_generic_worker(H=None, F=None, Z_events=self.Z_events, Z_fc=self.Z_fc[f_j_start(n):f_j_stop(n)], h_denom=self.h_denom, f_denom=self.f_denom, pipe_r=p2, f_start_index = f_j_start(n), len_fc = len(self.Z_fc))]
			workers += [ROC_generic_worker(H=None, F=None, Z_events=self.Z_events, Z_fc=self.Z_fc[f_j_start(n):f_j_stop(n)], h_denom=self.h_denom, f_denom=self.f_denom, pipe_r=p2, f_start = f_j_start(n), len_fc = len(self.Z_fc))]
			#
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
	@property
	def FH(self):
		return list(zip(self.F, self.H))
			
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
	#def __init__(self, H=None, F=None, pipe_r=None, f_start_index = None, len_fc=None, *args, **kwargs):
	def __init__(self, H=None, F=None, pipe_r=None, f_start = None, len_fc=None, *args, **kwargs):
		'''
		# worker for generic ROC mpp operations. for shared memory mpp.Array(), pass mpp.Array() objects (addresses) for H,F.
		# for piping approach, pass one end of a Pipe() into pipe_r. we'll check for pipe_r when we finish calculating and
		# send back if it exists.
		#
		# in this version of the worker, we can typically pass H,F as None values.
		'''
		#
		super(ROC_generic_worker, self).__init__(*args, **kwargs)
		mpp.Process.__init__(self)
		#self.H = (H or [])		# ... and we might change these defaults to be the length of the Z_fc array.
		#self.F = (F or [])	# can we set these like self.F = F[f_start:f_stop] ? will this assign by reference?
							# then, we can use f_start to augment f_denom and assign to F[j] instead of F[j+f_start].
		#
		# (this "if none" sytax might not work because H,F are arrays).
		# note: coming soon to python >3.5, X==None will be an element-wise comparison; X is None will compare the
		# list object itself. i'm not sure how the (X or X0) syntax will work; it has been problematic in the past.
		self.H = (H or numpy.zeros(len(self.Z_fc)))
		self.F = (F or numpy.zeros(len(self.Z_fc)))
		#self.f_start_index = (f_start_index or self.f_start)
		self.f_start = (f_start or self.f_start)
		self.len_fc = (len_fc or len(self.Z_fc))
		
		#print('worker fs: ', f_start_index, self.f_start_index)
							
		self.pipe_r = pipe_r
	def run(self):
		#self.roc_simple_sparse()
		self.roc_simple_sparse()
		if self.pipe_r!=None:
			# ... but actually, this will be a little more complicated. the f_start/f_stop parameters (indices) get a bit confused here.
			# 
			self.pipe_r.send(zip(self.F, self.H))
			self.pipe_r.close()
	#

def ROC_test1(N1=100, N2=10000, n_procs=None):
	# TODO: bench this for multiple processors (for both Array() and piped approaches) and also bench against a SPP ROC -- one from this module
	# and also from the optimizers/roc_tools.py module
	#
	# test the changes to this; i think there is a bias in the roc analysis...
	#
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
	
	print('finished 1 (shared): ', t1, " :: ", t1-t0)
	
	t2=time.time()
	B = ROC_generic_mpp_pipes(n_procs=n_procs, Z_events=Z_events, Z_fc=Z_fc)
	roc = B.calc_roc()
	t3=time.time()
	B.plot_FH(fignum=1)
	
	
	print("finished 2 (pipes): ", t2, " :: ", t3, " :: ", t3-t2)
	
	t4 = time.time()
	C = ROC_mpp(n_procs=n_procs+1, Z_events=Z_events, Z_fc=Z_fc)
	roc = C.calc_roc()
	t5 = time.time()
	C.plot_FH(fignum=2)
	
	print("finished 3 (npipes, n+1): ", t4, " :: ", t5, " :: ", t5-t4)
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
#
def roc_generic_testfig(N_ev=1000, N_fc=10000, n_cpu=1, fnum=0):
	R1=random.Random()
	R2=random.Random()
	#
	Z_ev = [R1.random() for _ in range(N_ev)]
	Z_fc = [R2.random() for _ in range(N_fc)]
	#
	kwargs = {'Z_events':Z_ev, 'Z_fc':Z_fc}
	#
	t0=time.time()
	X_spp = ROC_generic(**kwargs)
	FH_spp = X_spp.roc_simple_sparse()
	print('dt: ', time.time()-t0)
	#
	#X_mpp = ROC_generic_mpp_pipes(n_cpu=n_cpu, 
	t0=time.time()
	#X_mpp = ROC_generic_mpp_pipes(**kwargs)
	X_mpp = ROC_mpp(n_procs=n_cpu, **kwargs)		# short-hand for default mpp style.
	roc = X_mpp.calc_roc()
	print('dt: ', time.time()-t0)
	#
	fg=plt.figure(fnum)
	plt.clf()
	ax1=fg.add_axes([.1,.1,.4,.8])
	ax2=fg.add_axes([.5,.1,.4,.8])
	#
	ax1.plot(*zip(*FH_spp), marker='.')
	ax1.plot(range(2), range(2), lw='2', ls='-', marker='', color='r')
	
	ax2.plot(X_mpp.F, X_mpp.H, '.')
	ax2.plot(range(2), range(2), lw='2', ls='-', marker='', color='r')
	#
	return X_spp


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

