
'''
# bindex.py: an indexed binning data structure (and affiliated scripts and bits fo code).
# primary author: mark r. yoder, phd.
#
# code is available as is, with no gurantees of any kind.
'''

import numpy
import scipy
import datetime as dtm
import pytz
import matplotlib.dates as mpd
import multiprocessing as mpp
#
import random

def bindex2d_sample():
	#
	R=random.Random()
	dx=dy=.1
	x0=y0=0.
	prams = {'dx':dx, 'dy':dy, 'x0':x0, 'y0':y0, 'leaf_type':float}
	#
	# get a float type bindex:
	B=Bindex2D(**prams)
	#
	datas = [[x,y, R.random()] for x in numpy.arange(0., 2., .1) for y in numpy.arange(0.,3.,.1)]
	#
	B.add_list(datas)
	#
	b_datas = B.to_list()
	print "datas list: d1==d2: %d, len(d1)==len(d2): %d" % (b_datas==datas, len(b_datas)==len(datas))
	B2 = Bindex2D(**prams)
	B2.add_list(datas)
	#
	print "bindex_items==bindex_items: %d" % (B.items==B2.items)
	#
	B3 = Bindex2D(**prams)
	B3.leaf_type=list
	#
	B3.add_list([[x,y, range(j + j*i,j + j*i+5)] for i,x in enumerate(numpy.arange(0., 2., .1)) for j,y in enumerate(numpy.arange(0.,3.,.1))])
	
	#
	return B3
	

class Bindex2D(dict):
	
	# for now, let's just start with a 2D bindex and we'll generalize later
	def __init__(self, dx=1., dy=1., x0=0., y0=0., leaf_type=float):
		# leaf_type: should be a function that will evaluate as the default value for an object type for the leaf (end) nodes.
		#
		self.x0=x0
		self.y0=y0
		self.dx=dx
		self.dy=dy
		#
		self.bin_y_0 = 0
		self.bin_x_0 = 0
		#
		self.leaf_type=leaf_type
		#
		#self.items = {}
	#
	#def __iter__(self):
	#	return self.items.__iter__
	#def __sizeof__(self):
	#	return sum([len(rw) for rw in self.items.iteritems()])
	#
	def get_xbin_center(self, x, dx=None, x0=None, bin_0=0):
		'''
		# return a dict. with the bin index and center location closest to x.
		'''
		#
		j = self.get_xbin(x, dx=dx, x0=x0, bin_0=bin_0)
		xx = self.get_bin2x(j, dx=dx, x0=x0, bin_0=bin_0)
		#
		return {'index':j, 'center':xx}
	
	def get_ybin_center(self, y, dy=None, y0=None, bin_0=0):
		'''
		# return a dict. with the bin index and center location closest to x.
		'''
		#
		j = self.get_ybin(y, dy=dy, y0=y0, bin_0=bin_0)
		xx = self.get_bin2y(j, dy=dy, y0=y0, bin_0=bin_0)
		#
		return {'index':j, 'center':xx}
		
	def get_xbin(self, x, dx=None, x0=None, bin_0=0):
		dx = (dx or self.dx)
		x0 = (x0 or self.x0)
		#
		#return int((x-x0)/dx)
		return get_bin(x, dx=dx, x0=x0, bin_0=bin_0)
	#
	def get_bin2x(self, bin_num, dx=None, x0=None, bin_0=0):
		dx = (dx or self.dx)
		x0 = (x0 or self.x0)
		#
		return bin2x(bin_num, dx=dx, x0=x0, bin_0=bin_0)
	#
	def get_ybin(self, y, dy=None, y0=None, bin_0=0):
		dy = (dy or self.dy)
		y0 = (y0 or self.y0)
		return get_bin(y, dx=dy, x0=y0, bin_0=bin_0)
	#
	def get_bin2y(self, bin_num, dy=None, y0=None, bin_0=0):
		dy = (dy or self.dy)
		y0 = (y0 or self.y0)
		#
		return bin2x(bin_num, dx=dy, x0=y0, bin_0=bin_0)
	#
	def get_xybins(self, x, y, dx=None, dy=None,  bin_x_0=0, bin_y_0=0, r_type='dict'):
		dx = (dx or  self.dx)
		dy = (dy or self.dy)
		#
		x_bin = get_bin(x=x, dx=dx, x0=bin_x_0, bin_0=bin_x_0)
		y_bin = get_bin(x=y, dx=dy, x0=bin_y_0, bin_0=bin_y_0)
		#
		if r_type == 'tuple':
			return (x_bin,y_bin)
		elif r_type == 'list':
			return [x_bin, y_bin]
		elif r_type == 'dict':
			return {'x':x_bin, 'y':y_bin}
		else:
			return [x_bin,y_bin]
	#
	def get_bin_items(self, x=None, y=None):
		# get all the items in the bin that contains coordinates x,y
		x0, y0 = self.get_xybins(x,y, r_type='tuple')
		#
		#return self.items[x0][y0]
		return self[x0][y0]
	#
	#def add_item(self, x=None, y=None, z=None, unpack_z=False):
	#
	def add_list(self, xyz, j_x=0, j_y=1, j_z=2):
		[self.add_to_bin(rw[j_x], rw[j_y], rw[j_z]) for rw in xyz]
		#
		return None
	#
	def add_to_bin(self, x=None, y=None, z=None):
		# in other words, add z to the bin located at (x,y). (x,y) can be off-center; we'll find the appropriate bin.
		#
		# x,y: coordinates of item; z is the item.
		if z==None: z=self.leaf_type()
		#x0,y0 = self.get_xybins(x,y,r_type='tuple')
		x0 = get_bin(x, dx=self.dx, x0=self.x0, bin_0=0)
		y0 = get_bin(y, dx=self.dy, x0=self.y0, bin_0=0)
		#
		#if not self.items.has_key(x0): self.items[x0]={}
		#if not self.items[x0].has_key(y0):
		#	self.items[x0][y0]=self.leaf_type()
		
		if not self.has_key(x0): self[x0]={}
		if not self[x0].has_key(y0):
			self[x0][y0]=self.leaf_type()
			
		#
		#if unpack_z==True and hasattr(z, '__len__'):
		#	self.items[x0][y0]+=z
		#else:
		#	self.items[x0][y0]+=[z]
		
		#self.items[x0][y0]+=z
		self[x0][y0]+=z
		#
		return None
	#
	def to_list(self):
		# should probably figure out the proper way to override list conversion method(s).
		#
		r_list = []
		#for X, Y_rw in self.items.iteritems():
		#	x = bin2x(X, self.dx, self.x0, self.bin_0)
		#	#
		#	r_list += [[x, bin2x(Y, self.dy, self.y0, self.bin_y_0), z] for Y,z in Y_rw.iteritems()] 
		
		#
		#return [[bin2x(X, self.dx, self.x0, self.bin_x_0), bin2x(Y, self.dy, self.y0, self.bin_y_0), z] for X,Y_rw in self.items.iteritems()  for Y,z in Y_rw.iteritems()]
		return [[bin2x(X, self.dx, self.x0, self.bin_x_0), bin2x(Y, self.dy, self.y0, self.bin_y_0), z] for X,Y_rw in self.iteritems()  for Y,z in Y_rw.iteritems()]
	#
	def to_array(self):
		#
		# note: in this simplified format, this will only work for float or int type objects.
		#
		# numpy.core.records.fromarrays(zip(*best_fit_array), names = ['section_id', 'tau', 'beta', 'sigma_tau', 'sigma_beta', 'mean_chi_sqr'], formats = [type(x).__name__ for x in best_fit_array[0]])
		if not hasattr(self.leaf_type, '__len__'):
			z_type_name = type(leaf_type).__name__
			return numpy.core.records.fromarrays(zip(*self.to_list()), names = ['x', 'y', 'z'], formats = ['f8', 'f8', type(leaf_type).__name__])
		#
		# otherwise, we've got string, list, tuple, etc. types. find the max length of all z entries.

		max_z_len = max([len(x) for x in zip(*self.to_list())[2]])
		#
		if isinstance(self.leaf_type, str):
			return numpy.core.records.fromarrays(zip(*self.to_list()), names = ['x', 'y', 'z'], formats = ['f8', 'f8', '|S%d' % max_z_len])
		#
		# otherwise, it's some sort of list or tuple. for now, let's just extend each row. we could also convert to a string... or we could figure out how to
		# wrap list types into a recarray (which i think can be done, but there are likely issues with reference, scope, etc.... or we could use a PANDAS object,
		# but then we change our syntax... and anyway, we should just add a separate function call for PANDAS (that is dependent upon PANDAS being present).
		#
		# we will have to assume that all the values in the row are of the same type. for now, we'll just let that break upon exception.
		z_type = float
		do_break=False
		#for j,rw_j in enumerate(self.items.values()):
		for j,rw_j in enumerate(self.values()):
			#print "array-checking row: ", j
			for k,rw_k in enumerate(rw_j.values()):
				if len(rw_k)>0:
					z_type = type(rw_k[0])
					do_break=True
					break
			if do_break:
				#print "break index: ", j,k
				break
		#z_type = type(self.items.values()[0].values()[0][0])\
		z_type = type(self.values()[0].values()[0][0])
		#print "z_type: ", z_type
		z_col_names = ['z_%d' % j for j in xrange(max_z_len)]
		#
		return numpy.core.records.fromarrays(zip(*[rw[0:2] + rw[2] for rw in self.to_list()]), names = ['x', 'y'] + z_col_names, formats = ['f8', 'f8'] + [z_type.__name__ for j in z_col_names])
		

class Bindex(dict):
	# a single, 1D bindex component. stack these together for n-D bindex objects?
	# so an N-D bindex can be made from a bindex of bindices.
	#
	def __init__(self, x0=0., dx=0.):
		self.x0=x0
		self.dx=dx
	#
#
def get_bin(x, dx=1., x0=0., bin_0=0):
	# return bin_id/index along a single axis.
	return int(bin_0) + int(round((x-x0)/dx))

def bin2x(bin_num, dx=1., x0=0., bin_0=0):
	# return the position x of a bin along one axis.
	return (bin_num - bin_0)*dx + x0

	
