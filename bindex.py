
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
import matplotlib.dates as dtm
import multiprocessing as mpp


class Bindex2D(object):
	
	# for now, let's just start with a 2D bincex and we'll generalize later
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
		self.items = {}
	#
	def __iter__(self):
		return self.items.__iter__
	def __sizeof__(self):
		return sum([len(rw) for rw in self.items.iteritems()])
	#
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
	def get_xybins(self, x, y, dx=None, dy=None,  bin_x_0=0, bin_y_0=0, r_type='dict'):
		dx = (self.dx or  dx)
		dy = (self.dy or dy)
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
		return self.items[x0][y0]
	#
	#def add_item(self, x=None, y=None, z=None, unpack_z=False):
	#
	def add_list(self, xyz, j_x=0, j_y=1, j_z=2):
		[self.add_to_bin(rw[j_x], rw[j_y], rw[j_z]) for rw in xyz]
		#
		return None
	#
	def add_to_bin(self, x=None, y=None, z=None):
		# x,y: coordinates of item; z is the item.
		if z==None: z=self.leaf_type()
		#x0,y0 = self.get_xybins(x,y,r_type='tuple')
		x0 = get_bin(x, dx=self.dx, x0=self.x0, bin_0=0)
		y0 = get_bin(y, dx=self.dy, x0=self.y0, bin_0=0)
		#
		if not self.items.has_key(x0): self.items[x0]={}
		if not self.items[x0].has_key(y0):
			self.items[x0][y0]=self.leaf_type()
		#
		#if unpack_z==True and hasattr(z, '__len__'):
		#	self.items[x0][y0]+=z
		#else:
		#	self.items[x0][y0]+=[z]
		self.items[x0][y0]+=z
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
		return [[bin2x(X, self.dx, self.x0, self.bin_x_0), bin2x(Y, self.dy, self.y0, self.bin_y_0), z] for X,Y_rw in self.items.iteritems()  for Y,z in Y_rw.iteritems()]
	#
	def to_array(self):
		#
		# numpy.core.records.fromarrays(zip(*best_fit_array), names = ['section_id', 'tau', 'beta', 'sigma_tau', 'sigma_beta', 'mean_chi_sqr'], formats = [type(x).__name__ for x in best_fit_array[0]])
		return numpy.core.records.fromarrays(zip(*self.to_list()), names = ['x', 'y', 'z'], formats = ['f8', 'f8', 'f8'])

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
	return int(bin_0) + int(round((x-x0)/dx))

def bin2x(bin_num, dx=1., x0=0., bin_0=0):
	return (bin_num - bin_0)*dx + x0

	
