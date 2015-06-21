
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
	def __init__(self, x0=0., y0=0., dx=1., dy=1.):
		self.x0=x0
		self.y0=y0
		self.dx=dx
		self.dy=dy
		#
		self.items = {}
	#
	def get_xbin(self, x, dx=None, x0=None, bin_0=0):
		dx = (dx, self.dx)
		x0 = (x0, self.x0)
		#
		#return int((x-x0)/dx)
		return get_bin(x, dx=dx, x0=x0, bin_0=bin_0)
	#
	def get_ybin(self, y, dy=None, y0=None, bin_0=0):
		dy = (dy,self.dy)
		y0 = (y0, self.y0)
		return get_bin(y, dx=dy, x0=y0, bin_0=bin_0)
	#
	def get_xybins(self, x, y, dx=None, dy=None,  bin_x_0=0, bin_y_0=0, r_type='dict'):
		dx = (self.dx, dx)
		dy = (self.dy, dy)
		#
		x_bin = get_bin(x=x, dx=dx, x0=bin_x_0, bin_0=bin_x_0)
		y_bin = get_bin(x=x, dx=dy, x0=bint_y_0, bin_0=bin_y_0)
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
	def add_item(self, x=None, y=None, z=None, unpack_z=False):
		# x,y: coordinates of item; z is the item.
		x0,y0 = self.get_xybins(x,y,r_type='dict')
		#
		if not self.items.has_key(x0): self.items[x0]={}
		if not self.items[x0].has_key[y0]: self.items[x0][y0]=[]	# "leaf" objects are lists.
		#
		if unpack_z==True and hasattr(z, '__len__'):
			self.items[x0][y0]+=z
		else:
			self.items[x0][y0]+=[z]
		#
		return None
		
	

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
	return int(bin_0) + int((x-x0)/dx)

	
