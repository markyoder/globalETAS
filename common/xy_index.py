'''
# xy_index.py: xy binning index class and accessories.
# is it a proper "index"? maybe. (x,y) are presumably spatial coordinates,
# but they can be any (x,y). the basic structure is a dictionary of dictionaries containing lists
# of any objects. objects can be augmented like [x,y,object].
# or like {(x,y):object}
#
# my_xy ~ {0:{0:[(0,0) objects...], 1:[(0,1) objects], 2:[(0,2) objects]}, 1:{0:[(1,0) objects], 1:[(1,1) objects], 2:[(1,2) objects]}}
# 
# indices are defined by bin width: i_x = (x-x0)/dx
#
# note: i think i effectively rewrote this and called it bindex.py (binned incex); see globalETAS repository.
'''

class xy_index(dict):
	
	#
	def __init__(self, data_in=None, x0=0., y0=0., dx=.1, dy=.1):
		self.x0 = x0
		self.y0 = y0
		self.dx = dx
		self.dy = dy
		#
		self.data_in=data_in		# though eventually, we'll skip just adding this to populate the dict.
		#
		#self.catalog={}
		#
		#
	#
	def add_item(self, x, y, item=None):
		i_x = self.get_x_index(x)
		i_y = self.get_y_index(y)	# note that y0, dy (and x0, dx) will by default come from the class namespace.
		#
		if (i_x in self)==False:
			self[i_x]={}
		if (i_y in self[i_x])==False:
			self[i_x][i_y]={}
		#
	#
	#
	def get_x_index(self, x, x0=None, dx=None):
		if x0==None: x0=self.x0
		if dx==None: dx=self.dx
		#
		return (float(x)-x0)/float(dx)
	#
	def get_y_index(self, y, y0=None, dy=None):
		# note: this can be problematic with floating point errors. i think that can be fixed using round().
		if y0==None: y0=self.y0
		if dy==None: dy=self.dy
		#
		return (float(y)-y0)/float(dy)
	#
	
