'''
# earthquake_site.py
# mark yoder: mryoder@ucdavis.edu, mark.yoder@gmail.com
# dept. of physics, UC Davis and independent research.
#
# summary:
# site location objects for globalETAS and other applications. basic hierarchy:
#
# class site(): geo-spatial site base class. includes location, location translators (lat/lon, x/y, etc.),
#               distance calculators, and other geo-spatial bits.
# class site.grid_site(): specific to the characteristics of the lattice/grid being used. there might be
#                         more than one of these.
# class site.earthquake(): an earthquake is a site. also, of course, it includes: temporal coordinate (t), magnitude,
#                          L_r, A_r, theta_r, dip_r, ETAS parameters, etc.
'''
#
import math
deg_to_km = 111.1
#
class site(object):
	'''
	# base class for (geo-spatial) locations, including (ETAS) lattice sites, earthquakes, etc.
	# will have a "location" (center), extents, bounds, etc.
	#
	# and i think we can actually copy the corresponding/similar BASScast object for this.
	'''
	#
	def __init__(self, x=None, y=None, z=None, lat=None, lon=None, depth=None, f_to_lat_lon=None, f_from_lat_lon=None):
		# ... though these may change. basically, we want a universal locator (lat, lon) and some sort
		# of local/application specific coordinates (map transformations)
		# f_to_lat_lon, f_from_lat_lon are functions that define the default to/from lat-lon transformation.
		#
		self.x   = x
		self.y   = y
		self.z   = z
		self.lat = lat
		self.lon = lon
		self.depth = depth
		#
		if f_to_lat_lon   == None: f_to_lat_lon   = default_xy_to_lon_lat
		if f_from_lat_lon == None: f_from_lat_lon = default_lat_lon_to_xy
		#
		self.xy_from_lat_lon = f_from_lat_lon
		self.lat_lon_from_xy = f_to_lat_lon
		#
	#
	@property
	def loc(self):
		return [self.lon, self.lat]
	
	@property
	def loc_d(self):
		return {'lon': self.lon, 'lat':self.lat, 'depth':self.depth}
	

class grid_site(site):
	pass

class earthquake(site):
	pass

def default_xy_to_lon_lat(x=None, y=None):
	# a simple transformation from straight x,y to lat/lon.
	lat = None
	lon = None
	#
	if y! = None: lat = float(y)/deg_to_km
	if lat!= None and x!=None: lon = float(x)/(deg_to_km*math.cos(lat))
	#
	return (lon, lat)
#
def default_lon_lat_to_xy(lon=None, lat=None):
	x = None
	y = None
	#
	if lat!=None: y = float(lat)*deg_to_km
	if lon!=None and lat!=None: x = float(lon)*deg_to_km*math.cos(lat)
	#
	return (x,y)
	
