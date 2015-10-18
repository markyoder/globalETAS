import datetime as dtm
import matplotlib.dates as mpd
import pytz
tzutc = pytz.timezone('UTC')

import operator
import math
import random
import numpy
import scipy
#import scipy.optimize as spo
#import os
#from PIL import Image as ipp
import multiprocessing as mpp
#
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.mpl as mpl
#
#import shapely.geometry as sgp
#
from mpl_toolkits.basemap import Basemap as Basemap
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from geographiclib.geodesic import Geodesic as ggp
#
import ANSStools as atp

deg2km = 111.1


class ETAS_map(object):
	
	def __init__(self, lat_min=31., lat_max=42., lon_min=-125., lon_max=-115., mc=3.0, grid_dx=10., grid_dy=10., grid_units='km', date_min=dtm.datetime(1990,1,1, tzinfo=pytz.timezone('UTC')), date_max=None, N_max=None, cat_out_file=None):
		# make a lattice of grid sites (really points with some spacing; for now, assume rectalinearly equal spacing).
		# get a catalog of earthquakes; create an auxilary catalog as well.
		#
		self.lat_min=lat_min
		self.lat_max=lat_max
		self.lon_min=lon_min
		self.lon_max=lon_max
		self.mc=mc
		self.grid_dx = grid_dx
		self.grid_dy = grid_dy
		self.grid_units=grid_units
		self.date_min=date_min
		self.date_max=date_max		# note: these may differ from the values in the catalog, so we might want to write some code/commentary to be specific.
		self.N_max=N_max
		self.cat_out_file=cat_out_file
		#
		self.catalog = atp.catfromANSS(lon=[lon_min, lon_max], lat=[lat_min, lat_max], minMag=mc, dates0=[date_min, date_max], Nmax=N_max, fout=cat_out_file, rec_array=True)
		#
		
		
