import numpy
import scipy
import math
import datetime as dtm
import pytz

import pylab as plt
import matplotlib.dates as mpd
import matplotlib as mpl

import ANSStools as atp
import multiprocessing as mpp

def mpi_cat(lons=[-123., -114.], lats=[31.5, 43.], mc=2.5, date_range=None, cat_len=10.*365.25):
	#
	date_range = (date_range or [None, None])
	date_range[0] = (date_range[0] or dtm.datetime.now(pytz.timezone('UTC'))-dtm.timedelta(days=cat_len))
	date_range[1] = (date_range[1] or dtm.datetime.now(pytz.timezone('UTC')))
	#
	#
	cat0 = atp.catfromANSS(lon=lons, lat=lats, minMag=mc, dates0=date_range, Nmax=None, fout=None, rec_array=True)
	#
	# but for now, something simpler...
	n_cpu = mpp.cpu_count()
	ary = mpp.Array(list(range(10*n_cpu)))
	
	#
	return cat0
