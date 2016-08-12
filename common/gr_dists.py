import matplotlib
import numpy
import pylab as plt
import math
import random
import datetime as dtm
import pytz


from mpl_toolkits.mplot3d import Axes3D
#import matplotlib.colors as mcolor
#import matplotlib.colorbar as mcolorbar
#import matplotlib.lines as mlines
#import matplotlib.patches as mpatches

#import mpl_toolkits.basemap as bmp
#from mpl_toolkits.basemap import Basemap

import matplotlib as mpl
from matplotlib import cm
import itertools
plt.ion()
#
import math
import h5py
#
import numpy
import scipy
import scipy.optimize as spo
from scipy.spatial import KDTree
from scipy.spatial import cKDTree

import ANSStools

colors_ =  mpl.rcParams['axes.color_cycle']

def gr_thing(lats=[31., 43.], lons=[-126., -114.], mc=3.0, t0=dtm.datetime(1990,1,1, tzinfo=pytz.timezone('UTC')), t1=dtm.datetime.now(pytz.timezone('UTC')), m0=6.0, b0=1.0):
	N0=1000.
	cat = ANSStools.catfromANSS(lon=lons, lat=lats, minMag=mc, dates0=[t0, t1], Nmax=None, fout=None, rec_array=True)
	#
	# nominally, we'd get some sort of spatial subcats, but we won't bother just yet.
	# now, get subcats between each of the m>m_0 events; plot, fit, etc.
	#
	my_axes = []
	for fnum in range(4):
		plt.figure(fnum)
		plt.clf()
		my_axes += [plt.gca()]
		if fnum<2: my_axes[-1].set_yscale('log')
	b=b0		# though later, we might fit this
	
	#
	#for best speed and memory performance, we should just walk through this, but this will be easier to code:
	#
	m_index = [[j,m] for j,m in enumerate(cat['mag']) if m>=m0]
	abmag = []
	#
	for j,rw in enumerate(m_index[1:]):
		# color for plotting (optinal):
		this_color = colors_[j%len(colors_)]
		#
		these_mags = sorted(cat['mag'][m_index[j][0]+1:rw[0]].copy())
		#
		if len(these_mags)<2:
			print("insufficient sequence length: %d" % len(these_mags))
			continue
		#
		these_mags.reverse()
		Ns = numpy.arange(1., len(these_mags)+1, 1.)
		my_axes[0].plot(these_mags, Ns, 'o-')
		#
		# now, scale these for stacking. scale counts to m0.
		m_start = cat['mag'][m_index[j][0]]
		m_stop  = cat['mag'][rw[0]]
		m_mean = .5*(m_start+m_stop)
		#
		#Ns_prime = [n*(10.**(b*(m0-m))) for n,m in zip(Ns, these_mags)]
		#Ns_prime = [n*(10.**(b*(m0-m_mean))) for n,m in zip(Ns, these_mags)]
		Ns_prime = [n*Ns[0]/Ns[-1] for n in Ns]
		delta_log_ratio = abs(math.log10(Ns_prime[0]/Ns_prime[-1]))/math.log10(Ns[-1])
		ms_prime = [m*delta_log_ratio for m in these_mags]
		#my_axes[1].plot(these_mags, Ns_prime, '.-')
		my_axes[1].plot(ms_prime, Ns_prime, '.-')
		#
		# fit for b values:
		#print "len[%d]: (%d/%d:%f): %d " % (j, m_index[j][0], rw[0], rw[1], len(these_mags))
		#
		fits = numpy.linalg.lstsq(numpy.array([[m,1.0] for m in these_mags]), numpy.log10(Ns))[0]
		print("fits: ", fits[0])
		abmag += [[fits[1], fits[0], m_index[j][1], rw[1]]]
		#
		my_axes[2].plot([m_index[j][1], rw[1]], [fits[1], fits[1]], 'o-', color=this_color)
		#my_axes[2].plot([fits[1]], [rw[1]], 'o', color=this_color)
		#
		my_axes[3].plot([m_index[j][1], rw[1]], [fits[0], fits[0]], 'o-', color=this_color)
		#my_axes[3].plot([fits[0]], [rw[1]], 'o', color=this_color)
		
	#
	print("b-values:\n median: %f, mean: %f \\pm %f" % (numpy.median([rw[1] for rw in abmag]), numpy.mean([rw[1] for rw in abmag]), numpy.std([rw[1] for rw in abmag])))
	#
	plt.figure(2)
	#plt.clf()
	#plt.plot(*zip(*[[rw[0], rw[2]] for rw in abmag]), marker='o', ls='')
	#plt.plot(*zip(*[[rw[0], rw[3]] for rw in abmag]), marker='o', ls='')
	plt.ylabel('a-val')
	plt.xlabel('magnitude $m$')
	plt.figure(3)
	#plt.clf()
	#plt.plot(*zip(*[[rw[1], rw[2]] for rw in abmag]), marker='o', ls='')
	#plt.plot(*zip(*[[rw[1], rw[3]] for rw in abmag]), marker='o', ls='')
	plt.ylabel('b-val')
	plt.xlabel('magnitude $m$')
	ax_bcdf = my_axes[3].twiny()
	bs = [rw[1] for rw in abmag]
	bs.sort()
	ax_bcdf.plot(range(1,len(bs)+1), bs, 's-', lw=2, alpha=.7)
	ax_bcdf.set_ylabel('$N<b$')
	
