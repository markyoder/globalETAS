import math
from math import *	
#from scipy import *
import scipy
from pylab import *
from matplotlib import *
#import numpy.fft as nft
import scipy.optimize as spo
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.dates as mpd

# maping bits:
import matplotlib	# note that we've tome from ... import *. we should probably eventually get rid of that and use the matplotlib namespace.
matplotlib.use('Agg')
#from matplotlib.toolkits.basemap import Basemap
from mpl_toolkits.basemap import Basemap as Basemap
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
####
#
import string
import sys
#from matplotlib import *
#from pylab import *
import os
import random
import time
#
# gamma function lives here:
#import scipy.special
from scipy.special import gamma
#from scipy.optimize import leastsq
from matplotlib import axis as aa
from threading import Thread
#
#
import datetime
import calendar

import operator

#import yodapy as yp

#import omoribits as ob

###############################################
# rbIntervals.py                              #
#
# mark yoder
# UC Davis
# april 2009
#
# compute record breaking intervals between aftershocks
# namely, we will use events following the Parkfield Eq.
# catalog, etc. will be described inline.
###############################################

class intervalRecordBreaker:
	fullCat=[]
	shockCat=[]		# [[evDateTime, lat, lon, mag, a, b], [row1], ...]. where a,b are the eliptical x,y coordinates
	#catname='parkcat.cat'
	catname='cats/parkfieldfull10yrs.cat'
	eventDtTime=None
	
	# coordinate transformation variables:
	tLat=None # 35.9
	tLon=None # -120.5
	tTheta=None # 40.0		#47?	note: tTheta is the angle CCW of the x' (transformed) axis from the x axis.
	tA=None # .4		# ellipse axes
	tB=None # .15
	
	
	def __init__(self, catFname='cats/parkfieldfull10yrs.cat', theta=40.0, clat=35.9, clon=-120.5, ra=.4, rb=.15, eventDate=datetime.datetime(2004,9,28, 17,15,24), maxDate=datetime.datetime(2009,9,28, 17,15,24), skipSeconds=0):
		# don't know yet. we'll need a catalog and bits like that...
		self.initialize(catFname, theta, clat, clon, ra, rb)
	
	def initialize(self, catFname='cats/parkfieldfull10yrs.cat', theta=40.0, clat=35.9, clon=-120.5, ra=.4, rb=.15, eventDate=datetime.datetime(2004,9,28, 17,15,24), maxDate=datetime.datetime(2009,9,28, 17,15,24), skipSeconds=0):
		self.catname=catFname
		self.tLat=clat
		self.tLon=clon
		self.tTheta=theta
		self.tA=ra
		self.tB=rb
		
		if catFname==None: return None
		
		if theta!=None:
			print("set aftershock catalog")
			self.setAftershockCatalog(catFname, theta, clat, clon, ra, rb, eventDate, maxDate, skipSeconds )
		else:
			self.setNormalCat(catFname)	# default to socal...
	
	# 29 apr 2009 yoder: we want to do RB for a 2 year or so period in socal. create a normal active catalog.
	# (but then we abandoned this for a renewed approach; see recordBreaker.py)
	def setNormalCat(self, catFname=None, minlat=32, maxlat=36.5, minlon=-125, maxlon=-115, minDt=datetime.datetime(1984,0o1,0o1), maxDt=datetime.datetime(1985,12,31)):
		# set a normal (rectangular) catalog, no ellipse.
		if catFname!=None:
			self.setFullCat(catFname)
			
		tempcat=self.fullCat
		if minDt==None: minDt=tempcat[0][0]
		if maxDt==None: maxDt=tempcat[0][-1]
		
		for row in tempcat:
			# no rotations; just transfer the events in the space,time,mag space.
			if row[0]<=minDt: continue
			if row[0]>maxDt: break
		self.shockCat+=[[row[0], row[1], row[2], row[3], row[2], row[1]]]
		#
	#
	##################
	
	#def setAftershockCatalog(self, catFname=None, theta=tTheta, clat=35.9, clon=-120.5, ra=tA, rb=tB, eventDate=datetime.datetime(2004,9,28, 17,15,24), maxDate=datetime.datetime(2009,9,28, 17,15,24), skipNevents=0):
	def setAftershockCatalog(self, catFname=None, theta=tTheta, clat=35.9, clon=-120.5, ra=tA, rb=tB, eventDate=datetime.datetime(2004,9,28, 17,15,24), maxDate=datetime.datetime(2009,9,28, 17,15,24), skipSeconds=0):
		# keep all events in catFname defined by the following ellipse:
		# skipNevents: number of events after main-shock to skip (avoid messiness during mainshock)
		self.eventDtTime=eventDate
		self.catname=catFname
		self.shockCat=[]
		#print "event (start) date, catname: %s, %s, %s" % (eventDate, catFname, self.catname)
		#
		if catFname!=None:
			print(("setting catalog from [%s]" % catFname))
			self.setFullCat(catFname)
		#print "first cat row: %s" % self.fullCat[0]
		#
		tempCat=self.fullCat
		if eventDate==None:
			# we want the whole catalog (no min date):
			#print tempCat[-1][0]
			eventDate=tempCat[0][0]
		
		#nEventsSinceMS=0
		for row in tempCat:
			# rotate each element into our aftershock axis, is it in the ellipse?
			#print row[0]
			# 2009-06-04 yoder: skipping nEvents since mainshock... (note, at this time we're already excluding up to and including MS)
			#nEventsSinceMS+=1
			#if nEventsSinceMS<=skipNevents: continue
			# or some period of time-time:
			#if row[0]<=eventDate: continue
			if maxDate>=eventDate and row[0]<=eventDate+datetime.timedelta(seconds=skipSeconds): continue		# 1_day = 86400 seconds
			if maxDate<=eventDate and row[0]>=eventDate-datetime.timedelta(seconds=skipSeconds): continue		# 22 july 2009 yoder: facilitat reverse RB
			
			# remove the .01 day immediately after the event from the data-set:
			#if row[0]<=eventDate+datetime.timedelta(seconds=864): continue
			#if row[0]<=eventDate+datetime.timedelta(seconds=1600): continue
			if maxDate>=eventDate and row[0]>maxDate: break
			if maxDate<=eventDate and row[0]<maxDate: break		# 22 july 2009 yoder: facilitat reverse RB
			
			newVec=self.faultTransform(row[2], row[1], clat, clon, theta, ra, rb)
			#
			# is the rotated vector in our ellipse?
			if abs(newVec[0])>ra: continue
			Y=ellipseY(newVec[0], ra, rb)
			if abs(newVec[1])>Y: continue
			# dtm, lat, lon, mag, tX, tY 		(note this is like y,x, x`, y` for the space coordinates).
			self.shockCat+=[[row[0], row[1], row[2], row[3], newVec[0], newVec[1]]]	
		#
			
	def setFullCat(self, catFname=None, minmag=1.25):
		if catFname==None: catFname=self.catname
		self.fullCat=[]
		
		#print "fullcat catname: %s" % catFname
		f=open(catFname)
		nrows=0
		for row in f:
			if row=='\n' or row[0]=='#': continue
			#if nrows==0: print row
			nrows+=1
			delim=' '
			if '\t' in row: delim='\t'
			#
			rowl=row.split(delim)
			# 2008 12 2 -117.318 35.97 3.06
			mag=float(rowl[4])
			if mag<minmag: continue
			delim=rowl[0][4]
			thisDtm=datetimeFromStrings(rowl[0], rowl[1], delim)
			#
			# skip (probable) duplicate entries.
			if nrows>1:
				if thisDtm==self.fullCat[-1][0]: continue
				
			# dtm, lat, lon, mag
			self.fullCat+=[[thisDtm, float(rowl[2]), float(rowl[3]), mag]]
			#print 'catrow: %s' % self.fullCat[-1]
		print(("len of fullcat: %d, %s" % (len(self.fullCat), self.fullCat[0])))
	
	def getMainEvent(self, cat=None):
		# return catalog row of max magnitude (epicenter location (more or less)). note, by default we use ths shock-cat because it will
		# be faster than using the fullCat AND, the fullCat is likely to have other large earthquakes.
		if cat==None: cat=self.shockCat
		if len(cat)==0: cat=self.fullCat
		maxMag=cat[0][3]
		maxIndex=0
		for i in range(len(cat)):
			#print i, maxMag, maxIndex, cat[i][3]
			if cat[i][3]>maxMag:
				maxIndex=i
				maxMag=cat[i][3]
		#
		return cat[maxIndex] + [maxIndex]
			
	def fitAftershocksToOmori(self, minmag=1.2, p=None):
		# in any case, this current form seems to work pretty well so long as we guess the initial prams within an order of magnitude or so
		# for example: (50, .001, 1.5)
		#
		# p=scipy.array([0,1])
		# plsq=spo.leastsq(linRes, p, args=(scipy.array(y), scipy.array(x)), full_output=0, maxfev=20000)
		# omori's law: n(t) = k/(c+t)**p --> p[0]/(p[1]+t)**p[2]
		#p=[1.0/4000, 1/1000, 1]		# on the first day, we get 12 events; figure it takes a few hours to set up; this should give us reasonalbe starting prams.
		if p==None or len(p)<3:
			p=[50, .02, 1.5]	# in days...
			
		x=[]
		y=[]
		nthEvent=0
		for rw in self.shockCat:
			if rw[3]<minmag: continue
			#print "use %s" % x
			x+=[datetimeToFloat(rw[0])-datetimeToFloat(self.eventDtTime)]
			y+=[nthEvent]
			nthEvent+=1
		
		plsq=spo.leastsq(omoriCumRes, p, args=(scipy.array(y), scipy.array(x)), full_output=0, maxfev=200000)
		p=plsq[0]
		for i in range(1000):
			# loop pseudo-recursively to converge:
			plsq=spo.leastsq(omoriCumRes, p, args=(scipy.array(y), scipy.array(x)), full_output=0, maxfev=200000)
			p=plsq[0]
		#print plsq
		#return plsq[0]
		yfit=[]
		chiSqr=0
		ndof=-3
		
		for X in x:
			yfit+=[fomoriCum(X,p)]
			#
			chiSqr+=(yfit[ndof+3]-y[ndof+3])
			ndof+=1
		chiSqr/=ndof
		
		print((p[0], p[1], p[2], chiSqr))
		return [p.tolist()+[chiSqr], [x,y,yfit]]

	def plotfitAftershocksToOmori(self, minmag=1.2, p=None):
		prams=self.fitAftershocksToOmori(minmag, p)
		p=prams[0]
		x=prams[1][0]
		y=prams[1][1]
		yfit=prams[1][2]
		
		plt.figure(1)
		plt.plot(x,y,'.-', label='data')
		plt.plot(x,yfit, '.-', label='fit')
		plt.legend(loc='lower right')
		plt.show()
				
		return [p, [x,y,yfit]]
	
	def xyPlotShocks(self):
		#lats=map(operator.itemgetter(1), rbi.shockCat)
		plt.figure(0)
		plt.clf()
		plt.plot(list(map(operator.itemgetter(2), self.shockCat)), list(map(operator.itemgetter(1), self.shockCat)), '.')
		plt.show()
	
	def xyPlotFull(self):
		#lats=map(operator.itemgetter(1), rbi.shockCat)
		plt.figure(0)
		plt.clf()
		plt.plot(list(map(operator.itemgetter(2), self.fullCat)), list(map(operator.itemgetter(1), self.fullCat)), '.')
		plt.show()
	
	def getLatLonRange(self, cat=None, latloncols=[1,2]):
		if cat==None: cat=self.fullCat
		if latloncols==None: latloncols=[1,2]	# latitude, lon cols of catalog (order is lat, lon).
		#
		minLat=cat[0][latloncols[0]]
		maxLat=cat[0][latloncols[0]]
		minLon=cat[0][latloncols[1]]
		maxLon=cat[0][latloncols[1]]
		#
		for rw in cat:
			thisLat=rw[latloncols[0]]
			thisLon=rw[latloncols[1]]
			#
			if thisLat>maxLat: maxLat=thisLat
			if thisLat<minLat: minLat=thisLat
			if thisLon>maxLon: maxLon=thisLon
			if thisLon<minLon: minLon=thisLon
		#
		return [[minLat, minLon], [maxLat, maxLon]]
	
	def xyPlotCatsMap(self, doShow=True, doSave=False, saveName='catalogPlot.png', epicenter=None, legendLoc='upper left'):
		
		if epicenter==None: epicenter=[self.tLon, self.tLat]
		fcat=[]
		scat=[]
		for rw in self.shockCat:
			scat+=[rw[0:4]]
		for rw in self.fullCat:
			if rw not in scat: fcat+=[rw]
		#return [scat, fcat]
	
		f0=plt.figure(0)	
		plt.clf()
		#		
		#set up map:
		llr=self.getLatLonRange()	# latLonRange
		cntr=[float(llr[0][0])+(llr[1][0]-float(llr[0][0]))/2.0, float(llr[0][1])+(llr[1][1]-float(llr[0][1]))/2.0]
		catmap=Basemap(llcrnrlon=llr[0][1], llcrnrlat=llr[0][0], urcrnrlon=llr[1][1], urcrnrlat=llr[1][0], resolution ='l', projection='tmerc', lon_0=cntr[1], lat_0=cntr[0])
		canvas=FigureCanvas(f0)
		catmap.ax=f0.add_axes([0,0,1,1])
		f0.set_figsize_inches((8/catmap.aspect,8.))
		#
		catmap.drawcoastlines(color='gray')
		catmap.drawcountries(color='gray')
		catmap.fillcontinents(color='beige')
		xfull, yfull=catmap(list(map(operator.itemgetter(2), fcat)), list(map(operator.itemgetter(1), fcat)))
		xshock, yshock=catmap(list(map(operator.itemgetter(2), scat)), list(map(operator.itemgetter(1), scat)))
		epx, epy=catmap(epicenter[0], epicenter[1])
		catmap.plot(xfull, yfull, 'g+', label='Full Catalog')
		catmap.plot(xshock, yshock, 'b.', label='Aftershock zone')
		catmap.plot(epx, epy, 'ro')
		
		canvas.print_figure(saveName)
		
		#
		#ax=plt.gca()
		el = Ellipse((self.tLon, self.tLat), 2.0*self.tA, 2.0*self.tB, -self.tTheta, facecolor='b', alpha=0.4)
		#catmap.ax.add_artist(el)
		#ax.add_artist(el)
		#
		#plt.plot(map(operator.itemgetter(2), self.fullCat), map(operator.itemgetter(1), self.fullCat), '+')
		#plt.plot(map(operator.itemgetter(2), self.shockCat), map(operator.itemgetter(1), self.shockCat), '.')
		#plt.plot(map(operator.itemgetter(2), fcat), map(operator.itemgetter(1), fcat), '+', label='Full Catalog')
		#plt.plot(map(operator.itemgetter(2), scat), map(operator.itemgetter(1), scat), '.', label='Aftershock zone')
		#plt.plot([epicenter[0]], [epicenter[1]], 'ro', label='epicenter')
		plt.legend(loc=legendLoc, numpoints=1)
		if doSave: plt.savefig(saveName)
		if doShow: plt.show()

	def xyPlotCatalogs(self, doShow=True, doSave=False, saveName='catalogPlot.png', epicenter=None, legendLoc='upper left'):
		
		if epicenter==None: epicenter=[self.tLon, self.tLat]
		fcat=[]
		scat=[]
		for rw in self.shockCat:
			scat+=[rw[0:4]]
		for rw in self.fullCat:
			if rw not in scat: fcat+=[rw]
		#return [scat, fcat]
	
		plt.figure(0)	
		plt.clf()
		#
		ax=plt.gca()
		#el = Ellipse((-72.533,18.457), 1.25, .25, 15, facecolor='r', alpha=0.5)
		el = Ellipse((self.tLon, self.tLat), 2.0*self.tA, 2.0*self.tB, -self.tTheta, facecolor='b', alpha=0.4)
		ax.add_artist(el)
		#
		#plt.plot(map(operator.itemgetter(2), self.fullCat), map(operator.itemgetter(1), self.fullCat), '+')
		#plt.plot(map(operator.itemgetter(2), self.shockCat), map(operator.itemgetter(1), self.shockCat), '.')
		plt.plot(list(map(operator.itemgetter(2), fcat)), list(map(operator.itemgetter(1), fcat)), '+', label='Full Catalog')
		plt.plot(list(map(operator.itemgetter(2), scat)), list(map(operator.itemgetter(1), scat)), '.', label='Aftershock zone')
		plt.plot([epicenter[0]], [epicenter[1]], 'ro', label='epicenter')
		plt.legend(loc=legendLoc, numpoints=1)
		if doSave: plt.savefig(saveName)
		if doShow: plt.show()	
	
	def plotMagsIntervals(self, doShow=True, doSave=False, saveName='catalogMagsInts', pltTitle='Catalog Seismicity', nTime=False, minmag=1.5):
		fnameFull="%sFull.png"
		fnameShock="%sShock.png"
		
		cats=[self.fullCat, self.shockCat]
		fignums=list(range(len(cats)))
		#
		fsize=20
		#minmag=1.1	# we should get this from the data
		for fnum in fignums:
			currcat=cats[fnum]
			
			plt.figure(fnum)
			plt.clf()
			ax0=plt.axes([.1,.1,.85, .20])
			
			#
			#plt.title("Magnitudes")
			if nTime:
				plt.xlabel("Number of Events (n)", fontsize=fsize)
			else:
				plt.xlabel("date", fontsize=fsize)
			
			plt.ylabel("Mags", fontsize=fsize)
			#
			#ax1=plt.axes([.1, .4, .85, .50], sharex=ax0)
			ax1=plt.axes([.1, .55, .85, .35], sharex=ax0)
			#plt.title(pltTitle, fontsize=fsize)
			if currcat==self.fullCat: plt.title("Parkfield 4x4 catalog", fontsize=fsize)
			if currcat==self.shockCat: plt.title("Parkfield aftershock catalog", fontsize=fsize)
			plt.ylabel("intervals (days)", fontsize=fsize)
			plt.xlabel("")
			#
			intervals=[]
			intbars=[]
			magbars=[]
			for i in range(len(currcat)):
				thisx=currcat[i][0]
				magbars+=[[thisx,minmag], [thisx, currcat[i][3]], [thisx,minmag]]
				if i==0: continue
				intbars+=[[thisx,0], [thisx, date2num(currcat[i][0])-date2num(currcat[i-1][0])], [thisx,0]]
			
			#zeros=[0]
			#for i in xrange(1, len(currcat)):
			#	intervals+=[date2num(currcat[i][0])-date2num(currcat[i-1][0])]
			#	#zeros+=[0]
			#
			#ax0([datetime.datetime(2000,01,01), datetime.datetime(2010,12,31)])
			
			#ax0.plot_date(map(operator.itemgetter(0), self.fullCat), map(operator.itemgetter(3), self.fullCat), '-')
			#ax1.plot_date(map(operator.itemgetter(0), currcat[1:]), intervals, '-')
			
			if nTime==False:
				ax0.plot_date(list(map(operator.itemgetter(0), magbars)), list(map(operator.itemgetter(1), magbars)), '-')
				ax1.plot_date(list(map(operator.itemgetter(0), intbars)), list(map(operator.itemgetter(1), intbars)), '-')
				ax1.set_yscale('log')
				a = plt.gca()
				#a.set_xlim([currcat[0][0]-datetime.timedelta(days=20), currcat[-1][0]+datetime.timedelta(days=20)])
				
			if nTime==True:
				y1=list(map(operator.itemgetter(1), magbars))
				y2=list(map(operator.itemgetter(1), intbars))
				ax0.plot(list(range(len(y1))), y1, '-')
				ax1.plot(list(range(len(y2))), y2, '-')
				ax1.set_yscale('log')
				a = plt.gca()
				a.set_xlim([0, int(1.1*len(y2))])		
			#
			if doSave:
				if currcat==self.fullCat: catTag='-full'
				if currcat==self.shockCat: catTag='-shock'
				plt.savefig('%s-%s.png' % (saveName, catTag))
		if doShow: plt.show()
		
		
	def GRshock(self, doShow=True, fname='GRdist.png', plotTitle="Aftershock Region Magnitude Distribution"):
		# [[evDateTime, lat, lon, mag, a, b], [row1], ...]
		mags=list(map(operator.itemgetter(3), self.shockCat))
		return self.GRdist(mags, doShow, fname, plotTitle)
		
	def GRfullcat(self, doShow=True, fname='GRdist.png', plotTitle="Full Catalog Magnitude Distribution"):
		mags=list(map(operator.itemgetter(3), self.fullCat))
		return self.GRdist(mags, doShow, fname, plotTitle)
	
	def GRdist(self, mags, doShow=True, fname='GRdist.png', plotTitle="Magnitude Distribution"):
		# cat: a 1D array of magnitudes
		mags.sort()
		# get rid of biggest event (probably a large off-GR earthquake):
		mags.pop()
		#mags.reverse()
		#print mags
		#print len(mags)
		
		if doShow==True or fname!=None:
			# make a plot and show and/or save
			#Y=range(1, len(mags)+1)
			Y=list(map(float, list(range(1, len(mags)+1))))
			Y.reverse()
			#print Y
			#print len(Y)
			plt.figure(0)
			plt.clf()
			plt.semilogy(mags, Y, '.-')
			plt.xlabel("Magnitude, m")
			plt.ylabel("Number of Events, n")
			plt.title(plotTitle)
			if fname!=None: plt.savefig(fname)
			if doShow: plt.show()
		
		return mags
		#
		
	def fitOmoriRange(self, minmag=1.25, maxmag=4.5, dmag=.1, p=None):
		pramses=[[], [], [], [], []]

		fname="parkfieldOmoriFit.dat"
		f=open(fname,'w')
		f.write("#k\tc\tp\trChiSqr\n")
		f.close()
		
		while minmag<=maxmag:
			print(("fitting for mag %f." % minmag))
			prams=self.fitAftershocksToOmori(minmag,p)[0]
			print(prams)
			f=open(fname,'a')
			f.write("%f\t%f\t%f\t%f\t%f\n" % (minmag, prams[0], prams[1], prams[2], prams[3]))
			f.close()
			pramses[0]+=[minmag]
			pramses[1]+=[prams[0]]
			pramses[2]+=[prams[1]]
			pramses[3]+=[prams[2]]
			pramses[4]+=[prams[3]]
			
			minmag+=dmag
		return pramses
		
	def faultTransform(self, x, y, Lat=tLat, Lon=tLon, theta=tTheta, A=tA, B=tB):
		# x,y to transform via blah, blah.
		#
		theta=deg2rad(float(theta))
		xprime = (x-Lon)*cos(theta) - (y-Lat)*sin(theta)
		yprime = (x-Lon)*sin(theta) + (y-Lat)*cos(theta)
		
		return [xprime, yprime]
	#
	def getBigShocks(self, minmag=1.5, bigmag=5.0, cat=None):
		if cat==None: cat=self.shockCat
		rnum=0
		bigShocks=[]
		for rw in cat:
			if rw[3]<minmag: continue	# don't count rnum...
			if rw[3]>=bigmag: bigShocks+=[[rnum] + rw]
			rnum+=1
		return bigShocks
		
	def getIntervalRatios(self, minmag=3.0, windowLen=10, cat0=None, deltaipos=1, logZ=1.0):
		# deltaipos: event resolution; how many events to advance each step.
		if cat0==None: cat=self.shockCat
		if logZ==None:
			logZ = math.log10(float(windowLen))	
		logZinv=1.0/logZ	# this way, logZ is more intuitive; logZ=log(winlen)
		#
		# problem: if you give this function a catalog with magnitudes below the min-mag, the getRBintervals() function can break (it has a bit
		# that strips out sub-m0 events). SO, do that here.
		cat=[]
		#icount=0
		for rw in cat0:
			#print icount, rw[3], minmag, rw[3]>=minmag
			#icount+=1
			if rw[3]>=minmag: cat+=[rw]
		#
		ipos=0
		rbRatios=[]
		#print ipos, len(cat)-windowLen, minmag, len(cat), len(cat0)
		cat0=None
		while ipos<(len(cat)-windowLen):
			thisRatios=self.getRBintervals(minmag, cat[ipos:(ipos+windowLen)])	# note: if we want to do backwards RB, using cat[i:i+l].reverse() is a nice trick.
			#rbRatios+=[[ipos+windowLen, cat[(ipos+windowLen)][0], float(thisRatios[0][2][-1])/float(thisRatios[1][2][-1])]]
			r=(float(thisRatios[0][2][-1])/float(thisRatios[1][2][-1]))
			rbRatios+=[[ipos+windowLen, cat[(ipos+windowLen)][0], r**logZinv]]
			#ipos+=deltaipos
			ipos+=1
		# return [[n, dt, r]]
		thisRatios=None
		return rbRatios

	def plotIntervalRatiosAx(self, minmag=3.0, windowLen=10, cat0=None, hitThreshold=1.0, bigmag=5.0, thisAx=None, ratios=None, deltaipos=1, avlen=1, mainEV=None, logZ=1.0, rbLegLoc='best'):
		# avlen=10
		# eventually, this will probably be the sole version of this function....
		# this creates a plot of interval ratios in/on a specified axis. use this function by itself or to make complex figures.
		#
		if thisAx==None:
			f0=plt.figure()
			thisAx=f0.gca()
		#	
		legLoc=rbLegLoc
		eventName="Event RB ratios"
		if cat0==None: cat0=self.shockCat
		# getIntervalRatios(self, minmag=3.0, windowLen=10, cat0=None, deltaipos=1):
		#if ratios==None: ratios=self.getIntervalRatios(minmag, windowLen, cat0, deltaipos)
		if ratios==None: ratios=self.getIntervalRatios(minmag=minmag, windowLen=windowLen, cat0=cat0, deltaipos=deltaipos, logZ=logZ)
		fdts=[]
		for rw in ratios:
			 #fdts+=[rw[1].toordinal() + float(rw[1].hour)/24 + rw[1].minute/(24*60) + rw[1].second/(24*3600) + rw[1].microsecond/(24*3600000000)]
			 fdts+=[mpd.date2num(rw[1])]
		plaindts=list(map(operator.itemgetter(1), ratios))
		
		if mainEV==None: mainEV=self.getMainEvent(cat0)
		
		eventFloatDt=mpd.date2num(mainEV[0])
		#
		#f=plt.figure(fignum)
		#theseFigs+=[f]
		thisAx.cla()	
		#
		# ultimately, averaging the actual ratios is not meaningful; we need to average the logs.
		# so, <log(x)> = (1/N)(x1+x2+...xn) = log(Prod(x_i)**(1/N))
		ploty=logaverageOver(list(map(operator.itemgetter(2), ratios)), avlen)
		#
		thisAx.fill_between(plaindts, hitThreshold*scipy.ones(len(fdts),int), ploty[0], color='b', where=scipy.array([val>=hitThreshold for val in ploty[0]]))
		thisAx.fill_between(plaindts, hitThreshold*scipy.ones(len(fdts),int), ploty[0], color='r', where=scipy.array([val<=hitThreshold for val in ploty[0]]))
		#
		maxy=math.log10(max(ploty[0]))
		#thisAx.axvline(x=eventFloatDt, color='c', lw=3, label='mainshock' )
		mevmag=mainEV[3]*maxy
		thisAx.plot([eventFloatDt, eventFloatDt], [-mevmag, mevmag], color='c', lw=3, label='mainshock' )
		#miny=min(ploty[0])
		#print "maxy, miny: %d, %d" % (maxy, miny)
		
		# note: we don't set the y-log scale in these "fill()" commands. we can do that with axis.set_yscale('log') i think.
		# we achieve this by doing semilogy() plots below.
		fg=plt.gcf()
		#plt.title("%s rupture area, time-time, wLen=%d" % (eventName, windowLen))
		#plt.xlabel('time')
		#plt.ylabel('$r=N_{rb-long} / N_{rb-short}$')
		#thisAx.axvline(x=eventFloatDt)
		nbigshocks=0
		bigShocks=self.getBigShocks(minmag, bigmag, cat0)
		for rw in bigShocks:
			if nbigshocks==0:
				#thisAx.axvline(x=rw[1], color='g', label='m > %f' % bigmag)
				thisAx.plot([rw[1], rw[1]], [1.0/(rw[4]*maxy), rw[4]*maxy], 'g-', color='g', label='m > %f' % bigmag) 
				nbigshocks+=1
			else:
				thisAx.plot([rw[1], rw[1]], [1.0/(rw[4]*maxy), rw[4]*maxy], 'g-') 
				thisAx.plot([rw[1]], rw[4]*maxy, '*')
			
		#thisAx.semilogy([eventFloatDt], [1], 'r^', ms=10)
		thisAx.semilogy([mainEV[0]], [1], 'r^', ms=10)
		
		thisAx.axhline(y=1, color='k')
		thisAx.legend(loc=legLoc, numpoints=2)
		#ax=plt.gca()
		#fg=plt.gcf()
	#	thisAx.xaxis.set_major_formatter(dates.DateFormatter('%Y-%b-%d'))
		thisAx.set_ylim([.1,10])
	#	fg.autofmt_xdate()
		#plt.savefig('images/%sRuptureTimeTime-Wlen%d-mc%d.png' % (eventName, wLen, int(10*minMag)))
		#plt.show()		
	
	def plotIntervalRatios(self, minmag=3.0, windowLen=10, cat0=None, hitThreshold=1.0, bigmag=5.0, fignum=0, ratios=None, deltaipos=1, avlen=1, logZ=1.0, rbLegLoc='best'):
		# avlen=10
		# avtype: 0: normalize, then take mean, 1: take mean, then normalize the mean values
		if logZ==None:
			logZ=math.log10(windowLen)
		logZinv=1.0/logZ
		#
		legLoc=rbLegLoc
		eventName="Event RB ratios"
		if cat0==None: cat0=self.shockCat
		# getIntervalRatios(self, minmag=3.0, windowLen=10, cat0=None, deltaipos=1):
		if ratios==None:
			ratios=self.getIntervalRatios(minmag, windowLen, cat0, deltaipos, logZ=logZ)
		#
		fdts=[]
		for rw in ratios:
			 fdts+=[mpd.date2num(rw[1])]
			 #fdts+=[rw[1].toordinal() + float(rw[1].hour)/24 + rw[1].minute/(24*60) + rw[1].second/(24*3600) + rw[1].microsecond/(24*3600000000)]
		#
		# for avtype==0 or avtype==1, we take the average (always take the average here)
		# if avtype==0 or (equivalently), we've provided normalized ratios, we'll just take the average and skip the next step.
		ploty=logaverageOver(list(map(operator.itemgetter(2), ratios)), avlen)
		#
		mainEv=self.getMainEvent(cat0)
		#
		eventFloatDt = mpd.date2num(mainEv[0])
		#
		f=plt.figure(fignum)
		plt.clf()
		#	
		plt.axvline(x=eventFloatDt, color='c', lw=3, label='mainshock' )
		plt.fill_between(fdts, hitThreshold*scipy.ones(len(fdts),int), ploty[0], color='b', where=scipy.array([val>=hitThreshold for val in ploty[0]]))
		plt.fill_between(fdts, hitThreshold*scipy.ones(len(fdts),int), ploty[0], color='r', where=scipy.array([val<=hitThreshold for val in ploty[0]]))
		#
		# note: we don't set the y-log scale in these "fill()" commands. we can do that with axis.set_yscale('log') i think.
		# we achieve this by doing semilogy() plots below.
		plt.title("%s rupture area, time-time, wLen=%d" % (eventName, windowLen))
		plt.xlabel('time')
		plt.ylabel('$r=N_{rb-long} / N_{rb-short}$')
		plt.axvline(x=eventFloatDt)
		nbigshocks=0
		bigShocks=self.getBigShocks(minmag, bigmag, cat0)
		for rw in bigShocks:
			if nbigshocks==0:
				plt.axvline(x=mpd.date2num(rw[1]), color='g', label='m > %f' % bigmag)
				nbigshocks+=1
			else:
				plt.axvline(x=mpd.date2num(rw[1]), color='g')
				plt.plot([mpd.date2num(rw[1])], rw[4], '*')
			
		plt.semilogy([eventFloatDt], [1], 'r^', ms=10)
		
		plt.axhline(y=1, color='k')
		plt.legend(loc=legLoc, numpoints=2)
		ax=plt.gca()
		fg=plt.gcf()
		ax.xaxis.set_major_formatter(dates.DateFormatter('%Y-%b-%d'))
		ax.set_ylim([.1,10])
		fg.autofmt_xdate()
		#plt.savefig('images/%sRuptureTimeTime-Wlen%d-mc%d.png' % (eventName, wLen, int(10*minMag)))
		plt.show()
	
	def getMagSubset(self, mag=5.0, cat=None):
		# get a subset of a catalog with m>mmag
		if cat==None: cat=self.shockCat
		outCat=[]
		for rw in cat:
			if rw[3]>=mag: outCat+=[rw]
		return outCat
	
	def getLargeAftershocks(self, mag=5.0, cat=None):
		if cat==None: cat=self.shockCat
		mainShock=self.getMainEvent()
		outcat=[]
		for rw in cat:
			if rw[0]>mainShock[0] and rw[3]>=mag: outcat+=[rw]
		mainShock=None
		
		return outcat
		#
	def getEarthquakeRatioScore(self, ratios=None, earthquakes=None):
		# the idea here is to get the "current" value of r(t) when an earthquake occurs. was the event predicted?
		# so, for each earthquake in earthquakes[], what was the most recent value of r, in ratios[].
		# if ratios and earthquakes are not provided, use getIntervalRatios(self), and m5.0 events from self.shockCat, respectively.
		#
		if ratios==None: ratios=self.getIntervalRatios()	#[n, date, r]
		if earthquakes==None: earthquakes=self.getMagSubset(5.0, self.shockCat)		#[date, lat, lon, mag, a, b]
		returnQuakes=[]
		#
		# now, assign a ratio value to each earthquake:
		#
		#for eqrow in earthquakes:
		for i in range(len(earthquakes)):
			#if eqrow[0]<ratios[0][1]: continue	# there is a winLen lag; we dont' have a forecast yet.
			# add the r value of the closest ratios entry to each earthquake:
			for irat in range(1, len(ratios)):
				if earthquakes[i][0]==ratios[irat][1]:
				#if ratios[irat][1]>=earthquakes[i][0]:
					returnQuakes+=[[earthquakes[i]+[ratios[irat-1][2]]]]
					continue
		#		
		return returnQuakes
		
	# science:
	def getRBintervals(self, minmag=1.0, useCat=None):
		# default record-breaking. walk forward in the shockCat; look for the largest/smallest intervals between events:
		cat=[]
		# print "shockcat[0]: %s, %s, %s, %s" % (self.shockCat[0][0], self.shockCat[1][0], self.shockCat[2][0], self.shockCat[3][0])
		if useCat==None: useCat=self.shockCat
		#for row in self.shockCat:
		for row in useCat:
			#print "rw3, minmag: %s, %s" % (row[3], minmag)
			if row[3]>=minmag: cat+=[row]
		#
		#print "catlen: %d" % len(cat)
		biggest=abs(datetimeToFloat(cat[1][0])-datetimeToFloat(cat[0][0]))	# interval between second and first events.
		#print "biggest interval: %f" % biggest
		smallest=biggest
		nbigger=1
		nsmaller=1
		
		# these arrays will be returned plot-ready:
		biggers=[[cat[1][0]], [biggest], [nbigger], [1]]		# [[date], [bigInts], [NRB_big], [i]]
		smallers=[[cat[1][0]],[smallest], [nsmaller], [1]]
		
		#prevDtmBig=cat[0][0]
		#prevDtmSmall=cat[0][0]
		for i in range(1,len(cat)):
			#dT=(cat[i][0]-cat[i-1][0]).seconds
			dT=abs(datetimeToFloat(cat[i][0])-datetimeToFloat(cat[i-1][0]))		# aka, the interval between the current and previous event...
			# note: i guess this WILL work backwards. intrinsically, it will work in reverse - intervals will be negative. we could
			# use this as is, but it probably makes sense to use the absolute value. after all, we want the magnitude of the interval:
			#
			#if cat[i][0]<cat[i-1][0]:
			#	print cat[i]
			#	a=input("type something")
			
			#print cat[i][0], cat[i-1][0], cat[i][3], dT
			if dT>biggest:
				nbigger+=1
				biggers[0]+=[cat[i][0]]		# date record was broken (datetime object)
				biggers[1]+=[dT]				# interval since last event (RBID)
				biggers[2]+=[nbigger]		# number of records broken (NRB)
				biggers[3]+=[i]				# record broken on i'th earthquake since mainshock (natural time).	[i+1] ??
				biggest=dT
				#prevDtmBig=cat[i][0]
				#
			if dT<smallest:
				nsmaller+=1
				smallers[0]+=[cat[i][0]]
				smallers[1]+=[dT]
				smallers[2]+=[nsmaller]
				smallers[3]+=[i]				# record broken on i'th earthquake since mainshock.
				smallest=dT
				#prevDtmSmall=cat[i][0]
				#
			#
		#
		#print biggers[0]
		#print biggers[1]
		
		cat=None
		return[biggers, smallers]
	
	def plotRBintervalSet(self, minmag=1.0, maxmag=5.0, dmag=.1, outdir='images/parkfield/'):
		plt.clf()
		plt.cla()
		
		curmag=minmag
		intervalSet=[]
		mags=[]
		if outdir[-1]!='/': outdir="%s/" % outdir
		
		print(("begin plotRBintervalSet() {%f}" % curmag))
		
		while curmag<=maxmag:
			curRecords=self.getRBintervals(curmag)	# curr records -> [[bigs], [smalls]]
			#bigsX=curRecords[0][0]
			#bigsY=curRecorss[0][1]
			# big records from currRecords
			intervalSet+=[[curRecords[0][0], curRecords[0][1], curRecords[0][2], curRecords[0][3], curmag]]		# "biggest" records, [date occured, interval, Nth record, mag-bin]
			mags+=[curmag]
			curmag+=dmag
		
		print("intervalSet(s) assigned...")
			
		#print len(intervalSet)
		#print len(intervalSet[0][0])
		
		startTime=datetimeToFloat(self.eventDtTime)
		
		# for giggles, get all the intervals:
		Xev=[0]
		Yev=[0]
		XevLog=[1]
		NevLog=[1]
		YevLog=[1]
		for iev in range(1,len(self.shockCat)):
			if (datetimeToFloat(self.shockCat[iev][0]) - datetimeToFloat(self.shockCat[iev-1][0])==0) : continue
			#
			Xev+=[datetimeToFloat(self.shockCat[iev][0]) - datetimeToFloat(self.shockCat[0][0])]
			Yev+=[datetimeToFloat(self.shockCat[iev][0]) - datetimeToFloat(self.shockCat[iev-1][0])]
			#
			#if abs(log10(datetimeToFloat(self.shockCat[iev][0]) - datetimeToFloat(self.shockCat[iev-1][0])))>3: continue
			XevLog+=[log10(datetimeToFloat(self.shockCat[iev][0]) - datetimeToFloat(self.shockCat[0][0]))]
			NevLog+=[log10(iev)+1]
			YevLog+=[log10(datetimeToFloat(self.shockCat[iev][0]) - datetimeToFloat(self.shockCat[iev-1][0]))]
		
		fignum=0
		plt.figure(fignum)
		plt.clf()
		fignum+=1
		plt.loglog(Xev, Yev, '.')
		plt.title("All Intervals")
		plt.xlabel("days since earthquake")
		plt.ylabel("interval")
		
		plt.figure(fignum)
		plt.clf()
		fignum+=1
		plt.loglog(Xev, Yev, '.')
		for pset in intervalSet:
			#plt.plot_date(pset[0], pset[1], '.-', label=str(minmag))
			x=[]
			#x2=[]
			for elem in pset[0]:
				x+=[datetimeToFloat(elem)-startTime]
				#x2+=[datetimeToFloat(elem)]
			plt.loglog(x, pset[1], '-', label=str(pset[3]))
			#plt.loglog(x2, pset[1], '-')
			
		#
		plt.title("RB interval")
		plt.xlabel("days since earthquake")
		plt.ylabel("interval (days)")
		plt.savefig("%sRBintervals.pdf" % outdir)
		#plt.legend(loc='upper left')
		
		#print "saveFig: %sRBintervals.pdf" % outdir
		
		plt.figure(fignum)
		plt.clf()
		fignum+=1
		for pset in intervalSet:
			#print len(pset)
			#plt.plot_date(pset[0], pset[2], '.-', label=str(minmag))
			x=[]
			for elem in pset[0]:
				x+=[datetimeToFloat(elem)-startTime]
			plt.loglog(x, pset[2], '.-', label=str(pset[4]))
		
		plt.title("Number of New RB Intervals")
		plt.xlabel("days since earthquake")
		plt.ylabel("Number of Broken Records")
		plt.savefig("%snewRBintervals.pdf" % outdir)
		#plt.legend(loc='upper left')
		
		# and let's do some data fitting:
		# for now, fit t>=10**-1.5 in the Nrecords plot, maybe interval>10**-2 for RBinterval??
		# alternatively, always trim the first (maybe first 2) data points.
		#
		# we want the log/log fit, so we can eithe rconstruct a log-based error or we can take the log-log
		# of the data and do a linear fit.
		#
		fitPlotsInt=[]
		fitPlotsN=[]
		slopes=[]
		#print YevLog
		#print "lens: %d, %d" % (len(XevLog), len(YevLog))
		fbooga=open("%sXevLog.dat" % outdir, 'w')
		for booga in range(len(XevLog)):
			fbooga.write("%f\t%f\t%f\t%f\n" %(XevLog[booga], YevLog[booga], Xev[booga], Yev[booga]))
			#fbooga.write("%f\t%f\n" %(XevLog[booga], YevLog[booga]))
		fbooga.close()
		
		plt.figure(fignum)
		plt.clf()
		fignum+=1
		plt.plot(XevLog, YevLog, '.')
		for pset in intervalSet:
			x=[]
			y=[]
			yfit=[]
			for i in range(len(pset[0])):
				#if i<1: continue		# (skip first element)
				#print "dtm arg: %s" % (datetimeToFloat(pset[0][i])-startTime)
				if log10(datetimeToFloat(pset[0][i])-startTime)<-1.0: continue
				x+=[log10(datetimeToFloat(pset[0][i])-startTime)]
				y+=[log10(pset[1][i])]
				
			#
			# now, we have a log-log set of one dataset. fit it to a line...
			#print "x: %s" % str(x)
			#print "y: %s" % str(y)
			
			p=scipy.array([0,1])
			plsq=spo.leastsq(linRes, p, args=(scipy.array(y), scipy.array(x)), full_output=0, maxfev=20000)	# (function to minimize, initial prams, argument-arrays, max-iterations)
			slopes+=[plsq[0][1]]
			for X in x:
				fitval=plsq[0][0] + X*plsq[0][1]
				yfit+=[fitval]
			plt.plot(x,yfit, '.-', label=str(pset[3]))
		plt.title("Fit to RB Intervals")
		plt.xlabel("log10(days since mainshock)")
		plt.ylabel("log10(Days Since Prev eq m>=m0)")
		plt.savefig("%sRBintervalFits.pdf" % outdir)
		#plt.legend(loc='upper left')
		#
		plt.figure(fignum)
		plt.clf()
		fignum+=1
		plt.title("Interval Slopes")
		plt.xlabel("mag")
		plt.ylabel("slope")
		plt.plot(mags, slopes, '.-')
		plt.savefig("%sRBintervalSlopes.pdf" % outdir)
		
		#############
		fitPlotsInt=[]
		fitPlotsN=[]
		slopes=[]
		plt.figure(fignum)
		plt.clf()
		fignum+=1
		for pset in intervalSet:
			x=[]
			y=[]
			yfit=[]
			for i in range(len(pset[0])):
				#if i<1: continue		# (skip first element)
				if log10(datetimeToFloat(pset[0][i])-startTime)<-1.0: continue
				x+=[log10(datetimeToFloat(pset[0][i])-startTime)]
				y+=[log10(pset[2][i])]
				
			#
			# now, we have a log-log set of one dataset. fit it to a line...
			#print "x: %s" % str(x)
			#print "y: %s" % str(y)
			p=scipy.array([0,1])
			plsq=spo.leastsq(linRes, p, args=(scipy.array(y), scipy.array(x)), full_output=0, maxfev=50000)	# (function to minimize, initial prams, argument-arrays, max-iterations)
			slopes+=[plsq[0][1]]
			for X in x:
				#print X
				fitval=plsq[0][0] + X*plsq[0][1]
				yfit+=[fitval]
			plt.plot(x,yfit, '.-', label=str(pset[4]))
		plt.title("fit to Number of New Records")
		plt.xlabel("log10(days since event)")
		plt.ylabel("log10(Nrecord breaking Intervals)")
		plt.savefig("%sRBfitNewRecords.pdf" % outdir)
		#plt.legend(loc='upper left')
		
		plt.figure(fignum)
		plt.clf()
		fignum+=1
		plt.title("Slopes Nrecords")
		plt.xlabel("mag")
		plt.ylabel("slope")
		plt.plot(mags, slopes, '.-')
		plt.savefig("%sRBslopesNewRecords.pdf" % outdir)
		
		###########################################################
		###########################################################
		print("natural time plots...")
		# natrual time plots:
		plt.figure(fignum)
		plt.clf()
		fignum+=1
		plt.loglog(list(range(1,len(Yev)+1)), Yev, '.')	# all the events...
		#plt.semilogy(range(1,len(Yev)+1), Yev, '.')
		for pset in intervalSet:
			#print len(pset), pset[3]
			#print len(pset[2])
			#plt.plot_date(pset[0], pset[1], '.-', label=str(minmag))
			#
			plt.loglog(pset[3], pset[1], '+-', label=str(pset[4]))
			
			#plt.semilogy(pset[3], pset[1], '-', label=str(pset[4]))
		#
		plt.title("RB interval (natural time)")
		plt.xlabel("n events since mainshock")
		plt.ylabel("interval (days)")
		plt.savefig("%sRBintervalsNT.pdf" % outdir)
		#plt.legend(loc='upper left')
		#
		# the RB line does not appear to line up with the data-data. what are the values of the last elements?
		#print "elements: %d, %f, %f, %f" % (len(Yev)+1, Yev[-1], pset[3][-1], pset[1][-1])
		# output to text-file:
		fout=open('%sNRBNTdata.dat' % outdir, 'w')
		fout.write("#nthEvent\tinterval\n")
		for ii in range(len(Yev)):
			fout.write("%d\t%f\n" % (ii, Yev[ii]))
		fout.close()
		fout=open('%sNRBNTrecs.dat' % outdir, 'w')
		fout.write("#nEvents\tRBinterval\n")
		for ii in range(len(pset[3])):
			fout.write("%d\t%f\n" % (pset[3][ii], pset[1][ii]))
		fout.close()
		
		
		plt.figure(fignum)
		plt.clf()
		fignum+=1
		for pset in intervalSet:
			#print len(pset)
			#plt.plot_date(pset[0], pset[2], '.-', label=str(minmag))
			#
			plt.loglog(pset[3], pset[2], '.-', label=str(pset[4]))
			#plt.semilogy(pset[3], pset[2], '.-', label=str(pset[4]))
		
		plt.title("Number of RB Events (NT)")
		plt.xlabel("n events since mainshock")
		plt.ylabel("Number of Record Breaking Events")
		plt.legend(loc='lower right')
		plt.savefig("%snewRBintervalsNT.pdf" % outdir)
		
		# and let's do some data fitting:
		# for now, fit t>=10**-1.5 in the Nrecords plot, maybe interval>10**-2 for RBinterval??
		# alternatively, always trim the first (maybe first 2) data points.
		#
		# we want the log/log fit, so we can either construct a log-based error or we can take the log-log
		# of the data and do a linear fit.
		#
		fitPlotsInt=[]
		fitPlotsN=[]
		slopes=[]
		plt.figure(fignum)
		plt.clf()
		fignum+=1
		#plt.plot(NevLog, YevLog, '.')
		for pset in intervalSet:
			x=[]
			y=[]
			yfit=[]
			#print len(pset[0])
			for i in range(len(pset[0])):
				#if i<1: continue		# (skip first element)
				if log10(datetimeToFloat(pset[0][i])-startTime)<-1.0: continue
				x+=[log10(pset[3][i])]
				y+=[log10(pset[1][i])]
				
			#
			# now, we have a log-log set of one dataset. fit it to a line...
			#print "x: %s" % str(x)
			#print "y: %s" % str(y)
			p=scipy.array([0,1])
			plsq=spo.leastsq(linRes, p, args=(scipy.array(y), scipy.array(x)), full_output=0, maxfev=20000)	# (function to minimize, initial prams, argument-arrays, max-iterations)
			slopes+=[plsq[0][1]]
			for X in x:
				fitval=plsq[0][0] + X*plsq[0][1]
				yfit+=[fitval]
			#plt.plot(x,yfit, '.-', label=str(pset[4]))
		plt.title("Fit to RB Intervals (NT)")
		plt.xlabel("log10(nevents since mainshock)")
		plt.ylabel("log10(Days Since Prev eq m>=m0)")
		plt.savefig("%sRBintervalFitsNT.pdf" % outdir)
		#plt.legend(loc='upper left')
		#
		plt.figure(fignum)
		plt.clf()
		fignum+=1
		plt.title("Interval Slopes (NT)")
		plt.xlabel("mag")
		plt.ylabel("slope")
		plt.plot(mags, slopes, '.-')
		plt.savefig("%sRBintervalSlopesNT.pdf" % outdir)
		
		#############
		fitPlotsInt=[]
		fitPlotsN=[]
		slopes=[]
		plt.figure(fignum)
		plt.clf()
		fignum+=1
		for pset in intervalSet:
			x=[]
			y=[]
			yfit=[]
			for i in range(len(pset[0])):
				#if i<1: continue		# (skip first element)
				if log10(datetimeToFloat(pset[0][i])-startTime)<-1.0: continue
				x+=[log10(pset[3][i])]
				y+=[log10(pset[2][i])]
			#
			# now, we have a log-log set of one dataset. fit it to a line...
			#print "x: %s" % str(x)
			#print "y: %s" % str(y)
			p=scipy.array([0,1])
			plsq=spo.leastsq(linRes, p, args=(scipy.array(y), scipy.array(x)), full_output=0, maxfev=50000)	# (function to minimize, initial prams, argument-arrays, max-iterations)
			slopes+=[plsq[0][1]]
			for X in x:
				#print X
				fitval=plsq[0][0] + X*plsq[0][1]
				yfit+=[fitval]
			plt.plot(x,yfit, '.-', label=str(pset[4]))
		plt.title("fit to Number of New Records (NT)")
		plt.xlabel("log10(nevents since event)")
		plt.ylabel("log10(Nrecord breaking Intervals)")
		plt.savefig("%sRBfitNewRecordsNT.pdf" % outdir)
		#plt.legend(loc='upper left')
		
		plt.figure(fignum)
		plt.clf()
		fignum+=1
		plt.title("Slopes Nrecords (NT)")
		plt.xlabel("mag")
		plt.ylabel("slope")
		plt.plot(mags, slopes, '.-')
		plt.savefig("%sRBslopesNewRecordsNT.pdf" % outdir)
		#
		#plt.show()
	#####################	
	
	# utilities:
	def scatterShockCat(self):
		# create a scatter plot of the shock-catalog.
		vecs=[[],[],[],[]]
		for row in self.shockCat:
			vecs[0]+=[row[1]]
			vecs[1]+=[row[2]]
			vecs[2]+=[row[4]]
			vecs[3]+=[row[5]]
		#
		#print vecs[1]
		plt.figure(0)
		plt.plot(vecs[1], vecs[0], '.')
		
		plt.figure(1)
		plt.plot(vecs[2], vecs[3], '.')
		
		plt.show()
		return None
#
# utils:	
def getLogs(data=[], lbase=10):
	#output=[]
	#for rw in data:
	#	newrow=[]
	#	for elem in rw:
	#		newrow+=[math.log(elem, lbase)]
	#	output+=[newrow]
	output=[]
	for elem in data:
		if elem==0:
			output+=[None]
			continue
		output+=[math.log(elem, lbase)]
	return output
	
def logaverageOver(inData=[], n=1):
	# return average over n elements.
	# return 1:1 rows, note that the first n elements will be averaged over n`<n
	#import numpy
	outData=[[],[]]	# <x>, stdev
	#
	inlogs=getLogs(inData)
	N=1
	currVal=0
	for i in range(len(inData)):
		#outData[0]+=[sum(inData[i-N:i+1])/float(N)]
		#
		#outData[0]+=[numpy.mean(inData[i+1-N:i+1])]
		#outData[1]+=[numpy.std(inData[i+1-N:i+1])]
		#
		#outData[0]+=[(numpy.prod(inData[i+1-N:i+1]))**(1./N)]		# this is correct, but to get stdev (which has exactly what meaning in this case?), we need to get the logs anyway.
		outData[0]+=[10**numpy.mean(inlogs[i+1-N:i+1])]
		outData[1]+=[10**numpy.std(inlogs[i+1-N:i+1])]	# we want to return this in de-loged format. we'll plot in log-space.
		
		if N<n: N+=1
	return outData

def averageOver(inData=[], n=1):
	# return average over n elements.
	# return 1:1 rows, note that the first n elements will be averaged over n`<n
	#import numpy
	outData=[[],[]]	# <x>, stdev
	#
	N=1
	currVal=0
	for i in range(len(inData)):
		#outData[0]+=[sum(inData[i-N:i+1])/float(N)]
		outData[0]+=[numpy.mean(inData[i+1-N:i+1])]
		outData[1]+=[numpy.std(inData[i+1-N:i+1])]
		
		if N<n: N+=1
	return outData
#
def deg2rad(theta):
	return 2*pi*theta/360.0
	
def ellipseY(x, a, b):
	#print b, x, a
	return b*(1-x*x/(a*a))**.5

def datetimeFromStrings(strDt, strTm, dtdelim='/'):
	# looks like date[time].strptime(date_string, format) does this...
	if strTm=='': strTm='00:00:00.0'
	#
	ldt=strDt.split(dtdelim)
	ltm=strTm.split(':')
	
	lsecs=ltm[2].split('.')
	secs = int(lsecs[0])
	msecs=0
	if len(lsecs)>1:
		msecs = int(float("."+lsecs[1])*10**6)
	#
	#return datetime.datetime(long(ldt[0]), long(ldt[1]), long(ldt[2]), long(ltm[0]), long(ltm[1]), long(ltm[2]))
	return datetime.datetime(int(ldt[0]), int(ldt[1]), int(ldt[2]), int(ltm[0]), int(ltm[1]), secs, msecs)

def datetimeToFloat(dtm):
	# return number of days, including fractional bit.
	return float(dtm.toordinal()) + float(dtm.hour)/24.0 + float(dtm.minute)/(24.0*60.0) + float(dtm.second)/(24.0*60.0*60.0) + float(dtm.microsecond)/(24.0*60*60*10**6)

def floatToDateTime(fdate):
	# this might be off by a few hundred microseconds. in one example, i get 890000-> 889993
	datepart=datetime.date.fromordinal(int(fdate))
	hrs=24.0*(fdate-int(fdate))
	hr=int(hrs)
	mins=(hrs-hr)*60
	mn=int(mins)
	secs=(mins-mn)*60
	sec=int(secs)
	microsecs=int((secs-sec)*10**6)
	
	timepart=datetime.time(hr,mn,sec, microsecs)
	
	return datetime.datetime.combine(datepart, datetime.time(hr, mn, sec, microsecs))


