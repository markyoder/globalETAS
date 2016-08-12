from math import *	
from scipy import *
import scipy
from pylab import *
from matplotlib import *
import matplotlib.dates as mpd

import linefit

#
# maping bits:
import matplotlib	# note that we've tome from ... import *. we should probably eventually get rid of that and use the matplotlib namespace.
#matplotlib.use('Agg')
#from matplotlib.toolkits.basemap import Basemap
from mpl_toolkits.basemap import Basemap as Basemap
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
#
import pylab as plt
from matplotlib import rc
import numpy
#
import string
import sys
import os
import random
import time
#
#from scipy.optimize import leastsq
from matplotlib import axis as aa
#from threading import Thread
#
import rbIntervals as rbi
#
import datetime as dtm
import pytz
import calendar
import operator
import urllib.request, urllib.parse, urllib.error
#
from geographiclib.geodesic import Geodesic as ggp

#
Rearth = 6378.1	# km
deg2rad=2.0*math.pi/360.
#
tzutc=pytz.timezone('UTC')
#
class eqcatalog:
	# a simple catalog class. we'll need X,Y,wt,
	# call like: lf1=m.catalog()
	#
	#
	mc=None	# catalog threshold
	cat=[]	# let's use [[evDateTime, lat, lon, mag, a, b], [...]] and we'll plot with map(), etc.
	subcats=[]
	activefig=None
	catmap=None
	#
	# "transition" magnitude between large/medium earthquakes
	mt=7.6
	#
	#
	def __init__(self, inData=[]):
		#self.cat=[]
		#self.subcats=[]
		self.initialize(inData)
	
	def initialize(self, inData=[]):
		# what does the data look like?
		self.cat=[]
		self.subcats=[]
		self.mt=7.6
		#
		self.cat=inData
		self.sqlport=3306
		self.sqlhost='localhost'
		inData=None
		self.catmap=None
		self.__name__='eqcatalog'
		self.mapres='l'	# basemap map resolution. at some pont, add functions to sort out nonsensical values.
		#
		# massage data as necessary. we've recently added depth, so to be compatible with old data,
		# run through the data and insert a 0 depth if len(cat[i])<5:
		for i in range(len(self.cat)):
			if len(self.cat[i])>=5:
				# if one row has depth, they probably all do.
				break
			# missing depth. it goes on the end -- which is a little confusing because ANSS is depth, mag, but it's our
			# convention now...
			cat[i]+=[0.0]	# and just put them all on the surface.
		#
		self.checkdates()
	#
	def writeCatToFile(self, foutname=None, cat=None):
		if foutname==None: foutname='outcat.cat'
		if cat==None: cat=self.getcat(0)
		#
		fout=open(foutname, 'w')
		fout.write("#dtm, lat, lon, mag\n")
		for rw in cat:
			fout.write("%d/%d/%d %d:%d:%d.%d\t%f\t%f\t%f" % (rw[0].year, rw[0].month, rw[0].day, rw[0].hour, rw[0].minute, rw[0].second, rw[0].microsecond, rw[1], rw[2], rw[3]))
			if len(rw)>=5: fout.write('\t%f' % rw[4])
			fout.write('\n')
		fout.close()
			
	def loadCatFromFile(self, fname=None, minmag=None, cols=[0,1,2,3,4,5]):
		# standard file format: date \t time \t lat \t lon \t mag \t depth(maybe)
		#
		#print "fname: %s" % fname
		self.cat=[]
		fin=open(fname)
		for rw in fin:
			if rw[0] in ('#', '\t', '\n', ' '): continue
			#rws=rw.split('\t')
			rws=rw.split()
			#rint rws
			#if len(rws)<4: continue
			#print rw
			#print rws
			if len(rws)<5: continue
			if minmag!=None:
				#if float(rws[3])<minmag: continue
				if float(rws[4])<minmag: continue
			#self.cat+=[[datetimeFromString(rws[0], float(rws[1]), float(rws[2]), float(rws[3])]]
			#if '24:' in rws[0]: print rws[0]
			#print rws[1]
			#print rws
			#self.cat+=[[datetimeFromString(rws[0] + ' '+ rws[1]), float(rws[2]), float(rws[3]), float(rws[4])]]
			dpth=0.0
			if len(rws)>=6: dpth = float(rws[cols[5]])
			thisrw = [datetimeFromString(rws[cols[0]] + ' '+ rws[cols[1]]), float(rws[cols[2]]), float(rws[cols[3]]), float(rws[cols[4]]), dpth]
			self.cat+=[thisrw]
			#self.cat+=[[datetimeFromString(rws[cols[0]] + ' '+ rws[cols[1]]), float(rws[cols[2]]), float(rws[cols[3]]), float(rws[cols[4]])]]
		fin.close()
		#
		# add (UTC) timezone:
		self.checkdates()	# note, in the future we should check for TZ info from string. for now, assume utc.
		# and now, let's assume that we want this time-ordered:
		self.cat.sort(key = lambda x: x[0])
	#
	def checkdates(self, cat=None):
		if cat==None: cat=self.getcat(0)
		#for i in xrange(len(cat)):
		for i,rw in enumerate(cat):
			#
			if isinstance(rw[0], numpy.datetime64) or hasattr(rw[0], 'tzinfo')==False: 
				print('rw[0]', rw[0])
				cat[i][0] = mpd.num2date(mpd.datestr2num(str(rw[0])))
				#cat[i][0] = cat[i][0].tolist()
				print(cat[i][0], type(cat[i][0]))
			if cat[i][0].tzinfo==None:
				# no time-zone info. add UTC timezone.
				dt0=cat[i][0]
				thisdt=dtm.datetime(*dt0.timetuple()[:-2], tzinfo=pytz.timezone('UTC'))
				cat[i][0]=thisdt
			#
		#
		return None
	#
	def getMainEvent(self, thiscat=None):
		# return catalog row of max magnitude (epicenter location (more or less)) event. note, by default we use ths shock-cat because it will
		# be faster than using the fullCat AND, the fullCat is likely to have other large earthquakes.
		if thiscat==None: thiscat=self.cat
		maxMag=thiscat[0][3]
		maxIndex=0
		for i in range(len(thiscat)):
			#print i, maxMag, maxIndex, cat[i][3]
			if thiscat[i][3]>maxMag:
				maxIndex=i
				maxMag=thiscat[i][3]
		return thiscat[maxIndex] + [maxIndex]
	
	def getIndexDtm(self, mindt=None, cat=None, datecol=0):
		if cat==None or type(cat).__name__ not in ('list', 'tuple'): cat=self.cat
		if mindt==None or type(mindt).__name__!='datetime': mindt=self.getMainEvent()[0]
		#
		for rw in cat:
			if rw[datecol]>=mindt: return rw
		return None
		#
		return thiscat[maxIndex] + [maxIndex]
	
	def getSubCat(self, catindex=0):
		if len(self.subcats)==0: return None
		#
		return self.subcats[catindex][1]
	
	def getcat(self, catindex=0):
		# more general and simpler than getSubCat. 0 -> maincat, 1, etc. -> subcats:
		# it would probably be a good idea to also restructure how catalogs are stored inte class:
		# catalogs=[[maincat], [subcat1], [subcat2], ...]
		# self.cat -> catalogs[0]
		if catindex==0: return self.cat
		if catindex==-1: catindex=len(self.subcats)
		if len(self.subcats)<(catindex): return None	# this cat index does not exist
		return self.subcats[catindex-1][1]
	
	def getLatLonRange(self, cat=None, latloncols=[1,2]):
		if cat==None: cat=self.cat
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
	#
	def rotatexy(self, x, y, Lat, Lon, theta, degtype='deg'):
		# x,y to transform via blah, blah.
		#
		if degtype=='deg': theta=deg2rad * float(theta)
		xprime = (x-Lon)*cos(theta) - (y-Lat)*sin(theta)
		yprime = (x-Lon)*sin(theta) + (y-Lat)*cos(theta)
	
		return [xprime, yprime]
	# subcat (subcatalog) functions:
	def ellipseCat(self, fullcat=None, theta=0, clat=35.9, clon=-120.5, ra=1.0, rb=1.0):
		#
		
		#print "event (start) date, catname: %s, %s, %s" % (eventDate, catFname, self.catname)
		#
		if fullcat==None: fullcat=self.cat
		#self.subcats+=[[subcatname], []]
		tempcat=[]
		
		#nEventsSinceMS=0
		for row in fullcat:
			# rotate each element into our aftershock axis, is it in the ellipse?
			newVec=self.rotatexy(row[2], row[1], clat, clon, theta)
			#
			# is the rotated vector in our ellipse?
			if abs(newVec[0])>ra: continue
			Y=self.ellipseY(newVec[0], ra, rb)
			if abs(newVec[1])>Y: continue
			# dtm, lat, lon, mag, tX, tY 		(note this is like y,x, x`, y` for the space coordinates).
			#self.subcats[-1][1]+=[[row[0], row[1], row[2], row[3], newVec[0], newVec[1]]]
			#tempcat+=[[row[0], row[1], row[2], row[3], newVec[0], newVec[1]]]
			tempcat+=[row + [newVec[0], newVec[1]] ]
		return tempcat
	def ellipseY(self, x, a, b):
		#print b, x, a
		return b*(1.0-x*x/(a*a))**.5
		
	def polycat(self, cat=None, verts=None):
		# as per james' counsel, a "knot theory" approach is much simpler. independent of right/left handedness, the sum of above/below tests is >0 for
		# points inside, 0 for points outside (like EnM).
		# start by making verts -> vectors -> f(x)
		#
		# verts are like: [[x0,y0], [x1, y1], ..., [xn, yn], [x0, y0]]; last point is optional.
		if cat==None: cat=self.cat
		if verts==None or len(verts)<3:
			# don't know. if we don't have verts, what can we do?
			# also, we need at least 3 verts, or we just have a line.
			return None
		#
		if verts[-1]!=verts[0]: verts+=[verts[0]]
		#
		vecs=[]	# like [ [[x0,y0], [x1, y1]], [[x1, y1], [x2, y2]], ...]
		#vecdirs=[]	# vector directions; -1=left, 0=none, 1=right. this determines whether we want to be over or under. the first x-direction vector definds a right/left poly.
		# get lat, lon extrema and make vectors:
		extremeVerts=[verts[0][0], verts[0][0], verts[0][1], verts[0][1]]	# [minLon, maxLon, minLat, maxLat]
		for i in range(len(verts)-1):
			vecs+=[[verts[i], verts[i+1]]]
			if verts[i+1][0]>extremeVerts[1]: extremeVerts[1]=verts[i+1][0]
			if verts[i+1][0]<extremeVerts[0]: extremeVerts[0]=verts[i+1][0]
			if verts[i+1][1]>extremeVerts[3]: extremeVerts[3]=verts[i+1][1]
			if verts[i+1][1]<extremeVerts[2]: extremeVerts[2]=verts[i+1][1]
			#
			# and keep a list of vector directions (right,left; do we need up, down?):
			thisdir=0	# reserve for vertical elements.
			if verts[i+1][0]>verts[i][0]: thisdir=1 #CatMap
			if verts[i+1][0]<verts[i][0]: thisdir=-1
			#vecdirs+=[thisdir]
		#
		# we don't really need the center, but it might be useful later:
		center=scipy.array([extremeVerts[0] + (extremeVerts[1]-extremeVerts[0])/2.0, extremeVerts[2] + (extremeVerts[3]-extremeVerts[2])/2.0])
		#
		# and this way, we don't need the poly-direction.
		# now we can spin through the catalog. inout=sum(x^ * above/below). inout=0 means out; inout>0 means in.
		# where x^ is {-1, 1} for left, right; above/below is {-1, 1} for point is above, below. i don't think which one is -1 and which is 1 matters
		# so long as we are consistent. also, as per old-school gaussian integrals, the number of times we cross a boundary: odd-> in , even -> out
		# applies as well.
		polycat=[]
		for iev in range(len(cat)):
			event=cat[iev]
			x=event[2]
			y=event[1]
			# for speed, if we are outside the extreme vertices, move on:
			if (x<extremeVerts[0] or x>extremeVerts[1] or y<extremeVerts[2] or y>extremeVerts[3]):
					#print "extreme kill (%d, %d)" % (x, y)
					#keepEvent=0
					# and we're done...
					continue
			#
			#keepEvent=1	# start by assuming we keep the event.
			inPolyTracker=0	# running up/down score. by default, do not keep the event.
			#print "*#*#*#"
			for ivec in range(len(vecs)):
				vec=vecs[ivec]
				# make a line (if it's not vertical):
				if vec[1][0]-vec[0][0]==0: continue	# vertical segments do not contribute, and we'll get x/0 error.
				b=(vec[1][1]-vec[0][1])/(vec[1][0]-vec[0][0])
				a=vec[0][1]-b*vec[0][0]
				y0=a+b*x
				#
				# xrange:
				if vec[0][0]>vec[1][0]:
					bigX=vec[0][0]
					smallX=vec[1][0]
				if vec[0][0]<vec[1][0]:
					bigX=vec[1][0]
					smallX=vec[0][0]
				#
				# debug:
				#if iev<40:
				#	print vec[0][0], vec[1][0], x, lookUpDown, y, y0
				#	print (x>=smallX and x<=bigX), (lookUpDown==-1 and y>y0 ), (lookUpDown==1 and y<y0)
				#
				# are we in the current xrange?
				if (x<smallX or x>bigX): continue
				# if it's on the line, keep it:
				if y==y0:
					inPolyTracker=1
					continue
				# is it inside the polygon?
				if y>y0: isUp=1							# point is above
				if y<y0: isUp=-1							# point is below
				if vec[1][0]>vec[0][0]: vecDir=1		# to the right
				if vec[1][0]<vec[0][0]: vecDir=-1	# to the left
				inPolyTracker+=(vecDir*isUp)
				#
			#
			if inPolyTracker!=0: polycat+=[event]
			
		#print extremeVerts
		return polycat
		
		
	def polycat_cp(self, cat=None, verts=None):
		# my original version of polycat using cross products to determine the direction of the polygon. there is a faster way...
		# verts are like: [[x0,y0], [x1, y1], ..., [xn, yn], [x0, y0]]; last point is optional.
		if cat==None: cat=self.cat
		if verts==None or len(verts)<3:
			# don't know. if we don't have verts, what can we do?
			# also, we need at least 3 verts, or we just have a line.
			return None
		#
		if verts[-1]!=verts[0]: verts+=[verts[0]]
		#
		vecs=[]	# like [ [[x0,y0], [x1, y1]], [[x1, y1], [x2, y2]], ...]
		vecdirs=[]	# vector directions; -1=left, 0=none, 1=right. this determines whether we want to be over or under. the first x-direction vector definds a right/left poly.
		# get lat, lon extrema and make vectors:
		extremeVerts=[verts[0][0], verts[0][0], verts[0][1], verts[0][1]]	# [minLon, maxLon, minLat, maxLat]
		for i in range(len(verts)-1):
			vecs+=[[verts[i], verts[i+1]]]
			if verts[i+1][0]>extremeVerts[1]: extremeVerts[1]=verts[i+1][0]
			if verts[i+1][0]<extremeVerts[0]: extremeVerts[0]=verts[i+1][0]
			if verts[i+1][1]>extremeVerts[3]: extremeVerts[3]=verts[i+1][1]
			if verts[i+1][1]<extremeVerts[2]: extremeVerts[2]=verts[i+1][1]
			#
			# and keep a list of vector directions (right,left; do we need up, down?):
			thisdir=0	# reserve for vertical elements.
			if verts[i+1][0]>verts[i][0]: thisdir=1
			if verts[i+1][0]<verts[i][0]: thisdir=-1
			vecdirs+=[thisdir]
		#
		#print vecdirs
		# now, is the poly right or left handed? from the poly center, calculate the mean r x v.
		center=scipy.array([extremeVerts[0] + (extremeVerts[1]-extremeVerts[0])/2.0, extremeVerts[2] + (extremeVerts[3]-extremeVerts[2])/2.0])
		#print "verts: %s" % str(verts)
		#print "vecs: %s" % str(vecs)
		#print "center: %s" % str(center)
		polyDir=0
		for vec in vecs:
			# get the cross product r x vec:
			thisvec=scipy.array(vec[1])-scipy.array(vec[0])
			rvec=scipy.array(vec[0])-center
			cprod=numpy.cross(rvec, thisvec)
			#print vec, thisvec, rvec, cprod, type(cprod)
			# so i guess when your vectors are coplanar, scipy.array knows to retun just a scalar for the cross-product.
			polyDir+=cprod
		if polyDir>0: polyDir=1
		if polyDir<0: polyDir=-1
		#print "polyDir: %f" % polyDir
		#
		# now we can spin through the catalog to find elements above/below poly segments, depending on the direction of the segment and right/left
		# handedness of the poly.
		polycat=[]
		for iev in range(len(cat)):
			event=cat[iev]
			x=event[2]
			y=event[1]
			keepEvent=1	# start by assuming we keep the event.
			#print "*#*#*#"
			for ivec in range(len(vecs)):
				# test the event against each polygon segment. if it falls outside one or more, don't keep it...
				vec=vecs[ivec]
				# make a line:
				if vec[1][0]-vec[0][0]==0: continue
				#
				b=(vec[1][1]-vec[0][1])/(vec[1][0]-vec[0][0])
				a=vec[0][1]-b*vec[0][0]
				#
				lookUpDown=vecdirs[ivec]*polyDir
				y0=a+b*x
				#keep criteria:
				#if (x>=vec[0][0] and vec<=vec[1][0]) and ((lookUpdown==-1 and y<=y0 ) or (lookUpdown==1 and y>=0)):
				# so discard criteria is opposite (in y):
				if vec[0][0]>vec[1][0]:
					bigX=vec[0][0]
					smallX=vec[1][0]
				if vec[0][0]<vec[1][0]:
					bigX=vec[1][0]
					smallX=vec[0][0]
				#	
				if (x<extremeVerts[0] or x>extremeVerts[1] or y<extremeVerts[2] or y>extremeVerts[3]):
					print("extreme kill (%d, %d)" % (x, y))
					keepEvent=0
				if ((x>=smallX and x<=bigX) and ((lookUpDown==-1 and y>y0 ) or (lookUpDown==1 and y<y0))) :
					keepEvent=0
					print("f(x) kill (%d, %d)" % (x, y))
					# and for efficiency:
					continue
				#
			#
			if keepEvent==1: polycat+=[event]
		
			
		print(extremeVerts)
		return polycat
			
		
#[0,0], [2,0], [4,4], [2,6], [0,4]			
	
	def addEllipCat(self, subcatname='newcat', fullcat=None, theta=0, clat=35.9, clon=-120.5, ra=1.0, rb=1.0):
		#
		if fullcat==None: fullcat=self.cat
		
		newcat=self.ellipseCat(fullcat, theta, clat, clon, ra, rb)
		self.subcats+=[[subcatname, newcat]]
		
	def getMagSubcat(self, fullcat=None, minmag=2.5, magcol=3):
		newcat=[]
		for rw in fullcat:
			if rw[magcol]>=minmag: newcat+=[rw]
		return newcat
	def addMagSubcat(self, subcatname='magsubcat', fullcat=None, minmag=2.5, magcol=3):
		subcatname='%s-%s' % (subcatname, str(minmag))
		self.subcats+=[[subcatname, getMagSubcat(fullcat, minmag, magcol)]]
	#
	def getxytSubcat(self, fullcat=None, dts=[], lats=[], lons=[], llcols=[1,2]):
		if type(dts).__name__!='list': dts=[]
		if type(lats).__name__!='list': lats=[]
		if type(lons).__name__!='list': lons=[]
		while len(dts)<2: dts+=[None]
		#
		newcat=self.getTimeRangeCat(fullcat, dts[0], dts[1])
		newcat=self.getLatLonSubcat(newcat, lats, lons, llcols)
		#
		return newcat
	def addxytSubcat(self, subcatname='xytsubcat', fullcat=None, dts=[], lats=[], lons=[], llcols=[1,2]):
		self.subcats+=[[subcatname, self.getxytSubcat(fullcat, dts, lats, lons, llcols)]]
	
	def getTimeRangeCat(self, fullcat=None, dtFrom=None, dtTo=None):
		if fullcat==None: fullcat=self.cat
		if dtFrom==None: dtFrom=fullcat[0][0]
		if dtTo==None: dtTo=fullcat[-1][0]
		newcat=[]
		for rw in fullcat:
			if mpd.date2num(rw[0])>=mpd.date2num(dtFrom) and mpd.date2num(rw[0])<=mpd.date2num(dtTo): newcat+=[rw]
		#
		return newcat
	def addTimeRangeCat(self, subcatname='dtSubcat', fullcat=None, dtFrom=None, dtTo=None):
		self.subcats+=[[subcatname, self.getTimeRangeCat(fullcat, dtFrom, dtTo)]]
	
	def getLatLonSubcat(self, fullcat=None, lats=[], lons=[], llcols=[1,2]):
		# llcols: lat, lon
		llrange=None
		if fullcat==None: fullcat=self.getcat(0)
		#
		if lats in [[], None]:
			if llrange==None: llrange=self.getLatLonRange(fullcat, latloncols=[llcols[0], llcols[1]])
			deltaLats=llrange[1][0]-float(llrange[0][0])
			lats=[llrange[0][0]+deltaLats/2.0, llrange[1][0]-deltaLats/2.0]
		
		if lons in [[], None]:
			if llrange==None: llrange=self.getLatLonRange(fullcat, latloncols=[llcols[0], llcols[1]])
			deltaLons=llrange[1][1]-float(llrange[0][1])
			lons=[llrange[0][1]+deltaLons/2.0, llrange[1][1]-deltaLons/2.0]
			
			#lats={get min, max lat,lon from catalog}
		# and same for lons...
		#
		newcat=[]
		for rw in fullcat:
			if (rw[llcols[0]]>=lats[0] and rw[llcols[0]]<=lats[1]) and (rw[llcols[1]]>=lons[0] and rw[llcols[1]]<=lons[1]):
				newcat+=[rw]
			#
		#
		return newcat
	def addLatLonSubcat(self, subcatname='xysubcat', fullcat=None, lats=[], lons=[], llcols=[1,2]):
		self.subcats+=[[subcatname, self.getLatLonSubcat(fullcat, lats, lons, llcols)]] 
	
	def getxytmSubcat(self, fullcat=None, dts=[], lats=[], lons=[], minmag=2.5, llmcols=[1,2,3]):
		# just do the whole thing here, so it's fast:
		if type(dts).__name__!='list': dts=[]
		if type(lats).__name__!='list': lats=[]
		if type(lons).__name__!='list': lons=[]
		if type(llmcols).__name__!='list': llmcols=[1,2,3]
		while len(dts)<2: dts+=[None]
		while len(llmcols)<3: llmcols+=[llmcols[-1]+1]
		llrange=None
		#return llmcols
		if lats in [[], None] or len(lats)!=2:
			if llrange==None: llrange=self.getLatLonRange(fullcat, latloncols=[llmcols[0], llmcols[1]])
			deltaLats=llrange[1][0]-float(llrange[0][0])
			lats=[llrange[0][0]+deltaLats/2.0, llrange[1][0]-deltaLats/2.0]
		
		if lons in [[], None] or len(lons)!=2:
			if llrange==None: llrange=self.getLatLonRange(fullcat, latloncols=[llmcols[0], llmcols[1]])
			deltaLons=llrange[1][1]-float(llrange[0][1])
			lons=[llrange[0][1]+deltaLons/2.0, llrange[1][1]-deltaLons/2.0]
		#
		#newcat=self.getTimeRangeCat(fullcat, dts[0], dts[1])
		#newcat=self.getLatLonSubcat(newcat, lats, lons, llcols)
		newcat=[]
		print(lats, lons, dts)
		for rw in fullcat:
			if rw[llmcols[0]]>=lats[0] and rw[llmcols[0]]<=lats[1] and rw[llmcols[1]]>=lons[0] and rw[llmcols[1]]<=lons[1] and rw[llmcols[2]]>=minmag and rw[0]>=dts[0] and rw[0]<=dts[1]:
				newcat+=[rw]
		return newcat
		#
	def addxytmSubcat(self, subcatname='xytmsubcat', fullcat=None, dts=[], lats=[], lons=[], minmag=2.5, llmcols=[1,2,3]):
		#print llmcols
		self.subcats+=[[subcatname, self.getxytmSubcat(fullcat, dts, lats, lons, minmag, llmcols)]]
			
	def mapOverlay(self, catalog=None, fignum=0, dots='b.', doShow=False):
		# this does not quite work yet. the map does not rescale properly for the distinct catalogs with different lat/lon ranges.
		# it looks like a good approach might be to create a map-class, which can contain a catalog or vice-versa, or maybe one
		# could be a sub-class, but i dont' think that hierarchy is clear.
		# the basic idea: map: lat/lon range, lat/lon center, projection, etc., catalogOverlays [] (are a list of catalogs overlayed on the map. note
		# the lat/lon range will be (at least) max/min(lon/lat from any cat) overlayed onto the map). also, annotationOverlays (text, graphics, etc.),
		# other stuff too...
		if catalog==None: catalog=self.cat
		f0=plt.figure(fignum)
		#
		#set up map:
		llr=self.getLatLonRange(catalog)	# latLonRange
		llr[0][0]-=2.0
		llr[0][1]-=2.0
		llr[1][0]+=2.0
		llr[1][1]+=2.0
		
		cntr=[float(llr[0][0])+(llr[1][0]-float(llr[0][0]))/2.0, float(llr[0][1])+(llr[1][1]-float(llr[0][1]))/2.0]
		catmap=Basemap(llcrnrlon=llr[0][1], llcrnrlat=llr[0][0], urcrnrlon=llr[1][1], urcrnrlat=llr[1][0], resolution=self.mapres, projection='tmerc', lon_0=cntr[1], lat_0=cntr[0])
		canvas=FigureCanvas(f0)
		catmap.ax=f0.add_axes([0,0,1,1])
		f0.set_figsize_inches((8/catmap.aspect,8.))
		#
		catmap.drawcoastlines(color='gray')
		catmap.drawcountries(color='gray')
		catmap.fillcontinents(color='beige')
		xfull, yfull=catmap(list(map(operator.itemgetter(2), catalog)), list(map(operator.itemgetter(1), catalog)))
		#epx, epy=catmap(epicenter[0], epicenter[1])
		catmap.plot(xfull, yfull, dots, label='Full Catalog')
		#catmap.plot(epx, epy, 'ro')
		#canvas.print_figure(saveName)
		
		if doShow: plt.show()
		
		return None
		
	def plotCatsMap(self, catalogses=None, maincat=0, doShow=True, doSave=False, saveName='catalogPlot.png', epicenter=None, legendLoc='best', maincatname='full cat', fignum=0, doCLF=True, bigmag=6.0, padfactor=.25):
		# somehow, this is returning a skewed map - i think basically, the basemap object re-callibrates itself to the smaller catalog, so x,y=thisthing.basemapobject(lat, lon) returns something off by a bit.
		
		# same as plotCatMap, but multiple catalogs. we assume the lat/lon range comes from the first catalog.
		# maincat is the "main catalog", the subcat we care about most. we assume the primary catalog is the broadest; maincat contains the epicenter, etc.
		if catalogses==None: catalogses=[maincatname, self.cat] + self.subcats

		#catalogs=[self.cat] + map(operator.itemgetter(1), self.subcats)
		#catnames=[maincatname] + map(operator.itemgetter(0), self.subcats)
		catalogs=list(map(operator.itemgetter(1), catalogses))
		catnames=list(map(operator.itemgetter(0), catalogses))
		#return [catalogs, catnames]
		catalog=catalogs[0]
		
		if epicenter==None:
			#mainshock=self.getMainEvent(catalog)
			mainshock=self.getMainEvent(catalogs[maincat])
			epicenter=[mainshock[2], mainshock[1]]
		#
		f0=plt.figure(fignum)	
		if doCLF: plt.clf()
		#		
		#set up map:
		llr=self.getLatLonRange(catalog)	# latLonRange #return [[minLat, minLon], [maxLat, maxLon]]
		latpad=padfactor*(llr[1][0]-llr[0][0])
		lonpad=padfactor*(llr[1][1]-llr[0][1])
		llr[0][0]-= latpad	#.5
		#if llr[0][0]<90.: llr[0][0]=90.
		llr[0][1]-= lonpad	#.5
		#if llr[0][1]<-180.: llr[0][1]=-180.
		llr[1][0]+= latpad	#.5
		#if llr[1][0]>90.: llr[1][0]=90.
		llr[1][1]+= latpad	#.5
		#if llr[1][1]>180.: llr[1][1]=180.
		
		print("setting up map prams")
		
		cntr=[float(llr[0][0])+(llr[1][0]-float(llr[0][0]))/2.0, float(llr[0][1])+(llr[1][1]-float(llr[0][1]))/2.0]
		print("create basmap object.")
		catmap=Basemap(llcrnrlon=llr[0][1], llcrnrlat=llr[0][0], urcrnrlon=llr[1][1], urcrnrlat=llr[1][0], resolution =self.mapres, projection='tmerc', lon_0=cntr[1], lat_0=cntr[0])
		print("bm object created...")
		canvas=FigureCanvas(f0)
		catmap.ax=f0.add_axes([0,0,1,1])
		#f0.set_figsize_inches((8/catmap.aspect,8.))
		
		#f0.set_figsize_inches((10/catmap.aspect,10.))
		#f0.set_size_inches((10/catmap.aspect,10.))
		f0.set_size_inches((10.,15.))
		#
		print("draw stuff on map...")
		catmap.drawcoastlines(color='gray', zorder=0)
		catmap.drawcountries(color='gray', zorder=0)
		catmap.fillcontinents(color='beige', zorder=0)
		#catmap.drawrivers(color='b')
		catmap.drawstates()
		catmap.drawmeridians(list(range(int(llr[0][1]-2.0), int(llr[1][1]+2.0))), color='k', labels=[1,1,1,1])
		catmap.drawparallels(list(range(int(llr[0][0]-2.0), int(llr[1][0]+2.0))), color='k', labels=[1, 1, 1, 1])
		#
		'''
		catmap.llcrnrlon=llr[0][1]+2.0
		catmap.llcrnrlat=llr[0][0]+2.0
		catmap.urcrnrlon=llr[1][1]-2.0
		catmap.urcrnrlat=llr[1][0]-2.0
		'''
		print("plot catalogs...")
		icat=0
		for ct in catalogs:
			xfull, yfull=catmap(list(map(operator.itemgetter(2), ct)), list(map(operator.itemgetter(1), ct)))
			catmap.plot(xfull, yfull, '.', label='%s' % catnames[icat], ms=2, zorder=1, alpha=.5)
			icat+=1
		
		# now, plot all events m>m0 from the full catalog:
		#bigmag=5.0
		for rw in catalog:
			if rw[3]<bigmag: continue
			thisx, thisy=catmap(rw[2], rw[1])
			catmap.plot(thisx, thisy, '*', label='%s, %s\n (%s, %s)' % (str(rw[3]), str(rw[0]), str(rw[2]), str(rw[1])), ms=15, zorder=2)
			
		#epx, epy=catmap(epicenter[0], epicenter[1])
		#catmap.plot(epx, epy, 'ro', label='epicenter', zorder=1)
		#
		#################
		#################
		# this is how to draw an ellipse... obviously, this does not really belong in this part of the script;
		#  it was part of the learning process...
		###############
		#canvas.print_figure(saveName)
		#from matplotlib.patches import Ellipse
		#f=plt.figure(0)
		#
		#ax1=f0.gca()
		#el = Ellipse([-120.5, 35.9], .8, .3, -40, facecolor='b', alpha=0.4)
		#Xel, Yel = catmap(el.get_verts()[:,0],el.get_verts()[:,1])
		#catmap.plot(Xel, Yel, '-r', lw=2)
		#catmap.ax.fill(Xel, Yel, ec='r', fc='r', alpha=.4)
		###
		##################
		#ax1.add_artist(el)
		#catmap.ax.add_artist(el)
		#
		#ax=plt.gca()
		#el = Ellipse((self.tLon, self.tLat), 2.0*self.tA, 2.0*self.tB, -self.tTheta, facecolor='b', alpha=0.4)
		#catmap.ax.add_artist(el)
		#ax.add_artist(el)
		#
		#plt.plot(map(operator.itemgetter(2), self.fullCat), map(operator.itemgetter(1), self.fullCat), '+')
		#plt.plot(map(operator.itemgetter(2), self.shockCat), map(operator.itemgetter(1), self.shockCat), '.')
		#plt.plot(map(operator.itemgetter(2), fcat), map(operator.itemgetter(1), fcat), '+', label='Full Catalog')
		#plt.plot(map(operator.itemgetter(2), scat), map(operator.itemgetter(1), scat), '.', label='Aftershock zone')
		#plt.plot([epicenter[0]], [epicenter[1]], 'ro', label='epicenter')
		plt.legend(loc=legendLoc, numpoints=1)
		if doSave: plt.savefig('pltsave-%s' % saveName)
		
		canvas.print_figure(saveName)
		
		if doShow: plt.show()
		
		return catmap
	
	def plotCatMap(self, catalog=None, doShow=True, doSave=False, saveName='catalogPlot.png', epicenter=None, legendLoc='upper left', doCLF=True, eqicon='b,', myaxis=None, fignum=None, padfactor=.25, plotevents=True, meridian_labels=[0,0,0,1], parallel_labels=[0,1,0,0], mapLLlat=None, mapLLlon=None, mapURlat=None, mapURlon=None, deltaLat=1.0, deltaLon=1.0):
		# labels : [left, right, top, bottom]
		# temporary:
		padfactor=0.
		if catalog==None: catalog=self.cat
		
		if epicenter==None:
			mainshock=self.getMainEvent(catalog)
			epicenter=[mainshock[2], mainshock[1]]
		#
		if doShow>=1 and fignum==None: fnum=doShow
		if fignum!=None: fnum=fignum
		
		f0=plt.figure(int(doShow))	
		if doCLF: plt.clf()
		#		
		#set up map:
		llr=self.getLatLonRange(self.getcat(0))	# latLonRange
		#
		# permit custom lat/lon range.
		if mapLLlat!=None: 
			llr[0][0]=mapLLlat
			print("setting LL lat ", mapLLlat)
		if mapLLlon!=None: 
			llr[0][1]=mapLLlon
			print("setting LL lon ", mapLLlon)
		if mapURlat!=None: 
			llr[1][0]=mapURlat
			print("setting UR lat: ", mapURlat)
		if mapURlon!=None: 
			llr[1][1]=mapURlon
			print("setting UR lon: ", mapURlon)
		#
		#llr=self.getLatLonRange(catalog)	# latLonRange #return [[minLat, minLon], [maxLat, maxLon]]
		#llr = [[self.catmap.llcrnrlat, self.catmap.llcrnrlon],[self.catmap.urcrnrlat, self.catmap.urcrnrlon]]
		latpad=padfactor*(llr[1][0]-llr[0][0])
		lonpad=padfactor*(llr[1][1]-llr[0][1])
		llr[0][0]-= latpad	#.5
		#if llr[0][0]<90.: llr[0][0]=90.
		llr[0][1]-= lonpad	#.5
		#if llr[0][1]<-180.: llr[0][1]=-180.
		llr[1][0]+= latpad	#.5
		#if llr[1][0]>90.: llr[1][0]=90.
		llr[1][1]+= latpad	#.5
		#if llr[1][1]>180.: llr[1][1]=180.
		
		cntr=[float(llr[0][0])+(llr[1][0]-float(llr[0][0]))/2.0, float(llr[0][1])+(llr[1][1]-float(llr[0][1]))/2.0]
		if self.catmap==None or (mapLLlat!=None or mapLLlon!=None or mapURlat!=None or mapURlon!=None): self.catmap=Basemap(llcrnrlon=llr[0][1], llcrnrlat=llr[0][0], urcrnrlon=llr[1][1], urcrnrlat=llr[1][0], resolution=self.mapres, projection='tmerc', lon_0=cntr[1], lat_0=cntr[0])
		catmap=self.catmap
		
		canvas=FigureCanvas(f0)
		if myaxis==None: myaxis=f0.add_axes([0,0,1,1])
		#catmap.ax=f0.add_axes([0,0,1,1])
		catmap.ax=myaxis
		#f0.set_figsize_inches((8/catmap.aspect,8.))
		#
		catmap.drawcoastlines(color='gray', zorder=0)
		catmap.drawcountries(color='gray', zorder=0)
		catmap.drawstates(color='gray', zorder=0)
		catmap.drawrivers(color='gray', zorder=0)
		catmap.fillcontinents(color='beige', zorder=0)
		#
		# , deltaLat=1.0, deltaLon=1.0
		if deltaLat!=None:
			#catmap.drawparallels(range(int(llr[0][0]-2.0), int(llr[1][0]+2.0)), color='k', labels=parallel_labels)
			theseLats = [float(int(llr[0][0]-2.0*deltaLat))]
			while theseLats[-1]<=float(int(llr[1][0]+2.0*deltaLat)):
				theseLats+=[theseLats[-1]+deltaLat]
			catmap.drawparallels(theseLats, color='k', labels=parallel_labels)
		#	
		if deltaLon!=None:
			#catmap.drawmeridians(range(int(llr[0][1]-2.0), int(llr[1][1]+2.0)), color='k', labels=meridian_labels)
			theseLons = [float(int(llr[0][1]-2.0))]
			while theseLons[-1]<=float(int(llr[1][1]+2.0)):
				theseLons+=[theseLons[-1]+deltaLon]
			catmap.drawmeridians(theseLons, color='k', labels=meridian_labels)
		#
		if plotevents:
			xfull, yfull=catmap(list(map(operator.itemgetter(2), catalog)), list(map(operator.itemgetter(1), catalog)))
			epx, epy=catmap(epicenter[0], epicenter[1])
			#catmap.plot(xfull, yfull, 'b,', label='Full Catalog')
			catmap.plot(xfull, yfull, eqicon, label='earthquakes', alpha=.5, zorder=2)
			#catmap.plot(epx, epy, 'ro', zorder=2)
		
		# if we are inclned to save:
		if doSave and saveName!=None: canvas.print_figure(saveName)
		
		#
		#ax=plt.gca()
		#el = Ellipse((self.tLon, self.tLat), 2.0*self.tA, 2.0*self.tB, -self.tTheta, facecolor='b', alpha=0.4)
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
		
		return catmap
	#
	def getTargMag(self, m, mc=None, mt=7.6):
		if mc==None:
			mc=min(list(map(operator.itemgetter(3), self.getcat())))
		if mt==None: mt=self.mt
		if m<mt:
			# "small" earthquake
			winlen=10**(m-2.0-mc)	# where 2.0 is dmstar + dmprime
		if m>=mt:
			dms = (1.0*(mt-mc) + 1.5*(m-mt))/(m-mc)
			#winlen = 10**(1.0*(mt-mc) + 1.5*(targmag-mt-1.0) - dms)
			winlen = 10**(1.0*(mt-mc) + 1.5*(m-mt) - 2.0*dms)
		#
		return winlen

	def rbomoriQuadPlot(self, catnum=0, mc=None, winlen=None, targmag=None, rbthresh=1.0, bigmag=6.0, fignum=0, intlist=None, rbavelen=None, thislw=1.0, mainEV=None, plotevents=False, mapcatnum=None, rbLegLoc='best', logZ=1.0, mapLLlat=None, mapLLlon=None, mapURlat=None, mapURlon=None, deltaLat=1.0, deltaLon=1.0, weighted=False):
		# make an awesome quad plot: omori-times, rbRatios, mag-seismicity, map-catalog
		# catalog -> a yodapy.eqcatalog() object
		# thislw: linewidth
		# deltaLat, deltaLon: intervals between lat/lon lines.
		# for no lines, use None.
		# weighted: weight each value with 1/chi_sqr, where chi_sqr is from a local linear fit.
		#
		if mc==None:
			# let's just guess that the catalog was selected with a valid mc...
			mc=min(list(map(operator.itemgetter(3), self.cat, 3)))
		if targmag!=None:
			# estimate winlen from a target magnitude.
			mt=self.mt
			if targmag<mt:
				# "small" earthquake
				winlen=10**(targmag-2.0-mc)	# where 2.0 is dmstar + dmprime
			if targmag>=mt:
				dms = (1.0*(mt-mc) + 1.5*(targmag-mt))/(targmag-mc)
				#winlen = 10**(1.0*(mt-mc) + 1.5*(targmag-mt-1.0) - dms)
				winlen = 10**(1.0*(mt-mc) + 1.5*(targmag-mt-0.0) - 2.0*dms)
			#
			#winlen=int(10*round(winlen/10))
			#print "winlen0: %d" % winlen
			winlen=int(round(winlen,-1))
			#print "winlen0: %d" % winlen
			if winlen<1: winlen=1
			#print "winlen: %d" % winlen
		#
		# expected magnitude:
		#mofN = math.log10(winlen) + 2.0 + mc
		#if mofN>mt:
			# this is harder...
		#	mofN = 
		#
		#
		if rbavelen==None:
			rbavelen=int(winlen/10.)
			if rbavelen==0: rbavelen=1
		print("avelen: %d" % rbavelen)
		#
		if mapcatnum==None: mapcatnum=catnum
		if logZ==None:
			logZ = math.log10(float(winlen))
		#
		logZexp=1.0/logZ		# so logZ is more intuitive...
		#
		#rbavelen=1	# rb-averaging length (aka, <nrb>_thisnum
		#print "winlen: %d" % winlen
		#if catalog==None: catalog=getMexicaliCat()
		#intlist=[25, 256, 512]
		if intlist==None: intlist=[int(winlen/4), int(winlen/2), winlen, int(winlen*1.5)]
		catalog=self
		if len(catalog.getcat(catnum))<(winlen+rbavelen): return None
	
		plt.figure(fignum)
		plt.clf()
		#ax0=plt.axes([.1,.1,.85, .35])
		# define axis (subplots) boundaries:
		ydelim=.03
		xdelim=.05
		xTS=[0.05, .5]
		yTS0=0.05
		dyTS=.3
			
		myaxes=[]
		nax=0
		# mags:
		x0=xTS[0]
		y0=yTS0+dyTS*nax
		#myaxes+=[plt.axes([xTS[0], y0, xTS[1], dyTS])]
		myaxes+=[plt.axes([.1, .03, .45, .3])]
		nax+=1
		# intervals:
		#x0=xTS[0]
		y0=yTS0+dyTS*nax
		#myaxes+=[plt.axes([xTS[0], y0, xTS[1], dyTS], sharex=myaxes[0])]
		myaxes+=[plt.axes([.1, .37, .45, .27], sharex=myaxes[0])]
		nax+=1
		# ratios:
		#x0=xTS[0]
		y0=yTS0+dyTS*nax
		#myaxes+=[plt.axes([xTS[0], y0, xTS[1], dyTS], sharex=myaxes[0])]
		myaxes+=[plt.axes([.1, .68, .45, .3], sharex=myaxes[0])]
		#
		# map:
		nax+=1
		xs=[xTS[1]+xdelim, .95]
		ys=[yTS0, 1.0]
		#myaxes+=[plt.axes([xs[0], xs[1], ys[0], ys[1]])]
		#myaxes+=[plt.axes([.6, .05, .35, .90], sharex=myaxes[0])]
		myaxes+=[plt.axes([.6, .05, .35, .90])]
		#
		# get RB ratios:
		#try:
		#	catalog.rb
		#except:
		#	catalog.rb=rbi.intervalRecordBreaker(None)
		#ratios=self.getIntervalRatios(minmag, windowLen, cat0, deltaipos, avlen)
	#	ratios=catalog.rb.getIntervalRatios(mc, winlen, catalog.getcat(catnum), 1)
		#
		#plotIntervals(self, intervals=[10, 100, 1000], minmag=2.0, catalog=None, fignum=0, dtmlatlonmagCols=[0,1,2,3], plotDates=[None, None], thisAxes=None):
		#catalog.plotIntervals(intlist, mc, catalog.getcat(catnum), fignum, [0,1,2,3], [None, None], [myaxes[0], myaxes[1]])
		#format: plotInts(self, intervals=[10, 100, 1000], catalog=None, minmag=2.0, ax=None, dtmlatlonmagCols=[0,1,2,3], plotDates=[None, None]):
		epicen=None
		if mainEV!=None: epicen=[mainEV[2], mainEV[1]]
		#
		catalog.plotInts(intervals=intlist, catalog=catalog.getcat(catnum), minmag=mc, ax=myaxes[1], thislw=thislw, legendPos='upper left')
		#
		#myaxes[1].set_label('mean intervals $\\tau$')
		myaxes[1].set_ylabel('mean intervals $\\ < \\tau >$', size=14)
		#
		# format: #plotMags(self, catalog=None, minmag=2.0, ax=None, dtmlatlonmagCols=[0,1,2,3], plotDates=[None, None]):	
		catalog.plotMags(catalog.getcat(catnum), mc, myaxes[0])
		myaxes[0].set_ylabel('mag', size=14)
		#
		# plotIntervalRatiosAx(self, minmag=3.0, windowLen=10, cat0=None, hitThreshold=1.0, bigmag=5.0, thisAx=None, ratios=None, deltaipos=1, avlen=1, mainEV=None)
		# let's pull this out of the recordbreaking module (and get rid of the rb module...).
		
		#catalog.rb.plotIntervalRatiosAx(minmag=mc, windowLen=winlen, cat0=catalog.getcat(catnum), hitThreshold=rbthresh, bigmag=bigmag, thisAx=myaxes[2], ratios=None, deltaipos=1, avlen=rbavelen, mainEV=mainEV, rbLegLoc=rbLegLoc, logZ=logZ)
		myaxes[2].set_ylabel('RB ratio $r(N=%d)$' % winlen, size=14)
		#
		# plotIntervalRatiosAx(winlen=10, cat=None, hitThreshold=1.0, bigmag=None, thisAx=None, ratios=None, delta_t=1, avlen=1, mainEv=None, logZ=None, rbLegLoc='best', reverse=False):
		if weighted == True:
			catalog.plotWeightedIntervalRatiosAx(winlen=winlen, cat=catalog.getcat(catnum), hitThreshold=rbthresh, bigmag=bigmag, thisAx=myaxes[2], ratios=None, delta_t=1, avlen=rbavelen, mainEv=mainEV, logZ=logZ, rbLegLoc=rbLegLoc, reverse=False)
		else:
			catalog.plotIntervalRatiosAx(winlen=winlen, cat=catalog.getcat(catnum), hitThreshold=rbthresh, bigmag=bigmag, thisAx=myaxes[2], ratios=None, delta_t=1, avlen=rbavelen, mainEv=mainEV, logZ=logZ, rbLegLoc=rbLegLoc, reverse=False)
	
		#plt.figure(fignum)
		myfsize=12
		myaxes[1].text(.1, .1, '$<\\tau>$', None, rotation='vertical', size=myfsize)
		myaxes[0].text(.1, .1, 'mags', None, rotation='vertical', size=myfsize)
		myaxes[2].text(.1, .1, '$r(1) = frac{n_{rb-large}}{n{rb-small}}$', None, rotation='vertical', size=myfsize)
	
		#
		# plotCatMap(self, catalog=None, doShow=True, doSave=False, saveName='catalogPlot.png', epicenter=None, legendLoc='upper left', doCLF=True, eqicon='b,', myaxis=None, fignum=None, padfactor=.25, plotevents=True)
		X=catalog.plotCatMap(catalog=catalog.getcat(mapcatnum), doShow=True, doSave=False, saveName=None, epicenter=epicen, legendLoc='best', doCLF=False, eqicon='b,', myaxis=myaxes[3], fignum=0, padfactor=.15, plotevents=plotevents, mapLLlat=mapLLlat, mapLLlon=mapLLlon, mapURlat=mapURlat, mapURlon=mapURlon, deltaLat=deltaLat, deltaLon=deltaLon)
		# and large events:
		bigEq=[]
		bigeqindex=0
		#print catalog.catmap
		for rw in catalog.getcat(catnum):
			if rw[3]>=bigmag and rw[0]>=catalog.getMainEvent(catalog.getcat(catnum))[0]:
				eqx,eqy = catalog.catmap(rw[2], rw[1])
				#catalog.catmap.plot(eqx, eqy, '*', label='m=%f, %s' % (rw[3], str(rw[0])))
				#catalog.catmap.plot(eqx, eqy, '*', label='m=%.2f, %d' % (rw[3], bigeqindex), ms=15)
				catalog.catmap.plot(eqx, eqy, '*', label='m=%.2f' % (rw[3]), ms=15)
				bigeqindex+=1
				#print "eq: %d, %f, %s" % (bigeqindex, rw[3], str(rw[0]))
		catalog.catmap.ax.legend(loc='best', numpoints=1)
	#	#if len(bigEqs)>0:
		#	#XX=catalog.plotCatMap(bigEqs, True, False, None, 'best', False, '*', myaxes[3], fignum)
		#	XX=catalog.plotCatMap(bigEqs, True, False, None, None, 'best', False, '*', myaxes[3], fignum)
		# sharex=ax0
		#
		#
		plt.show()
	
		# problem: the x-axis of RB-ratios is in float form, formatted as datetimes (basically using "plot" type functions). the intervals/magnitudes are date_plot
		# functions. in short, they have different x-axis variable types, so we can't share the rb-ratios x-axis with the other two. it probably makes sense to convert
		# the interval-ratio plots to the rb-ratio style (because the rb-rations uses the "fillbetween" function).
	
		catalog.axlist=myaxes
		return catalog

	def plotIntervalRatiosAx(self, winlen=10, cat=None, hitThreshold=1.0, bigmag=None, thisAx=None, ratios=None, delta_t=1, avlen=1, mainEv=None, logZ=None, rbLegLoc='best', reverse=False):
		# getNRBratios()
		# def getNRBratios(self, intervals=None, winlen=10, delta_t=1, reverse=False)
		# def rb.plotIntervalRatiosAx(self, minmag=3.0, windowLen=10, cat0=None, hitThreshold=1.0, bigmag=5.0, thisAx=None, ratios=None, deltaipos=1, avlen=1, mainEV=None, logZ=1.0, rbLegLoc='best')
		#
		if ratios==None:
			intervals = self.getIntervals(interval_length=1, catList=cat)
			ratios = self.getNRBratios(intervals=intervals, winlen=winlen, delta_t=delta_t, reverse=reverse)
		if thisAx==None:
			thisAx=plt.gca()
		if thisAx==None:
			plt.figure()
			thisAx=plt.gca()
		#
		#if len(ratios[0])<6: ratios[0]+=[ratios[0][4]]
		for i in range(1,len(ratios)+1):
			# "i" is the index of the next entry, like subset = fullset[(i-len):i], returning up to the (i-1)th entry.
			if len(ratios[i-1])>=6:
				#ratios[i-1]=ratios[i-1][:-1]	# remove mean value entry and recalculate (just in case we've changed something).
				continue	# this is sloppy, but we've already averaged.
							# probably, the right thing to do is to do a pre-loop to strip out the averaged values before hand, or
							# create a proper class object.
			#
			i0=max(0, i-avlen)	# early entries we average over what we've got so far.
			#theseRs=map(math.log10, map(operator.itemgetter(4), ratios[i0:i]))
			theseRs = list(map(operator.itemgetter(4), ratios[i0:i]))		# note: returning through the (i-1)th entry, hence i -> len()+1.
			#ratios[i-1]+=[scipy.mean(theseRs)]
			ratios[i-1]+=[scipy.prod(theseRs)**(1.0/len(theseRs))]	# ... and this mean value corresponds to the (i-1)th entry.
			#
		#thisAx.set_yscale('log')
		X = list(map(operator.itemgetter(1), ratios))
		Y = list(map(operator.itemgetter(-1), ratios))
		X2,Y2 = self.zeroFillInts(X,Y)
		Ythresh = hitThreshold*scipy.ones(len(Y2))
		Ygt = list(map(scipy.greater_equal, Y2, Ythresh))	# scipy.greater_equal(), scipy.greter()?
		Ylt = list(map(scipy.less_equal,  Y2, Ythresh))
		#
		thisAx.set_yscale('log')
		#print "xylen: %d, %d" % (len(X), len(Y))
		#thisAx.plot(X, Y, 'k-', lw=1)
		thisAx.plot([X[0], X[-1]], [1., 1.], 'k-', zorder=0)
		thisAx.plot([X[0], X[-1]], [hitThreshold, hitThreshold], 'k--', zorder=0)
		#
		thisAx.fill_between(X2,Y2,y2=Ythresh, where=Ylt, color='r', alpha=.9)
		thisAx.fill_between(X2,Y2, y2=Ythresh, where=Ygt, color='b', alpha=.9)
		#
		return ratios
	#
	def plotWeightedIntervalRatiosAx(self, winlen=10, cat=None, hitThreshold=1.0, bigmag=None, thisAx=None, ratios=None, delta_t=1, avlen=1, mainEv=None, logZ=None, rbLegLoc='best', reverse=False):
		# plot weighted intervals:
		# fit each interval sub-set to a line; weight by the chi-sqr.
		# note this also suggests using linear interpolated segments rather than just average value for smoothing... something to do later.
		#
		# getNRBratios()
		# def getNRBratios(self, intervals=None, winlen=10, delta_t=1, reverse=False)
		# def rb.plotIntervalRatiosAx(self, minmag=3.0, windowLen=10, cat0=None, hitThreshold=1.0, bigmag=5.0, thisAx=None, ratios=None, deltaipos=1, avlen=1, mainEV=None, logZ=1.0, rbLegLoc='best')
		#
		if ratios==None:
			intervals = self.getIntervals(interval_length=1, catList=cat)
			ratios = self.getNRBratios(intervals=intervals, winlen=winlen, delta_t=delta_t, reverse=reverse)
		if thisAx==None:
			thisAx=plt.gca()
		if thisAx==None:
			plt.figure()
			thisAx=plt.gca()
		#
		#if len(ratios[0])<6: ratios[0]+=[ratios[0][4]]
		for i in range(1,len(ratios)+1):
			# "i" is the index of the next entry, like subset = fullset[(i-len):i], returning up to the (i-1)th entry.
			if len(ratios[i-1])>5:
				#ratios[i-1]=ratios[i-1][:-1]	# remove mean value entry and recalculate (just in case we've changed something).
				continue	# this is sloppy, but we've already averaged.
							# probably, the right thing to do is to do a pre-loop to strip out the averaged values before hand, or
							# create a proper class object.
			#
			i0=max(0, i-avlen)	# early entries we average over what we've got so far.
			#theseRs=map(math.log10, map(operator.itemgetter(4), ratios[i0:i]))
			theseRs = list(map(operator.itemgetter(4), ratios[i0:i]))		# note: returning through the (i-1)th entry, hence i -> len()+1.
			#ratios[i-1]+=[scipy.mean(theseRs)]
			ratios[i-1]+=[scipy.prod(theseRs)**(1.0/len(theseRs))]	# ... and this (geometric) mean value corresponds to the (i-1)th entry.
			#
		#thisAx.set_yscale('log')
		#
		# the weighted fit will be log(r)/chi_sqr
		#
		X = list(map(operator.itemgetter(1), ratios))
		Y0 = list(map(operator.itemgetter(-1), ratios))
		chi_sqrs=self.get_ratio_fits(ratios=ratios, fitlen=avlen, x_col=1, y_col=5)
		#
		# and let's use the std-dev, not the var so numbers don't explode so badly (should not make a huge difference).
		chi_sqrs = [x**.5 for x in chi_sqrs]
		#
		# ... and chi_sqrs can be really small, so Y --> really big. too big for computer-numbers, so let's quasi-normalize this;
		# what we need is some measure of relative weight. the absolute weight is not important.
		mean_chi_sqr = numpy.mean(chi_sqrs)
		chi_sqrs = [x/mean_chi_sqr for x in chi_sqrs]		# and a bit of normalization to keep the numbers under control...
		#
		#Y = numpy.array(map(math.log10, Y0)[-len(chi_sqrs):])/numpy.array(chi_sqrs)
		Y=map(math.log10, Y0)[-len(chi_sqrs):]
		Y=[Y[i]/x for i, x in enumerate(chi_sqrs)]
		#Y = Y.tolist()
		Y = [math.pow(10., x) for x in Y]
		X = X[-len(Y):]
		#
		#return [X, Y, Y0, chi_sqrs]
		X2,Y2 = self.zeroFillInts(X,Y, dolog=False)
		Ythresh = hitThreshold*scipy.ones(len(Y2))
		Ygt = list(map(scipy.greater_equal, Y2, Ythresh))	# scipy.greater_equal(), scipy.greter()?
		Ylt = list(map(scipy.less_equal,  Y2, Ythresh))
		#
		thisAx.set_yscale('log')
		#print "xylen: %d, %d" % (len(X), len(Y))
		#thisAx.plot(X, Y, 'k-', lw=1)
		thisAx.plot([X[0], X[-1]], [1., 1.], 'k-', zorder=0)
		thisAx.plot([X[0], X[-1]], [hitThreshold, hitThreshold], 'k--', zorder=0)
		#
		thisAx.fill_between(X2,Y2,y2=Ythresh, where=Ylt, color='r', alpha=.9)
		thisAx.fill_between(X2,Y2, y2=Ythresh, where=Ygt, color='b', alpha=.9)
		#
		return ratios
	
	def zeroFillInts(self, X0,Y0, thresh=1.0, dolog=True):
		X=scipy.array(X0).copy().tolist()
		Y=scipy.array(Y0).copy().tolist()
		if type(X0[0])!=float:
			# this is a little bit sloppy, but it will do for now.
			X=list(map(mpd.date2num, X0))
		#
		i=1
		while i<len(X):
			# we want to stick a y=0 where y crosses thresh. note, this is default parameterized for log data.
			#if ((Y[i]>=thresh)!=(Y[i-1]>=thresh)):
			if ((Y[i]>thresh and Y[i-1]<thresh) or (Y[i]<thresh and Y[i-1]>thresh)):
				#we've crossed the line...
				y2=Y[i]
				y1=Y[i-1]
				x2=X[i]
				x1=X[i-1]
				if dolog==True:
					b = math.log10(y2/y1)/(x2-x1)
					deltax = math.log10(thresh/y1)/b
				if dolog==False:
					b = (y2-y1)/(x2-x1)
					deltax = (thresh-y1)/b
				#
				newx = x1+deltax
				X.insert(i, newx)
				Y.insert(i, thresh)
				
				#
			i+=1
			#
		#
		X=mpd.num2date(X)
		#
		return (X,Y)
	#
	def get_ratio_fits(self, ratios=None, fitlen=1, x_col=1, y_col=5):
		# ratios: return from rbratios (or any time series)
		# fitlen: sequence length over which to fit.
		# ratiocol: data column
		#
		# simplified: return [chi_sqrs]
		#
		while len(ratios[0])<(y_col+1): y_col-=1
		#
		chi_sqrs=[]	# and note that this will return a set shorter than the input.
		#
		r_vals = list(map(operator.itemgetter(y_col), ratios))
		X_vals = list(map(operator.itemgetter(x_col), ratios))
		#
		if type(X_vals[0])==type(dtm.datetime.now()):
			X_vals = list(map(mpd.date2num, X_vals))
		#
		for i in range(3,len(ratios)):
			# ... and we fit the logs of the inputs (linearized values)...
			these_r = list(map(math.log10, r_vals[max(0, (i-fitlen)):i]))
			#if type(X_vals[0])==type(dtm.datetime.now()):
			#	these_x = map(mpd.date2num, X_vals[(i-n):i])
			#else:
			#	these_x = X_vals[(i-n):i]
			#
			these_x = X_vals[max(0, (i-fitlen)):i]
			#
			these_p = (0., 0.)	
			#
			# now, get a line-fit for this segment...
			#
			#print "range: ", max(0, (i-fitlen)), ", ", i
			#print "these_x: ", these_x
			#print "these_y: ", these_r
			lf  = linefit.linefit([these_x, these_r])
			lf.doFit()
			#
			#fitsets[n]['means']          += [meanval]
			##fitsets[n]['fit_prams']      += [list(fit_vals) + [cov]]
			##fitsets[n]['mean_fit_prams'] += [list(mean_fit_vals) + [chisqr]]
			#fitsets[n]['fit_prams'] += [[lf.a, lf.b, lf.meanVar()]]
			#fitsets[n]['mean_fit_prams'] += [[lf2.a, lf2.b, lf2.meanVar()]]
			#
			chi_sqrs += [lf.meanVar()]
			#
			#print 'fitses: ', fitsets[n]['fit_prams'][-1]
		#
		#
		return chi_sqrs
	#
	#
	def testMap(self):
		import pickle
		import time
		#
		fig=plt.figure()
		#
		t1 = time.clock()
		#m = Basemap(width=920000, height=1100000, resolution='f', projection='tmerc', lon_0=-4.2, lat_0=54.6)
		#m = Basemap(llcrnrlon=-11.0, llcrnrlat=45.0, urcrnrlon=3.0, urcrnrlat=59.0, resolution='f', projection='tmerc', lon_0=-4.2, lat_0=54.6)
		
		lllon=-115
		lllat=32
		urlon=-105
		urlat=42
		lon0=lllon + (urlon-lllon)/2.0
		lat0=lllat + (urlat-urlat)/2.0
		print("center: %f, %f" % (lon0, lat0))
		m = Basemap(llcrnrlon=lllon, llcrnrlat=lllat, urcrnrlon=urlon, urcrnrlat=urlat, resolution=self.mapres, projection='tmerc', lon_0=lon0, lat_0=lat0)
		m.drawcountries()
		m.drawrivers()
		print(time.clock()-t1,' secs to create original Basemap instance')

		# cPickle the class instance.
		pickle.dump(m,open('map.pickle','wb'),-1)

		# clear the figure
		plt.clf()
		# read cPickle back in and plot it again (should be much faster).
		t1 = time.clock()
		m2 = pickle.load(open('map.pickle','rb'))
		# draw coastlines and fill continents.
		m.drawcoastlines()
		# fill continents and lakes
		m.fillcontinents(color='coral',lake_color='aqua')
		# draw political boundaries.
		m.drawcountries(linewidth=1)
		# fill map projection region light blue (this will
		# paint ocean areas same color as lakes).
		m.drawmapboundary(fill_color='aqua')
		# draw major rivers.
		m.drawrivers(color='b')
		print(time.clock()-t1,' secs to plot using using a pickled Basemap instance')
		# draw parallels
		circles = np.arange(48,65,2).tolist()
		m.drawparallels(circles,labels=[1,1,0,0])
		# draw meridians
		meridians = np.arange(-12,13,2)
		m.drawmeridians(meridians,labels=[0,0,1,1])
		plt.title("High-Res British Isles",y=1.04)
		plt.show()
	
	def setSpecialCatSQL(self, catname='parkfield'):
		#reload(yp)
		if catname in ['parkfield', 'pf', 'PF', 'park']:
			#theta=40.0, clat=35.9, clon=-120.5, ra=.4, rb=.15
			#self.setCatFromSQL(dtm.datetime(1969,1,1, tzinfo=tzutc), dtm.datetime.now(tzutc), [34.9, 36.9], [-121.5, -119.5], 1.5, "Earthquakes", 523, 'asc')
			self.setCatFromSQL(dtm.datetime(1972,1,1, tzinfo=tzutc), dtm.datetime.now(tzutc), [34.4, 37.4], [-121.11, -119.4], 1.5, "Earthquakes", 523, 'asc')
			self.addEllipCat('PFshock (.8 x .15)', self.cat, 40.0, 35.9, -120.5, 0.8, 0.15)
			self.addEllipCat('PFshock (.4 x .15)', self.cat, 40.0, 35.9, -120.5, 0.4, 0.15)
		if catname == 'PF5yr':
			#theta=40.0, clat=35.9, clon=-120.5, ra=.4, rb=.15
			#self.setCatFromSQL(dtm.datetime(1969,1,1, tzinfo=tzutc), dtm.datetime.now(tzutc), [34.9, 36.9], [-121.5, -119.5], 1.5, "Earthquakes", 523, 'asc')
			self.setCatFromSQL(dtm.datetime(1999, 9, 28, tzinfo=tzutc), dtm.datetime(2009,9,29, tzinfo=tzutc), [34.4, 37.4], [-121.11, -119.4], 1.5, "Earthquakes", 523, 'asc')
			self.addEllipCat('PFshock (.8 x .15)', self.cat, 40.0, 35.9, -120.5, 0.8, 0.15)
			self.addEllipCat('PFshock (.4 x .15)', self.cat, 40.0, 35.9, -120.5, 0.4, 0.15)
		if catname in ['taiwan']:
			self.setCatFromSQL(dtm.datetime(1980,1,1, tzinfo=tzutc), dtm.datetime(2010,6,1, tzinfo=tzutc), [-90, 90], [-180, 180], 2.0, 'Earthquakes', 21)
			
	def setCatFromSQL(self, startDate=dtm.datetime(1999,9,28, 17,15,24, tzinfo=tzutc), endDate=dtm.datetime(2009,9,28, 17,15,24, tzinfo=tzutc), lats=[32.0, 37.0], lons=[-125.0, -115.0], minmag=3.0, catalogName='Earthquakes', catalogID=523, ordering='asc'):
		self.cat=self.getCatFromSQL(startDate, endDate, lats, lons, minmag, catalogName, catalogID, ordering)
		return None
	
	def getCatFromSQL(self, startDate=dtm.datetime(1999,9,28, 17,15,24, tzinfo=tzutc), endDate=dtm.datetime(2009,9,28, 17,15,24, tzinfo=tzutc), lats=[32.0, 37.0], lons=[-125.0, -115.0], minmag=2.0, catalogName='Earthquakes', catalogID=523, ordering='asc'):
		# return a catalog:
		if lats[0]>lats[1]: lats.reverse()
		if lons[0]>lons[1]: lons.reverse()
		#
		if startDate>endDate:
			middledate=startDate
			startDate=endDate
			endDate=middledate
			middledate=None
		if ordering not in ['asc', 'desc']: ordering='desc'
		
		import _mysql
		import MySQLdb
		#
		#sqlHost = 'localhost'
		sqlHost = self.sqlhost
		sqlUser = 'myoder'
		sqlPassword = 'yoda'
		sqlPort = self.sqlport
		sqlDB = 'QuakeData'
		con=MySQLdb.connect(host=sqlHost, user=sqlUser, passwd=sqlPassword, port=sqlPort, db=sqlDB)
		c1=con.cursor()
		sqlstr='select eventDateTime, lat, lon, mag from %s where catalogID=%d and lat between %f and %f and lon between %f and %f and mag>=%f and eventDateTime between \'%s\' and \'%s\' order by eventDateTime %s' % (catalogName, catalogID, lats[0], lats[1], lons[0], lons[1], minmag, str(startDate), str(endDate), ordering)
		catList=[]
		#print sqlstr
		#
		c1.execute(sqlstr)
		rw=c1.fetchone()
		while rw!=None:
		#	# spin through the cursor; write a catalog. note formatting choices...
			catList+=[[rw[0], float(rw[1]), float(rw[2]), float(rw[3])]]
			rw=c1.fetchone()
		#catList=self.fetchall()
		c1.close()
		con.close()
		# now we have a catalog of the parkfield area (note, it is partially defined by our "parkfieldquakes" MySQL view.
		#
		#makeShockCat(incat, outcat)
		#makeShockCat(fullcatout, shockcatout)
		return catList
	
	def plotGRdistsFromTo(self, frmDt=None, toDt=None, catlist=None, fignum=0):
		# plot GRdist between two dates for a catalog.
		# assume standard catalog format.
		if catlist==None: catlist=self.cat
		if frmDt==None: frmDt=catlist[0][0]
		if toDt==None: toDt=catlist[-1][0]
		#
		mymags=[]
		for rw in catlist:
			if rw[0]>=frmDt and rw[0]<=toDt: mymags+=[rw[3]]
			# if rw[0]>toDt: break	# unless the cat is out of sequence. computers are fast; just spin through for now.
		#
		return self.plotGRdist(mags=mymags, fignum=fignum)
		
	def plotGRdist(self, mags=None, doShow=True, fname='GRdist.png', plotTitle="Magnitude Distribution", fignum=0):
		# mags: a 1D array of magnitudes
		if mags==None: mags=list(map(operator.itemgetter(3), self.cat))
		# if mags rows are not scalar, assume a full standard type catalog has been passed.
		try:
			if len(mags[0])>=3: mags=list(map(operator.itemgetter(3), mags))
		except TypeError:
			# a list of scalars will throw a "can't get len." error. we should be able to skip without doing anything at all.
			# maybe a better approach is to test the type of mags[0] for list or tuple...
			dummyvar=None	# place-holder
		#
		mags.sort()
		# get rid of biggest event (probably a large off-GR earthquake):
		#mags.pop()
		#mags.reverse()
		#print mags
		#print len(mags)
		
		if doShow==True or fname!=None:
			# make a plot and show and/or save
			Y=list(range(1, len(mags)+1))
			Y.reverse()
			#Y=frange(1, len(mags)+1, -1)
			#print Y
			#print len(Y)
			plt.figure(fignum)
			plt.clf()
			plt.semilogy(mags, Y, '.-')
			plt.xlabel("Magnitude, m")
			plt.ylabel("Number of Events, n")
			plt.title(plotTitle)
			if fname!=None: plt.savefig(fname)
			if doShow: plt.show()
		
		return mags	

	def getIntervals(self, catList=None, interval_length=1):
		# returns raw intervals in days (aka, <dt(winLen=10)> ~ 10*<dt(winLen=1)> )
		if catList==None: catList = self.getcat(0)
		#if interval_length==None: interval_length = 10
		#
		# let's get fancy and allow for the list to be reverse-sorted.
		# and in fact, let's sacrifice a bit of speed to be sure we're getting a sorted list.
		# we only need the dates...
		event_dates=list(map(operator.itemgetter(0), catList))
		mylist=list(map(mpd.date2num, event_dates))	# float type in units of "days"
		mylist.sort()	# now in ascending order...
		mylist=scipy.array(mylist)
		intervals=mylist[interval_length:]-mylist[0:-interval_length]
		#print intervals[0:5]
		#
		return_intervals = scipy.array([event_dates[interval_length:], intervals])
		#return_intervals = return_intervals.transpose()
		#ary_shape=return_intervals.shape
		#print "shape: ", ary_shape, len(event_dates[interval_length:]), len(intervals)
		#return_intervals.shape=(ary_shape[1], ary_shape[0])
		#print "shape: ", return_intervals.shape, len(event_dates[interval_length:]), len(intervals)
		#		
		return return_intervals.transpose()
	
	def getNRBratios(self, intervals=None, winlen=10, delta_t=1, reverse=False, catnum=0):
		# slicker replacement for the current record-breaker engine. we'll do our own plots too.
		# intervals should be like [ [t_end, dt], [t_end, dt], ... ]
		# 
		# return same as current module -- something like: [ [i, t, r**1/log(Z)], ...]
		#
		# inervals are in ascending time order.
		#
		if intervals==None: intervals=self.getIntervals(catList=self.getcat(catnum), interval_length=1)
		# intervals likeL [[date, interval], ...] in ascending order.
		#
		rbratios=[]
		#i=winlen
		#imax=len(intervals)
		#while i<imax:
		#print intervals[0:10]
		for i in range(winlen, len(intervals)):
			#print "***", intervals[(i-winlen): i]
			thisX=list(map(operator.itemgetter(1), intervals[i-winlen:i]))
			thisrbdata=self.getnrbs(thisX, reverse=reverse)
			# like: [xmax, xmin, nbigger, nsmaller]
			r=(thisrbdata[0]/float(thisrbdata[1]))
			rbratios+=[[i-1, intervals[i-1][0], thisrbdata[0], thisrbdata[1], r]]
		#
		return rbratios
	#
	def getnrbs(self, X, reverse=False):
		# (this still needs to be verified).
		#
		# record-breaking stats from a subset X (start to finish)
		# in "reverse", go back to front.
		#
		nbigger, nsmaller = 1, 1
		if reverse==False:
			xmax, xmin = X[0], X[0]
			for x in X[1:]:
				# the first element is everybody's record-breaking event...
				if x>xmax:
					xmax=x
					nbigger+=1
				if x<xmin:
					xmin=x
					nsmaller+=1
			#
		if reverse==True:
			xmax, xmin = X[-1], X[-1]
			for i in range(2, len(X)+1):
				x=X[-i]
				# the first element is everybody's record-breaking event...
				if x>xmax:
					xmax=x
					nbigger+=1
				if x<xmin:
					xmin=x
					nsmaller+=1
		#
		# there's really not much point in returning the largest/smallest. if we become interested
		# in magnitude, it is more meaningful to look at the sequence of rb events.
		#return [xmax, xmin, nbigger, nsmaller]
		return [nbigger, nsmaller]
	#	
	def getIntervals_depricated(self, catList=None, winLen=None):
		if catList==None: catList = self.cat
		if winLen==None: winLen = 10
		#
		catLen=len(catList)
		i=(catLen-1-winLen)	# start winLen positions from the end.
		thisInterval=0
		#N=1
		intervals=[]	# [[eventDateTime, totalInterval]]
		while i>=0:
			#
			thisInterval=mpd.date2num(catList[i+winLen][0])-mpd.date2num(catList[i][0])
			intervals+=[[catList[i+winLen][0], thisInterval]]
			i-=1
		#
		#return [intervals, catList]
		return intervals
	#	
	def getDistances(self, catList=None, center=None):
		if catList==None: catList=self.cat
		if center==None:
			ctr = self.getMainEvent()
			#print 'ctr: ', ctr
			center = [ctr[2], ctr[1]]	# lon, lat (currently in degrees)
		#
		#phis = center[0]*deg2rad
		#lambs = center[1]*deg2rad		# "source" prams now in rads
		phis = center[1]*deg2rad
		lambs = center[0]*deg2rad		# "source" prams now in rads
		#
		phis0 = center[0]					
		lambs0 = center[1]				#phis0, lambs0 still in degrees...
		#
		# can we just use basemap for this, or do we have to do our own spherical geometry?
		#bm = Basemap(width=width,height=width,projection='aeqd', lat_0=center[1],lon_0=center[0])
		#		#
		dists = []	# [dtm, dx, dy]	# one day, maybe dz as well?
		for rw in catList:
			#print rw
			#thisx = rw[2]
			#thisy = rw[1]
			#thisx, thisy = bm(catList[2], catList[1])
			#phif = rw[2]*deg2rad
			#lambf = rw[1]*deg2rad		# catalog coords in degrees, now converted to rads...
			phif = rw[1]*deg2rad
			lambf = rw[2]*deg2rad		# catalog coords in degrees, now converted to rads...
			#
			#print 'phif: ', phif
			#print 'lambf: ', lambf
			#
			dphi = (phif - phis)
			dlambda = (lambf - lambs)		# in rads...
			#
			thisx=Rearth*dphi
			thisy=Rearth*dlambda
			#
			# now, a bunch of approximations for sighat, the angle between the two points.
			# but has problems with small angles
			#R0 = Rearth*math.sqrt( (dlambda)**2.0 + (math.cos(lambf)*phif-math.cos(lambs)*phis)**2.0 )
			#x1=math.cos(lambf)*phif
			x1=math.cos(lambs)*phif
			x2=math.cos(lambs)*phis	#it's a right-triangle, so we want cos along one coord (probably the mean lambda would be better).
			R0 = Rearth*math.sqrt( (dlambda)**2.0 + (x1-x2)**2.0 )
			#
			# now, various approximations from spherical geometry...  
			#print phis, phif, dlambda
			sighat1 = math.acos(math.sin(phis)*math.sin(phif) + math.cos(phis)*math.cos(phif)*math.cos(dlambda))
			R1 = Rearth * sighat1
			#
			sighat2 = 2.0*math.asin(math.sqrt((math.sin(dphi/2.))**2 + math.cos(phis)*math.cos(phif)*(math.sin(dlambda/2.0))**2.0 ))
			R2 = Rearth * sighat2
			#
			# this one is supposed to be bulletproof:
			sighat3 = math.atan( math.sqrt((math.cos(phif)*math.sin(dlambda))**2.0 + (math.cos(phis)*math.sin(phif) - math.sin(phis)*math.cos(phif)*math.cos(dlambda))**2.0 ) / (math.sin(phis)*math.sin(phif) + math.cos(phis)*math.cos(phif)*math.cos(dlambda))  )
			R3 = Rearth * sighat3
			#
			# and now, use geographiclib to calculate distance:
			g1=ggp.WGS84.Inverse(lambs0, phis0, rw[1], rw[2])
			#g1=ggp.WGS84.Inverse(phis0, lambs0, rw[1], rw[2])
			# this returns a geodesic dictionary
			R4=g1['s12']/1000.0	# and converting from m to km.
			#
			# R0: rect-lin distance, R1: approx1, R2: approx2 (better than 1), R3: bomb-proof spherical, R4: from geographicLib,
			dists += [[rw[0], thisx, thisy, R0, R1, R2, R3, R4, rw[4]]]	# and R3 is the only distance we really need; rw[4]->depth.
		#
		return dists
		
	#
	def plotIntervals(self, intervals=[10, 100, 1000], minmag=2.0, catalog=None, fignum=0, dtmlatlonmagCols=[0,1,2,3], plotDates=[None, None], thisAxes=None):
		if type(plotDates).__name__!='list': plotDates=[None, None]
		while len(plotDates)<2: plotDates+=[None]
		#
		if catalog==None: catalog=self.cat
		#X = plotIntervals(intervals, minmag, catalog, fignum, dtmlatlonmagCols)
		#return X
		#zonedat=[35.9, -120.5, .4, .15, 40.0]
		cols=dtmlatlonmagCols	# for efficient notation
		#zonedat=[35.9, -120.5, .4, .05, 40.0]	# this will be done in advance of this function call, when the catalog is made.
		#minmag=2.0
		#dts=['1950-01-01', str(dtm.datetime.now(tzutc))]
		#sqlcat="Earthquakes"
		#catid=523
	
		#
		plt.figure(fignum)
		if thisAxes==None:
			#plt.figure(fignum)
			plt.clf()
			ax0=plt.axes([.1,.1,.85, .35])
			plt.xlabel("time")
			plt.ylabel("mags")
			ax1=plt.axes([.1, .55, .85, .35], sharex=ax0)
			plt.ylabel("mean interval")
			plt.xlabel("")
			plt.title("Mean intervals, $m_c=%s$" % str(minmag))
		else:
			ax0=thisAxes[0]
			ax1=thisAxes[1]
	
		#dtms=map(operator.itemgetter(cols[0]), catalog)
		#lats=map(operator.itemgetter(cols[1]), catalog)
		#lons=map(operator.itemgetter(cols[2]), catalog)
		#mags=map(operator.itemgetter(cols[3]), catalog)
		mags=[]
		activecat=[]
		for rw in catalog:
			if rw[cols[3]]<minmag: continue
			mags+=[[rw[cols[0]], rw[cols[3]]]]
			activecat+=[rw]
		mags=vlinePadList(mags, minmag-abs(minmag)*.1)	# return the mags data padded for vertical line style plotting. this is just a trick to get width=1 histograms.
		#
		ax0.plot_date(list(map(operator.itemgetter(0), mags)), list(map(operator.itemgetter(1), mags)), '-')
		shockints=[]
		
		#print "plotdates: %s" % str(plotDates)
		for wlen in intervals:
			#print "plotting for wlen: %d" % wlen
			##
			#shockints+=[getIntervals(catalog, wlen)]
			shockints+=[self.getIntervals(activecat, wlen)]
			#
			# trim off max/min date ends for prettier plots:
			#print "mindt: %s" % str(shockints[-1][0])
			while (plotDates[1]!=None and plotDates[1]<shockints[-1][0][0]): a=shockints[-1].pop(0)
			while plotDates[0]!=None and plotDates[0]>shockints[-1][-1][0]: a=shockints[-1].pop()
			#
			#plt.plot(map(operator.itemgetter(0), shockints[-1]), scipy.array(map(operator.itemgetter(1), shockints[-1]))/float(wlen), '-', label='winLen=%d' % wlen)
			#
			X=list(map(operator.itemgetter(0), shockints[-1]))
			# pylab.date2num(dtm)
			#XX=date2num(X)
			ax1.plot(X, scipy.array(list(map(operator.itemgetter(1), shockints[-1])))/float(wlen), '-', label='$N=%d$' % wlen, lw=1.0)
			#ax1.semilogy(map(operator.itemgetter(0), shockints[-1]), scipy.array(map(operator.itemgetter(1), shockints[-1]))/float(wlen), '-', label='winLen=%d' % wlen)
			# fg.autofmt_xdate()
			
		
			#plt.plot(map(operator.itemgetter(0), shockints[-1]), scipy.array(map(operator.itemgetter(1), shockints[-1])), '-', label='winLen=%d' % wlen)
		#
		plt.legend(loc='upper left')
			
		#plt.legend(loc='lower left')
	
	
		plt.show()
	
		return shockints
		
	def plotInts(self, intervals=[10, 100, 1000], catalog=None, minmag=None, ax=None, dtmlatlonmagCols=[0,1,2,3], plotDates=[None, None], thislw=1.0, legendPos='best'):
		# plot mean intervals only. having figured out how to put a bunch of independently generated plots onto a single canvas, we
		# split up some of these plots...
		#
		if type(plotDates).__name__!='list': plotDates=[None, None]
		while len(plotDates)<2: plotDates+=[None]
		#
		if catalog==None: catalog=self.cat
		try:
			if minmag==None: minmag=self.minmag
		except:
			if minmag==None: minmag=min(list(map(operator.itemgetter(3), catalog)))
		#
		cols=dtmlatlonmagCols	# for efficient notation
		activecat=[]
		for rw in catalog:
			if rw[cols[3]]>=minmag: activecat+=[rw]	# build active catalog of events m>mc
		#
		shockints=[]
		for wlen in intervals:
			shockints+=[self.getIntervals(activecat, wlen)]
			# trim off catalog elements outside min/max date range:
			while (plotDates[1]!=None and plotDates[1]<shockints[-1][0][0]): a=shockints[-1].pop(0)
			while plotDates[0]!=None and plotDates[0]>shockints[-1][-1][0]: a=shockints[-1].pop()
			#
			ax.plot(list(map(operator.itemgetter(0), shockints[-1])), scipy.array(list(map(operator.itemgetter(1), shockints[-1])))/float(wlen), '-', label='$N=%d$' % wlen, lw=thislw)
			#ax1.semilogy(map(operator.itemgetter(0), shockints[-1]), scipy.array(map(operator.itemgetter(1), shockints[-1]))/float(wlen), '-', label='winLen=%d' % wlen)
			
		
			#plt.plot(map(operator.itemgetter(0), shockints[-1]), scipy.array(map(operator.itemgetter(1), shockints[-1])), '-', label='winLen=%d' % wlen)
		#
		ax.legend(loc=legendPos)
			
		return shockints
		
	def plotMags(self, catalog=None, minmag=2.0, ax=None, dtmlatlonmagCols=[0,1,2,3], plotDates=[None, None]):
		if type(plotDates).__name__!='list': plotDates=[None, None]
		while len(plotDates)<2: plotDates+=[None]
		if ax==None: ax=plt.gca()
		#
		if catalog==None: catalog=self.cat
		cols=dtmlatlonmagCols	# for efficient notation
		#
		mags=[]
		activecat=[]
		for rw in catalog:
			if rw[cols[3]]>=minmag:
				mags+=[[rw[cols[0]], rw[cols[3]]]]
				ax.plot([rw[cols[0]], rw[cols[0]]], [minmag-.25, rw[cols[3]]], 'b-')
		#mags=vlinePadList(mags, minmag-abs(minmag)*.1)	# return the mags data padded for vertical line style plotting. this is just a trick to get width=1 histograms.
		#
		#ax.plot(map(operator.itemgetter(0), mags), map(operator.itemgetter(1), mags), '-')
		ax.legend(loc='best')
		return mags

	
def datetimeFromString(t):
	return mpd.num2date(mpd.datestr2num(t))
