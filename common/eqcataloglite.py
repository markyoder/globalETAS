import math
#import scipy
# matplotlib might give us a problem (RuntimeError: Failed to create /var/www/.matplotlib; consider setting MPLCONFIGDIR to a writable directory for matplotlib configuration data)

import string
import sys


#import scipy.optimize as spo
#import numpy

import os
import random
import time
#
#from threading import Thread
#
#
import datetime as dtm
import pytz
import calendar
import operator
import urllib.request, urllib.parse, urllib.error


'''
# maping bits:
#import matplotlib	
import matplotlib.pyplot as mplt
matplotlib.use('Agg')
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from mpl_toolkits.basemap import Basemap as Basemap
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
'''
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
	#
	def __init__(self, inData=[]):
		#self.cat=[]
		#self.subcats=[]
		self.initialize(inData)
	
	def initialize(self, inData=[]):
		# what does the data look like?
		self.cat=[]
		self.subcats=[]
		#
		self.cat=inData
		#self.sqlport=3306
		#self.sqlhost='localhost'
		inData=None
		self.catmap=None
		self.__name__='eqcatalog'
		self.mapres='l'	# basemap map resolution. at some pont, add functions to sort out nonsensical values.
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
			fout.write("%d/%d/%d %d:%d:%d.%d\t%f\t%f\t%f\n" % (rw[0].year, rw[0].month, rw[0].day, rw[0].hour, rw[0].minute, rw[0].second, rw[0].microsecond, rw[1], rw[2], rw[3]))
		fout.close()
			
	def loadCatFromFile(self, fname=None, minmag=None):
		# standard file format: date \s time \t lat \t lon \t mag
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
			if len(rws)<5: continue
			if minmag!=None:
				#if float(rws[3])<minmag: continue
				if float(rws[4])<minmag: continue
			#self.cat+=[[datetimeFromString(rws[0], float(rws[1]), float(rws[2]), float(rws[3])]]
			#if '24:' in rws[0]: print rws[0]
			#print rws[1]
			#print rws
			thisdt=datetimeFromString(rws[0] + ' '+ rws[1])
			thisdt=dtm.datetime(*thisdt.timetuple()[:-2], tzinfo=pytz.timezone('UTC'))
			self.cat+=[[thisdt, float(rws[2]), float(rws[3]), float(rws[4])]]
		fin.close()
		# add (UTC) timezone:
		self.checkdates()	# note, in the future we should check for TZ info from string. for now, assume utc.
		# and now, let's assume that we want this time-ordered:
		self.cat.sort(key = lambda x: x[0])
	#
	def checkdates(self, cat=None):
		if cat==None: cat=self.getcat(0)
		for i in range(len(cat)):
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
			newVec=rotatexy(row[2], row[1], clat, clon, theta)
			#
			# is the rotated vector in our ellipse?
			if abs(newVec[0])>ra: continue
			Y=ellipseY(newVec[0], ra, rb)
			if abs(newVec[1])>Y: continue
			# dtm, lat, lon, mag, tX, tY 		(note this is like y,x, x`, y` for the space coordinates).
			#self.subcats[-1][1]+=[[row[0], row[1], row[2], row[3], newVec[0], newVec[1]]]
			tempcat+=[[row[0], row[1], row[2], row[3], newVec[0], newVec[1]]]
		return tempcat
	
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
			if rw[0]>=dtFrom and rw[0]<=dtTo: newcat+=[rw]
		#
		return newcat
	def addTimeRangeCat(self, subcatname='dtSubcat', fullcat=None, dtFrom=None, dtTo=None):
		self.subcats+=[[subcatname, self.getTimeRangeCat(fullcat, dtFrom, dtTo)]]
	
	def getLatLonSubcat(self, fullcat, lats=[], lons=[], llcols=[1,2]):
		# llcols: lat, lon
		llrange=None
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
	
	def plotCatMap(self, catalog=None, doShow=True, doSave=False, saveName='catalogPlot.png', epicenter=None, legendLoc='upper left', doCLF=True, eqicon='b,', myaxis=None, fignum=None, padfactor=.25, plotevents=True):
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
		#llr=self.getLatLonRange(catalog)	# latLonRange
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
		
		cntr=[float(llr[0][0])+(llr[1][0]-float(llr[0][0]))/2.0, float(llr[0][1])+(llr[1][1]-float(llr[0][1]))/2.0]
		if self.catmap==None: self.catmap=Basemap(llcrnrlon=llr[0][1], llcrnrlat=llr[0][0], urcrnrlon=llr[1][1], urcrnrlat=llr[1][0], resolution=self.mapres, projection='tmerc', lon_0=cntr[1], lat_0=cntr[0])
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
		
		catmap.drawmeridians(list(range(int(llr[0][1]-2.0), int(llr[1][1]+2.0))), color='k', labels=[1,1,1,1])
		catmap.drawparallels(list(range(int(llr[0][0]-2.0), int(llr[1][0]+2.0))), color='k', labels=[1, 1, 1, 1])
		
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

	def rbomoriQuadPlot(self, catnum=0, mc=2.5, winlen=501, rbthresh=1.0, bigmag=5.0, fignum=0, intlist=None, rbavelen=1, thislw=1.0, mainEV=None, plotevents=False, mapcatnum=None, rbLegLoc='best', logZ=None):
		# make an awesome quad plot: omor-times, rbRatios, mag-seismicity, map-catalog
		# catalog -> a yodapy.eqcatalog() object
		# thislw: linewidth
		#
		if rbavelen==None:
			rbavelen=int(winlen/10.)
			if rbavelen==0: rbavelen=1
		print("avelen: %d" % rbavelen)
		#
		if mapcatnum==None: mapcatnum=catnum
		#
		# normalization factor fo rb. sequences.
		if logZ==None:
			logZ = math.log10(float(winlen))
		#
		logZexp=1.0/logZ		# so logZ is more intuitive...
		#
		#rbavelen=1	# rb-averaging length (aka, <nrb>_thisnum
		#print "winlen: %d" % winlen
		#if catalog==None: catalog=getMexicaliCat()
		#intlist=[25, 256, 512]
		if intlist==None: intlist=[int(winlen/2), winlen, winlen*2]
		catalog=self
	
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
		try:
			catalog.rb
		except:
			catalog.rb=rbi.intervalRecordBreaker(None)
		#
		#format: plotInts(self, intervals=[10, 100, 1000], catalog=None, minmag=2.0, ax=None, dtmlatlonmagCols=[0,1,2,3], plotDates=[None, None]):
		epicen=None
		if mainEV!=None: epicen=[mainEV[2], mainEV[1]]
		#
		catalog.plotInts(intervals=intlist, catalog=catalog.getcat(catnum), minmag=mc, ax=myaxes[1], thislw=thislw, legendPos='upper left')
		#myaxes[1].set_label('mean intervals $\\tau$')
		myaxes[1].set_ylabel('mean intervals $\\ < \\tau >$', size=14)
		#
		# format: #plotMags(self, catalog=None, minmag=2.0, ax=None, dtmlatlonmagCols=[0,1,2,3], plotDates=[None, None]):	
		catalog.plotMags(catalog.getcat(catnum), mc, myaxes[0])
		myaxes[0].set_ylabel('mag', size=14)
		#
		# plotIntervalRatiosAx(self, minmag=3.0, windowLen=10, cat0=None, hitThreshold=1.0, bigmag=5.0, thisAx=None, ratios=None, deltaipos=1, avlen=1, mainEV=None)
		catalog.rb.plotIntervalRatiosAx(minmag=mc, windowLen=winlen, cat0=catalog.getcat(catnum), hitThreshold=rbthresh, bigmag=bigmag, thisAx=myaxes[2], ratios=None, deltaipos=1, avlen=rbavelen, mainEV=mainEV, rbLegLoc=rbLegLoc, logZ=logZ)
		myaxes[2].set_ylabel('RB ratio $r(N=%d)$' % winlen, size=14)
	
		#plt.figure(fignum)
		myfsize=12
		myaxes[1].text(.1, .1, '$<\\tau>$', None, rotation='vertical', size=myfsize)
		myaxes[0].text(.1, .1, 'mags', None, rotation='vertical', size=myfsize)
		myaxes[2].text(.1, .1, '$r(1) = frac{n_{rb-large}}{n{rb-small}}$', None, rotation='vertical', size=myfsize)
	
		#
		# plotCatMap(self, catalog=None, doShow=True, doSave=False, saveName='catalogPlot.png', epicenter=None, legendLoc='upper left', doCLF=True, eqicon='b,', myaxis=None, fignum=None, padfactor=.25, plotevents=True)
		X=catalog.plotCatMap(catalog=catalog.getcat(mapcatnum), doShow=True, doSave=False, saveName=None, epicenter=epicen, legendLoc='best', doCLF=False, eqicon='b,', myaxis=myaxes[3], fignum=0, padfactor=.15, plotevents=plotevents)
		# and large events:
		bigEq=[]
		bigeqindex=0
		#print catalog.catmap
		for rw in catalog.getcat(catnum):
			if rw[3]>=bigmag and rw[0]>=catalog.getMainEvent(catalog.getcat(catnum))[0]:
				eqx,eqy = catalog.catmap(rw[2], rw[1])
				catalog.catmap.plot(eqx, eqy, '*', label='m=%.2f' % (rw[3]), ms=15)
				bigeqindex+=1
				print("eq: %d, %f, %s" % (bigeqindex, rw[3], str(rw[0])))
		catalog.catmap.ax.legend(loc='best', numpoints=1)
		#
		#
		plt.show()
	
		# problem: the x-axis of RB-ratios is in float form, formatted as datetimes (basically using "plot" type functions). the intervals/magnitudes are date_plot
		# functions. in short, they have different x-axis variable types, so we can't share the rb-ratios x-axis with the other two. it probably makes sense to convert
		# the interval-ratio plots to the rb-ratio style (because the rb-rations uses the "fillbetween" function).
	
		catalog.axlist=myaxes
		return catalog

	
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
			#self.setCatFromSQL(dtm.datetime(1969,1,1), dtm.datetime.now(), [34.9, 36.9], [-121.5, -119.5], 1.5, "Earthquakes", 523, 'asc')
			self.setCatFromSQL(dtm.datetime(1972,1,1, tzinfo=pytz.timezone('UTC')), dtm.datetime.now(pytz.timezone('UTC')), [34.4, 37.4], [-121.11, -119.4], 1.5, "Earthquakes", 523, 'asc')
			self.addEllipCat('PFshock (.8 x .15)', self.cat, 40.0, 35.9, -120.5, 0.8, 0.15)
			self.addEllipCat('PFshock (.4 x .15)', self.cat, 40.0, 35.9, -120.5, 0.4, 0.15)
		if catname == 'PF5yr':
			#theta=40.0, clat=35.9, clon=-120.5, ra=.4, rb=.15
			#self.setCatFromSQL(dtm.datetime(1969,1,1), dtm.datetime.now(), [34.9, 36.9], [-121.5, -119.5], 1.5, "Earthquakes", 523, 'asc')
			self.setCatFromSQL(dtm.datetime(1999, 9, 28, tzinfo=pytz.timezone('UTC')), dtm.datetime(2009,9,29,tzinfo=pytz.timezone('UTC')), [34.4, 37.4], [-121.11, -119.4], 1.5, "Earthquakes", 523, 'asc')
			self.addEllipCat('PFshock (.8 x .15)', self.cat, 40.0, 35.9, -120.5, 0.8, 0.15)
			self.addEllipCat('PFshock (.4 x .15)', self.cat, 40.0, 35.9, -120.5, 0.4, 0.15)
		if catname in ['taiwan']:
			self.setCatFromSQL(dtm.datetime(1980,1,1, tzinfo=pytz.timezone('UTC')), dtm.datetime(2010,6,1, tzinfo=pytz.timezone('UTC')), [-90, 90], [-180, 180], 2.0, 'Earthquakes', 21)
			
	def setCatFromSQL(self, startDate=dtm.datetime(1999,9,28, 17,15,24, tzinfo=pytz.timezone('UTC')), endDate=dtm.datetime(2009,9,28, 17,15,24, tzinfo=pytz.timezone('UTC')), lats=[32.0, 37.0], lons=[-125.0, -115.0], minmag=3.0, catalogName='Earthquakes', catalogID=523, ordering='asc'):
		self.cat=self.getCatFromSQL(startDate, endDate, lats, lons, minmag, catalogName, catalogID, ordering)
		return None
	
	def getCatFromSQL(self, startDate=dtm.datetime(1999,9,28, 17,15,24, tzinfo=pytz.timezone('UTC')), endDate=dtm.datetime(2009,9,28, 17,15,24, tzinfo=pytz.timezone('UTC')), lats=[32.0, 37.0], lons=[-125.0, -115.0], minmag=2.0, catalogName='Earthquakes', catalogID=523, ordering='asc'):
		# return a catalog:
		if lats[0]>lats[1]: lats.reverse()
		if lons[0]>lons[1]: lons.reverse()
		#if yp.datetimeToFloat(startDate)>yp.datetimeToFloat(endDate):
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
			#Y=range(1, len(mags)+1)
			Y=frange(1, len(mags), -1)
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

	def getIntervals(self, catList=None, winLen=1):
		if catList==None: catList=self.getcat(0)
		#
		catLen=len(catList)
		i=(catLen-1-winLen)	# start winLen positions from the end.
		thisInterval=0
		#N=1
		intervals=[]	# [[eventDateTime, totalInterval]]
		while i>=0:
			#
			thisInterval=datetimeToFloat(catList[i+winLen][0])-datetimeToFloat(catList[i][0])
			intervals+=[[catList[i+winLen][0], thisInterval]]
			i-=1
		
		#
		#return [intervals, catList]
		return intervals
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
		#dts=['1950-01-01', str(dtm.datetime.now(pytz.timezone('UTC')))]
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
		
	def plotInts(self, intervals=[10, 100, 1000], catalog=None, minmag=2.0, ax=None, dtmlatlonmagCols=[0,1,2,3], plotDates=[None, None], thislw=1.0, legendPos='best'):
		# plot mean intervals only. having figured out how to put a bunch of independently generated plots onto a single canvas, we
		# split up some of these plots...
		#
		if type(plotDates).__name__!='list': plotDates=[None, None]
		while len(plotDates)<2: plotDates+=[None]
		#
		if catalog==None: catalog=self.cat
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
		#
		if catalog==None: catalog=self.cat
		cols=dtmlatlonmagCols	# for efficient notation
		#
		mags=[]
		activecat=[]
		for rw in catalog:
			if rw[cols[3]]>=minmag: mags+=[[rw[cols[0]], rw[cols[3]]]]
		mags=vlinePadList(mags, minmag-abs(minmag)*.1)	# return the mags data padded for vertical line style plotting. this is just a trick to get width=1 histograms.
		#
		ax.plot(list(map(operator.itemgetter(0), mags)), list(map(operator.itemgetter(1), mags)), '-')
		ax.legend(loc='best')
		return mags
		
	def getANSStoFilehandler(self, lon=None, lat=None, minMag=4.92, dates0=[dtm.date(2001,0o1,0o1), dtm.date(2010, 12, 31)], keywds='', etype='E', depths=[0, 200], Nmax=999999):
		# fetch data from ANSS; return a file handler.
		#
		# use urllib in "post" mode. an example from http://www.python.org/doc/current/library/urllib.html#urllib.FancyURLopener)
		# using "get" (aka, query-string method; note the ?%s string at the end of the URL, this is a single pram call to .urlopen):
		#
		#>>> import urllib
		#>>> params = urllib.urlencode({'spam': 1, 'eggs': 2, 'bacon': 0})
		#>>> f = urllib.urlopen("http://www.musi-cal.com/cgi-bin/query?%s" % params)
		#>>> print f.read()
		#
		# using "post" (note this is a 2 pram call):
		#>>> import urllib
		#>>> params = urllib.urlencode({'spam': 1, 'eggs': 2, 'bacon': 0})
		#>>> f = urllib.urlopen("http://www.musi-cal.com/cgi-bin/query", params)
		#>>> print f.read()
		#
		# make ANSS prams dictionary (thank james for the bash-template):
		# ANSSquery has day-resolution:
		dates=dates0
		if type(dates[0]).__name__=='datetime':
			dates=[dtm.date(dates0[0].year, dates0[0].month, dates0[0].day), dtm.date(dates0[1].year, dates0[1].month, dates0[1].day)]
		
		#anssPrams={'format':'cnss', 'output':'readable', 'mintime':str(dates[0]).replace('-', '/'), 'maxtime':str(dates[1]).replace('-', '/'), 'minmag':str(minMag), 'minlat':lat[0], 'maxlat':lat[1], 'minlon':lon[0], 'maxlon':lon[1], 'etype':'E', 'searchlimit':Nmax}
		anssPrams={'format':'cnss', 'output':'readable', 'mintime':str(dates[0]).replace('-', '/'), 'maxtime':str(dates[1]).replace('-', '/'), 'etype':etype, 'minmag':str(minMag), 'mindepth':depths[0], 'maxdepth':depths[1], 'searchlimit':Nmax}
		#  'minlat':lat[0], 'maxlat':lat[1], 'minlon':lon[0], 'maxlon':lon[1],
		if lat[0]!=None and lat[0]!='' and lat[1]!=None and lat[1]!='':
			anssPrams['minlat']=lat[0]
			anssPrams['maxlat']=lat[1]
		if lon[0]!=None and lon[0]!='' and lon[1]!=None and lon[1]!='':
			anssPrams['minlon']=lon[0]
			anssPrams['maxlon']=lon[1]
		if keywds!=None and keywds!='':
			#anssPrams+={"keywds":keywds}
			anssPrams['keywds']=keywds
		# 
		f = urllib.request.urlopen('http://www.ncedc.org/cgi-bin/catalog-search2.pl', urllib.parse.urlencode(anssPrams))
		#
		# we might return f, a string of f, or maybe a list of lines from f. we'll work that out shortly...
		return f

		#def getANSSlist(lon=[-125, -115], lat=[32, 45], minMag=4.92, dates0=[dtm.date(2001,01,01), dtm.date(2010, 12, 31)], Nmax=999999, fin=None):
	def getANSSlist(self, lon=[-125, -115], lat=[32, 45], minMag=4.92, dates0=[dtm.datetime(2001,0o1,0o1, tzinfo=pytz.timezone('UTC')), dtm.datetime(2010, 12, 31, tzinfo=pytz.timezone('UTC'))], keywds='', etype='E', depths=[0, 200], Nmax=999999, fin=None):
		#
		# note: this appears to be a bad idea for global downloads. a full catalog is ~4GB, which kills my computer.
		#
		# note: this may be repeated exactly in ygmapbits.py
		# fetch new ANSS data; return a python list object of the data.
		# fin: data file handler. if this is None, then get one from ANSS.
		#dates=[dtm.date(dates0[0].year, dates0[0].month, dates0[0].day), dtm.date(dates0[1].year, dates0[1].month, dates0[1].day)]
		dates=dates0
		if fin==None:
			#print "get data from ANSS...(%s, %s, %s, %s, %s)" % (lon, lat, minMag, dates, Nmax)
			#fin = getANSStoFilehandler(lon, lat, minMag, dates, Nmax)
			fin = self.getANSStoFilehandler(lon, lat, minMag, dates, keywds, etype, depths, Nmax)
			#fin = getANSStoFilehandler([-180, 180], [-90, 90], 0, [datetime.date(1910,01,01), datetime.date(2010, 01, 16)], 9999999)
			#print "data handle fetched..."
		return self.getANSSlistFile(fin)
	
	def getANSSlistFile(self, fin):
		anssList=[]
		for rw in fin:
			if rw[0] in ["#", "<"] or rw[0:4] in ["Date", "date", "DATE", "----"]:
				#print "skip a row... %s " % rw[0:10]
				continue
			#anssList+=[rw[:-1]]
			# data are fixed width delimited
			# return date-time, lat, lon, depth, mag, magType, nst, gap, clo, rms, src, catEventID (because those are all the available bits)
			#print "skip a row... %s " % rw
			rwEvdt=rw[0:22].strip()
			rwLat=rw[23:31].strip()
			if rwLat=='' or isnumeric(str(rwLat))==False or rwLat==None:
				continue
				#rwLat=0.0
			else:
				rwLat=float(rwLat)
			rwLon=rw[32:41].strip()
			if rwLon=='' or isnumeric(str(rwLon))==False or rwLon==None:
				#rwLon=0.0
				continue
			else:
				rwLon=float(rwLon)
			rwDepth=rw[42:48].strip()
			if rwDepth=='' or isnumeric(str(rwDepth))==False or rwDepth==None or str(rwDepth).upper() in ['NONE', 'NULL']:
				#rwDepth=0.0
				rwDepth=None
				continue
			else:
				rwDepth=float(rwDepth)
			rwMag=rw[49:54].strip()
			if rwMag=='' or isnumeric(str(rwMag))==False or rwMag==None:
				#rwMag=0.0
				continue
			else:
				rwMag=float(rwMag)
			#rwMagType=rw[55:59].strip()
			#rwNst=rw[60:64].strip()
			#if rwNst=='':
			#	rwNst=0.0
			#else:
			#	rwNst=float(rwNst)
			#rwGap=rw[65:68].strip()
			#rwClo=rw[69:73].strip()
			#rwrms=rw[74:78].strip()
			#if rwrms=='':
			#	rwrms=0.0
			#else:
			#	rwrms=float(rwrms)		
			#rwsrc=rw[79:83].strip()
			#rwCatEventId=rw[84:96].strip()
		
			#anssList+=[[rw[0:22].strip(), float(rw[23:31].strip()), float(rw[32:41].strip()), float(rw[42:48].strip()), float(rw[49:54].strip()), rw[55:59].strip(), float(rw[60:64].strip()), rw[65:68].strip(), rw[69:73].strip(), float(rw[74:78].strip()), rw[79:83].strip(), rw[84:96].strip()]]
			#anssList+=[[rwEvdt, rwLat, rwLon, rwDepth, rwMag, rwMagType, rwNst, rwGap, rwClo, rwrms, rwsrc, rwCatEventId]]
			# ['2009/10/16 10:03:39.83',
			dtmEvDt=dtm.datetime.strptime(rwEvdt, '%Y/%m/%d %H:%M:%S.%f')
			dtmEvDt=dtm.datetime(*dtmEvDt.timetuple()[:-2], tzinfo=pytz.timezone('UTC'))
			anssList+=[[dtmEvDt, rwLat, rwLon, rwMag, rwDepth]]
		return anssList
	
	def getANSSlistFileFullrow(self, fin):
		anssList=[]
		for rw in fin:
			if rw[0] in ["#", "<"] or rw[0:4] in ["Date", "date", "DATE", "----"]:
				#print "skip a row... %s " % rw[0:10]
				continue
			#anssList+=[rw[:-1]]
			# data are fixed width delimited
			# return date-time, lat, lon, depth, mag, magType, nst, gap, clo, rms, src, catEventID (because those are all the available bits)
			#print "skip a row... %s " % rw
			rwEvdt=rw[0:22].strip()
			rwLat=rw[23:31].strip()
			if rwLat=='' or isnumeric(str(rwLat))==False or rwLat==None:
				continue
				#rwLat=0.0
			else:
				rwLat=float(rwLat)
			rwLon=rw[32:41].strip()
			if rwLon=='' or isnumeric(str(rwLon))==False or rwLon==None:
				#rwLon=0.0
				continue
			else:
				rwLon=float(rwLon)
			rwDepth=rw[42:48].strip()
			if rwDepth=='' or isnumeric(str(rwDepth))==False or rwDepth==None or str(rwDepth).upper() in ['NONE', 'NULL']:
				#rwDepth=0.0
				rwDepth=None
				continue
			else:
				rwDepth=float(rwDepth)
			rwMag=rw[49:54].strip()
			if rwMag=='' or isnumeric(str(rwMag))==False or rwMag==None:
				#rwMag=0.0
				continue
			else:
				rwMag=float(rwMag)
			rwMagType=rw[55:59].strip()
			rwNst=rw[60:64].strip()
			if rwNst=='':
				rwNst=0.0
			else:
				rwNst=float(rwNst)
			rwGap=rw[65:68].strip()
			rwClo=rw[69:73].strip()
			rwrms=rw[74:78].strip()
			if rwrms=='':
				rwrms=0.0
			else:
				rwrms=float(rwrms)		
			rwsrc=rw[79:83].strip()
			rwCatEventId=rw[84:96].strip()
		
			#anssList+=[[rw[0:22].strip(), float(rw[23:31].strip()), float(rw[32:41].strip()), float(rw[42:48].strip()), float(rw[49:54].strip()), rw[55:59].strip(), float(rw[60:64].strip()), rw[65:68].strip(), rw[69:73].strip(), float(rw[74:78].strip()), rw[79:83].strip(), rw[84:96].strip()]]
			anssList+=[[rwEvdt, rwLat, rwLon, rwMag, rwDepth, rwMagType, rwNst, rwGap, rwClo, rwrms, rwsrc, rwCatEventId]]
		return anssList
			#####
			# end ANSS query functions
		#

			
def helloWorld(mystr='hello world.'):
	# use this as a module-load diagnostic.
	return mystr	

def isnumeric(value):
  return str(value).replace(".", "").replace("-", "").isdigit()
	
