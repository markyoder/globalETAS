from math import *	
from scipy import *
import scipy
from pylab import *
from matplotlib import *
import numpy.fft as nft
import scipy.optimize as spo
#from matplotlib import pyplot as plt
import pylab as plt
from matplotlib import rc
import numpy


import string
import sys
#from matplotlib import *
#from pylab import *
import os
import random
#
# gamma function lives here:
#import scipy.special
from scipy.special import gamma
#from scipy.optimize import leastsq
from matplotlib import axis as aa
from threading import Thread
#
import datetime as dtm
import time
import pytz
import calendar
import operator
import urllib.request, urllib.parse, urllib.error
import MySQLdb

import rbIntervals as rbi

from eqcatalog import *
# 20120719 yoder: (this might have impacts on apps. that auto-fetch from ANSS)
from ANSStools import *

# maping bits:
import matplotlib	# note that we've tome from ... import *. we should probably eventually get rid of that and use the matplotlib namespace.
#matplotlib.use('Agg')
#from matplotlib.toolkits.basemap import Basemap
from mpl_toolkits.basemap import Basemap as Basemap
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
#
tzutc = pytz.timezone('UTC')

################################################
# Utilities
################################################

class linefit:
	# a simple line-fit class. we'll need X,Y,wt,
	# call like: lf1=m.linefit(dataSet)
	#
	#
	datas=[]		# [X, Y, Wt]
	totalVar=0			# sum( (x_i - f(x_i)**2 )
	errs=[]
	Ndof=0
	a=None
	b=None
	activeRes=None
	AIC=None
	nPrams=None	# len(fitData)=Ndof+nPrams
	
	#
	def __init__(self, inData=[]):
		self.initialize(inData)
	
	def initialize(self, inData=[]):
		# what does the data look like?
		self.activeRes=self.linRes
		dLen=len(inData)		# if dLen==2,3, then we probably have [[X], [Y], [wt]]
		#
		if dLen<=3:
			# we should probably sort the arrays...
			self.datas=inData
		if dLen>3:
			# data is probably like [[x,y,w], [x,y,w]...]
			if len(inData[0])==1:
				# 1D array; assume a plain ole sequence for X
				self.datas=[]
				self.datas+=[list(range(len(inData)))]
				self.datas+=[list(map(operator.itemgetter(0), inData))]
			if len(inData[0])>=2:
				# assume the data are ordered pairs, so we can sort them on the x coordinate.
				inData.sort(key=operator.itemgetter(0))
				self.datas=[]
				self.datas+=[list(map(operator.itemgetter(0), inData))]
				self.datas+=[list(map(operator.itemgetter(1), inData))]
			if len(inData[0])>=3: self.datas+=[list(map(operator.itemgetter(2), inData))]
		if len(self.datas)==2:
			# add even weight.
			self.datas+=[[]]
			for x in self.datas[0]:
				self.datas[2]+=[1]
		#
	#
	def meanVar(self):
		# aka, rChiSqr
		return self.totalVar/self.Ndof
		
	#def doOmoriFit(self, p0=[None, None, None], xmin=None, xmax=None, r0=0):
	def doOmoriFit(self, p0=[None, None, None], xmin=None, xmax=None, r0=0, fitres=None):
		# the pram r0 is for the OmoriIntRes2() residual, which incorporates backtround seismicity.
		# do a linear fit:
		# a,b parameters are first guesses.
		# if they are None, guess starting prams:
		if fitres==None: fitres=self.omoriIntRes
		if p0==None: p0=[None, None, None]
		if p0[1]==None:
			p0[1]=(float(self.datas[1][-1])-float(self.datas[1][0]))/(float(self.datas[0][-1])-float(self.datas[0][0]))
		if p0[0]==None:
			p0[0]=self.datas[1][0]-p0[1]*self.datas[0][0]
		if p0[2]==None: p0[2]=1.0
		#self.a=a
		#self.b=b
		p=scipy.array(p0)
		#
		if xmin==None: xmin=self.datas[0][0]
		if xmax==None: xmax=self.datas[0][-1]
		# 
		# get X,Y,W:
		X=[]
		Y=[]
		W=[]
		for i in range(len(self.datas[0])):
			if self.datas[0][i]<xmin: continue
			#
			X+=[self.datas[0][i]]
			Y+=[self.datas[1][i]]
			W+=[self.datas[2][i]]
			if X[-1]>=xmax: break

		#
		# now, fit data:
		
		#print "do the fit..."
		# note: args are (y, x, wt)
		#plsq=spo.leastsq(self.linRes, p, args=(scipy.array(self.datas[1]), scipy.array(self.datas[0]), scipy.array(self.datas[2])), full_output=1)
		print("prams: %s" % str(p))
		#plsq=spo.leastsq(fitres, p, args=(scipy.array(Y), scipy.array(X), scipy.array(W), r0), full_output=1)
		plsq=spo.leastsq(fitres, p, args=(scipy.array(Y), scipy.array(X), scipy.array(W)), full_output=1)
		#print "fit done. sum error..."
		for sig in self.errs:
			self.totalVar+=sig*sig
		#self.Ndof=len(self.datas[0])-len(p)
		self.Ndof=len(X)-len(p)
		self.nPrams=len(p)
		self.AIC=self.totalVar+2*nPrams
		self.a=plsq[0][0]
		self.b=plsq[0][1]
		'''
		plsq=spo.leastsq(linRes, p, args=(ymax,x), full_output=0, maxfev=200000)
		amax=plsq[0][0]
		bmax=plsq[0][1]
		'''
		return plsq
	
	def doLogFit(self, lbase=10.0, a=None, b=None, xmin=None, xmax=None, thisdatas=None):
		# fit the log-log representation.
		if thisdatas==None: thisdatas=self.datas
		#print "datalen: %d" % len(thisdatas)
		logdatas=[[], [], []]
		#
		# get logarithms of data:
		for i in range(len(thisdatas[0])):
			logdatas[0]+=[math.log(thisdatas[0][i], lbase)]
			logdatas[1]+=[math.log(thisdatas[1][i], lbase)]
			wt=1
			if len(thisdatas)>=3: wt=thisdatas[2][i]
			logdatas[2]+=[wt]				
		#
		#return logdatas
		
		return self.doFit(a, b, xmin, xmax, logdatas)
		
	def doFit(self, a=None, b=None, xmin=None, xmax=None, thisdatas=None, fop=1):
		if thisdatas==None: thisdatas=self.datas
		# do a linear fit:
		# a,b parameters are first guesses.
		# if they are None, guess starting prams:
		self.errs=[]
		self.totalVar=0
		if b==None:
			b=(float(thisdatas[1][-1])-float(thisdatas[1][0]))/(float(thisdatas[0][-1])-float(thisdatas[0][0]))
		if a==None:
			a=thisdatas[1][0]-b*thisdatas[0][0]
		self.a=a
		self.b=b
		p=scipy.array([a,b])
		
		if xmin==None: xmin=thisdatas[0][0]
		if xmax==None: xmax=thisdatas[0][-1]
		# 
		# get X,Y,W:
		X=[]
		Y=[]
		W=[]
		for i in range(len(thisdatas[0])):
			if thisdatas[0][i]<xmin: continue
			#
			X+=[thisdatas[0][i]]
			Y+=[thisdatas[1][i]]
			W+=[thisdatas[2][i]]
			if X[-1]>=xmax: break
		#
		# now, fit data:
		
		#print "do the fit..."
		# note: args are (y, x, wt)
		#plsq=spo.leastsq(self.linRes, p, args=(scipy.array(self.datas[1]), scipy.array(self.datas[0]), scipy.array(self.datas[2])), full_output=1)
		plsq=spo.leastsq(self.activeRes, p, args=(scipy.array(Y), scipy.array(X), scipy.array(W)), full_output=fop)
		#print "fit done. sum error..."
		for sig in self.errs:
			self.totalVar+=sig*sig
		#self.Ndof=len(thisdatas[0])-len(p)
		self.AIC=self.totalVar+2.*len(p)
		self.Ndof=len(X)-len(p)
		self.dataLen=len(X)
		self.a=plsq[0][0]
		self.b=plsq[0][1]
		'''
		plsq=spo.leastsq(linRes, p, args=(ymax,x), full_output=0, maxfev=200000)
		amax=plsq[0][0]
		bmax=plsq[0][1]
		'''
		#
		# now, get error for a,b.
		# start by calculating Delta and DeltaPrime (error-partition for w!=1, w=1):
		self.delta=0
		self.deltaPrime=0
		aX=scipy.array(X)
		aXX=aX**2
		aY=scipy.array(Y)
		aYY=aY**2
		aW=scipy.array(W)
		aWW=aW**2
		
		#sigsSq=scipy.array(self.errs)**2	# i think this does not get squared (we're using linRes() )... it would not be a bad idea to just calc them here...
		#self.delta=sum
		# delta first (this might not be quite right except when w=1/sig^2, so be careful...)
		# note: we assume W=1/sigma**2
		#self.delta=sum(aW)*sum(aXX*aW) - (sum(aX*aW))**2	# check this formulation for order of operations...
		self.delta=sum(aWW)*sum(aXX*aWW) - (sum(aX*aWW))**2.	# check this formulation for order of operations...
		self.deltaPrime=float(len(X))*sum(aXX) - sum(aX)**2.
		#
		# weighted pram errors (variance):
		#thisSig=self.totalVar**.5
		#
		self.vara=(1.0/self.delta)*sum(aXX*aWW)
		self.varb=(1.0/self.delta)*sum(aWW)	# note that w=1/sig^2 in most texts, as per gaussian stats. this is more general than that, and it might make a big fat
		#										# mess when applied outside gaussian distributions.
		# w_i=1 case (note, this can be generallized to w_i=w; 1/delta -> 1/(w*delta):
		self.varaprime=(self.meanVar()/self.deltaPrime)*sum(aXX)
		self.varbprime=float(len(X))*self.meanVar()/self.deltaPrime
			
		
		return plsq		

	def tofile(self, fname='data/lfdata.dat', lfheader='#data from linefit object\n'):
		fout=open(fname, 'w')
		fout.write(lfheader)
		for i in range(len(self.datas[0])):
			fout.write('%f\t%f\t%f\n' % (self.datas[0][i], self.datas[1][i], self.datas[2][i]))
		fout.close()
	
	def fLin(self, x,p):
		return p[0] + x*p[1]
		
	def fPL(self, x,p):
		return (10**p[0])*(x**p[1])
	
		
	def linRes(self, p,y,x,w):
		err=y-(p[0]+x*p[1])
		self.errs=w*err
		#return w*err
		return self.errs
	
	def omoriRateRes(self, p, y, x, w):
		err=y-(1.0/(p[0]*(1+(x/p[1]))**p[2]))
		self.errs=w*err
		#return w*err
		return self.errs
	
	def omoriIntRes(self, p, y, x, w):
		err=y-(p[0]*(1+(x/p[1]))**p[2])
		werr=w*err
		self.errs=w*err
		#return werr
		return self.errs
	
	def omoriIntRes2(self, p, y, x, w, r0):
		# this is an omori like function that includes a "background rate" r0.
		#err=w*(y-(p[0]*(1+(x/p[1]))**p[2]))
		#r0=0.2222
		#err=w*(y- (1/(r0 + 1/(p[0]*(1+(x/p[1]))**p[2]) )) )
		err=w*(y- 1/(r0 + 1/(p[0]*(1+(x/p[1]))**p[2])) )
		self.errs=err
		return err
	
	def getFitPlotAry(self):
		print("from getFitPlotAry(): %s, %s, %s, %s" % (self.datas[0][0], self.datas[0][-1], self.a, self.b))
		Fx=[self.datas[0][0], self.datas[0][-1]]
		Fy=[self.datas[0][0]*self.b + self.a, self.datas[0][-1]*self.b + self.a]
		return [Fx, Fy]
		
	def quickPlot(self, toFig=True, colorNum=0):
		fitXY=self.getFitPlotAry()
		# toFig: do a figure here or assume that there is an active figure in the ether. this way, we can put many quickPlots onto a single canvas.
		
		if toFig:
			plt.figure(0)
			plt.clf()
		plt.plot(self.datas[0], self.datas[1], color=pycolor(colorNum))
		plt.plot(fitXY[0], fitXY[1], color=pycolor(colorNum), label='a=%f, b=%f, rChi^2=%f' % (self.a, self.b, self.meanVar()**.5))
		if toFig:
			plt.legend(loc='upper right')
			plt.show()
		return None
			
# utils:
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
		
def pycolor(num=0):
	if num==None: num=0
	num=int(num%7)
	clrs=['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
	
	return clrs[num]
def pyicon(num=0):
	if num==None: num=0
	num=int(num%20)
	# note: skip ',', the "pixel" marker
	icns=['.', 'o', 'v', '^', '<', '>', '1', '2', '3', '4', 's', 'p', '*', 'h', 'H', '+', 'x', 'D', 'd', '|', '_']
	return icns[num]
	
def getMonthName(monthNum=1):
	if monthNum==1:
		return "January"
	if monthNum==2:
		return "February"
	if monthNum==3:
		return 'March'
	if monthNum==4:
		return 'April'
	if monthNum==5:
		return 'May'
	if monthNum==6:
		return 'June'
	if monthNum==7:
		return 'July'
	if monthNum==8:
		return 'August'
	if monthNum==9:
		return 'September'
	if monthNum==10:
		return 'October'
	if monthNum==11:
		return 'November'
	if monthNum==12:
		return 'December'

def plotPolygons(polygons=None):
	import ygmapbits as yg
	# getPIsquareRays
	if polygons==None: polygons=yg.getReducedPolys()
	#if polygons==None: polygons=getPIsquareRays()
	#
	plt.figure(0)
	for ply in polygons:
		#if len(ply)<5: continue
		plt.fill(list(map(operator.itemgetter(1), ply)), list(map(operator.itemgetter(0),ply)), '.-')
	
	plt.show()
	
	return polygons
	
def printPolyLens(polygons=None):
	import ygmapbits as yg
	if polygons==None: polygons=yg.getReducedPolys()
	#
	i=0
	for ply in polygons:
		print("poly(%d): %d" % (i, len(ply)))
		i+=1
	return None

##################
# write an all-python method to fetch ANSS data and insert into MySQL
# this will eventually replace {NEICANSS2sql.py}.fetchAndInsertANSS()

#
def replaceANSS2SQL(anssList=None, catID=523, tempCatID=None):
	# just a familiar-name wrapper:
	return ANSSlist2SQL(anssList, catID, tempCatID)
	
def ANSSlist2SQL(anssList=None, catID=523, tempCatID=None):
	# this is the end-product; call it and update the WHOLE CATALOG. BUT, this replaces the whole catalog, so it takes forever.
	# we need an update system...
	#
	# insert an ANSSlist into MySQL as catalogID={catID}. use tempCatID to temporarily dump existing data to a safe place,
	# should something go wrong. reserve an option to skip this step for brevity.
	#
	sqlHost = 'localhost'
	sqlUser = 'myoder'
	sqlPassword = 'yoda'
	sqlPort = 3306
	sqlDB = 'QuakeData'
	myConn = MySQLdb.connect(host=sqlHost, user=sqlUser, passwd=sqlPassword, port=sqlPort, db=sqlDB)
	myConn2 = MySQLdb.connect(host=sqlHost, user=sqlUser, passwd=sqlPassword, port=sqlPort, db=sqlDB)
	strDel=""
	# if no list has been provided, get a current, complete world catalog:
	#if anssList==None: anssList=getANSSlist([-180, 180], [-90, 90], 0, [dtm.date(1932,01,01), dtm.date.fromordinal(dtm.datetime.now().toordinal())], 9999999)
	if anssList==None: anssList=getANSSlist([-180., 180.], [-90., 90.], 0, [dtm.datetime(1932,0o1,0o1, tzinfo=pytz.timezone('UTC')), dtm.datetime.now(pytz.timezone('UTC'))], 9999999)
	print("ANSS list retrieved, len: %d" % len(anssList))
	#
	if tempCatID==None:
		# get a save catID; use that one.
		#strGetMaxID="select min(catalogID) from Earthquakes where catalogID>%d" % catID
		strGetMaxID="select max(catalogID) from Earthquakes"
		myConn.query(strGetMaxID)
		rID=myConn.store_result()
		# there should only be one value. for now, assume this is true (which is, in general, really bad DB programming form):
		tempCatID=int(rID.fetch_row()[0][0])+1
		#rID=None
		strDel='delete from Earthquakes where catalogID=%d; update Earthquakes set catalogID=%d where catalogID=%d' % (tempCatID, tempCatID, catID)
		
	elif tempCatID==0:
		# skip this step; just delete the catalog and hope all goes well.
		strDel='delete from Earthquakes where catalogID=%d' % catID		
	#
	else:
		# do we have a real backup catalogID? update the "real" id with the tempCatID value.
		strDel='update Earthquakes set catalogID=%d where catalogID=%d' % (tempCatID, catID)
		
	myConn.query(strDel)
	
	#
	print("begin sql insert loop...")
	for rw in anssList:
		# write an insert string or command or whatever...
		# one string, or a bunch? i don't think it matters for inserts.
		# ... and this whole thing hast to be error handled (None, Null, quotes, etc.)...
		#print "rw: %s" % rw
		#insStr = "insert into QuakeData.Earthquakes (catalogID, eventDateTime, lat, lon, depth, mag, magType, nst, gap, clo, rms, src, catEventID) values (%d, '%s',%s, %s,%s, %s,%s, %s,%s, '%s',%s, '%s', '%s') " % (catID, str(rw[0]), str(rw[1]), str(rw[2]), str(rw[3]), str(rw[4]), str(rw[5]), str(rw[6]), str(rw[7]), str(rw[8]), str(rw[9]), str(rw[10]), str(rw[11]) )
		insStr = "insert into QuakeData.Earthquakes (catalogID, eventDateTime, lat, lon, depth, mag, magType, src, catEventID) values (%d, '%s', %f, %f, %f, %f, '%s', '%s', '%s')" % (catID, str(rw[0]), float(rw[1]), float(rw[2]), float(rw[3]), float(rw[4]), str(rw[5]), str(rw[10]), str(rw[11]) )
		#print insStr
		# it might be faster to bundle these statements. also, it might not be a bad idea to remove unique constraints from the Earthquakes table to facilitate faster
		# inserts. nominally, actual inserts could be done on threads to make it super speedy.
		myConn2.query(insStr)
	
	myConn.close()
	myConn2.close()


def updateANSS2SQL(catID=523):
	# this is the end-product. call this to update MySQL.QuakeData.Earthquakes with the most recent events; it will call other functions as necessary.
	# this function assumes all existing data are correct and complete.
	# here, we get the date of the most recent event and select from ANSS all events after that.
	#
	# also, ANSS allows only day-level queries, so we have to delete everything date>date0, then replace with all events date>date0
	#
	sqlHost = 'localhost'
	sqlUser = 'myoder'
	sqlPassword = 'yoda'
	sqlPort = 3306
	sqlDB = 'QuakeData'
	myConn = MySQLdb.connect(host=sqlHost, user=sqlUser, passwd=sqlPassword, port=sqlPort, db=sqlDB)
	myConn2 = MySQLdb.connect(host=sqlHost, user=sqlUser, passwd=sqlPassword, port=sqlPort, db=sqlDB)
	strDel=""
	# if no list has been provided, get a current, complete world catalog:
	#if anssList==None: anssList=getANSSlist([-180, 180], [-90, 90], 0, [dtm.date(1932,01,01), datetime.date.fromordinal(datetime.datetime.now(tzutc).toordinal())], 9999999)
	anssList=[]
	# get a maximum date:
	#sqlMaxDate="select max(eventDateTime) from Earthquakes where catalogID=%d" % catID
	maxdates=myConn.cursor()
	maxdates.execute("select max(eventDateTime) from Earthquakes where catalogID=%d" % catID)
	maxDt=maxdates.fetchone()[0]
	print("getList prams: (%s, %s, %s, %s, %s)" % ([-180, 180], [-90, 90], 0, [maxDt, dtm.date.fromordinal(dtm.datetime.now(tzutc).toordinal())], 9999999))
	anssList=getANSSlist([-180, 180], [-90, 90], 0, [maxDt, dtm.datetime.now(tzutc)], 9999999)
	
	# now, delete most recent days events (they will be replaced):
	strDel = "delete from Earthquakes where eventDateTime>='%s' and catalogID=%d" % (str(dtm.date(maxDt.year, maxDt.month, maxDt.day)), catID)
	myConn.query(strDel)
	for rw in anssList:
		# write an insert string or command or whatever...
		# one string, or a bunch? i don't think it matters for inserts.
		# ... and this whole thing hast to be error handled (None, Null, quotes, etc.)...
		#print "rw: %s" % rw
		#insStr = "insert into QuakeData.Earthquakes (catalogID, eventDateTime, lat, lon, depth, mag, magType, nst, gap, clo, rms, src, catEventID) values (%d, '%s',%s, %s,%s, %s,%s, %s,%s, '%s',%s, '%s', '%s') " % (catID, str(rw[0]), str(rw[1]), str(rw[2]), str(rw[3]), str(rw[4]), str(rw[5]), str(rw[6]), str(rw[7]), str(rw[8]), str(rw[9]), str(rw[10]), str(rw[11]) )
		insStr = "insert into QuakeData.Earthquakes (catalogID, eventDateTime, lat, lon, depth, mag, magType, src, catEventID) values (%d, '%s', %f, %f, %f, %f, '%s', '%s', '%s')" % (catID, str(rw[0]), float(rw[1]), float(rw[2]), float(rw[3]), float(rw[4]), str(rw[5]), str(rw[10]), str(rw[11]) )
		#print insStr
		# it might be faster to bundle these statements. also, it might not be a bad idea to remove unique constraints from the Earthquakes table to facilitate faster
		# inserts. nominally, actual inserts could be done on threads to make it super speedy.
		myConn2.query(insStr)
	#
	myConn.close()
	myConn2.close()
	
	return anssList

###############

def isnumeric(value):
  return str(value).replace(".", "").replace("-", "").isdigit()

def getValsAbove(inList, aboveVal=0):
	# input a 1D list. return: the value if it's above aboveVal, otherwise aboveVal:
	retList=[]
	for x in inList:
		if x>=aboveVal:
			retList+=[x]
		else:
			retList+=[aboveVal]
	return retList

def getValsBelow(inList, belowVal=0):
	# input a 1D list. return: the value if it's above aboveVal, otherwise aboveVal:
	retList=[]
	for x in inList:
		if x<=belowVal:
			retList+=[x]
		else:
			retList+=[belowVal]
	return retList

def frange(xmax=10.0, x0=0.0, dx=1.0):
	if dx==0: dx=1.0
	if (xmax<x0 and dx>0) or (xmax>x0 and dx<0): dx=-dx
	X=[float(x0)]
	while (X[-1]<xmax and xmax>x0) or (X[-1]>xmax and xmax<x0):
		X+=[float(X[-1]+dx)]
	return X

def datetimeToFloat(dtIn=dtm.datetime.now(tzutc)):
	# returns a float version in units of days.
	# and this too can be replace by matplotlib.dates.date2num()
	fdt= dtIn.toordinal() + float(dtIn.hour)/24.0 + dtIn.minute/(24.0*60.0) + dtIn.second/(24.0*3600.0) + dtIn.microsecond/(24.0*3600000000.0)
	return fdt

def timeDeltaToFloat(tdelta=None, timeunit='day'):
	if tdelta==None: return None	# or maybe some default val.
	# time deltas are (days, secs., microsecs.)
	# let's, by default, convert to days. when we add other time units, we'll just multiply or divide accordingly.
	#
	if timeunit.lower() in ['day', 'days', 'dy', 'dys']:
		#
		rval=tdelta.days + (tdelta.seconds + (float(tdelta.microseconds)/10**6)/(3600*24.0) )
	else:
		rval=timeDeltaToFloat(tdelta, 'day')
	#
	return rval

#def datetimeFromFloat(fin=1000.0):
#	dt=dtm.datetime.fromordinal(int(fin))
#	dayfrac=fin%1
#	hrs=int(round(dayfrac*24))
#	mins=int(round(dayfrac*24*60/24))
#	secs=int(round(dayfrac*24*3600/(24*60)))
#	msecs=int(round(dayfrac*24*3600%1)*10**6)
#	print hrs, mins, secs, msecs
#	#
#	return dtm.datetime(dt.year, dt.month, dt.day, hrs, mins, secs, msecs)

def deg2rad(theta):
	return 2.0*pi*theta/360.0
	
def ellipseY(x, a, b):
	#print b, x, a
	return b*(1.0-x*x/(a*a))**.5
	
def rotatexy(x, y, Lat, Lon, theta):
	# x,y to transform via blah, blah.
	#
	theta=deg2rad(float(theta))
	xprime = (x-Lon)*cos(theta) - (y-Lat)*sin(theta)
	yprime = (x-Lon)*sin(theta) + (y-Lat)*cos(theta)
	
	return [xprime, yprime]

def datetimeFromString(strDtin=dtm.datetime.now(tzutc).isoformat(' '), delim='-'):
	# note: this can probably be replaced by a matplotlib.dates function
	possibleDelims=['/', '-']
	if delim not in strDtin:
		for dlm in possibleDelims:
			if dlm in strDtin: delim=dlm
			
	#print strDtin
	strDt=strDtin.split(' ')[0]
	#print strDt
	if ' ' in strDtin and ':' in strDtin:
		strTm=strDtin.split(' ')[1]
	else:
		strTm="0:0:0.00"
	if '.' not in strTm: strTm=strTm+'.00'
	#
	dts=strDt.split(delim)
	tms=strTm.split(':')
	
	secFract='.'+tms[2].split('.')[1]
	if secFract=='':
		secFract=0.0
	else:
		secFract=float(secFract)
	
	yr=int(dts[0])
	mnth=int(dts[1])
	dy=int(dts[2])
	#
	hrs=int(tms[0])
	mns=int(tms[1])
	secs=tms[2].split('.')[0]
	if secs=='': secs=0
	secs=int(secs)
	msecs=int(secFract*10**6)
	#
	# and we keep seeing seconds>60...
	if secs>=60:
		secs=secs%60
		mns+=secs/60
	if mns>=60:
		mns=mns%60
		hrs+=mns/60
	retDate=dtm.datetime(yr, mnth, dy, hrs, mns, secs, msecs, tzutc)
	if hrs>=24:
		hrs=hrs%24
		retDate+=dtm.timedelta(days=hrs/24)
	#
	
	
	#print yr, mnth, dy, hrs, mns, secs, msecs
	
	#return None
	#return dtm.datetime(yr, mnth, dy, hrs, mns, secs, msecs, tzutc)
	return retDate
	


def greaterof(a,b):
	if a>=b: return a
	if b>a: return b

def lesserof(a,b):
	if a<=b: return a
	if b<a: return b	

def vlinePadList(lst, minVal=0):
	# pad list vals so they plot as vertical spikes.
	newLst=[]
	for rw in lst:
		# assume: [[x,y]...]
		if rw[1]<minVal: continue
		newLst+=[[rw[0], minVal], [rw[0], rw[1]], [rw[0], minVal]]
	#
	return newLst

def getIntervals(catList, winLen):
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
	
# eqctalog moved to separate module (included at top)

def yodaecho(something=None):
	# this is for loadFileToH/VList() so we can have a null conversion (return whatever is there).
	return something

def loadFileToHlist(fname=None, castFunct=yodaecho):
	# in a later version, consider castFunct=[int, int, float, str, int...] so we can make a mask to convert
	if fname==None: return None
	try:
		x=castFunct(42.42)
	except:
		castFunct=yodaecho
		
	#
	f=open(fname, 'r')
	X=[]
	for rw in f:
		if rw[0] in ['#', ' ', '\t']: continue
		rws=rw.split('\t')
		X+=[[]]
		for elem in rws:
			X[-1]+=[castFunct(elem)]
			
		#
	#
	f.close()
	#
	return X

def floatYear(thisdt):
	yr=thisdt.year
	Ndays=dtm.datetime(yr,12,31).timetuple()[7]
	T=thisdt.timetuple()
	nday=T[7]-1
	dayfrac=(T[3]*60.*60.*10**6+T[4]*60*10**6 + T[5]*10**6 + thisdt.microsecond)/((24.*60.*60. + 60*60 + 60)*10**6)
	#
	return yr + (nday+dayfrac)/Ndays
floatyear=floatYear

def decistring(floatval, ndeciplaces):
	#return a string of a float with a fixed number of decimals
	strval=str(round(floatval, ndeciplaces))
	strvals=strval.split('.')
	return strvals[0]+'.'+strvals[1][0:ndeciplaces]

'''
def logSpacedLogPoints(XY, logbase=10.0):
	# we want [X,Y]. we'll figure out how to reshape [[x,y]...] arrays later.
	# this version for power-law (or exponential) data (aka, we'd plot it logX).
	#
	outsies=[[],[]]	#output array.
	X=XY[0]
	Y=XY[1]
	logsies=[int(math.log(X[0], logbase))]
	outsies[0]+=[X[0]]
	outsies[1]+=[Y[0]]
	for i in xrange(1, len(X)):
		logx=int(math.log(X[i], logbase))
		#if logx==math.log(outsies[0][-1], logbase): continue
		if logx==logsies[-1]: continue
		logsies+=[logx]
		outsies[0]+=[X[i]]
		outsies[1]+=[Y[i]]
	#
	return outsies
'''
fitMarkerShortList=['o', '^', 's', 'p', '*', 'h', '+', 'H', 'D', 'x']

def integerSpacedPoints(XY, intFactor=10):
	# linear data come in (presumably the log of PL/exp data), so the data density
	# are probably higher at one end of the distribution.
	# intFactor: basically, how many points between base-10 integers (effectively, base for log-binning).
	outsies=[[],[]]	#output array.
	X=XY[0]
	Y=XY[1]
	intsies=[int(X[0]/intFactor)]
	outsies[0]+=[X[0]]
	outsies[1]+=[Y[0]]
	#
	for i in range(1, len(X)):
		intx=int(X[i]/intFactor)
		if intx==intsies[-1]: continue
		intsies+=[intx]
		outsies[0]+=[X[i]]
		outsies[1]+=[Y[i]]
	return outsies
'''
# see ANSStools.py
def catfromANSS(lon=[135., 150.], lat=[30., 41.5], minMag=4.0, dates0=[dtm.date(2005,01,01), None], Nmax=999999, fout='cats/mycat.cat'):
	# oops. this is a little bit of a disaster. this version of catfromANSS returns the whole catlist, in which
	# [..., dept, mag, ...]
	# another version in ANSStools (about to become standard) returns a spacialized list [...,mag,depth].
	# going forward, it will be necessary to be very very careful.
	# perhaps we'll load ANSS tools as a named namespace so legacy apps will break, not do something random.
	# that said, most of the legacy apps only use this function to write files, not return catalog lists...
	#
	# get a basic catalog. then, we'll do a poly-subcat. we need a consistent catalog.
	# eventually, cut up "japancatfromANSS()", etc. to call this base function and move to yodapy.
	if dates0[1]==None:
		# i think this needs a "date" object, and datetime breaks.
		# so, make a Now() for date.
		nowdtm=dtm.datetime.now()
		dates0[1]=dtm.date(nowdtm.year, nowdtm.month, nowdtm.day)
	#	
	catlist=getANSSlist(lon, lat, minMag, dates0, Nmax, None)
	f=open(fout, 'w')
	f.write("#anss catalog\n")
	f.write("#lon=%s\tlat=%s\tm0=%f\tdates=%s\n" % (str(lon), str(lat), minMag, str(dates0)))
	f.write("#dtm, lat, lon, mag, depth\n")
	for rw in catlist:
		# simple, right? except that ANSS has a habit of writind useless date-times like "2001/10/08 24:00:07.62" (hour=24), or
		# where minute=60. we could toss these. for now, assume 2001/10/8 24:00:00 -> 2001/10/9/00:00:00. change by proper time-arithmetic.
		# first, parse the date-string:
		strDt, strTm=rw[0].split()[0], rw[0].split()[1]
		if '/' in strDt: delim='/'
		if '-' in strDt: delim='-'
		strDts=strDt.split(delim)
		strTms=strTm.split(':')
		yr=int(strDts[0])
		mnth=int(strDts[1])
		dy=int(strDts[2])
		hr=int(strTms[0])
		mn=int(strTms[1])
		sc=float(strTms[2])
		microsecs=(10**6)*sc%1.
		# one approach is to start with year, month and add all the subsequent quantities using datetime.timedelta objects, which we have to
		# do once we get into callendar addition anyway...
		#so let's assume the date part is correct:
		myDt=dtm.datetime(yr, mnth, dy)
		#mytimedelta=dtm.timedelta(hours=hr)
		myDt+=dtm.timedelta(hours=hr)
		myDt+=dtm.timedelta(minutes=mn)
		myDt+=dtm.timedelta(seconds=sc)
		myDt+=dtm.timedelta(microseconds=microsecs)
		#
		myDtStr='%d/%d/%d %d:%d:%d.%d' % (myDt.year, myDt.month, myDt.day, myDt.hour, myDt.minute, myDt.second, myDt.microsecond)	
		#
		#f.write('%s\t%s\t%s\t%s\n' % (rw[0], rw[1], rw[2], rw[4]))
		#f.write('%s\t%s\t%s\t%s\n' % (myDtStr, rw[1], rw[2], rw[4]))
		f.write('%s\t%s\t%s\t%s\t%s\n' % (myDtStr, rw[1], rw[2], rw[4], rw[3]))
	f.close()
	
	 
	return catlist
'''

'''
# see ANSStools.py
def getANSStoFilehandler(lon=[-125, -115], lat=[32, 45], minMag=4.92, dates0=[dtm.date(2001,01,01), dtm.date(2010, 12, 31)], Nmax=999999, keywds=''):
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
	dates=[dtm.date(dates0[0].year, dates0[0].month, dates0[0].day), dtm.date(dates0[1].year, dates0[1].month, dates0[1].day)]
	anssPrams={'format':'cnss', 'output':'readable', 'mintime':str(dates[0]).replace('-', '/'), 'maxtime':str(dates[1]).replace('-', '/'), 'minmag':str(minMag), 'minlat':lat[0], 'maxlat':lat[1], 'minlon':lon[0], 'maxlon':lon[1], 'etype':'E', 'searchlimit':Nmax, 'keywds':keywds}
	f = urllib.urlopen('http://www.ncedc.org/cgi-bin/catalog-search2.pl', urllib.urlencode(anssPrams))
	#
	# we might return f, a string of f, or maybe a list of lines from f. we'll work that out shortly...
	return f

def getANSSlist(lon=[-125., -115.], lat=[32., 45.], minMag=4.92, dates0=[dtm.date(2001,01,01), dtm.date(2010, 12, 31)], Nmax=999999, fin=None, keywds=''):
	#
	# note: this appears to be a bad idea for global downloads. a full catalog is ~4GB, which kills my computer.
	#
	# note: this may be repeated exactly in ygmapbits.py
	# fetch new ANSS data; return a python list object of the data.
	# fin: data file handler. if this is None, then get one from ANSS.
	dates=[dtm.date(dates0[0].year, dates0[0].month, dates0[0].day), dtm.date(dates0[1].year, dates0[1].month, dates0[1].day)]
	anssList=[]
	if fin==None:
		#print "get data from ANSS...(%s, %s, %s, %s, %s)" % (lon, lat, minMag, dates, Nmax)
		fin = getANSStoFilehandler(lon=lon, lat=lat, minMag=minMag, dates0=dates, Nmax=Nmax, keywds=keywds)
		#fin = getANSStoFilehandler([-180, 180], [-90, 90], 0, [datetime.date(1910,01,01), datetime.date(2010, 01, 16)], 9999999)

		print "data handle fetched..."
		
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
		anssList+=[[rwEvdt, rwLat, rwLon, rwDepth, rwMag, rwMagType, rwNst, rwGap, rwClo, rwrms, rwsrc, rwCatEventId]]
	return anssList
'''
		
