import math
import scipy
#import pylab
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as dtm
#import numpy.fft as nft
import scipy.optimize as spo
#from matplotlib import pyplot as plt
#import pylab as plt
#from matplotlib import rc
import numpy


#import string
#import sys
#from matplotlib import *
#from pylab import *
#import os
import random
#
# gamma function lives here:
#import scipy.special
#from scipy.special import gamma
#from scipy.optimize import leastsq
from matplotlib import axis as aa
#from threading import Thread
#
import datetime as dtm
#import time
import pytz
import calendar
import operator
import urllib.request, urllib.parse, urllib.error

import rbIntervals as rbi

#from eqcatalog import *
# 20120719 yoder: (this might have impacts on apps. that auto-fetch from ANSS)
from ANSStools import *

# maping bits:
#import matplotlib	# note that we've tome from ... import *. we should probably eventually get rid of that and use the matplotlib namespace.
#matplotlib.use('Agg')
#from matplotlib.toolkits.basemap import Basemap
#from mpl_toolkits.basemap import Basemap as Basemap
#from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
#

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
				self.datas[2]+=[1.0]
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
		print(("prams: %s" % str(p)))
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
		
		if xmin==None:
			# xmin=thisdatas[0][0]
			xmin = min(thisdatas[0])
		if xmax==None:
			# xmax=thisdatas[0][-1]
			xmax = max(thisdatas[0])
		# 
		# get X,Y,W:
		X=[]
		Y=[]
		W=[]
		#
		for i in range(len(thisdatas[0])):
			#if thisdatas[0][i]<xmin: continue
			#
			X+=[thisdatas[0][i]]
			Y+=[thisdatas[1][i]]
			W+=[thisdatas[2][i]]
			if X[-1]>=xmax: continue
		
		#
		# now, fit data:
		#print scipy.array(Y), scipy.array(X), scipy.array(W)
		#print len(scipy.array(Y)), len(scipy.array(X)), len(scipy.array(W))
		#print self.activeRes
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
		self.varb=(1.0/self.delta)*sum(aWW)	# note that w=1/sig^2 in most texts, as per gaussian stats.
														# this is more general than that, and it might make a big fat
		#												# mess when applied outside gaussian distributions.
		#
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
		print(("from getFitPlotAry(): %s, %s, %s, %s" % (self.datas[0][0], self.datas[0][-1], self.a, self.b)))
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
	
