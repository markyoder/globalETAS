#! /usr/bin/python

# import pylab as plt
import operator
import urllib.request, urllib.parse, urllib.error
#import httplib	# apparently, urllib is a happy wrapper that uses http lib; httplib is rarely used directly.
import datetime as dtm
import glob
import os

# output like:
# new google.maps.LatLng(33.1598,-115.537), new google.maps.LatLng(32.4545,-118.1633), new google.maps.LatLng(32.3067, -115.2278)
# print '<script type="text/javascript">\n'

def manualQuakeInfo():
	mystr=''
	mystr+= 'pyquakesLL = new Array('\
		'new google.maps.LatLng(33.1598,-115.537), new google.maps.LatLng(32.4545,-118.1633), '\
    	'new google.maps.LatLng(32.3067, -115.2278), new google.maps.LatLng(32.0, -115.0), '\
		'new google.maps.LatLng(43.0, -118.0)'\
    	');\n'
	mystr+= 'pyInfo = new Array('\
		'"5.11 2005/09/02", "4.99 2005/10/16", "5.37 2006/05/24", "11 2009/10/10", "dunno..."' \
		');'
		
	return mystr

def sayHello(strHello='billyjoebob'):
	return "pyHello %s!" % strHello

def sayHello2(strHello='billyjoebob2'):
	return "pyHello %s!" % strHello


def fullModpyDir(fulldir, pathToRoot='../../../../../'):
	rstr= pathToRoot + fulldir
	rstr.replace('//', '/')	# assume it can only happen once...
	return rstr
#
def getPIsquareRays(fname='relm_pi.dat', L=.1):
	# get a list of lists of poly-vertexes from a file of PI centers.
	# L = square side length.
	f=open(fname)
	squares=[]
	for rw in f:
		if rw[0]=='#': continue
		while rw.find('  ')>=0: rw.replace('  ', ' ')	# remove double-spaces
		#
		rws=rw.split(' ')
		lat=float(rws[1])
		lon=float(rws[0])
		#pi=rws[2]	#we're not using this just yet...
		squares+=[[[lat-L/2.0, lon-L/2.0], [lat+L/2.0, lon-L/2.0], [lat+L/2.0, lon+L/2.0], [lat-L/2.0, lon+L/2.0], [lat-L/2.0, lon-L/2.0]]]
		
	return squares

def getPIsquareRaysData(data, L=.1):
	# get a list of lists of poly-vertexes from a file of PI centers.
	# L = square side length.
	squares=[]
	for rws in data:
		lat=float(rws[1])
		lon=float(rws[0])
		#pi=rws[2]	#we're not using this just yet...
		squares+=[[[lat-L/2.0, lon-L/2.0], [lat+L/2.0, lon-L/2.0], [lat+L/2.0, lon+L/2.0], [lat-L/2.0, lon+L/2.0], [lat-L/2.0, lon-L/2.0]]]
		
	return squares
#
def getGmapPolyAry(polyData, sColor="#FF0000", sOpacity=0.8, sWeight=2, fColor="#AA2500", fOpacity=0.35):
	'''
	mypoly1 = new google.maps.Polygon({paths: myPolyRay1,
	strokeColor: "#FF0000",
	strokeOpacity: 0.8,
	strokeWeight: 2,
	fillColor: "#AA2500",
	fillOpacity: 0.35
	});
	
	new google.maps.LatLng(34.0, -120.5)
	'''
	# polyData are PI polygon vertices.
	#
	strPolyRay=' new Array('
	for sqr in polyData:
		strPolyRay+='new google.maps.Polygon({paths: new Array('
		for vert in sqr:
			strPolyRay+='new google.maps.LatLng(%f, %f),' % (vert[0], vert[1])
		strPolyRay=strPolyRay[:-1]		# get rid of the last comma
		strPolyRay+='), strokeColor: "%s", strokeOpacity: %f, strokeWeight: %d, fillColor: "%s", fillOpacity: %f}),' % (sColor, sOpacity, sWeight, fColor, fOpacity)
	strPolyRay=strPolyRay[:-1]+')'
	#
	
	return strPolyRay
#
def getPolyAryString(fname='relm_pi.dat', L=.1, polyData=None, sColor="#FF0000", sOpacity=0.8, sWeight=2, fColor="AA2500", fOpacity=0.35):
	# a little wrapper so we don't always have to read a file to get kml-polygons from a data-set.
	#pdata=getPIsquareRays(fname, L)
	pdata=getReducedPolys(fname, L)	# same type of data as getPIsquareRays(), but reduced to compact polygons.
	return getGmapPolyAry(pdata, sColor, sOpacity, sWeight, fColor, fOpacity)
	
#
def getScorecardFile(fname='scripts/scorecardevents.dat'):
	# get earthquake list from scorecard data file.
	strQuakes='pyquakesLL = new Array('
	strInfo='pyInfo = new Array('
	
	f=open(fname)
	datamode=None
	for rw in f:
		if rw=="#!CA earthquakes\n":
			datamode='equakes'
		if rw=="#!CA earthquake labels\n":
			datamode='labels'
		#
		if rw[0] in [' ', '\t', '#']: continue
		#
		if datamode=='equakes':
			rws=rw.split('\t')
			if len(rws)<4: continue
			rws[3]=rws[3][:-1]
			#
			#lon=float(rws[0])
			#lat=float(rws[1])
			#mag=float(rws[2])
			#dt=rws(5)
			
			strQuakes+='new google.maps.LatLng(%s, %s),' % (rws[1], rws[0])
			#strInfo+='"%s %s",' % (rws[2], rws[5])
			strInfo+='"%s %s",' % (rws[2], rws[3])

		if datamode=='labels': continue
	#
	strQuakes=strQuakes[:-1] + ');\n'
	strInfo=strInfo[:-1] + ');\n'
	
	#print strQuakes
	#print strInfo
	#
	return strQuakes + strInfo

def getReducedPoints(indata, L=.1):
	# reduces a cluster of contiguous points by making squares around center points
	# and then removing vertices that are evenly redundant (aka, n=2)
	# at this point, this function is mostly an academic exercise. a similar approach with vectors (taken from the square sets perhaps)
	# could be useful. basically, given an array of vectors, look for vectors with tip/tail relationships and in the same direction. replace
	# that with one longer vector. then, remove or rotate other subsequent vectors appropriately (aka, rotate the "down" vector on the right side of a square
	# about the tip to oppose the direction of the longer vector on the opposite side of the square.
	sqrs = getPIsquareRaysData(indata, L)
	rsqrs=[]
	for rw in sqrs:
		# each square is a closed polygon; x[0]=x[-1]
		rsqrs+=[rw[:-1]]
	sqrs=None
	#
	# make a plain list of unique indices and count (or skip this step and only admit odd counts):
	#rverts=[[],[]]		#[[vertices], [counts]]
	rverts=[]
	for rw in rsqrs:
		for elem in rw:
			# for each vertex, is it in rverts? if so, remove it.
			# if not, add it. this will allow only odd counts to persist.
			if elem in rverts:
				rverts.remove(elem)
				continue
			else:
				rverts+=[elem]
	#
	return rverts

def nearlyIn(vec, ary, precis=4.0):
	# we seem to have a floating point error problem identifying elements that are "in" an array.
	# for now, vec is a [tail, tip] pair; ary is a list of tt-pairs.
	# ...though this does not seem to be the problem...
	
	rval=False
	#
	tailx=float(vec[0][0])
	taily=float(vec[0][1])
	tipx=float(vec[1][0])
	tipy=float(vec[1][1])
	#
	i=0
	for vc in ary:
		thistailx=float(vc[0][0])
		thistaily=float(vc[0][1])
		thistipx=float(vc[1][0])
		thistipy=float(vc[1][1])
		#
		# are they close enough?
		# in this case, we want to check the reverse vector as well.
		#
		#taildx=(1-(thistailx/tailx))**2
		#taildy=(1-(thistaily/taily))**2
		#tipdx=(1-(thistipx/tipx))**2
		#tipdy=(1-(thistipy/tipy))**2
		taildx=abs(1-(thistailx/tailx))
		taildy=abs(1-(thistaily/taily))
		tipdx=abs(1-(thistipx/tipx))
		tipdy=abs(1-(thistipy/tipy))
		#
		if taildx+taildy+tipdx+tipdy<(10**(-precis)):
			#rval=True
			print("nearlyIn: %s, %s" % (vec, vc))
			rval=i+1
			break
		i+=1
	#
	return rval	
	
def getReducedPolys(fname='relm_pi.dat', L=.1):
	# this replaces (with a reduced set) the base call: getPIsquareRays()
	#
	# this is an experiment. reduce a bunch of squares (from say a bunch of PI data) to polynomials.
	squares = getPIsquareRays(fname, L)	# array of square polynomial vertices. [[sq1], [sq2],...]
	# now, convert the polynomials to vectors like: [[[x1, y1], [x2, y2]] ]
	allVecs=[]
	for rw in squares:
		for i in range(1,len(rw)):
			#allVecs+=[[ [rw[i-1][0], rw[i-1][1]], [rw[i][0], rw[i][1]] ]]
			#allVecs+=[[ [float(rw[i-1][0]), float(rw[i-1][1])], [float(rw[i][0]), float(rw[i][1])] ]]
			# somewhere in the float handling processes, we get funny inequalities, aka, 32.650000000000006!=32.65, which seems obvious.
			# the point is that python will turn 32.65 to 32.650000000000006 as the process of some float conversion and consequently, we cannot
			# find repeated vectors in our lists. note also that: 
			#  32.650000000000006==32.650000000000002 : False
			#  32.650000000000006==32.650000000000003 : True
			allVecs+=[[ [round(float(rw[i-1][0]),6), round(float(rw[i-1][1]),6)], [round(float(rw[i][0]),6), round(float(rw[i][1]),6)] ]]
		#
	#
	# now, remove all the redundant vectors, aka vectors that are identical or opposite (x1->x2 = -(x2->x1)).
	
	rvecs=[]
	for vec in allVecs:
		revVec=[vec[1], vec[0]]
		if vec in rvecs:
			while vec in rvecs:
				rvecs.remove(vec)
			#while vec in allVecs:
			#	allVecs.remove(vec)
			continue
		elif revVec in rvecs:
			while revVec in rvecs:
				rvecs.remove(revVec)
			#while revVec in allVecs:
			#	allVecs.remove(revVec)
			continue
		else:
			rvecs+=[vec]
	#return rvecs
	
	# get a count:
	rvecstemp=[]
	rvecscount=[]
	for vec in rvecs:
		revVec=[vec[1], vec[0]]
		doLoop=1
		if vec in rvecstemp:
			indx=rvecstemp.index(vec)
			rvecscount[indx]+=1
		elif revVec in rvecstemp:
			indx=rvecstemp.index(revVec)
			rvecscount[indx]+=1
		else:
			rvecstemp+=[vec]
			rvecscount+=[1]
	# return [rvecstemp, rvecscount]
	newRvecs=[]
	for i in range(len(rvecstemp)):
		if rvecscount[i]==1: newRvecs+=[rvecstemp[i]]
		if rvecscount[i]>1: print("%d, %s" % (rvecscount[i], rvecstemp[i]))
	
	# print "lens: %d, %d" % (len(rvecs), len(newRvecs))
	rvecs=newRvecs
	
			
	allVecs=None
	#
	polys=[]
	polygons=[]
	while len(rvecs)>0:
		currVec=rvecs.pop()
		##polys+=[[rvecs.pop()]]
		polys+=[[currVec]]
		polygons+=[[currVec[0], currVec[1]]]	# start with the tail
		
		#thisTail=polys[-1][-1][0]
		#thisTip=polys[-1][-1][1]
		thisTail=currVec[0]
		thisTip=currVec[1]
		#
		#print "len(rvecs): %d " % len(rvecs)
		i=0
		while i<len(rvecs):
			if i>=len(rvecs): break
			#
			# check for tip/tails:
			#print 'rvecs[%d]: %s' % (i, rvecs[i])
			newTail=rvecs[i][0]
			newTip=rvecs[i][1]
			if newTail==thisTip:
				#thisTail=newTail
				thisTip=newTip
				newVec=rvecs[i]
				polys[-1]+=[newVec]
				polygons[-1]+=[newVec[1]]
				rvecs.remove(newVec)
				i=0
				continue
			if len(rvecs)==0: break
			'''
			if newTip==thisTail:
				thisTail=newTail
				newVec=rvecs[i]
				polys[-1].insert(0, newVec)
				polygons[-1].insert(0, newVec[0])
				rvecs.remove(newVec)
				i=0
				continue
			if len(rvecs)==0: break
			'''
			i+=1				
			#
		#
		#if polygons[-1][-1]!=polygons[-1][0]: polygons[-1]+=[polygons[-1][0]]
	#
	# now, compress the polygons:
	for ii in range(len(polygons)):
		if len(polygons[ii])<=5: continue		# just a square...
		#for i in xrange(1,len(polygons[ii])):
		i=1
		while i<(len(polygons[ii])-2):
			if polygons[ii][i-1][0]==polygons[ii][i][0] and polygons[ii][i-1][0]==polygons[ii][i+1][0]:
				# 3 consecutive equal x-values: a vertical line. the middle point is redundant.
				polygons[ii].pop(i)
				continue
			if polygons[ii][i-1][1]==polygons[ii][i][1] and polygons[ii][i-1][1]==polygons[ii][i+1][1]:
				# horizontal line:
				# 3 consecutive equal x-values: a vertical line. the middle point is redundant.
				polygons[ii].pop(i)
				continue
			i+=1
			if i>= (len(polygons[ii])-1): break
	
	#return polys
	#return [rvecstemp, rvecscount]
	return polygons

####
# a set of functions to get ANSS data at at differend levels:
# getANSShtml(): returns a file object. read with f.read(); this includes the HTML tags, etc.
# getANSSlist(): calls getANSShtml, then parses to a list of datetime, lat, lon, depth, mag, magt, Nst, Gap, Clo, RMS, Src, EventID
# get ANSSregion

def getANSSregion(region='socal', minMag=5.0, dates=[dtm.date(2001,0o1,0o1), dtm.date(2010, 12, 31)], Nmax=999999):
	# some pre-defined regions like socal, norcal, cal, parkfield, global, blah, blah.
	# since we can't overload the declaration, first check the 'region' for variable type.
	# assemble lat, lon limits from the following lists. otherwise, if it's a list or tuple or whatever, use those vals.
	#
	# set parameters from region...
	#
	# return getANSSdate(lon, lat, minMag, dates)
	return None
	
def getANSShtml(lon=[-125, -115], lat=[32, 45], minMag=4.92, dates=[dtm.date(2001,0o1,0o1), dtm.date(2010, 12, 31)], Nmax=999999):
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
	anssPrams={'format':'cnss', 'output':'readable', 'mintime':str(dates[0]).replace('-', '/'), 'maxtime':str(dates[1]).replace('-', '/'), 'minmag':str(minMag), 'minlat':lat[0], 'maxlat':lat[1], 'minlon':lon[0], 'maxlon':lon[1], 'etype':'E', 'searchlimit':Nmax}
	f = urllib.request.urlopen('http://www.ncedc.org/cgi-bin/catalog-search2.pl', urllib.parse.urlencode(anssPrams))
	#
	# we might return f, a string of f, or maybe a list of lines from f. we'll work that out shortly...
	return f

def getANSSlist(lon=[-125, -115], lat=[32, 45], minMag=4.92, dates=[dtm.date(2001,0o1,0o1), dtm.date(2010, 12, 31)], Nmax=999999):
	anssList=[]
	fin=getANSShtml(lon, lat, minMag, dates, Nmax)
	for rw in fin:
		if rw[0] in ["#", "<"] or rw[0:4] in ["Date", "date", "DATE", "----"]: continue
		#anssList+=[rw[:-1]]
		# data are fixed width delimited
		# return date-time, lat, lon, depth, mag, magType, nst, gap, clo, rms, src, catEventID (because those are all the available bits)
		anssList+=[[rw[0:22].strip(), float(rw[23:31].strip()), float(rw[32:41].strip()), float(rw[42:48].strip()), float(rw[49:54].strip()), rw[55:59].strip(), float(rw[60:64].strip()), rw[65:68].strip(), rw[69:73].strip(), float(rw[74:78].strip()), rw[79:83].strip(), rw[84:96].strip()]]
	return anssList
			
			

		####################3
#		strValues += str(catID).strip() + ", '"
#			strValues += myLine[0:22].strip() + "', "								# date-time
#			strValues += myLine[23:31].strip() + ", " 							# lat
#			strValues += myLine[32:41].strip() + ", "								# lon
#			strValues += strNullIfBlank(myLine[42:48].strip()) + ", "		# depth
#			strValues += myLine[49:54].strip() + ", '"							# mag
#			strValues += myLine[55:59].strip() + "', "							# magType
#			strValues += strNullIfBlank(myLine[60:64].strip()) + ", "		# nst
#			strValues += strNullIfBlank(myLine[65:68].strip()) + ", "		# gap
#			strValues += strNullIfBlank(myLine[69:73].strip()) + ", "		# clo
#			strValues += strNullIfBlank(myLine[74:78].strip()) + ", '"		# rms
#			strValues += myLine[79:83].strip() + "', '"							# src
#			strValues += myLine[84:96].strip() + "'"								# catEventID
		##################



def updateEventslist(lon=[-125, -115], lat=[32, 45], minMag=4.92, dates=[dtm.date(2001,0o1,0o1), dtm.date(2010, 12, 31)], Nmax=999999, foutname='scorecardevents.dat', expTimeDelta=dtm.timedelta(days=1)):
	import datetime as dtm
	# check currency of the events list. if it is expired, create a new data file from ANSS. return 0,1,2 (error, current, updated) as appropriate (?dunno?)
	#
	# is there an earthquakes datafile?
	isAfile=True
	getNewANSSdata=True	# we start with true default in case the file has no mod_date.
	#
	fouts=glob.glob(foutname)
	if len(fouts)<0:
		isAfile=False
		getNewANSSdata=True
	#
	modDate=dtm.datetime(1900, 1,1, 0,0,0)	# just in case there is no mod_date...
	# is the file current?
	# we could use something like fileInfo=os.stat(fname); modDate=fileInfo.st_mtime
	# which returns the mod-date in seconds maybe? my experience is that this is a relatively unreliable way to determine whether
	# or not a file is current (funny things can happen to change mod dates and all that good stuff). it is much more reliable
	# to evaluate information in the file (or a database query is even better). also, different platforms might use different
	# time formats (milliseconts, days, etc.), so .st_mtime might not compare to datetime.datetime.now() properly.
	if isAfile:
		# is it current:
		modDate=dtm.datetime(1900, 1,1, 0,0,0)	# just in case there is no mod_date...
		fList=[]
		f=open(foutname, 'r')
		for rw in f:
			# look for info:
			if rw[0:11] == '#!mod_date:':
				# this is the modification date...
				rws=rw.split('\t')
				modDtm=rws[1][:-1].split(' ')	#date on in [0], time in [1]
				delim='-'
				if '/' in modDtm[0]: delim=['/']
				lstDate=modDtm[0].split(delim)
				lstTime=[0,0,0]
				if len(modDtm)>1:
					lstTime=modDtm[1].split(':')
				while len(lstTime)<3:
					lstTime+=[0]
				for i in range(len(lstDate)): lstDate[i]=int(lstDate[i])
				for i in range(len(lstTime)): lstTime[i]=int(float(lstTime[i]))
				modDate=dtm.datetime(lstDate[0], lstDate[1], lstDate[2], lstTime[0], lstTime[1], int(lstTime[2]))	# note, we skip fractional seconds.
				freshDate=dtm.datetime.today()-expTimeDelta
				#print modDate, freshDate, modDate<freshDate
				if modDate<freshDate:
					getNewANSSdata=True
				else:
					getNewANSSdata=False
				break
				#
			#
		f.close()
		#
	# at this point, we know if we have an events file and if it is current. if the file does not exist
	# or is not current, fetch new data from ANSS. if not, simply return
	#
	if getNewANSSdata==False: return 1
	#
	#print "getting fresh data"
	# otherwise, retrieve data from ANSS and write a new data file:
	# for speed, just code this here; write direcly to the output file:
	fin=getANSShtml(lon, lat, minMag, dates, Nmax)
	# tempFname='scripts/temp/tempEvents.dat'
	tempFname=fullModpyDir('/var/www/quakemaps/scripts/temp/tempEvents.dat')
	tempFile=open(tempFname, 'w')
	tempFile.write("#!mod_date:\t%s\n" % str(dtm.datetime.today()))	# date the file
	tempFile.write("#!CA earthquakes\n")
	for rw in fin:
		if rw[0] in ["#", "<"] or rw[0:4] in ["Date", "date", "DATE", "----"]: continue
		#anssList+=[rw[:-1]]
		# data are fixed width delimited
		# return date-time, lat, lon, depth, mag, magType, nst, gap, clo, rms, src, catEventID (because those are all the available bits)
		#anssList+=[[rw[0:22].strip(), float(rw[23:31].strip()), float(rw[32:41].strip()), float(rw[42:48].strip()), float(rw[49:54].strip()), rw[55:59].strip(), float(rw[60:64].strip()), rw[65:68].strip(), rw[69:73].strip(), float(rw[74:78].strip()), rw[79:83].strip(), rw[84:96].strip()]]
		# we want tab delimited file with lon, lat, mag, date
		dtm=rw[0:22].strip()
		dt=dtm.split(' ')[0].replace('-', '/')	# i think we want / delimiters...
		tempFile.write("%f\t%f\t%s\t%s\n" % (float(rw[32:41]), float(rw[23:31]), float(rw[49:54]), dtm.split(' ')[0] ))
		#
	tempFile.close()
	#
	# now, we've written a temp file. swap the old one out for the new...
	#
	# write a backup (overwrites daily) and copy the temp file into production place.
	bkpFname=''
	fnames=foutname.split('.')
	fnames[-2]+='-bkp'
	for bit in fnames:
		bkpFname+=bit+'.'
	bkpFname=bkpFname[:-1]
	os.system('cp %s %s' % (foutname, bkpFname))
	#
	os.system('cp %s %s' % (tempFname, foutname))
	
####################
###################
#####
def pi2kml(inFilename = 'data/deltaPdither.xyz', outFilename = 'data/myPI1.kml', delim='\t', xy=[0,1]):
	#from kmlparser import *
	import kmlperser as kmp

	#inFilename = 'data/deltaPdither.xyz'
	#outFilename = 'data/myPI1.kml'
	
	#newfile = 'deltaPsmall.xyz'
	#
	myfilein = open(inFilename, 'r')
	myfileout = open(outFilename, 'a')
	#myNewFile = open(newfile, 'w')
	#
	myfileout.write('<?xml version="1.0" encoding="UTF-8"?>\n<kml xmlns="http://www.opengis.net/kml/2.2">\n\n')

	# clearly, this is a stupid way to add styles. eventually, make a standard set of KML to include...
	myfileout.write(kmp.add_kmlstyle('redsies'))
	myfileout.write("\n")

	myfileout.write(kmp.add_kmlstyle('transBluePoly'))
	myfileout.write("\n")

	myfileout.write(add_kmlstyle('orangsies'))
	myfileout.write("\n")


	ict=0
	for lines in myfilein:
		#print lines
		#print len(lines)
		linestuple = lines.split(delim)		# assuming a bit about the input file of course...
		#linestuple[2]=linestuple[2][:-2]
		linestuple[-1]=linestuple[-1][:-1]	#CRLF at the end of the line...
		#ict = ict+1
		if abs(eval(linestuple[2]))>=.1:		#PI value...
		#if eval(linestuple[2])!=0:
			thisStyle='orangsies'
			if abs(eval(linestuple[2]))>=.5:
				thisStyle='redsies'
			myfileout.write(kmp.add_square(linestuple[xy[1]], linestuple[xy[0]], .1, 0, '', '', thisStyle))
			myfileout.write("\n")
			#print linestuple
			#myNewFile.write(lines)

	#print ict
	#myfileout = open(outFilename, 'a')
	myfileout.write("</kml>")
	myfileout.close
	#myNewFile.close

	

	


