from Scientific.Geometry import *
from math import *
import string
#from matplotlib import *
from pylab import *
import copy

def sortXYtuples(data, sortIndex1, sortIndex2, S1):
	# note: this was a fun exercise before finding the map(operator.itemgetter()) built in way to sort
	# d>1 tuples. that said, can we use this pseudo-built-in method to sort on two indeces? i think you can do this
	# by creating a function to sort on (aka, create an index S[i][x].S[i][y], where a.b indicates some sort of integer
	# digit transformation so that 1,1 -> 11, 1,2, -> 12, 2,3-> 23, etc. (here, x_max=10).
	#
	# data is clearly your data - nested lists/tuples, whatever.
	# sortIndex1,2 is the index on which to sort (aka, ix, iy)
	# S1,2 are the dimensions of the system. aka, an x then y sort
	# on a 100, 200 grid.
	#
	# strategy: transform x,y into a 1D index: x,y -> y*X + x
	#
	dd = {}
	retuple = []
	#
	for rows in data:
		thisIndex = rows[sortIndex1]*S1+rows[sortIndex2]
		dd[thisIndex] = rows
	#
	listIndex = list(dd.keys())
	listIndex.sort()
	#x
	for myindex in listIndex:
		retuple = retuple + [dd[myindex]]
	#similarly, i think, though there is some proving to be done, that if we subtract a - scrambled(a), 
	return retuple
		
def getDataTuplesFromFile(filename, delim=None):
	myfile = open(filename, 'r')
	outdata = []
	for lines in myfile:
		if lines[0]=='#': continue
		#
		if delim==None:
			# guess a delimiter (obviously, this is not a sophisticated guess):
			delim=' '
			if '\t' in lines: delim='\t'
		#
		linestuple = lines.split(delim)
		for i in range(len(linestuple)):
			linestuple[i] = float(linestuple[i])
		#linestuple[2] = linestuple[2][:-2]
		outdata = outdata + [linestuple]
	#
	return outdata

def boxyContour2(data, z0=0, ix=0, iy=1, iz=2, dx=.1, dy=.1):
	# so i think this basically takes a bunch of contiguous boxes and makes a single polygon, which facilitates much more compact KML.
	# it might start with just the box center data points.
	
	# get max x,y values:
	maxX = 0
	maxY = 0
	for rows in data:
		#print rows
		if rows[ix]>maxX:
			maxX = rows[ix]
		if rows[iy]>maxY:
			maxY = rows[iy]
	#print maxX
	#print maxY
	# first, sort y, then x (scan in x):
	data1 = sortXYtuples(data, iy, ix, maxX)
	# scan the sorted data. when we get a state change (aka, step onto or off of a z0 site), set a point.
	cdata = []		# contour data set
	dstate = 0
	lastPos = [data1[0][ix]-dx, data1[0][iy]]
	for rows in data1:
		#  and rows[ix]=lastx+dx and rows[iy]=lasty
		# new region:
		if rows[iz]>=z0 and dstate==0:
			# we've hit a new "high" spot
			dstate=1
			if [rows[ix], rows[iy]] not in cdata:	
				#cdata = cdata + [rows[ix], rows[iy], rows[iz]]
				cdata = cdata + [[rows[ix], rows[iy]]]
		elif (dstate==1 and rows[iz]<z0):
			# we fall off the contour onto a data element
			dstate=0
			if lastPos not in cdata:
				# but we want to keep the previous element in the contour set...
				cdata = cdata + [lastPos]
		elif (dstate==1 and rows[iz]>=z0 and rows[ix]!=lastPos[0]+dx or rows[iy]!=lastPos[1] ):
			# we're on a contour and skip over some empty data elements onto a new contour.
			if lastPos not in cdata:
				cdata = cdata + [lastPos]
			if [rows[ix], rows[iy]] not in cdata:	
				cdata = cdata + [[rows[ix], rows[iy]]]
		#
		lastPos = [rows[ix], rows[iy]]

	#
	# now, sort by x then y, repeat with vertical scan, add to cdata:
	data1 = sortXYtuples(data, ix, iy, maxY)
	dstate = 0	
	lastPos = [data1[0][ix], data1[0][iy]-dy]
	# scan the sorted data. when we get a state change (aka, step onto or off of a z0 site), set a point.
	for rows in data1:
		#  and rows[ix]=lastx+dx and rows[iy]=lasty
		# new region:
		if rows[iz]>=z0 and dstate==0:
			# we've hit a new "high" spot
			dstate=1
			if [rows[ix], rows[iy]] not in cdata:	
				#cdata = cdata + [rows[ix], rows[iy], rows[iz]]
				cdata = cdata + [[rows[ix], rows[iy]]]
		# we fall off the contour onto a data element
		elif (dstate==1 and rows[iz]<z0):
			dstate=0
			if lastPos not in cdata:
				cdata = cdata + [lastPos]
		# we're on a contour and skip over some empty data elements onto a new contour.
		elif (dstate==1 and rows[iz]>=z0 and rows[iy]!=lastPos[1]+dy or rows[ix]!=lastPos[0] ):
			if lastPos not in cdata:
				cdata = cdata + [lastPos]
			if [rows[ix], rows[iy]] not in cdata:	
				cdata = cdata + [[rows[ix], rows[iy]]]
		#
	#
	# now we have a contour set. we need to connect these to make a polygon.
	# the way we've done this, lines can connect to nn or nnn lattice sites.
	# start with a random point in cdata (maybe the first point?).
	# it does not matter which direction we go; we just can't backtrack.
	# --> this basic method, but instead employing an index (ie, the xy index)
	# would probably be faster. the yx scan should probably calculate and use
	# the xy index as well...
	#
	#just in case (this might take a long time, so should maybe be skipped)
	#print "starting set"
	#cdata = set(cdata)
	#print  "finishing set"
	#
	return cdata

def findInTuples(mytuple, theData, startIndex=0):
	# returns the index of the found data.
	tLen = len(mytuple)
	retVal=-1
	rowth=0

	for rows in theData:
		isIn=[]
		for i in range(tLen):
			if mytuple[i]!=rows[i+startIndex]:
				# at least one pair of elements are not equal
				isIn=isIn+[False]
				break
			else:
				isIn=isIn+[True]	# this is actually unnecessary, but intuitive and instructional... an i guess it protects against empty tuples?
		# we've spun through a tuple. is it a hit?
		# if isIn contains no "false" entries, the (left nTuples of) theData record match mytuple
		if False not in isIn and True in isIn:
			# we've found a row that equals our tuple.
			retVal=rowth
			break
			#
		#
		rowth = rowth+1
	#
	#
	return retVal

def findHeadTailInTuples(mytuple, theData):
	# returns the index of the found data.
	# let's customize this. we'll always use the first 4 of the tuples.
	#tLen = len(mytuple)	# this should always be 4...
	retVal=-1
	rowth=0
	headTail=-1
	#
	for rows in theData:
		#print 'from findHeadTail: ' + str(mytuple) + " :: " + str(rows)
		if [mytuple[2], mytuple[3]] == [rows[0], rows[1]]:
			# head-to-tail
			retVal = rowth
			headTail=0
			break
			#
		if [mytuple[2], mytuple[3]] == [rows[2], rows[3]]:
			# head-to-head
			retVal = rowth
			headTail=1
			break
			#
		rowth = rowth+1
	#
	return [headTail, retVal]


def simpleContour2(data, z0=0, ix=0, iy=1, iz=2, dx=.1, dy=.1):
	# get max x,y values:
	maxX = 0
	maxY = 0
#	for rows in data:
#		#print rows
#		if rows[ix]>maxX:
#			maxX = rows[ix]
#		if rows[iy]>maxY:
#			maxY = rows[iy]
	# we need to make a copy of just the x,y parts of the data:
	d2file = open('data/data2.dat', 'w')
	data2=[]
	for rows in data:
		if abs(rows[iz])>=z0:
			#print rows[iz]
			data2 = data2 + [[int(round(rows[ix]/dx)), int(round(rows[iy]/dy))]]
			d2file.write(str(int(round(rows[ix]/dx))) + '\t' + str(int(round(rows[iy]/dy))) + '\n')
	d2file.close()
	#
	# now, hollow out data2 (in place? be careful of in place; index handling gets messy.)
	cdata = []
	for rows in data2:
		if [rows[ix]+1, rows[iy]] not in data2 or [rows[ix]-1, rows[iy]] not in data2 or [rows[ix], rows[iy]+1] not in data2 or [rows[ix], rows[iy]-1] not in data2 :
			print(str([rows[ix], rows[iy]]) + ' :: ' + str([rows[ix]-1, rows[iy]]))
			if not [[rows[ix], rows[iy]]] in cdata:
				cdata = cdata + [[rows[ix], rows[iy]]]
#		if findInTuples([rows[ix]+1, rows[iy]], data2)<0 or findInTuples([rows[ix]-1, rows[iy]], data2)<0:
#			print str([rows[ix], rows[iy]]) + ' :: ' + str([rows[ix]+1, rows[iy]])
#			if not [[rows[ix], rows[iy]]] in cdata:
#				cdata = cdata + [[rows[ix], rows[iy]]]	
	#
	return cdata

def boxyContour(data, z0=0, ix=0, iy=1, iz=2, dx=.1, dy=.1):
	# get max x,y values:
	maxX = 0
	maxY = 0
	minX = 9999999999
	minY = 9999999999
	for rows in data:
		#print rows
		if rows[ix]>=maxX:
			maxX = rows[ix]
		if rows[iy]>=maxY:
			maxY = rows[iy]
		if rows[ix]<=minX:
			minX = rows[ix]
		if rows[iy]<=minY:
			minY = rows[iy]
	# we need to make a copy of just the x,y parts of the data:
	d2file = open('data/data2.dat', 'w')
	data2=[]
	for rows in data:
		#if abs(rows[iz])>=z0:
		if rows[iz]>=z0:
			#print rows[iz]
			# subtract bias (starting position). this might not be QUITE right, as it starts the binning (if there is binning)
			# from the lowest value. it may be necessary to specify minY, minX.
			rows[ix] = rows[ix]-minX
			rows[iy] = rows[iy]-minY
			data2 = data2 + [[int(round(rows[ix]/dx)), int(round(rows[iy]/dy))]]
			d2file.write(str(int(round(rows[ix]/dx))) + '\t' + str(int(round(rows[iy]/dy))) + '\n')
			#d2file.write(str([[int(round(rows[ix]/dx)), int(round(rows[iy]/dy))]]) + '\n')
	d2file.close()
	data=None
	#
	# note: data2 values are in terms of integer-dx,dy values
	#
	# now, hollow out data2 (in place? be careful of in place; index handling gets messy.)
	#cdata = []
	rightsies = []
	leftsies = []
	topsies = []
	bottomsies = []
	fcont = open('data/contpoints.dat', 'w')
	for rows in data2:
		if [rows[ix]+1, rows[iy]] not in data2 or [rows[ix]-1, rows[iy]] not in data2 or [rows[ix], rows[iy]+1] not in data2 or [rows[ix], rows[iy]-1] not in data2 :
			fcont.write(str(rows[0]) + '\t' + str(rows[1]) + '\n')
			if [rows[ix]+1, rows[iy]] not in data2:
				# no right neighbor:
				rightsies.append(rows)
				#cdata = cdata + [[rows[ix]+.5, rows[iy]]]
			if [rows[ix]-1, rows[iy]] not in data2:
				# no left neighbor:
				leftsies.append(rows)
				#cdata = cdata + [[rows[ix]-.5, rows[iy]]]
			if [rows[ix], rows[iy]+1] not in data2:
				topsies.append(rows)
				#cdata = cdata + [[rows[ix], rows[iy]+.5]]
			if [rows[ix], rows[iy]-1] not in data2 :
				bottomsies.append(rows)
				#cdata = cdata + [[rows[ix]+.5, rows[iy]-.5]]
			#
			#if not [[rows[ix], rows[iy]]] in cdata:
			#	cdata = cdata + [[rows[ix], rows[iy]]]
		#
		# now, compact these to longer segments (aka, string together contiguous bits)
		#
		# sort each:
		# by x then y:
	# clean up a bit:
	data2=None		# finished with this array; save some memory for later...
	fcont.close()
	#print 'rightsies: ' + str(len(rightsies))
	#print 'leftsies: ' + str(len(leftsies))
	#print 'topsies: ' + str(len(topsies))
	#print 'bottomsies: ' + str(len(bottomsies))

	#vTail = [rightsies[0][0]+.5, rightsies[0][1]-.5]
	#thisHead = [rightsies[0][0]+.5, rightsies[0][1]+.5]
	rightSegs = []
	rightsies = sortXYtuples(rightsies, 0, 1, maxY)
	rightsies.append(rightsies[0])		# this tricks the subsequent code into closing the last vector (note we could alternatively do an
	vTail = rightsies[0]						# "append" operation after the loop
	vHead = vTail
	for row in rightsies:
		if row[0]==vTail[0] and row[1] in [vHead[1]+1, vHead[1]]:	# a vector can be its own head/tail
			vHead = row
		else:
			# rightSegs.append(vTail + vHead)
			#print vTail, vHead
			rightSegs.append([vTail[0]+.5, vTail[1]-.5, vHead[0]+.5, vHead[1]+.5])
			vTail = row
			vHead = row
	#rightSegs.append([vTail[0]+.5, vTail[1]-.5, vHead[0]+.5, vHead[1]+.5])
	#
	leftSegs = []
	leftsies = sortXYtuples(leftsies,0,1, maxY)
	leftsies.append(leftsies[0])
	vTail = leftsies[0]
	vHead = vTail
	for row in leftsies:
		if row[0]==vTail[0] and row[1] in [vHead[1]+1, vHead[1]]:	# a vector can be its own head/tail
			vHead = row
		else:
			#leftSegs.append(vTail + vHead)
			leftSegs.append([vTail[0]-.5, vTail[1]-.5, vHead[0]-.5, vHead[1]+.5])
			vTail = row
			vHead = row
	#leftSegs.append([vTail[0]-.5, vTail[1]-.5, vHead[0]-.5, vHead[1]+.5])
	#
	topSegs = []
	topsies = sortXYtuples(topsies,1,0, maxX)
	topsies.append(topsies[0])
	vTail = topsies[0]
	vHead = vTail
	for row in topsies:
		if row[1]==vTail[1] and row[0] in [vHead[0]+1, vHead[0]]:	# a vector can be its own head/tail
			vHead = row
		else:
			#topSegs.append(vTail + vHead)
			topSegs.append([vTail[0]-.5, vTail[1]+.5, vHead[0]+.5, vHead[1]+.5])
			vTail = row
			vHead = row
	#topSegs.append([vTail[0]-.5, vTail[1]+.5, vHead[0]+.5, vHead[1]+.5])
	#
	bottomSegs = []
	bottomsies = sortXYtuples(bottomsies,1,0, maxX)
	bottomsies.append(bottomsies[0])
	vTail = bottomsies[0]
	vHead = vTail
	vlen=0
	for row in bottomsies:
		if row[1]==vTail[1] and row[0] in [vHead[0]+1, vHead[0]]:	# a vector can be its own head/tail; we allow 2 conditions to jump-start the process.
			vHead = row
		else:
			# bottomSegs.append(vTail + vHead)
			bottomSegs.append([vTail[0]-.5, vTail[1]-.5, vHead[0]+.5, vHead[1]-.5])
			vTail = row
			vHead = row
	#bottomSegs.append([vTail[0]-.5, vTail[1]-.5, vHead[0]+.5, vHead[1]-.5])
	#
	# clean up a bit:
	rightsies = None
	leftsies = None
	topsies = None
	bottomsies = None
	#
	# and now, we have (4) sets of vectors outlining shapes. we need to arrange them head to tail.
	# we use each vector once, so we can "pop()" them out of the stack as we go.
	# we could combine these four sets of tuples, then for each row, search the whole stack.
	# arguably, we actually save a little bit of time by separating them. we know a LEFT vector
	# will not connect to another LEFT or RIGHT vector, so we can skip that stack. create a tuple of polygon-tuples
	polygons = []
	# these segment vectors are godless and without direction or morality, so we look for a fellow segment
	# with an overlapping head OR tail. we select our initial direction arbitrarily to be this vecor's head.
	# we loop until we connect to its tail.
	#

	while len(rightSegs + leftSegs + topSegs + bottomSegs)>0:
		#print 'remaining segments: ' + str(len(rightSegs + leftSegs + topSegs + bottomSegs))
		# start a new polygon:
		thisPoly=[]
		# start with any segment:
		if len(rightSegs)>0:
			thisSeg = rightSegs.pop(0)		# remove rightSegs[0] into thisSeg
		elif len(leftSegs)>0:
			thisSeg = leftSegs.pop(0)
		elif len(topSegs)>0:
			thisSeg = topSegs.pop(0)
		elif len(bottomSegs)>0:
			thisSeg = bottomSegs.pop(0)
		#
		#thisPoly.append(thisSeg)
		thisPoly.append([thisSeg[0]*dx + minX, thisSeg[1]*dy + minY])
		thisPoly.append([thisSeg[2]*dx + minX, thisSeg[3]*dy + minY])
		nextSeg = []
		polyTail = [thisSeg[0], thisSeg[1]]
		while [thisSeg[2], thisSeg[3]] != polyTail:		# note, this breaks for length zero polygons. could do a head & tail condition.
			# use only horizontal->vertical for vertical->horizontal vectors:
			# there is some clumsy, not super intuitive code here to optimize search time...
			# i should probably use the "ctypes" library to make pointers to these arrays
			# but i think for now, we'll just write more code.
			#
			#print thisSeg
			if thisSeg[0]==thisSeg[2]:	
				# horizontal vector; check for vertical next-vector:
				segIndex=0
				foundOne = findHeadTailInTuples(thisSeg, topSegs)
				# we get back a tuple like: [{0 if match tail, 1 if match head}, {index if found, else -1}]
				if foundOne[1]==-1:
					segIndex=1
					foundOne = findHeadTailInTuples(thisSeg, bottomSegs)
			elif thisSeg[1]==thisSeg[3]:
				# a vertical vector. we should probably throw in an error trap if both fail.
				segIndex=2
				foundOne = findHeadTailInTuples(thisSeg, leftSegs)
				if foundOne[1]==-1:
					segIndex=3
					foundOne = findHeadTailInTuples(thisSeg, rightSegs)
			#print 'foundone' + str(foundOne)
			#
			if foundOne[1]>=0:	# at this point, it better be...
				# fetch the segment from the appropriate tuple:
				if segIndex==0:
					nextSeg=topSegs.pop(foundOne[1])
				if segIndex==1:
					nextSeg=bottomSegs.pop(foundOne[1])
				if segIndex==2:
					nextSeg=leftSegs.pop(foundOne[1])
				if segIndex==3:
					nextSeg=rightSegs.pop(foundOne[1])
				#
				if foundOne[0]==1:
					# we matched head-to-head; switch order.
					nextSeg=[nextSeg[2], nextSeg[3], nextSeg[0], nextSeg[1]]
				#
				# finally, we have the next segment. append to thisPoly:
				#thisPoly.append(nextSeg)
				thisPoly.append([(nextSeg[2]*dx)+minX, (nextSeg[3]*dy)+minY])
				thisSeg=nextSeg
			else:
				print('failed to find continuing vector.')
				thisPoly.append([(polyTail[0]*dx)+minX, (polyTail[1]*dy)+minY])
				break
				break
				#
			#
		#
		polygons.append(thisPoly)
		#
	# function root level
				

	# 
	#return rightsies + leftsies + topsies + bottomsies
	#return cdata
	#
	#return rightSegs+leftSegs+topSegs+bottomSegs
	#
	return polygons
	#
	# end function

def simpleContour(data, z0=0, ix=0, iy=1, iz=2, dx=.1, dy=.1):
	# get max x,y values:
	maxX = 0
	maxY = 0
	for rows in data:
		#print rows
		if rows[ix]>maxX:
			maxX = rows[ix]
		if rows[iy]>maxY:
			maxY = rows[iy]
	#print maxX
	#print maxY
	# first, sort y, then x (scan in x):
	data1 = sortXYtuples(data, iy, ix, maxX)
	# scan the sorted data. when we get a state change (aka, step onto or off of a z0 site), set a point.
	cdata = []		# contour data set
	dstate = 0
	lastPos = [data1[0][ix]-dx, data1[0][iy]]
	for rows in data1:
		#  and rows[ix]=lastx+dx and rows[iy]=lasty
		# new region:
		if rows[iz]>=z0 and dstate==0:
			# we've hit a new "high" spot
			dstate=1
			if [rows[ix], rows[iy]] not in cdata:	
				#cdata = cdata + [rows[ix], rows[iy], rows[iz]]
				cdata = cdata + [[rows[ix], rows[iy]]]
		# we fall off the contour onto a data element
		elif (dstate==1 and rows[iz]<z0):
			dstate=0
			if lastPos not in cdata:
				cdata = cdata + [lastPos]
		# we're on a contour and skip over some empty data elements onto a new contour.
		elif (dstate==1 and rows[iz]>=z0 and rows[ix]!=lastPos[0]+dx or rows[iy]!=lastPos[1] ):
			if lastPos not in cdata:
				cdata = cdata + [lastPos]
			if [rows[ix], rows[iy]] not in cdata:	
				cdata = cdata + [[rows[ix], rows[iy]]]
		#
		lastPos = [rows[ix], rows[iy]]

	#
	# now, sort by x then y, repeat with vertical scan, add to cdata:
	data1 = sortXYtuples(data, ix, iy, maxY)
	dstate = 0	
	lastPos = [data1[0][ix], data1[0][iy]-dy]
	# scan the sorted data. when we get a state change (aka, step onto or off of a z0 site), set a point.
	for rows in data1:
		#  and rows[ix]=lastx+dx and rows[iy]=lasty
		# new region:
		if rows[iz]>=z0 and dstate==0:
			# we've hit a new "high" spot
			dstate=1
			if [rows[ix], rows[iy]] not in cdata:	
				#cdata = cdata + [rows[ix], rows[iy], rows[iz]]
				cdata = cdata + [[rows[ix], rows[iy]]]
		# we fall off the contour onto a data element
		elif (dstate==1 and rows[iz]<z0):
			dstate=0
			if lastPos not in cdata:
				cdata = cdata + [lastPos]
		# we're on a contour and skip over some empty data elements onto a new contour.
		elif (dstate==1 and rows[iz]>=z0 and rows[iy]!=lastPos[1]+dy or rows[ix]!=lastPos[0] ):
			if lastPos not in cdata:
				cdata = cdata + [lastPos]
			if [rows[ix], rows[iy]] not in cdata:	
				cdata = cdata + [[rows[ix], rows[iy]]]
		#
	#
	# now we have a contour set. we need to connect these to make a polygon.
	# the way we've done this, lines can connect to nn or nnn lattice sites.
	# start with a random point in cdata (maybe the first point?).
	# it does not matter which direction we go; we just can't backtrack.
	# --> this basic method, but instead employing an index (ie, the xy index)
	# would probably be faster. the yx scan should probably calculate and use
	# the xy index as well...
	#
	#just in case (this might take a long time, so should maybe be skipped)
	#print "starting set"
	#cdata = set(cdata)
	#print  "finishing set"
	#
	return cdata

	#
	# oops... this method will not always outline a polygon...
	'''
	polyPairs = []
	currentPoly = []		# current polygon. eventually, we'll reconnect to the first polygon, then we'll delete these records from cdata.
	rows = cdata[0]
	checkPoint = []
	fromPoint = []
	for rowsis in cdata:
		# look for nn or nnn point:
		foundPoint = 0
		#
		#print rows[0]
		#print rows[0]+dx
		checkPoint = [rows[0]findHeadTailInTuples(thisSeg, rightSegs)[1]>=0+dx, rows[1]]
		if checkPoint in cdata and foundPoint==0 and checkPoint!=fromPoint:
			foundPoint=1
			nextPoint = checkPoint
		#
		checkPoint = [rows[0]+dx, rows[1]+dy]
		if checkPoint in cdata and foundPoint==0 and checkPoint!=fromPoint:
			foundPoint=1
			nextPoint = checkPoint
		#
		checkPoint = [rows[0], rows[1]+dy]
		if checkPoint in cdata and foundPoint==0 and checkPoint!=fromPoint:
			foundPoint=1
			nextPoint = checkPoint
		#
		checkPoint = [rows[0]-dx, rows[1]+dy]
		if checkPoint in cdata and foundPoint==0 and checkPoint!=fromPoint:
			foundPoint=1
			nextPoint = checkPoint
		#
		checkPoint = [rows[0]-dx, rows[1]]
		if checkPoint in cdata and foundPoint==0 and checkPoint!=fromPoint:
			foundPoint=1
			nextPoint = checkPoint
		#
		checkPoint = [rows[0]-dx, rows[1]-dy]
		if checkPoint in cdata and foundPoint==0 and checkPoint!=fromPoint:
			foundPoint=1
			nextPoint = checkPoint
		#
		checkPoint = [rows[0], rows[1]-dy]
		if checkPoint in cdata and foundPoint==0 and checkPoint!=fromPoint:
			foundPoint=1
			nextPoint = checkPoint
		#
		checkPoint = [rows[0]+dx, rows[1]-dy]
		if checkPoint in cdata and foundPoint==0 and checkPoint!=fromPoint:
			foundPoint=1
			nextPoint = checkPoint
		#
		polyPairs = polyPairs + [[rows[0], rows[1], nextPoint[0], nextPoint[1]]]
		rows = [nextPoint[0], nextPoint[1]]
		fromPoint = [rows[0], rows[1]]
	#
	return polyPairs
	'''

#

