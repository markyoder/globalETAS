import math
import numpy

def findClustersSquare(gridlist=[], L=None):
	# find clusters and their size (for square lattice).
	# two passes total (i think) through the data:
	# - walk gridk "left" to "right" (aka, i -> i+1). L/R adjacent elements are elements of a cluster; record clusters in clustList=[[j0, [indeces]], [j1, [ind.]] ]
	# - if there is an element "below" (aka z[i-L]):
	#		- z[i] -> that cluster
	#		- members of z[i]'s cluster are tagged for renaming (add an entry to clustMap -> [j, j']. use new clust number going forward.
	# - after walking the grid, update clustList indeces from top to bottom according to clustMap, also top to bottom.
	#
	mygridlist=scipy.array(gridlist)
	clusterIDs=[]	# will be same dim. as gridlist, but will be populated with clusterID indices.
	clustList=[[0, []]]	#[[j0, [indices for j0]], [j1, [indeces for j1]...]	# and we might not need the index name j0, j1, ... we can probably just use order as index.
	clustMap=[]		# when we join two clusters: [[j0, j0'], [j1, j1'], ...]
	#
	# 1D or  2D array?
	#if type(gridlist[0]).__name__ in ('float', 'int')
	# more robustly?
	'''
	if type(gridlist[0])==type(5.0) or type(gridlist[0])==(5):
		# in case type names change in later pythons? aka, float64, intlong, or whatever. of course, this dos not guarantee that the flavor
		# of float/int we provide is the same.
		#
		# anyway, this is a 1D array. infer 2D from L.
		# let's just cast this into a 2D array (it's got to be one or the other).
		mygridlist.shape=(len(gridlist)/L, L)	# assuming config. [x0,y0; x1, y0; x2,y0... we should also trap for not integer division, etc.
	#
	# otherwise, it's a 2D array (or something that will make a big mess and break or produce nonsense).
	'''
	# but i seem to be happier with 1D array management, so let's do it this way:
	if type(gridlist[0]).__name__ in ('tuple', 'list', 'ndarray'):
		# elemtens are list-like, aka a 2+D array. assume 2D.
		L=len(gridlist[0])
		mygridlist.shape=(len(mygridlsit))
	#
	gridlen=len(gridlist)
	if L==None:
		# before giving up, guess a square:
		squareL=((float(gridlen))**.5)
		if squareL%1==0:
			# it could be a square...
			L=squareL
		else:
			# we don't know how to parse it.
			print("array width not defined.")
			return None
	#
	clustIndex=0
	#
	i=0
	print("L=%d" % L)
	while i<gridlen:
		if mygridlist[i]==None or mygridlist[i]<=0:
			clusterIDs+=[None]
			i+=1
			continue	# though we might just use NONE and to allow for 0-value perimeter padding.
		#
		# at this point, we've found an occupied site.
		#
		# is it a (potentially) new cluster on the left edge?
		#if i%L==0:
		# for starters, assume it is:
		myclustID=clustList[-1][0]+1
		if i%L>0:
			# not on left edge
			if clusterIDs[i-1]!=None:
				# join cluster to left
				myclustID=clusterIDs[i-1]
				
		if i>=L:		# a more efficient (aka, not always checking i>=L) approach might be to process the first row separately.
			# we've passed the first row. is it part of a prev-row cluster?
			if clusterIDs[i-L]!=None:
				myclustID=clusterIDs[i-L]	# join the cluster below.
				# was it part of a cluster to the left?
				if clusterIDs[i-1]!=None and clusterIDs[i-1]!=myclustID:
					# i suppose we could change this guy while we're here, but it is not necessary. write the label change:
					#clustMap+=[[myclustID, clusterIDs[i-1]]] #[[i_from, i_to]] as executed from the top of the list.
					clustMap+=[[clusterIDs[i-1], myclustID]] #[[i_from, i_to]] as executed from the top of the list.
				
				#
			#
		#
		# write to cluterIDS:
		clusterIDs+=[myclustID]
		#
		# add this element to its cluter list.
		if len(clustList)<(myclustID+1):
			# we've added a new cluster.
			clustList+=[[myclustID, []]]
		clustList[myclustID][1]+=[i]
		#
		i+=1
	# now, stitch together the clusters:
	newclusts={}
	for i in range(len(clustList)):
		newclusts[i]=i
		
	for rw in clustMap:
		srcindex=newclusts[rw[0]]
		targetindex=newclusts[rw[1]]
		#
		#clustList[rw[1]][1]+=clustList[rw[0]][1]
		#clustList[rw[0]][1]=[]
		clustList[targetindex][1]+=clustList[srcindex][1]
		clustList[srcindex][1]=[]
		newclusts[srcindex] = targetindex
	
	#return [clustList, clustMap, clusterIDs]
	# return clusters.
	rlist=[]
	for rw in clustList:
		if len(rw[1])==0: continue
		rlist+=[rw[1]]
	return rlist

def getNeighbors(i, nxyz=[], geom=0, bounds=0):
	# nxyz: [width, length, depth, etc.]
	if geom==0 or geom=='square' or geom==None:
		# square geom.
		X=nxyz[0]
		Y=nxya[1]
		L=X*Y
		#
		dns=[-1, -X, 1, X]
		#
		neighbs=[]
		for dn in dns:
			# quasi-heliacal BC. note, to do proper periodic BC, we first calc. the x,y coords.
			newi=i+dn
			if newi<0: newi+=L
			neighbs+=[newi%L]
			#
		#
		
	return neighbs
		
		
		
		
		
		
