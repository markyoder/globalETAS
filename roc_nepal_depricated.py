# depricated code from the globalETAS/nepal paper.


class ROC_base(object):
	# DEPRICATION: most or all of the ROC stuff from this module will eventually be removed. use optimizers.roc_tools.py instead.
	#  note that when we do ROC correctly, we don't get any value from mpp, so while we might write some parallel script to run
	#  multiple ROC analyses, we won't be splitting ROC to multiple processors.
	#
	# a base class object for MPP ROC stuff. we'll inherit this for both the worker and handerl mpp classes.
	# note: this is a ROC tool specific to these map forecast analyses. we should separate the get_site() part and then use
	# the roc_tools() optimizers.roc_tools bits.
	#
	def __init__(self, fc_xyz=None, test_catalog=None, from_dt=None, to_dt=None, dx=None, dy=None, cat_len=120., mc=5.0):
		#
		if isinstance(fc_xyz, str):
			# if we're given a filename...
			with open(fc_xyz, 'r') as froc:
				fc_xyz= [[float(x) for x in rw.split()] for rw in froc if rw[0] not in('#', ' ', '\t', '\n')]
		#
		if to_dt   == None: to_dt = dtm.datetime.now(pytz.timezone('UTC'))
		if from_dt == None: from_dt = to_dt-dtm.timedelta(days=120)			# but this really doesn't make much sense for this version of ROC...
		#
		#return fc_xyz
		if not hasattr(fc_xyz, 'dtype'):
			fc_xyz = numpy.core.records.fromarrays(zip(*fc_xyz), names=('x','y','z'), formats=['>f8', '>f8', '>f8'])
		#
		lats = [min(fc_xyz['y']), max(fc_xyz['y'])]
		lons = [min(fc_xyz['x']), max(fc_xyz['x'])]
		#
		#mc   = mc_roc
		X_set = sorted(list(set(fc_xyz['x'])))
		Y_set = sorted(list(set(fc_xyz['y'])))
		d_lon = (dx or abs(X_set[1] - X_set[0]))
		d_lat = (dy or abs(Y_set[1] - Y_set[0]))
		nx = len(X_set)
		ny = len(Y_set)
		#
		# (for this application, we can also just get nasty and to a loop-loop with geodetic distancing).
		# (maybe ought to set this as a regular class function?)
		#get_site = lambda x,y: int(round((x-lons[0]+.5*d_lon)/d_lon)) + int(round((y-lats[0]+.5*d_lat)/d_lat))*nx
		#
		print("get cataog: ", lons, lats, mc, from_dt, to_dt)
		if test_catalog==None: test_catalog = atp.catfromANSS(lon=lons, lat=lats, minMag=mc, dates0=[from_dt, to_dt])
		print("catlen: ", len(test_catalog))
		#
		# note: this will reduce the number of z values a bit, but what we really need to do is bin them into N bins... or maybe
		# we do explicitly want to remove each site, one at a time to get a proper ROC...
		Zs = sorted(list(fc_xyz['z'].copy()))
		ROCs = [[0,0,0,0]]
		# nominally, we'd make a copy of the whole catalog, but for memory conservation, just index the lsit.
		#eq_ks = [[get_site(eq['lon'], eq['lat']), eq] for eq in test_catalog]
		#eq_site_indices = [get_site(eq['lon'], eq['lat']) for eq in test_catalog]
		#
		self.__dict__.update(locals())
		#
		# ... but running this here won't work properly for mpp instances... well, it will work, but we'll have to set eq_z_vals twice.
		# the first time, when it runs as super(), it will use the local get_site(), which we need to override in the inherited versions.
		self.eq_z_vals = None
		#self.eq_z_vals       = [fc_xyz[self.get_site(eq['lon'], eq['lat'])] for eq in test_catalog]	# ... but not this won't scale up
		#
	#
	def get_site(self, x,y):
		# TODO: this binning needs to be double-checked. do lons[0], lats[0] define the edge of a bin or the boundary in which we
		# select earthquakes (in which case, their outer edges are +/- dx/2.
		# ... but i think this is correct.
		return int(round((x-self.lons[0]+.5*self.d_lon)/self.d_lon)) + int(round((y-self.lats[0]+.5*self.d_lat)/self.d_lat))*self.nx

	def calc_ROCs(self, H_denom=None, F_denom=None, m_c=None):
		#
		if m_c == None:
			m_c = min(self.catalog['mag'])
		#
		if self.eq_z_vals == None or m_c!=None:
			print("setting default eq_z_vals")
			self.eq_z_vals       = [self.fc_xyz['z'][self.get_site(eq['lon'], eq['lat'])] for eq in self.test_catalog if eq['mag']>=m_c]	
			#
		#
		#Zs = sorted(list(fc_xyz['z'].copy()))
		#ROCs = [[0,0,0,0]]
		#eq_site_indices = [get_site(eq['lon'], eq['lat']) for eq in test_catalog]
		Zs   = sorted(self.Zs)
		#ROCs = self.ROCs
		#eq_site_indices = self.eq_site_indices
		fc_xyz = self.fc_xyz

		# might want to specify the denominators for H and F for MPP calculations.
		H_denom = (H_denom or len(self.test_catalog))
		F_denom = (F_denom or len(self.fc_xyz))
		#
		ROCs = [[0,0,0,0]]		# hits, falsies, misses, correct-non-occurrence
		for j_z, z0 in enumerate(Zs):
			# z0 is the threshold z for predicted=True/False
			#
			#for j_eq, eq in enumerate(self.test_catalog):
			#for k in eq_site_indices:
			for k,z in enumerate(self.eq_z_vals):
				# TODO: and parse it up for MPP.
				#
				#k = eq_site_indices[j_eq]			# we can probably optimize this process by just looping through eq_site_indices directly.
				#									# ...might even be able to do an all numpy mod as well. 
				#print('site: ', k)
				#z_val = fc_xyz['z'][k]
				#
				#if z_val>=z0:
				# this whole process could (probably) be sped up by taking a different approach with vectors:
				# hits_mask = [(fc_xyz['z'][k]>=z0) for k in eq_site_indices]
				# hits_sum = sum(hits_mask)
				# false_sum = n_sites - hits_sum
				# miss_sum = n_quakes - hits_sum
				# correct_nonoccurrences = ... dunno... or something like that.
				#if fc_xyz['z'][k]>=z0:
				#print("**debug: ", z, z0)
				if z>=z0:
					# predicted!
					ROCs[-1][0]+=1
					#
					# ... and subtract from falsies; in the end, we'll assume all sites>z were false alarms:
					ROCs[-1][1]-=1
				else:
					# missed it
					ROCs[-1][2]+=1
					ROCs[-1][3]-=1		# an earthquake occurred in this site. it did not correctly predict non-occurrence.
				#
			#
			#n_gt = float(len([z for z in fc_xyz['z'] if z>=z0]))
			#n_lt = float(len([z for z in fc_xyz['z'] if z<z0]))
			
			#n_gt, n_lt= 0, 0
			# this should be a little bit faster (but it needs to be modified to add in place)...
			#[operator.add(ROCs[-1][1],1) if z>=z0 else operator.add(ROCs[-1][3],1) for z in fc_xyz['z']]
			for z in fc_xyz['z']:
				ROCs[-1][1]+=(z>=z0)
				ROCs[-1][3]+=(z<z0)
			#	n_gt += (z>=z0)
			#	n_lt += (z<z0)
			##
			#ROCs[-1][1]+=n_gt
			#ROCs[-1][3]+=n_lt
			#
			ROCs += [[0,0,0,0]]
		#
		self.ROCs=ROCs
		#
		FH = [[rw[1]/F_denom, float(rw[0])/H_denom] for rw in ROCs[:-1]]
		#
		self.FH = FH
	#
	def plot_HF(self, fignum=0, do_clf=True):
		plt.figure(fignum)
		if do_clf: plt.clf()
		#
		if 'FH' not in self.__dict__.keys():
			self.calc_ROCs()
		#
		plt.plot(*zip(*self.FH), ls='-', label='ROC_approx.', lw=2., alpha=.8)
		#plt.plot(Fs2, Hs2, '-', label='ROC', lw=2., alpha=.8)
		plt.plot(range(2), range(2), 'r--', lw=2.5, alpha=.6)



class ROC_mpp_handler(ROC_base):
	# you know, what we should start doing is write up these "handler" classes as mpp.Pool() subclasses...
	# ... or write up like:
	# def __init__(self, n_procs, *args, **kwargs):
	#	super().__init__(*args, **kwargs):
	#def __init__(self, fc_xyz=None, test_catalog=None, from_dt=None, to_dt=None, dx=None, dy=None, cat_len=120., mc=5.0, n_procs=None):
	def __init__(self, n_procs=None, *args, **kwargs):
		# note that the call signature using *args (ordered parameters) won't be identical to non-mpp subclasses; kwargs or key-word calling syntax won't change.
		# call signature will be like:
		# ROC_mpp_handler(n_procs, fc_xyz=None, test_catalog=None, from_dt=None, to_dt=None, dx=None, dy=None, cat_len=120., mc=5.0)
		# (all prams after n_procs are from ROC_base)
		#
		super(ROC_mpp_handler,self).__init__(*args,**kwargs)
		# ... also might need to call mpp.Process.__init__() explicitly.
		#
		if n_procs==None: n_procs = min(1, mpp.cpu_count()-1)		# we should also catch the case n_procs==1 and defer to an spp solution.
		self.n_procs = n_procs
		#
		# hit and false-alarm denominators. we'll be distributing ROC to multiple processors, so we need to know these in advance.
		self.h_denom = float(len(self.test_catalog))
		self.f_denom = float(len(self.fc_xyz))
		#
		# for now, set the z_indices here and parse them to the workers. we'll sort out indexing to do this on the worker level later
		# (so we'll pass less data and distribute more work, but for now just make it work.)
		self.eq_z_vals = None
	#
	def set_eq_z_vals(self, m_c=None):
		if m_c==None:
			m_c = min(self.test_catalog['mag'])
		self.eq_z_vals = [self.fc_xyz['z'][self.get_site(eq['lon'], eq['lat'])] for eq in self.test_catalog if eq['mag']>=m_c]	
		#
	#
	def calc_ROCs(self, n_procs=None, m_c=None):
		n_procs = (n_procs or self.n_procs)
		# ((this is still totally ****** in mpp). n>0 nodes are not getting proper indexing.
		#
		# now, split up fc_xyz into n_procs, start each proc and pipe back the results.
		# (this logic tree needs some work. we want to be able to change the eq_z_values, but maybe not all the time?)
		# this approach (if m_c!=None) potentially produces some undesirable behavior; if m_c hasn't changed, we re-calc eq_z_vals, which
		# might take some time
		if self.eq_z_vals == None or m_c!=None:
			print('re-calculating eq_z_vals for m_c=%f' % m_c)
			self.set_eq_z_vals(m_c=m_c)
		#
		self.FH=[]
		fc_len = int(numpy.ceil(len(self.fc_xyz)/n_procs))
		roc_workers = []
		roc_pipes = []
		Zs_fc = self.fc_xyz['z']
		Zs_fc.sort()
		print("calculating ROCs: %d processors" % (int(n_procs)))
		for j in range(n_procs):
			# make an ROC_worker instance for each process, fc_xyz = fc_xyz[j*fc_len:(j+1)*fc_len]
			p1,p2 = mpp.Pipe()
			#
			roc_workers += [ROC_worker_simple(Z_fc=Zs_fc[j*fc_len:(j+1)*fc_len], Z_events=self.eq_z_vals, h_denom=len(self.eq_z_vals), f_denom=len(Zs_fc), f_start = j*fc_len, pipe_r=p1)]
			
			#
			roc_pipes += [p2]
			roc_workers[-1].start()
			#
		#
		FH = []
		Hs=[]
		Fs=[]
		plt.clf()
		for j,p in enumerate(roc_pipes):
			# ... and i always forget if there's still any value in doing a join() loop, noting that (as i recall), the join()s must follow the recv()s.
			# ... and we should probalby rewrite this usint Pool(), where we can still use ROC_worker() objects, but i think it handles the sequence
			# in which workers start/finish more efficiently.
			#ary = p.recv()
			#
			###
			#fh = list(roc_pipes[j].recv())
			#plt.plot(*zip(*fh), marker='.', ls='')
			#FH += fh
			###
			# this looks right, but it needs to be verified.			
			Hs += list(roc_pipes[j].recv())
			roc_pipes[j].close()
		#
		l_F = len(Zs_fc)
		l_h = len(self.test_catalog)
		#
		# let's handle the two cases: 1) we return just H; 2) we return FH
		if not hasattr(Hs[0], '__len__'):
			#FH = list(zip([(l_F-float(j+1))/l_F for j in range(l_F)], Hs))
			FH = list(zip([(l_F-float(j)-Hs[j]*l_h)/l_F for j in range(l_F)], Hs))
		else:
			FH.sort(key=lambda x: x[0])

		self.FH = FH
	#
class ROC_worker(ROC_base, mpp.Process):
	# worker class for ROC calculations.
	#def __init__(self, fc_xyz, test_catalog=None, from_dt=None, to_dt=None, dx=None, dy=None, cat_len=120., mc=5.0, pipe_r=None, H_denom=None, F_denom=None):
	def __init__(self, pipe_r=None, H_denom=None, F_denom=None, eq_z_vals=None, *args, **kwargs):
		# note: the __init__ should be satisfied with fc_xyz, test_catalog. the mpp_handler or some other autonomouse form
		# should automate catalog fetching and etas running.
		#super(ROC_worker,self).__init__(**{key:val for key,val in locals().items() if key not in ('self', '__class__', 'pipe_r', 'H_denom', 'F_denom')})
		super(ROC_worker,self).__init__(*args, **kwargs)
		mpp.Process.__init__(self)
		self.pipe_r=pipe_r
		self.H_denom = H_denom
		self.F_denom = F_denom
		#self.z_start_index = z_start_index
		#
		#self.set_eq_z_vals(eq_z_vals)
		self.eq_z_vals = eq_z_vals
	#
	#def get_site(self, x,y,z_start_index=None):
	#	# overriding base class definition to include z_start_index parameter, so indexing will translate correctly to n>0 mpp nodes.
	#	z_start_index = int((z_start_index or self.z_start_index) or 0)
	#	#print("z_start_index = %d " % z_start_index)
	#	#
	#	return int(round((x-self.lons[0]+.5*self.d_lon)/self.d_lon)) + int(round((y-self.lats[0]+.5*self.d_lat)/self.d_lat))*self.nx - z_start_index
	def set_eq_z_vals(self, vals=None):	
		if vals==None:
			self.eq_z_vals = [self.fc_xyz[self.get_site(eq['lon'], eq['lat'])] for eq in self.test_catalog]
		else:
			self.eq_z_vals = vals
		#
	def run(self):
		# we can probably speed this up a bit by eliminating some of the copies...
		self.calc_ROCs(H_denom=self.H_denom, F_denom=self.F_denom)
		#self.calc_HF()
		#
		self.pipe_r.send(self.FH)
		self.pipe_r.close()
#
class ROC_worker_simple(mpp.Process):
	# super lean ROC worker. i think we don't need to inherit from other ROC classes, in fact we might re-write to base from this class.
	# next time, let's start over and just write an ROC analyzer... which we may have done in vc_parser... but either suck it up and 
	# do the guts in C++/extension or use shared memory arrays.
	def __init__(self, Z_fc=None, Z_events=None, h_denom=None, f_denom=None, f_start = 0., pipe_r=None):
		super(ROC_worker_simple,self).__init__()
		self.__dict__.update({key:val for key,val in locals().items() if not key in ('self', 'class')})
		self.F=[]
		self.H=[]
		self.pipe_r=pipe_r
		
		#
	#def roc_simple_sparse(Z_fc=None, Z_events=None, h_denom=None, f_denom=None, f_start = 0.):
	def roc_simple_sparse(self, h_denom=None, f_denom=None, f_start=0):
		# simple ROC calculator.
		h_denom = (h_denom or self.h_denom)
		f_denom = (f_denom or self.f_denom)
		Z_events = self.Z_events
		Z_fc = self.Z_fc		# these just make a reference, but it still writes a copy of a variable, so it might be a bit faster to just reff self.x, etc.
		#
		h_denom = (h_denom or len(Z_events))
		f_denom = (f_denom or len(Z_fc))
		#
		# ... though we might modify this, or write a different version that just returns the hits, which would be faster for piping back in mpp.
		# ... actually, this isn't quite right; it assumes that all z values are unique (which they are in some cases, but it's not an accurate assumption).
		#
		# more precuse, but not quite as fast:
		r_XY = []
		for j,z_fc in enumerate(Z_fc):
			n_f = sum([(z_ev>=z_fc) for z_ev in Z_events])
			r_XY += [[(f_start + j - n_f)/f_denom, n_f/h_denom]]
		#
		return r_XY
		#
		# a minor approximation when not sparse:
		#return [[(f_start + j)/f_denom, sum([(z_ev>=z_fc) for z_ev in Z_events])/h_denom] for j,z_fc in enumerate(Z_fc)]
	#
	def roc_hits_sparse(self, h_denom=None):
		h_denom = (h_denom or self.h_denom)
		# return just the hits. for mpp operations, this might be faster, since the Fs array is just [j/f_denom for j,x in enumerate(z_fc)] (or equivalent).
		# it'll taka longer to pipe it back than to calc it.
		# ... and note, this isn't quite right; it assumes that all z values are unique... which they may not be. that said, it is nominally a reasonable 
		# approximation and it should be much faster than if we actually check z_fc like f = sum([z>=z0 for z in Z_fc]) and easier than if we have to 
		# 'manually' check for unequal values.
		return [sum([float(z_ev>=z_fc) for z_ev in self.Z_events])/h_denom for j,z_fc in enumerate(self.Z_fc)]
	#
	def run(self):
		# note for different roc calcs (sparse, exact, hits only, etc. we can either pass a pram or create special subclasses).
		self.pipe_r.send(self.roc_hits_sparse())
		#self.pipe_r.send(self.roc_simple_sparse(h_denom=self.h_denom, f_denom=self.f_denom))
		self.pipe_r.close()

#
def roc_class_test_1(n_cpus=None):
	'''
	# for mpp ROC calcs, make this script work. also, check spp scripts by uncommenting the SPP
	# bits.
	'''
	# test roc etas classes.
	# for now, just use the nepal ETAS classes:
	#
	if n_cpus==None: n_cpus = max(2,min(1, mpp.cpu_count()-1))
	mc_roc = 5.0
	#

	etas_fc = get_nepal_etas_fc(n_procs=n_cpus, cat_len=100.)
	#etas_fc = get_nepal_etas_test()
	# if we want to load from file, we'll need some of these data:
	# {'bigmag': 7.0,'bigquakes': None, 'catdir': 'kml', 'catlen': 1825.0, 'cmfnum': 0, 'contour_intervals': None, 'contres': 5, 'doplot': False, 'eqeps': None,'eqtheta': None, 'fignum': 1, 'fitfactor': 5.0, 'fnameroot': 'nepal', 'gridsize': 0.1, 'kmldir': 'kml', 'lats': [23.175, 33.175], 'lons': [79.698, 89.698],'mc': 3.5, 't_now': None}
	#
	fc_date = max(etas_fc.catalog['event_date']).tolist()
	fc_date = dtm.datetime(2015,5,7,tzinfo=pytz.timezone('UTC'))
	to_dt = fc_date + dtm.timedelta(days=120)
	test_catalog = atp.catfromANSS(lon=etas_fc.lons, lat=etas_fc.lats, minMag=mc_roc, dates0=[fc_date, to_dt])
	#
	xyz = etas_fc.ETAS_array.copy()	# might not need to copy, but i don't think this is usually a significant performance issue.
	#
	h_denom = float(len(test_catalog))
	f_denom = float(len(xyz))
	'''
	#roc = ROC_base(fc_xyz=xyz, test_catalog=test_catalog, from_dt=fc_date, to_dt=to_dt, dx=None, dy=None, cat_len=120., mc=5.0)
	roc = ROC_base(fc_xyz=xyz, test_catalog=test_catalog)
	#
	fh_obj = roc.calc_ROCs(H_denom=h_denom, F_denom=f_denom)
	roc.plot_HF()
	#
	#ZZ = roc_normal_from_xyz(fc_xyz=xyz, test_catalog=test_catalog, from_dt=fc_date, to_dt=to_dt, dx=None, dy=None, cat_len=120., mc=5.0, fignum=2, do_clf=True)
	fh_fxyz = roc_normal_from_xyz(fc_xyz=xyz, test_catalog=test_catalog, fignum=2, do_clf=True)
	#
	fh_fetas = roc_normal(etas_fc=etas_fc, test_catalog=test_catalog, mc_roc=5.0, fignum=4, do_clf=True)
	#
	'''
	#
	roc_mpp = ROC_mpp_handler(fc_xyz=xyz, test_catalog=test_catalog, n_procs = n_cpus)
	#roc_mpp = ROC_mpp_handler(fc_xyz='data/nepal_etas_xyz.csv', test_catalog=test_catalog, n_procs = n_cpus)
	#fh_mpp = roc_mpp.calc_ROCs(H_denom=h_denom, F_denom=f_denom)
	fh_mpp = roc_mpp.calc_ROCs()
	roc_mpp.plot_HF(fignum=6)
	#
	FH_spp = roc_normal_from_xyz(fc_xyz=etas_fc.ETAS_array, test_catalog=test_catalog)
	#
	plt.figure(1)
	plt.clf()
	plt.plot(*zip(*roc_mpp.FH), ls='--', marker='.')
	plt.plot(*zip(*FH_spp), ls='-', marker='o')
	
	
	#return (roc, fh_obj, fh_fxyz, fh_fetas, fh_mpp)
	return roc_mpp, fh_mpp
#
# #########################
# Comopnents and helper scripts, in various states of repair...


def global_roc3(fc_xyz='global/global_xyz_20151129.xyz', n_cpu=None, fnum=0, m_cs=[4.0, 5.0, 6.0, 6.5], test_catalog=None, fc_start_date=None, fc_end_date=None, cat_len=120, fc_frac=1.0, fout='global_roc3.png'):
	# ok, so all of this is a mess. the roc bits need to be consolidated, cleaned up, and better modularized. we'll do some of that here, at lest by example.
	# @fc_frac: fraction of fc sites to use (aka, skip the first (1-frac)*N (low z) sites.
	# fetch the fc_xyz, generate a test catalog, then use generic_roc tools (mpp implementations) to do the roc.
	# with the test catalog: 1) get the test catalog for the region, time domain (maybe write  script to do this)
	# then, get the fc_z-values, but keep those as pairs like [[z, mag], ...]. this way, we can quickly make z_events lists.
	#
	# note: use optimizers.roc_tools ROC calculator(s). also, see the new notebook script for nepal ROC.
	#
	n_cpu = (n_cpu or mpp.cpu_count())
	#
	# for now, with this script, we're assuming that we are using this specific file, but we might pre-load it.
	if isinstance(fc_xyz, str):
		with open(fc_xyz, 'r') as froc:
			fc_xyz= [[float(x) for x in rw.split()] for rw in froc if rw[0] not in('#', ' ', '\t', '\n')]
		#
	#
	#return fc_xyz
	fc_start_date = (fc_start_date or dtm.datetime(2015,11,30, tzinfo=pytz.timezone('UTC')))
	fc_end_date   = (fc_end_date or fc_start_date + dtm.timedelta(days=cat_len))
	#
	if not hasattr(fc_xyz, 'dtype'):
		fc_xyz = numpy.core.records.fromarrays(zip(*fc_xyz), names=('x','y','z'), formats=['>f8', '>f8', '>f8'])
	#
	#lons=[-180., 180.]
	#lats=[-89., 89.]
	mc=min(m_cs)
	#
	X_set = sorted(list(set(fc_xyz['x'])))
	Y_set = sorted(list(set(fc_xyz['y'])))
	nx = len(X_set)
	ny = len(Y_set)
	lons = [min(X_set), max(X_set)]
	lats = [min(Y_set), max(Y_set)]
	d_lon = abs(X_set[1] - X_set[0])
	d_lat = abs(Y_set[1] - Y_set[0])
	get_site = lambda x,y: int(round((x-lons[0]+.5*d_lon)/d_lon)) + int(round((y-lats[0]+.5*d_lat)/d_lat))*nx
	#	
	print("get cataog: ", lons, lats, mc, fc_start_date, fc_end_date)
	if test_catalog==None:
		test_catalog = atp.catfromANSS(lon=lons, lat=lats, minMag=min(m_cs), dates0=[fc_start_date, fc_end_date])
	print("catlen: ", len(test_catalog))
	#
	# forecast z-values.
	#Zs = sorted(list(fc_xyz['z'].copy()))
	Zs = sorted(list(fc_xyz['z']))
	#
	# now, get both the eq magnitudes and site_z values, so we can change the mc threshold later.
	eq_site_zs = [[fc_xyz['z'][get_site(eq['lon'], eq['lat'])], eq['mag']] for eq in test_catalog]
	#eq_site_zs.sort(key=lambda x: x[1])
	#
	plt.figure(fnum)
	plt.clf()
	plt.plot(range(2), range(2), ls='--', color='m', lw=3., alpha=.75, zorder=2)
	FHs = {}		# we'll use mc as a key, FH as a val: {mc:[FH]...}
	#				# it won't key well because of floating point error (aka, FHs[5.5] will not be reliable. but it will make a decent container.
	#
	for j,mc in enumerate(m_cs):
		print('doing ROC for mc=%f' % mc)
		#
		# n_procs,Z_events, Z_fc, h_denom=None, f_denom=None, f_start=0., f_stop=None
		roc = roc_generic.ROC_mpp(n_procs=n_cpu, Z_events=[z for z,m in eq_site_zs if m>=mc], Z_fc=Zs, h_denom=None, f_denom=None, f_start=0, f_stop=None)
		a=roc.calc_roc()		# no parameters, and in fact no return value, but it's never a bad idea to leave a place-holder for one.
		#
		clr = colors_[j%len(colors_)]
		plt.plot(roc.F, roc.H, ls='-', color=clr, marker='', lw=2.5, label='$m_c=%.2f$' % mc)
		#FHs[mc]=[[f,h] for f,h in zip(roc.F, roc.H)]
		#
		plt.show()	# just in case...
	bins, mins, maxes = roc_random(n_events=100, n_fc=10000, n_rocs=100, n_cpus=None, ax=plt.gca(), n_bins=100, line_color='m', shade_color='m')
	#plt.fill_between(bins, mins, maxes, color='m', alpha=.3)
	#plt.plot(bins, mins, '-', lw=1.5, alpha=.8)
	#plt.plot(bins, maxes, '-', lw=1.5, alpha=.8)
	plt.legend(loc=0, numpoints=1)
	plt.title('Global ROC', size=18)
	plt.xlabel('False Alarm Rate $F$', size=18)
	plt.ylabel('Hit Rate $H$', size=18)
	plt.savefig(fout)
	#
	#
	return FHs
#
def global_roc_comparison(fc_xyz='global/global_xyz_20151129.xyz', n_cpu=None, fnum=0, mc=6.0, roc_fracs=[1.0, .8, .5], test_catalog=None, fc_start_date=None, fc_end_date=None, fc_frac=1.0):
	#
	# (turns out that i don't think we really need this script. its apparent value appears to have resulted from an error
	# in the global_roc scripts.
	n_cpu = (n_cpu or mpp.cpu_count())
	#
	# for now, with this script, we're assuming that we are using this specific file, but we might pre-load it.
	if isinstance(fc_xyz, str):
		with open(fc_xyz, 'r') as froc:
			fc_xyz= [[float(x) for x in rw.split()] for rw in froc if rw[0] not in('#', ' ', '\t', '\n')]
		#
	#
	fc_start_date = (fc_start_date or dtm.datetime(2015,11,30, tzinfo=pytz.timezone('UTC')))
	fc_end_date   = (fc_end_date or fc_start_date + dtm.timedelta(days=120))
	#
	if not hasattr(fc_xyz, 'dtype'):
		fc_xyz = numpy.core.records.fromarrays(zip(*fc_xyz), names=('x','y','z'), formats=['>f8', '>f8', '>f8'])
	#
	X_set = sorted(list(set(fc_xyz['x'])))
	Y_set = sorted(list(set(fc_xyz['y'])))
	nx = len(X_set)
	ny = len(Y_set)
	lons = [min(X_set), max(X_set)]
	lats = [min(Y_set), max(Y_set)]
	d_lon = abs(X_set[1] - X_set[0])
	d_lat = abs(Y_set[1] - Y_set[0])
	get_site = lambda x,y: int(round((x-lons[0]+.5*d_lon)/d_lon)) + int(round((y-lats[0]+.5*d_lat)/d_lat))*nx
	#	
	print("get cataog: ", lons, lats, mc, fc_start_date, fc_end_date)
	if test_catalog==None:
		test_catalog = atp.catfromANSS(lon=lons, lat=lats, minMag=mc, dates0=[fc_start_date, fc_end_date])
	print("catlen: ", len(test_catalog))
	#
	# forecast z-values.
	#Zs = sorted(list(fc_xyz['z'].copy()))
	Zs = sorted(list(fc_xyz['z']))
	#
	# now, get both the eq magnitudes and site_z values, so we can change the mc threshold later.
	#eq_site_zs = [[Zs[get_site(eq['lon'], eq['lat'])], eq['mag']] for eq in test_catalog if eq['mag']>=mc]
	eq_site_zs = [fc_xyz['z'][get_site(eq['lon'], eq['lat'])] for eq in test_catalog if eq['mag']>=mc]
	#eq_site_zs.sort(key=lambda x: x[1])
	#
	plt.figure(fnum)
	plt.clf()
	plt.plot(range(2), range(2), ls='--', color='r', lw=3., alpha=.75, zorder=2)
	#
	for j,frac in enumerate(roc_fracs):
		#
		print('calcing roc for mc=%f, frac=%f' % (mc, frac))
		j_roc = int(len(fc_xyz)*(1.-frac))
		Zs = sorted(list(fc_xyz['z'][j_roc:].copy()))
		#
		# n_procs,Z_events, Z_fc, h_denom=None, f_denom=None, f_start=0., f_stop=None
		roc = roc_generic.ROC_mpp(n_procs=n_cpu, Z_events=eq_site_zs, Z_fc=Zs, h_denom=None, f_denom=None, f_start=0, f_stop=None)
		a=roc.calc_roc()		# no parameters, and in fact no return value, but it's never a bad idea to leave a place-holder for one.
		#
		plt.plot(roc.F, roc.H, ls='-', label='$m_c=%f, frac=%f' % (mc, frac))
	#
	plt.title('ROC Comparison')
	plt.xlabel('False Alarm Rate $F$')
	plt.ylabel('Hit Rate $H$')
#

def nepal_roc_script_a():
	#
	# i think this script works correctly, but slowly. see the notebook script (which will probably be moved to nepal_figs.py before too long)
	# this script will be moved to the deprication, discard holding module...
	#
	
	etas_fc = etas_analyzer.get_nepal_etas_fc()
	#nepal_etas_test = get_nepal_etas_test()
	Xs = sorted(list(set(etas_fc.ETAS_array['x'])))
	Ys = sorted(list(set(etas_fc.ETAS_array['y'])))
	get_site = lambda x,y: int(numpy.floor((x-lons[0])/d_lon)) + int(numpy.floor((y-lats[0])/d_lat))*nx
	
	print('do real nepal roc')
	A=etas_analyzer.roc_normalses(etas_fc, test_catalog=None, to_dt=None, cat_len=120., mc_rocs=[4.0, 5.0, 6.0, 7.0], fignum=1, do_clf=True, roc_ls='-')
	#
	
	# now, get roc for a 1/r map (i think there's a script for that)
	etas_toy = etas_analyzer.Toy_etas_invr(etas_in=etas_fc, mainshock=nepal_mainshock)
	r0 = 10.**(.5*7.8-1.76)
	x0=nepal_epi_lon
	y0=nepal_epi_lat
	#for j,(x,y,z) in enumerate(etas_fc.ETAS_array): ETAS_array['z'][j]=1/(dist_to(x,y,x0,y0) + r0)
	etas_toy.d_lat=etas_fc.d_lat
	etas_toy.d_lon=etas_fc.d_lon
	B=etas_analyzer.roc_normalses(etas_toy, test_catalog=etas_fc.catalog, to_dt=None, cat_len=120., mc_rocs=[4.0, 5.0, 6.0, 7.0], fignum=1, do_clf=False, roc_ls='--') 
	ax=plt.gca()
	#
	#
	# and draw in roc for random...
	print('do toy, 1/r roc')
	bins, mins, maxes = roc_random(n_events=100, n_fc=10000, n_rocs=100, n_cpus=None, ax=ax, n_bins=100, line_color='m', shade_color='m')
	plt.draw()
#

