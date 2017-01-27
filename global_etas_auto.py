'''
# an implementation of etas_auto() for 1) globalETAS, and 2) in python3
#
'''

import globalETAS as gep
import yodiipy.ANSStools as atp
import datetime as dtm
import pytz
import math
import matplotlib as mpl

lon2km=111.31
colors_ =  mpl.rcParams['axes.color_cycle']
tzutc = pytz.timezone('UTC')
#
days2secs = 60.*60.*24.
year2secs = 60.*60.*24.*365.
deg2km=111.31
deg2rad = 2.0*math.pi/360.
alpha_0 = 2.28
#
# gep.ETAS_mpp_handler_xyz(**np_prams)
def auto_etas(lon_center=None, lat_center=None, d_lat_0=.25, d_lon_0=.5, dt_0=10, Lr_factor=4.0, mc=2.5, mc_0=4.5, d_lat=.1, d_lon=.1, to_dt=None, fnameroot='etas_auto', catlen=5.0*365.0, doplot=False, kmldir='kml_auto', **kwargs):
	my_prams = locals()
	#
	etas = auto_etas_global_params(**my_prams)
	#
	return etas
	
def auto_etas_global_params(lon_center=None, lat_center=None, d_lat_0=.25, d_lon_0=.5, dt_0=10,  mc=2.5, mc_0=4.5, d_lat=.1, d_lon=.1, to_dt=None, fnameroot='etas_auto', catlen=5.0*365.0, doplot=False, kmldir='kml_auto', etas_range_factor=10., etas_range_padding=.5, etas_fit_factor=2.0, transform_type='equal_area',transform_ratio_max=2.5, calc_etas=True, n_contours=15, etas_cat_range=None, etas_xyz_range=None, p_cat=1.1, q_cat=1.5, p_etas=None, Lr_factor=3., **kwargs):
	'''
	# this is probably a major staple of future work in ETAS; when in doubt, start here and point it at a big earthquake.
	#
	# a starter script to auto-select some parameters for an ETAS run. in practice, give this script a center location (probably a mainshock epicenter).
	# the script will find the largest earthquake in that region and scale up an ETAS parameter set accordingly.
		# d_lat_0, d_lon_0, dt_0 are the starting catalog parameters (largest earthquake in d_lat_0 x d_lon_0 x dt_0 cube).
	'''
	#
	#if to_dt == None: to_dt = dtm.datetime.now(pytz.timezone('UTC'))
	to_dt = (to_dt or dtm.datetime.now(pytz.timezone('UTC')))
	mc_0  = (mc_0 or mc)
	#
	if lon_center==None and lat_center==None:
		# let's look for any large earthquake in the world. assume for this, mc
		mc_0=6.0
		lat_center = 0.
		lon_center = 0.
		d_lat_0 = 88.
		d_lon_0 = 180.
	#
	# get a preliminary catalog:
	cat_0 = atp.catfromANSS(lon=[lon_center-d_lon_0, lon_center+d_lon_0], lat=[lat_center - d_lat_0, lat_center+d_lat_0], minMag=mc_0, dates0=[to_dt-dtm.timedelta(days=dt_0), to_dt], fout=None, rec_array=True)
	#
	#biggest_earthquake = filter(lambda x: x['mag']==max(cat_0['mag']), cat_0)[0]
	mainshock = {cat_0.dtype.names[j]:x for j,x in enumerate(list(filter(lambda x: x['mag']==max(cat_0['mag']), cat_0))[0])}
	#
	# now, get new map domain based on rupture length, etc.
	L_r = .5*10.0**(.5*mainshock['mag'] - 1.76)
	delta_lat = Lr_factor*L_r/lon2km
	delta_lon = Lr_factor*L_r/(lon2km*math.cos(deg2rad*mainshock['lat']))
	print("mainshock data: ", mainshock, L_r, delta_lat, delta_lon)
	#
	working_cat = atp.catfromANSS(lon=[mainshock['lon']-delta_lon, mainshock['lon']+delta_lon], lat=[mainshock['lat']-delta_lat, mainshock['lat']+delta_lat], minMag=mc, dates0=[to_dt-dtm.timedelta(days=catlen), to_dt], fout=None, rec_array=True)
	#
	print("biggest event(s): ", [rw for rw in working_cat if rw['mag']==max(working_cat['mag'])])
	#
	# now, do some ETAS:
	# skip working_cat above, but parse lon, lat, etc. parameters similarly. pass those (and other) params to make_etas_fcfiles()
	# looks like this: make_etas_fcfiles(root_prams=nepal_ETAS_prams, **kwargs), and:
	# nepal_ETAS_prams = {'todt':None, 'gridsize':.1, 'contres':5, 'mc':4.5, 'kmldir':kmldir, 'catdir':kmldir, 'fnameroot':'nepal', 'catlen':5.0*365.0, 'doplot':False, 'lons':[nepal_epi_lon-nepal_dlon, nepal_epi_lon+nepal_dlon], 'lats':[nepal_epi_lat-nepal_dlat, nepal_epi_lat+nepal_dlat], 'bigquakes':None, 'bigmag':7.00, 'eqtheta':None, 'eqeps':None, 'fitfactor':5.0, 'cmfnum':0, 'fignum':1, 'contour_intervals':None}
	root_prams = {'t_now':None, 'd_lat':d_lat, 'd_lon':d_lon, 'contres':10, 'mc':mc, 'kmldir':kmldir, 'catdir':kmldir, 'fnameroot':fnameroot, 'catlen':catlen, 'doplot':False, 'lons':[mainshock['lon']-delta_lon, mainshock['lon']+delta_lon], 'lats':[mainshock['lat']-delta_lat, mainshock['lat']+delta_lat], 'etas_range_factor':etas_range_factor, 'etas_range_padding':etas_range_padding, 'etas_fit_factor':etas_fit_factor, 'transform_type':transform_type, 'transform_ratio_max':transform_ratio_max, 'calc_etas':calc_etas, 'n_contours':n_contours, 'etas_cat_range':etas_cat_range, 'etas_xyz_range':etas_xyz_range, 'p_cat':p_cat, 'q_cat':q_cat, 'p_etas':p_etas}
	#
	etas = gep.ETAS_mpp_handler_xyz(**root_prams)
	
	# 'bigquakes':None, 'bigmag':mainshock['mag']-1.5, 'eqtheta':None, 'eqeps':None, 'fitfactor':5.0, 'cmfnum':0, 'fignum':1, 'contour_intervals':None}
	#root_prams = {}
	print("now execute with root_prams: ", root_prams)
	#my_kwargs = {}
	#
	return etas

