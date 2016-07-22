
import optimizers
import optimizers.roc_tools as ortp
import numpy
import pylab as plt
import globalETAS as gep
from yodiipy import ANSStools as atp
import datetime as dtm
import pytz
#
#
class ETAS_roc(ortp.ROC_xyz_handler):
    #
    # class to construct an optimizers.ROC_xyz_handler() from etas and etas related inputs.
    # as needed, load or compute ETAS, load catalog.
    #def __init__(self, fc_xyz, events_xyz=None, dx=None, dy=None, fignum=0, do_clf=True, z_event_min=None, z_events_as_dicts=False):
    def __init__(self, etas_xyz=None, events_xyz=None, lats=None, lons=None, z_event_min=None,
                 fc_date=None, cat_len=5*365, fc_len=120, d_lat=.1, d_lon=.1, n_procs=None, mc_etas=2.5,
                 mc_catalog=None, etas_range_factor=25, etas_range_padding=1.0, transform_type='equal_area',
                 transform_ratio_max=2.5, etas_fit_factor=1.0, *args, **kwargs):
        #
        mc_catalog = (mc_catalog or mc_etas)
        #
        # first handle the ETAS (forecast) bits.
        if etas_xyz is None or (hasattr(etas_xyz,'__len__') and not isinstance(etas_xyz,str) and len(etas_xyz[0])>3): 
            # get compute ETAS from input parameters.
            # ... eventually.
            # globalETAS call signature: def __init__(self, catalog=None, lats=None, lons=None, mc=2.5, mc_etas=None,
            #d_lon=.1, d_lat=.1, bin_lon0=0., bin_lat0=0., etas_range_factor=10.0, etas_range_padding=.25, etas_fit_factor=1.0,
            #t_0=dtm.datetime(1990,1,1, tzinfo=tz_utc), t_now=dtm.datetime.now(tzutc), transform_type='equal_area', 
            #transform_ratio_max=2.5, cat_len=5.*365., 
            #calc_etas=True, n_contours=15, etas_cat_range=None, etas_xyz_range=None, p_cat=1.1, q_cat=1.5, p_etas=None,**kwargs)
            #
            # can we figure out if a list/array passed is a catalog? maybe if len(X[0]>3, then it's probably a catalog?)
            #
            if not etas_xyz is None: catalog = etas_xyz

            etas = ETAS_mpp_handler(n_processes=n_procs, catalog=catalog, lats=lats, lons=lons, mc=mc_etas, d_lat=d_lat,
                                    d_lon=d_lon, etas_range_factor=etas_range_factor, etas_range_padding=etas_range_padding,
                                    etas_fit_factor=etas_fit_factor, t_now=fc_date, cat_len=cat_len, 
                                    transform_type=transform_type, transform_ratio_max=transform_ratio_max)
            #etas_xyz = etas.ETAS_array.tolist()
            etas_xyz = etas.ETAS_array
            #
        if isinstance(etas_xyz, str):
            # it's a file; load it. assume x,y,z format.
            with open(etas_xyz, 'r') as f:
                etas_xyz = [[float(x) for x in rw.split()] 
                            for rw in f if rw[0] not in (chr(9), chr(32), chr(10), chr(13), '#')]
            #
        #
        # now, let's handle events. gather some lats, lons from the ETAS input. we'll need to go over the default
        # behaviors a few times until we get it right.
        #
        #lons = (lons or [min(x for x,y,z in etas_xyz), max(x for x,y,z in etas_xyz)])
        #lats = (lats or [min(y for x,y,z in etas_xyz), max(y for x,y,z in etas_xyz)])
        Xs = sorted(list(set(x for x,y,z in etas_xyz)))
        Ys = sorted(list(set(y for x,y,z in etas_xyz)))
        lons = (lons or [min(Xs), max(Xs)])
        lats = (lats or [min(Ys), max(Ys)])
        d_lon = abs(Xs[1] - Xs[0])
        d_lat = abs(Ys[1] - Ys[0])
        #
        if events_xyz is None:
            # ... and maybe other criteria? load a catalog from ANSS
            anss_cat = atp.catfromANSS(lon=[lons[0]-d_lon/2., lons[1]+d_lon/2.], lat=[lats[0]-d_lat/2., lats[1]+d_lat/2.], minMag=mc_etas, 
                                   dates0=[fc_date, fc_date + dtm.timedelta(days=cat_len)],
                                   Nmax=None, fout=None, rec_array=True)
            self.anss_cat = anss_cat
        if isinstance(events_xyz,str):
            with open(events_xyz, 'r') as f:
                events_xyz = [[float(x) for x in rw.split()] 
                              for rw in f if rw[0] not in (chr(9), chr(32), chr(10), chr(13), '#')]
            #
        if not hasattr(events_xyz, 'dtype'):
            events_xyz = numpy.core.records.fromarrays((anss_cat['lon'], anss_cat['lat'], anss_cat['mag']), 
                                                  dtype=[('x', 'double'), ('y', 'double'), ('z', 'double')])
            
    
        super(ETAS_roc, self).__init__(etas_xyz, events_xyz, z_event_min=z_event_min)

def roc_tools_test_nepal_global():
	# run test ROC on Nepal and global sets.
	#
	etas_input = 'data/nepal_etas_xyz.csv'   # this is the old etas, with the 1/area instead of 1/L error... and it 
	#                                        # happens to produce a slightly better ROC... so we should look into
	#                                        # steeper spatial scaling exponents.
	fcdt = dtm.datetime(2015,5,7, tzinfo=pytz.timezone('UTC'))

	z_event_min=None
	z_event_min = 4.0
	etas_roc = ETAS_roc(etas_input, fc_date=fcdt, z_event_min=z_event_min)

	etas_roc2 = ETAS_roc('data/nepal_etas_20150507_b.xyz', fc_date=fcdt, z_event_min=z_event_min)

	FH_molch = etas_roc.calc_molchan()
	plt.figure()
	plt.plot(*zip(*FH_molch), marker='', ls='-')
	plt.plot(*zip(*etas_roc.calc_roc()), marker='', ls='-')
	#
	plt.plot(*zip(*etas_roc2.calc_roc()), marker='', ls='-')
	plt.plot(*zip(*etas_roc2.calc_molchan()), marker='', ls='-')
	#
	plt.plot(range(2), range(2), ls='-', lw=2., color='r')
	#
	etas_input_glob = 'global/global_xyz_20151129.xyz'   # this is the old etas, with the 1/area instead of 1/L error... and it 
	#                                        # happens to produce a slightly better ROC... so we should look into
	#                                        # steeper spatial scaling exponents.
	fcdt = dtm.datetime(2015,11,29, tzinfo=pytz.timezone('UTC'))
	#
	z_event_min=None
	z_event_min = 4.0
	etas_roc_glob = ETAS_roc(etas_input_glob, fc_date=fcdt, z_event_min=z_event_min)
	#
	#etas_roc2 = ETAS_roc('data/nepal_etas_20150507_b.xyz', fc_date=fcdt, z_event_min=z_event_min)
	#
	plt.figure()
	plt.plot(*zip(*etas_roc_glob.calc_molchan()), marker='', ls='-')
	plt.plot(*zip(*etas_roc_glob.calc_roc()), marker='', ls='-')
	#
	#plt.plot(*zip(*etas_roc2.calc_roc()), marker='', ls='-')
	#plt.plot(*zip(*etas_roc2.calc_molchan()), marker='', ls='-')
	#
	plt.plot(range(2), range(2), ls='-', lw=2., color='r')



