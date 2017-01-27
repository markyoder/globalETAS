from globalETAS import *
from etas_analyzer import *
#
def cali_rocs(to_dt=None, from_dt=None, cat_len=120):
	# rocs for california +nevata (rectangular catalog).
	#
	to_dt   = (to_dt or dtm.datetime.now(pytz.timezone('UTC')))
	from_dt = (from_dt or to_dt = dtm.timedelta(days=120)
	#
	lats = [31., 43.]
	lons = [-125., -114.]
	#
	# note: calc etas to our roc_from_dt; roc will be from from_dt -> to_dt.
	etas_cali = ETAS_mpp_handler_xyz(n_processes = 6, lats=lats, lons=lons, to_dt=from_dt
	
