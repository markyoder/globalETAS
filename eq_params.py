import math
import datetime as dtm
import pytz
kmldir='kml'
catdir='kml'

lon2km=111.1
deg2rad = math.pi*2.0/360.
tz_utc = pytz.timezone('UTC')
#
# some parameter sets for important earthquakes:
#
tohoku_ETAS_prams = {'todt':dtm.datetime.now(pytz.timezone('UTC')), 'gridsize':.1, 'contres':3, 'mc':4.5, 'kmldir':kmldir, 'catdir':kmldir, 'fnameroot':'tohoku', 'catlen':5.0*365.0, 'doplot':False, 'lons':[135., 146.], 'lats':[30., 41.5], 'bigquakes':None, 'bigmag':7.50, 'eqtheta':None, 'eqeps':None, 'fitfactor':5.0, 'cmfnum':0, 'fignum':1, 'contour_intervals':None}
#
chengdu_ETAS_prams = {'todt':dtm.datetime.now(pytz.timezone('UTC')), 'gridsize':.1, 'contres':5, 'mc':4.5, 'kmldir':kmldir, 'catdir':kmldir, 'fnameroot':'chengdu', 'catlen':5.0*365.0, 'doplot':False, 'lons':[100.367, 106.367], 'lats':[31.06-3., 31.06+3.], 'bigquakes':None, 'bigmag':6.7, 'eqtheta':None, 'eqeps':None, 'fitfactor':5.0, 'cmfnum':0, 'fignum':1, 'contour_intervals':None}

#'lat_center':31.021, 'lon_center':103.367
#
nepal_epi_lon = 84.698
nepal_epi_lat = 28.175
nepal_dlon = 5.
nepal_dlat = 5.
nepal_ETAS_prams = {'todt':None, 'gridsize':.1, 'contres':5, 'mc':3.5, 'kmldir':kmldir, 'catdir':kmldir, 'fnameroot':'nepal', 'catlen':5.0*365.0, 'doplot':False, 'lons':[nepal_epi_lon-nepal_dlon, nepal_epi_lon+nepal_dlon], 'lats':[nepal_epi_lat-nepal_dlat, nepal_epi_lat+nepal_dlat], 'bigquakes':None, 'bigmag':7.00, 'eqtheta':None, 'eqeps':None, 'fitfactor':5.0, 'cmfnum':0, 'fignum':1, 'contour_intervals':None}#
