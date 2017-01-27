import datetime as dtm
import matplotlib.dates as mpd
import pytz
import calendar
import operator
import sys

import urllib
try:
	# should work with python 3.x
	import urllib.request, urllib.parse, urllib.error
except:
	print("failed while loading: urllib.request, urllib.parse, urllib.error.\n probably Python 2.x?")
	#urllib.request.urlopen = urllib.urlopen

try:
	import ullib2
except:
	print("failed while loading urllib and/or urllib. maybe python 3.x?")
	
#import urllib.request, urllib.error, urllib.parse
import requests
import numpy
#
# note on datetimes:
# timezone awareness is confusing but not as bad as it looked a minute ago.
# datetime.utcnow() gives a UTC time, but without tz awareness.
# probably the better approach is to use datetime.now(pytz.timezone('UTC'))
tzutc=pytz.timezone('UTC')

def anssDateStr(x=dtm.datetime.now(pytz.timezone('UTC')) ):
	# yoder, 13 july 2015: ANSS seems to have made some changes. these date formats are breaking. probalby a matter of leading 0's in dates; might be fractional seconds.
	#yr=x.year
	#mo=x.month
	#dy=x.day
	#hr=x.hour
	#mn=x.minute
	#sc=x.second
	#ms=x.microsecond
	#fsecs=float(sc) + float(ms)*(10**(-6.0))
	#
	yr = str(x.year)
	mo = ('00' + str(x.month))[-2:]
	dy = ('00' + str(x.day))[-2:]
	hr = ('00' + str(x.hour))[-2:]
	mn = ('00' + str(x.minute))[-2:]
	sc = ('00' + str(x.second))[-2:]
	#
	# ANSS seems to be complaining about fractional seconds, so skip this and return integer seconds.
	'''
	ms=x.microsecond
	fsecs=float(sc) + float(ms)*(10**(-6.0))
	#
	# trim extra zeros:
	fsecs_str = str(fsecs)
	while ('.' in fsecs_str and len(fsecs_str)>3 and fsecs_str[-1]=='0'):
		fsecs_str = fsecs_str[:-1]
	'''
	#
	#return '%s/%s/%s,%s:%s:%f' % (yr, mo, dy, hr, mn, fsecs)
	return '%s/%s/%s,%s:%s:%s' % (yr, mo, dy, hr, mn, sc)
	
def getANSStoFilehandler(lon=[-125, -115], lat=[32, 45], minMag=4.92, dates0=[dtm.datetime(2001,1,1, tzinfo=tzutc), dtm.datetime(2010, 12, 31, tzinfo=tzutc)], Nmax=999999):
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
	# revision: ANSS has time resolution, but you have to replace "-" -> "/" and the " " (space) -> ","
	#dates=[dtm.date(dates0[0].year, dates0[0].month, dates0[0].day), dtm.date(dates0[1].year, dates0[1].month, dates0[1].day)]
	dates=dates0
	datestr1 = anssDateStr(dates[0])
	datestr2 = anssDateStr(dates[1])
	#print datestr1, datestr2
	#
	#anssPrams={'format':'cnss', 'output':'readable', 'mintime':str(dates[0]).replace('-', '/'), 'maxtime':str(dates[1]).replace('-', '/'), 'minmag':str(minMag), 'minlat':lat[0], 'maxlat':lat[1], 'minlon':lon[0], 'maxlon':lon[1], 'etype':'E', 'searchlimit':Nmax}
	# so this is better, but i think it is still limited to 1 second resolution.
	#
	anssPrams={'format':'cnss', 'output':'readable', 'mintime':datestr1, 'maxtime':datestr2, 'minmag':str(minMag), 'minlat':lat[0], 'maxlat':lat[1], 'minlon':lon[0], 'maxlon':lon[1], 'etype':b'E', 'searchlimit':Nmax}
	#anssPrams={'format':b'cnss', 'output':b'readable', 'mintime':bytearray(datestr1, 'utf-8'), 'maxtime':bytearray(datestr2, 'utf-8'), 'minmag':bytearray(str(minMag), 'utf-8'), 'minlat':lat[0], 'maxlat':lat[1], 'minlon':lon[0], 'maxlon':lon[1], 'etype':b'E', 'searchlimit':Nmax}
	#print "debug: ", anssPrams
	#post_data = urllib.parse.urlencode(anssPrams)
	#binary_post_data = post_data.encode('ascii')
	#
	# now, let's support some backwards compatibility, at least for a little while:
	if sys.version_info.major == 2:
		# old python...
		f = urllib.urlopen('http://www.ncedc.org/cgi-bin/catalog-search2.pl', urllib.urlencode(anssPrams))
	else:
		binary_post_data = urllib.parse.urlencode(anssPrams).encode('ascii')
		f = urllib.request.urlopen('http://www.ncedc.org/cgi-bin/catalog-search2.pl', binary_post_data )
	#
	# we might return f, a string of f, or maybe a list of lines from f. we'll work that out shortly...
	return f
#
def catfromANSS(lon=[135., 150.], lat=[30., 41.5], minMag=4.0, dates0=[dtm.datetime(2005,1,1, tzinfo=tzutc), None], Nmax=None, fout=None, rec_array=True):
	# get a basic catalog. then, we'll do a poly-subcat. we need a consistent catalog.
	# eventually, cut up "japancatfromANSS()", etc. to call this base function and move to yodapy.
	#
	# note: there may be a version inconsisency. older versions of this function may have returned catlist raw, in which
	# [..., depth, mag], where we regurn [..., mag, depth] here.
	#
	if Nmax==None: Nmax=999999
	#
	if dates0[1]==None:
		# i think this needs a "date" object, and datetime breaks.
		# so, make a Now() for date.
		#nowdtm=dtm.datetime.now()
		#dates0[1]=dtm.date(nowdtm.year, nowdtm.month, nowdtm.day)
		dates0[1]=dtm.datetime.now(tzutc)
	#	
	catlist=getANSSlist(lon, lat, minMag, dates0, Nmax, None)
	if fout==None: print(" no file.")
	
	if fout!=None:
		f=open(fout, 'w')
		f.write("#anss catalog\n")
		f.write("#query_date(UTC): %s" % str(dtm.datetime.now(pytz.timezone('UTC'))))
		f.write("#lon=%s\tlat=%s\tm0=%f\tdates=%s\n" % (str(lon), str(lat), minMag, str(dates0)))
	
	rlist=[]
	for rw in catlist:
		# simple, right? except that ANSS has a habit of writing useless date-times like "2001/10/08 24:00:07.62" (hour=24), or
		# where minute=60. we could toss these. for now, assume 2001/10/8 24:00:00 -> 2001/10/9/00:00:00. change by proper time-arithmetic.
		#
		# it might be worth checking to see if numpy.datetime64() handles these exceptions.
		#
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
		myDt=dtm.datetime(yr, mnth, dy, tzinfo=tzutc)
		#mytimedelta=dtm.timedelta(hours=hr)
		myDt+=dtm.timedelta(hours=hr)
		myDt+=dtm.timedelta(minutes=mn)
		myDt+=dtm.timedelta(seconds=sc)
		myDt+=dtm.timedelta(microseconds=microsecs)
		#
		# note: we switch the order of depth, mag here. 
		#"list" gives [dt, lat, lon, depth, mag]; "cat" gives [dt, lat, lon, mag, depth?]
		# if we add a float-date value to the end of this, does it screw up any conventions? might screw up BASScasts. so we can add it here -- but
		# it has to be at the end. so let's give it a go (sometimes facilitates easier date-handling), but be prepared to drop it.
		# 
		#rlist +=[[myDt, float(rw[1]), float(rw[2]), float(rw[4]), float(rw[3])]]
		rlist +=[[myDt, float(rw[1]), float(rw[2]), float(rw[4]), float(rw[3]), mpd.date2num(myDt)]]
		if fout!=None:
			myDtStr='%d/%d/%d %d:%d:%d.%d' % (myDt.year, myDt.month, myDt.day, myDt.hour, myDt.minute, myDt.second, myDt.microsecond)	
			#
			#f.write('%s\t%s\t%s\t%s\n' % (rw[0], rw[1], rw[2], rw[4]))
			#f.write('%s\t%s\t%s\t%s\n' % (myDtStr, rw[1], rw[2], rw[4]))
			
			#f.write('%s\n' % '\t'.join([str(x) for x in [myDtStr] + rw[1:]]))
			f.write('%s\n' % '\t'.join([str(x) for x in [myDtStr, float(rw[1]), float(rw[2]), float(rw[4]), float(rw[3]), mpd.date2num(myDt)]]))
			
			#f.write('%s\t%s\t%s\t%s\t%s\n' % (myDtStr, rw[1], rw[2], rw[4], rw[3]))	# indlude depth...
	if fout!=None:
		f.close()
	# 
	#return catlist
	# to do:
	# re-cast rlist as a recarray. here's probaby the best way to do this (there are lots of ways to
	# cast recarrays; this appears to be the most direct:
	# rlist=numpy.rec.array(rlist, dtype=[('event_date', 'M8[us]'), ('lat','f'), ('lon','f'), ('mag','f'), ('depth','f')])	# note: numpy.rec.array() also has "names=" and "formats=" keywords.
	# (but we'll want to test existing programs to be sure this doesn't break).
	#
	# yoder: cast as recarray:
	if rec_array:
		#rlist=numpy.rec.array(rlist, dtype=[('event_date', 'M8[us]'), ('lat','f'), ('lon','f'), ('mag','f'), ('depth','f')])
		# note: specify f8 (or greater), or float dates might be truncated).
		rlist=numpy.rec.array((rlist if len(rlist)>0 else [[]]), dtype=[('event_date', 'M8[us]'), ('lat','f8'), ('lon','f8'), ('mag','f8'), ('depth','f8'), ('event_date_float', 'f8')])
	
	return rlist
#
def cat_from_usgs(duration='week', mc=2.5, rec_array=True):
	# use: http://earthquake.usgs.gov/earthquakes/feed/v1.0/summary/2.5_day.csv
	# or better (maybe), use the geojson format:
	#    http://earthquake.usgs.gov/earthquakes/feed/v1.0/summary/2.5_day.geojson
	#
	# for one week:
	# http://earthquake.usgs.gov/earthquakes/feed/v1.0/summary/2.5_week.csv
	# http://earthquake.usgs.gov/earthquakes/feed/v1.0/summary/2.5_week.geojson
	#
	# the geojson format is pretty awesome, but we just want a catalog, so for now, let's keep it old-school simple and use the .csv
	# incidentally, to read the geojson, use something like:
	# with open(gj_file, 'r') as f:
	#	data = json.load(f)
	# quakes = data['features']		# it will have (presently) 4 primary entries, ['type', 'features', 'bbox', 'metadata']
	# geojson returns the event date as miliseconds since 1-1-1970 (i think). note that matplotlib.dates.date2num() returns
	# (fractional) number of days since 0001-01-01 00:00:00 UTC, plus one, where the "plus one" is a "historical artifact" of some sort.
	#
	# one way to fetch as a file handle from the url...
	# f = urllib.urlopen('http://www.ncedc.org/cgi-bin/catalog-search2.pl', urllib.urlencode(anssPrams))
	# also, a = requests.get(url)	returns an iterable.
	#
	if duration not in ('day', 'week', 'month'): duration = 'week'
	if isinstance(mc, int): mc=float(mc)
	if isinstance(mc, float): mc=str(mc)
	if mc not in ('2.5', '4.5', 'significant', 'all'): mc='2.5'
	#
	#mag_string = '%s_%s' % (mc, duration)
	#
	# for now, stick with the csv:
	cat_out = []
	print("url_str: %s" % ('http://earthquake.usgs.gov/earthquakes/feed/v1.0/summary/%s_%s.csv' % (mc, duration)))
	#url_data = requests.get('http://earthquake.usgs.gov/earthquakes/feed/v1.0/summary/%s_%s.csv' % (mc, duration))
	#with urllib.urlopen('http://earthquake.usgs.gov/earthquakes/feed/v1.0/summary/%s_%s.csv' % (mc, duration)) as furl:
	# 3.x likes:
	#furl = urllib.request.urlopen('http://earthquake.usgs.gov/earthquakes/feed/v1.0/summary/%s_%s.csv' % (mc, duration))
	# 2.x likes:
	# furl = urllib.urlopen('http://earthquake.usgs.gov/earthquakes/feed/v1.0/summary/%s_%s.csv' % (mc, duration))
	if True:
	#for url_rw in url_data:
		if sys.version_info.major == 2:
			furl = urllib.urlopen('http://earthquake.usgs.gov/earthquakes/feed/v1.0/summary/%s_%s.csv' % (mc, duration))
			cols = furl.readline().replace('\n', '').split(',')
		else:
			furl = urllib.request.urlopen('http://earthquake.usgs.gov/earthquakes/feed/v1.0/summary/%s_%s.csv' % (mc, duration))
			cols = (furl.readline()).decode('utf-8').replace('\n', '').split(',')
		#
		# index of cols we care about (... later):
		#my_cols = ['latitude'
		#
		#my_col_names = ['event_date', 'lat', 'lon', 'mag', 'depth']
		#my_col_types = ['M8[us]', 'float', 'float', 'float', 'float']
		#
		for rw_0 in furl:
			if sys.version_info.major == 2:
				rw = rw_0
			else:
				rw=rw_0.decode('utf-8')
			#
			if rw[0] in (' ', '\n', '\r', '\t', '#'): continue
			#
			rw=rw.replace('\n', '')
			rws = rw.split(',')
			dt_str = rws[0]
			dt_str, tm_str = dt_str.split('T')
			yr, mnth, dy = [int(x) for x in dt_str.split('-')]
			hr, mn, secs = [x for x in tm_str.split(':')]
			secs, msecs = [int(x) for x in (secs[:-1] if secs[-1]=='Z' else secs).split('.') ]
			hr, mn = [int(x) for x in [hr, mn]]
			#
			this_dt = dtm.datetime(yr, mnth, dy, hr, mn, secs, msecs, tzinfo=pytz.timezone('UTC'))
			#
			# and put together a row in the standard order/format (which will give way to more descriptive formats soon enough...):
			# anssList+=[[rwEvdt, rwLat, rwLon, rwDepth, rwMag, rwMagType, rwNst, rwGap, rwClo, rwrms, rwsrc, rwCatEventId]]
			#
			#out_rw  = [this_dt] + [float(x) for x in [rws[1:3] + rws[4] + rws[3]]] + rws[5:7]
			# gap, rwClo + rms
			#out_rw += [float(x) for x in rws[rw[7:10]]]
			# for these, the "id" field is just the datetime, and anyway we only really care about a few of these fields. so let's wrap this up in a recarray,
			#
			#
			cat_out += [[this_dt] + [float(x) for x in rws[1:3] + [rws[4], rws[3]]]]
		#
		if rec_array:
			cat_out=numpy.rec.array(cat_out, dtype=[('event_date', 'M8[us]'), ('lat','f'), ('lon','f'), ('mag','f'), ('depth','f')])	
		#
		return cat_out
#
def dictfromANSS(lons=[135., 150.], lats=[30., 41.5], mc=4.0, date_range=[dtm.datetime(2005,1,1, tzinfo=tzutc), None], Nmax=999999, fout='cats/mycat.cat'):
	#
	# get a dictionary type catalog (aka, a list of dicts[{}, {}, ...]
	# note the modified syntax to the newer standards (lon-->lons, minMag --> mc, dates0-->date_range) for
	# consistency with newer modules and codes.
	#
	if date_range[1]==None:
		# in the past, this field requred a DATE type object; DATETIME would break.
		# that problem appears to be fixed now.
		date_range[1]=dtm.datetime.now(tzutc)
	#	
	catlist=getANSSlist(lons, lats, mc, date_range, Nmax, None)
	if fout==None: print(" no file.")
	
	if fout!=None:
		f=open(fout, 'w')
		f.write("#anss catalog\n")
		f.write("#lons=%s\tlats=%s\tm0=%f\tdates=%s\n" % (str(lons), str(lats), mc, str(date_range)))
	
	rlist=[]
	for rw in catlist:
		# simple, right? except that ANSS has a habit of writing useless date-times like "2001/10/08 24:00:07.62" (hour=24), or
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
		myDt=dtm.datetime(yr, mnth, dy, tzinfo=tzutc)
		#mytimedelta=dtm.timedelta(hours=hr)
		myDt+=dtm.timedelta(hours=hr)
		myDt+=dtm.timedelta(minutes=mn)
		myDt+=dtm.timedelta(seconds=sc)
		myDt+=dtm.timedelta(microseconds=microsecs)
		#
		# note: we switch the order of depth, mag here. 
		#"list" gives [dt, lats, lons, depth, mag]; "cat" gives [dt, lats, lons, mag, depth?]
		#
		#rlist +=[[myDt, float(rw[1]), float(rw[2]), float(rw[4]), float(rw[3])]]
		rlist +=[{'event_date':myDt, 'lats':float(rw[1]), 'lons':float(rw[2]), 'mag':float(rw[4]), 'depth':float(rw[3])}]
		#
		if fout!=None:
			# if we wanted to get nutty, we could output this as JSON...
			#
			myDtStr='%d/%d/%d %d:%d:%d.%d' % (myDt.year, myDt.month, myDt.day, myDt.hour, myDt.minute, myDt.second, myDt.microsecond)	
			#
			#f.write('%s\t%s\t%s\t%s\n' % (rw[0], rw[1], rw[2], rw[4]))
			#f.write('%s\t%s\t%s\t%s\n' % (myDtStr, rw[1], rw[2], rw[4]))
			f.write('%s\t%s\t%s\t%s\t%s\n' % (myDtStr, rw[1], rw[2], rw[4], rw[3]))	# indlude depth...
	if fout!=None:
		f.close()
	
	 
	#return catlist
	return rlist
#
def getANSSlist(lon=[-125, -115], lat=[32, 45], minMag=4.92, dates0=[dtm.datetime(2001,1,1, tzinfo=tzutc), dtm.datetime(2010, 12, 31, tzinfo=tzutc)], Nmax=999999, fin=None):
	#
	# this is typically a preliminary function call. it returns a list object-catalog. the date will be in string format.
	# typicall, use catfromANSS() for a more useful list. also see the (new) dictfromANSS() for a dict. type catalog.
	#
	# note: this appears to be a bad idea for global downloads. a full catalog is ~4GB, which kills my computer.
	#
	# note: this may be repeated exactly in ygmapbits.py
	# fetch new ANSS data; return a python list object of the data.
	# fin: data file handler. if this is None, then get one from ANSS.
	#dates=[dtm.date(dates0[0].year, dates0[0].month, dates0[0].day), dtm.date(dates0[1].year, dates0[1].month, dates0[1].day)]
	dates=dates0	# date/datetime issue is fixed.(ish)
	anssList=[]
	if fin==None:
		#print "get data from ANSS...(%s, %s, %s, %s, %s)" % (lon, lat, minMag, dates, Nmax)
		fin = getANSStoFilehandler(lon, lat, minMag, dates, Nmax)
		#fin = getANSStoFilehandler([-180, 180], [-90, 90], 0, [datetime.date(1910,1,1), datetime.date(2010, 01, 16)], 9999999)

		print("data handle fetched...")
		
	for rw_0 in fin:
		if sys.version_info.major ==  2:
			rw = rw_0
		else:
			rw=rw_0.decode('utf-8')
		#
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
#
	
def auto_cat_params(lon_center=None, lat_center=None, d_lat_0=.25, d_lon_0=.5, dt_0=10, mc_0=4.5, to_dt=None, range_factor=5., **kwargs):
	'''
	# (auto_cat functions still need to be tested.)
	# auto_cat: for a given input area of interest, find the largest earthquake, then get a new catalog around that earthquake,
	# based on rupture length scaling, etc. 
	#
	# d_lat/d_lon: lat/lon spread for preliminary catalog.
	# mc_0: mc for preliminary catalog. we're looking for big earthquakes, so we can save some compute cycles and
	# make this fairly large.
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
		d_lon_0 = -180.
	#
	# get a preliminary catalog:
	#print('stuff: ', lat_center, lon_center, d_lon_0, d_lat_0, mc_0, to_dt)
	cat_0 = catfromANSS(lon=[lon_center-d_lon_0, lon_center+d_lon_0], lat=[lat_center - d_lat_0, lat_center+d_lat_0], minMag=mc_0, dates0=[to_dt-dtm.timedelta(days=dt_0), to_dt], fout=None, rec_array=True)
	#
	#biggest_earthquake = filter(lambda x: x['mag']==max(cat_0['mag']), cat_0)[0]
	mainshock = {cat_0.dtype.names[j]:x for j,x in enumerate(list(filter(lambda x: x['mag']==max(cat_0['mag']), cat_0))[0])}
	#
	# now, get new map domain based on rupture length, etc.
	L_r = .5*10.0**(.5*mainshock['mag'] - 1.76)
	delta_lat = range_factor*L_r/lat2km
	delta_lon = range_factor*L_r/(lat2km*math.cos(deg2rad*mainshock['lat']))
	#print("mainshock data: ", mainshock, L_r, delta_lat, delta_lon)
	#
	return {'lon':[mainshock['lon']-delta_lon, mainshock['lon']+delta_lon], 'lat':[mainshock['lat']-delta_lat, mainshock['lat']+delta_lat], 'mainshock_date':mainshock['event_date'], 'mainshock_lat':mainshock['lat'], 'mainshock_lon':mainshock['lon']}
	
def auto_cat(lon_center=None, lat_center=None, d_lat_0=.25, d_lon_0=.5, dt_0=10,  mc=2.5, mc_0=4.5, to_dt=None, catlen_before=5.0*365.0, catlen_after=5.0*365., range_factor=5., rec_array=True, **kwargs):
	#
	cat_params = auto_cat_params(lon_center=lon_center, lat_center=lat_center, d_lat_0=d_lat_0, d_lon_0=d_lon_0, dt_0=dt_0, mc_0=mc_0, to_dt=to_dt, range_factor=range_factor)
	#
	#return atp.catfromANSS(lon=[mainshock['lon']-delta_lon, mainshock['lon']+delta_lon], lat=[mainshock['lat']-delta_lat, mainshock['lat']+delta_lat], minMag=mc, dates0=[to_dt-dtm.timedelta(days=catlen), to_dt], fout=None, rec_array=True)
	return catfromANSS(lat=cat_params['lat'], lon=cat_params['lon'], minMag=mc, dates0=[to_dt-dtm.timedelta(days=catlen), to_dt], fout=None, rec_array=True)


#
def isnumeric(value):
  return str(value).replace(".", "").replace("-", "").isdigit()
 
def numpy_date_to_datetime(numpy_date, tz='UTC'):
	#
	if isinstance(numpy_date, dtm.datetime): return numpy.date
	#
	if isinstance(numpy_date,float): return mpd.num2date(numpy_date)
	#
	return dtm.datetime(*list(numpy_date.tolist().timetuple())[:6] + [numpy_date.tolist().microsecond], tzinfo=pytz.timezone(tz))
#

