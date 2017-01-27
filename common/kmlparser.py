# eventually, this will be dynamically called from the web. for now, let's use this to parse data files/queries and
# output to files which we can upload to google maps.

import contours

#def add_kmlstyle(filepath, stylename):
def add_kmlstyle(stylename="redsies"):
	strStyle=""
	
	if stylename=='transBluePoly':
		strStyle= \
		'<Style id="transBluePoly">\n'\
      ' <LineStyle>\n'\
      '    <width>1</width>\n'\
      '    <color>5dff0000</color>\n'\
      '  </LineStyle>\n'\
      '  <PolyStyle>\n'\
      '    <color>7dff0000</color>\n'\
      '  </PolyStyle>\n'\
      '</Style>\n'

	elif stylename=='redsies':
		print("doing redsies.")
		strStyle= \
		'<Style id="redsies">\n'\
      ' <LineStyle>\n'\
      '    <width>1</width>\n'\
      '    <color>ff0000ff</color>\n'\
      '  </LineStyle>\n'\
      '  <PolyStyle>\n'\
      '    <color>dd0000ff</color>\n'\
      '  </PolyStyle>\n'\
      '</Style>\n'

	elif stylename=='redbg':
		print("doing redbg")
		strStyle= \
		'<Style id="redbg">\n'\
      ' <LineStyle>\n'\
      '    <width>1</width>\n'\
      '    <color>3d0000ff</color>\n'\
      '  </LineStyle>\n'\
      '  <PolyStyle>\n'\
      '    <color>4d0000ff</color>\n'\
      '  </PolyStyle>\n'\
      '</Style>\n'

	elif stylename=='yellowsies':
		strStyle= \
		'<Style id="yellowsies">\n'\
      ' <LineStyle>\n'\
      '    <width>1</width>\n'\
      '    <color>5d007d7d</color>\n'\
      '  </LineStyle>\n'\
      '  <PolyStyle>\n'\
      '    <color>7d007d7d</color>\n'\
      '  </PolyStyle>\n'\
      '</Style>\n'

	elif stylename=='orangsies':
		strStyle= \
		'<Style id="orangsies">\n'\
      ' <LineStyle>\n'\
      '    <width>1</width>\n'\
      '    <color>ff25708d</color>\n' \
      '  </LineStyle>\n'\
      '  <PolyStyle>\n'\
      '    <color>7d2525ee</color>\n'\
      '  </PolyStyle>\n'\
      '</Style>\n'
	#else: strStyle=add_kmlstyle()
	if strStyle=="": strStyle=add_kmlstyle()
	
	return strStyle

def add_placemarker(latitude, longitude, altitude = 0.0, description = " ", name = " ", range = 6000, tilt = 45, heading = 0, styleString='redsies'):
	#"adds the point to a kml file"
	#file = open(filepath,"a")
	# style: normalPlaceMarker
	kmlstr = \
	"<Placemark>\n"\
	" <description>" + description + "</description>\n"\
	" <name>" + name + "</name>\n"\
	" <styleUrl>#normalPlaceMarker</styleUrl>" + \
	" <LookAt>\n"\
	"   <longitude>" + str(longitude) + "</longitude>\n"\
	"   <latitude>" + str(latitude) + "</latitude>\n"\
	"   <range>" + str(range) + "</range>\n"\
	"   <tilt>" + str(tilt) + "</tilt>\n"\
	"   <heading>" + str(heading) + "</heading>\n"\
	" </LookAt>\n"\
	" <visibility>0</visibility>\n"\
	" <Point>\n"\
	"   <extrude>1</extrude>\n"\
	"   <altitudeMode>relativeToGround</altitudeMode>\n"\
	"   <coordinates>" + str(longitude) + "," + str(latitude) +", " + str(altitude) + "</coordinates>\n"\
	" </Point>\n"\
	" </Placemark>\n"
	#file.write(kmlstr)
	#file.close()
	#
	return kmlstr

def add_polygon2ring(outcoords, incoords, altitude = 0.0, description = " ", name = " ", range = 6000, tilt = 45, heading = 0):
	#file = open(filepath,"a")
	#file.write(
	kmlStr = \
	"<Placemark>\n"\
	" <description>" + description + "</description>\n"\
	" <name>" + name + "</name>\n"\
	" <styleUrl>#transBluePoly</styleUrl>\n"\
	" <visibility>0</visibility>\n"\
	"<Polygon>\n"\
   "   <extrude>1</extrude>\n"\
   "   <altitudeMode>relativeToGround</altitudeMode>\n"\
   "   <outerBoundaryIs>\n"\
   "     <LinearRing>\n"\
   "       <coordinates>" + outcoords + "</coordinates>\n"\
   "     </LinearRing>\n"\
   "   </outerBoundaryIs>\n"\
   "   <innerBoundaryIs>\n"\
   "     <LinearRing>\n"\
   "       <coordinates>" + inoutcoords + "</coordinates>\n"\
   "     </LinearRing>\n"\
   "   </innerBoundaryIs>\n"\
   " </Polygon>\n"\
	 " </Placemark>\n"
	
	#file.close()
	return kmlStr
	#	if we want to plat a placemark somewhere in the polygon?
	#	 " <styleUrl>#normalPlaceMarker</styleUrl>" + 
	#	 " <LookAt>\n"\
	#	 "   <longitude>" + str(longitude) + "</longitude>\n"\
	#	 "   <latitude>" + str(latitude) + "</latitude>\n"\
	#	 "   <range>" + str(range) + "</range>\n"\
	#	 "   <tilt>" + str(tilt) + "</tilt>\n"\
	#	 "   <heading>" + str(heading) + "</heading>\n"\
	#	 " </LookAt>\n"\

def getCoordsString(coords):
	# should be like [[lat, lon], [lat, lon]...[]]
	strCoords='\n'
	for rw in coords:
		strCoords+='%f, %f\n' % (rw[1], rw[0])
	return strCoords

def add_simplePoly(coords, altitude = 0.0, description = " ", name = " ", range = 6000, tilt = 45, heading = 0, styleStr='redsies'):
	#file = open(filepath,"a")
	#file.write(
	kmlStr = \
	"<Placemark>\n"\
	" <description>" + description + "</description>\n"\
	" <name>" + name + "</name>\n"\
	" <styleUrl>#" + styleStr + "</styleUrl>\n"\
	" <visibility>0</visibility>\n"\
	"<Polygon>\n"\
   "   <extrude>1</extrude>\n"\
   "   <altitudeMode>relativeToGround</altitudeMode>\n"\
   "   <outerBoundaryIs>\n"\
   "     <LinearRing>\n"\
   "       <coordinates>" + getCoordsString(coords) + "</coordinates>\n"\
   "     </LinearRing>\n"\
   "   </outerBoundaryIs>\n"\
   " </Polygon>\n"\
	 " </Placemark>\n"
	
	#file.close()
	return kmlStr

def add_square(lat0, lon0, sidelen, altitude = 0.0, description = " ", name = " ", style="redsies"):
#def add_square(lat1, lon1, sidelen):
	# provide square centers.
	#"adds the point to a kml file"
	lat1 = eval(lat0)-sidelen/2
	lon1 = eval(lon0)-sidelen/2
	#print str(5)
	#
	strCoords = str(lon1) + ", " + str(lat1) + "\n" +\
	str(lon1+sidelen) + ", " + str(lat1) + "\n" +\
	str(lon1+sidelen) + ", " + str(lat1+sidelen) + "\n" +\
	str(lon1) + ", " + str(lat1+sidelen) + "\n" +\
	str(lon1) + ", " + str(lat1) + "\n"
	#print strCoords
	#file = open(filepath,"a")
	#file.write(
	strKml = \
	"<Placemark>\n"\
	" <description>" + description + "</description>\n"\
	" <name>" + name + "</name>\n"\
	" <styleUrl>#" + style + "</styleUrl>\n"\
	" <visibility>0</visibility>\n"\
	"<Polygon>\n"\
   "   <extrude>1</extrude>\n"\
   "   <altitudeMode>relativeToGround</altitudeMode>\n"\
   "   <outerBoundaryIs>\n"\
   "     <LinearRing>\n"\
   "       <coordinates>" + strCoords + "</coordinates>\n"\
   "     </LinearRing>\n"\
   "   </outerBoundaryIs>\n"\
   " </Polygon>\n"\
	 " </Placemark>\n"
		
	#file.close()
	#
	return strKml

def PI2KMLpolys(data, z0=0, ix=0, iy=1, iz=2, dx=.1, dy=.1, strStyle='orangsies', altitude = 0.0, description = " ", name = " ", range = 6000, tilt = 45, heading = 0):
	# take in a wad of PI (pattern informatics) data points, spit out KML polygons (huge string-o-text)
	#
	# get a tuple of poly-tuples.
	mypolys = contours.boxyContour(data, z0, ix, iy, iz, dx, dy)
	#
	# each entry is a polygon (with several point-elements)
	kmlStr=''
	for poly in mypolys:
		kmlStr = kmlStr +\
		"<Placemark>\n"\
		" <description>" + description + "</description>\n"\
		" <name>" + name + "</name>\n"\
		" <styleUrl>#" + strStyle + "</styleUrl>\n"\
		" <visibility>0</visibility>\n"\
		"<Polygon>\n"\
		"   <extrude>1</extrude>\n"\
		"   <altitudeMode>relativeToGround</altitudeMode>\n"\
		"   <outerBoundaryIs>\n"\
		"     <LinearRing>\n"\
		"       <coordinates>"
		for corner in poly:
			#
			kmlStr = kmlStr + str(corner[0]) + ', ' + str(corner[1]) + '\n'
			#
		#
		kmlStr = kmlStr + \
		"</coordinates>\n"\
		"     </LinearRing>\n"\
		"   </outerBoundaryIs>\n"\
		" </Polygon>\n"\
		 " </Placemark>\n"
	#
	#
	return kmlStr

			
	


