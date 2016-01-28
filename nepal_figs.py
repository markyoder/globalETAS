import pylab as plt
import numpy
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap

import globalETAS
import roc_generic
import contours2kml


class Map_drawer(object):
	# container class to draw maps and stuff like that.
	def __init__(self, xyz='data/global_xyz_20151129.xyz'):
		#self.__dict__.update(locals())
		#self.fignum=fignum
		if isinstance(xyz, str):
			with open(xyz, 'r') as f:
				xyz = [[float(x) for x in rw.split()] for rw in f if rw[0]!='#']
			#
		#
		self.XYZ = xyz
		#
		self.lonses=sorted(list(set([rw[0] for rw in self.XYZ])))
		self.latses=sorted(list(set([rw[1] for rw in self.XYZ])))
		self.lons = (min(self.lonses), max(self.lonses))
		self.lats = (min(self.latses), max(self.latses))
		#
		
			
	def draw_map(self, fignum=0, fig_size=(6.,6.), map_resolution='i', map_projection='cyl', d_lon_range=1., d_lat_range=1.):
		'''
		# plot contours over a map.
		'''
		#
		# first, get contours:
		#etas_contours = self.calc_etas_contours(n_contours=n_contours, fignum=fignum, contour_fig_file=contour_fig_file, contour_kml_file=contour_kml_file, kml_contours_bottom=kml_contours_bottom, kml_contours_top=kml_contours_top, alpha_kml=alpha_kml, refresh_etas=refresh_etas)
		#
		# now, clear away the figure and set up the basemap...
		#		
		plt.figure(fignum, fig_size)
		plt.clf()
		#
		lons, lats = self.lons, self.lats
		cntr = [numpy.mean(lons), numpy.mean(lats)]
		cm = Basemap(llcrnrlon=self.lons[0], llcrnrlat=self.lats[0], urcrnrlon=self.lons[1], urcrnrlat=self.lats[1], resolution=map_resolution, projection=map_projection, lon_0=cntr[0], lat_0=cntr[1])
		#
		#cm.drawlsmask(land_color='0.8', ocean_color='b', resolution=map_resolution)
		cm.drawcoastlines(color='gray', zorder=1)
		cm.drawcountries(color='black', zorder=1)
		cm.drawstates(color='black', zorder=1)
		#cm.drawrivers(color='blue', zorder=1)
		cm.fillcontinents(color='beige', lake_color='blue', zorder=0)
		# drawlsmask(land_color='0.8', ocean_color='w', lsmask=None, lsmask_lons=None, lsmask_lats=None, lakes=True, resolution='l', grid=5, **kwargs)
		#cm.drawlsmask(land_color='0.8', ocean_color='c', lsmask=None, lsmask_lons=None, lsmask_lats=None, lakes=True, resolution=self.mapres, grid=5)
		#
		#
		cm.drawmeridians(numpy.arange(int(lons[0]), int(lons[1]), d_lon_range), color='k', labels=[0,0,1,1])
		cm.drawparallels(numpy.arange(int(lats[0]), int(lats[1]), d_lat_range), color='k', labels=[1, 1, 0, 0])
		#
		return cm
	def make_etas_contour_map(self, n_contours=None, fignum=0, fig_size=(6.,6.), contour_fig_file=None, contour_kml_file=None, kml_contours_bottom=0., kml_contours_top=1.0, alpha=.6, alpha_kml=.5, refresh_etas=False, map_resolution='i', map_projection='cyl', map_cmap='spectral'):
		n_contours = (n_contours or self.n_contours)
		#
		cm = self.draw_map(fignum=fignum, fig_size=fig_size, map_resolution=map_resolution, map_projection=map_projection)
		#X,Y = cm(numpy.array(self.lonses), numpy.array(self.latses))
		X,Y = sorted(list(set([rw[0] for rw in self.XYZ]))), sorted(list(set([rw[1] for rw in self.XYZ])))
		Zs = numpy.array(numpy.log10([rw[2] for rw in self.XYZ]))
		#
		Zs.shape = (len(Y), len(X))
		#
		etas_contours = plt.contourf(X,Y, Zs, n_contours, zorder=8, alpha=alpha)
		plt.colorbar()
		#
		self.cm=cm
		#
		return cm
		#

def draw_global_etas_contours(xyz='data/global_xyz_20151129.xyz', fignum=0, n_conts=12):
	mm = Map_drawer(xyz=xyz)
	mm.draw_map(d_lat_range=10., d_lon_range=10., fignum=0)
	#return mm
	#
	lns, lts = mm.lonses, mm.latses
	Zs = numpy.log10([rw[2] for rw in mm.XYZ])
	Zs.shape=(len(lts), len(lns))
	#
	plt.figure(fignum)
	plt.contourf(lns, lts, Zs, n_conts, alpha=.5, zorder=11, cmap=plt.cm.coolwarm)
	
	
	
