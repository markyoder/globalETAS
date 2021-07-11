'''
# mapping functions for nETAS / globalETAS
'''
#
import datetime as dtm
import matplotlib.dates as mpd
import pytz
tzutc = pytz.timezone('UTC')
#
import math
import random
import numpy
import scipy
import scipy.optimize as spo
import itertools
import sys
import os
import multiprocessing as mpp
#
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy
#
# let's see if we can compile some of thes functions.
import numba
#
from yodiipy import contours2kml as c2kml
#
import geopy
from geopy.distance import great_circle
from geographiclib.geodesic import Geodesic as ggp
#
class nETAS_mapper():
    def __init__(self, nETAS_input=None, **kwargs):
        # TODO: everything. we need a viable inheritance model... I think this can inherit from Global_ETAS_model(),
        #  which we might rename to nETAS_model or use nETAS_compute or something). MPP handler(s) can also inherit,
        #  but we'll want to come up with a streamlined way for the compute, plotting, and mapping classes to interac.
        #
        # this can stand alone(ish) if we pass it an ETAS object as an input. Alternatively,
        #  we'll create a composite class like globalETAS(nETAS_compute, nETAS_mapper)
        #
        self.__dict__.update({ky:vl for ky,vl in kwargs.items() if not ky in ('self', '__class__')})
        #self.__dict__.update({ky:vl for ky,vl in locals().items() if not ky in ('self', '__class__')})
        #if not nETAS_input is None:
        #    self.__dict__.update({ky:vl} for ky,vl in nETAS_input.__dict__.items() if not ky in ('self', '__class__'))
        
    #
    def draw_map(self, fignum=0, fig_size=(6.,6.), map_resolution='i', map_projection=None, lon_line_spacing=None, lat_line_spacing=None, lats_map=None, lons_map=None, ax=None, do_states=True, do_rivers=True, lake_color='blue', lat_label_indices=[1,1,0,0], lon_label_indices=[0,0,1,1], **kwargs):
        '''
        # TODO: we end up matching up a bunch of procedural calls, which is a big pain. we should write an ETAS_Map() class
        # which includes the contour,etc. figures... but we can keep the variables, like lon_label_indices, etc.
        # in one place...
        #
        # @map_projection: a cartopy.crs.{Projection} instance. TODO: create a string->object
        #  dict to pass as short strings; for now, I think it has to be a class instance.
        # @ax: Axis to plot onto. NOTE: the projection for the axis MUST BE SET in advance. So,
        #  if using this function to plot into an existing axis, set the prjection when the axis
        #  is created.
        '''
        #
        if lons_map is None: lons_map = self.lons
        if lats_map is None: lats_map = self.lats
        cntr = [numpy.mean(lons_map), numpy.mean(lats_map)]
        #self.lon0 = cntr[0]
        self.lon0 = 0.0
        #
        lon_line_spacing = (lon_line_spacing or 1.)
        lat_line_spacing = (lat_line_spacing or 1.)
        #
        if ax==None:
            #
            if map_projection is None or isinstance(map_projection, str):
                map_projection=cartopy.crs.PlateCarree(central_longitude = self.lon0)
            #
            plt.figure(fignum, fig_size)
            plt.clf()
            ax=plt.axes(projection=map_projection)
        #
        self.map_projection = map_projection
        #
        ax.set_extent(numpy.ravel([lons_map, lats_map]) )
        ax.coastlines(color='gray', zorder=1)
        
        #cm.drawcountries(color='black', zorder=1)
        #if do_states: cm.drawstates(color='black', zorder=1)
        #if do_rivers: cm.drawrivers(color='blue', zorder=1)
        #cm.fillcontinents(color='beige', lake_color=lake_color, zorder=0)
        #
        #
        ax.gridlines(label=True, color='black', alpha=.7, linestyle='--', xlocs=numpy.arange(int(lons_map[0]/lon_line_spacing)*lon_line_spacing),\
            ylocs=numpy.arange(int(lats_map[0]/lat_line_spacing)*lat_line_spacing)\
            )

        #
        return ax
    #
    def make_etas_contour_map(self, n_contours=None, fignum=0, fig_size=(6.,6.), contour_fig_file=None, contour_kml_file=None, kml_contours_bottom=0., kml_contours_top=1.0, alpha=.5, alpha_kml=.5, refresh_etas=False, map_resolution='i', map_projection=None, map_cmap=None, lat_interval=None, lon_interval=None, lats_map=None, lons_map=None, ax=None, do_colorbar=True, do_states=True, do_rivers=True, lake_color='blue', Z=None , **kwargs):
        #
        ax = self.draw_map(**{k:v for k,v in locals().items() if not k in ('self', '__class')})
        if map_cmap is None: map_cmap = self.cmap_contours
        #
        n_contours = (n_contours or self.n_contours)
        #
        fg=plt.gcf()
        #
        # X,Y = cm(numpy.array(self.lonses), numpy.array(self.latses))
        X = self.lonses - self.lon0
        Y = self.latses
        if Z is None: Z = numpy.log10(self.lattice_sites)
        #
        etas_contours = ax.contourf(X,Y, Z, n_contours, zorder=8, alpha=alpha,\
            transform=ax.projection, cmap=map_cmap)
        # ax.colorbar() ??
        if do_colorbar:
            try:
                plt.colorbar(etas_contours, cax=None, ax=ax, cmap=map_cmap)
            except:
                print('DEBUG: error creating colorbar() in globalETAS.make_etas_contourmap()')
        #
        self.etas_contours = etas_contours
        #
        return ax
        #
    #
    def make_etas_boxy_map(self, n_contours=None, fignum=0, fig_size=(6.,6.), contour_fig_file=None, contour_kml_file=None, kml_contours_bottom=0., kml_contours_top=1.0, alpha=.6, alpha_kml=.5, refresh_etas=False, map_resolution='i', map_projection=None, map_cmap='jet', ax=None):
        #
        ax = self.draw_map(**{k:v for k,v in locals().items() if not k in ('self', '__class')})
        #
        c_map = plt.get_cmap(map_cmap)
        zs = numpy.log10(self.ETAS_array['z'])
        cNorm = mpl.colors.Normalize(vmin=min(zs), vmax=numpy.nanmax(zs))
        scalarMap = mpl.cm.ScalarMappable(norm=cNorm, cmap=c_map)
        #
        for rw in self.ETAS_array:
            #print('** rw: {}, {}, {}'.format(rw['x'], rw['y'], rw['z']) )
            #ax.plot([rw['x'], rw['y']], marker='.', transform=self.map_projection)
            ax.add_patch(mpl.patches.Rectangle(xy=[rw['x']-self.d_lon-self.lon0, rw['y']-self.d_lat], width=self.d_lon, height=self.d_lat,
                                    facecolor=scalarMap.to_rgba(numpy.log10(rw['z'])),
                                    alpha=0.8,
                                    transform=ax.projection)
                                    )
#            ax.fill_between(x=[rw['x']-self.d_lon/2., rw['x']+self.d_lon/2.], y1=[rw['y']-.5*self.d_lat, rw['y']-.5*self.d_lat], y2=[rw['y']+.5*self.d_lat, rw['y']+.5*self.d_lat], color=scalarMap.to_rgba(numpy.log10(rw['z'])), transform=ax.projection, alpha=alpha, zorder=11 )
        #plt.colorbar()        # not sure how to make this work for non-contour plot...
        #
        return ax
        #
    def plot_mainshock_and_aftershocks(self, m0=6.0, n_contours=25, mainshock=None, fignum=0, ax=None):
        #
        #map_etas = self.make_etas_contour_map(n_contours=n_contours, fignum=fignum, ax=ax)
        if mainshock is None:
            mainshock = self.catalog[0]
            for rw in self.catalog:
                if rw['mag']>mainshock['mag']: mainshock=rw
        ms=mainshock
        #
        ax = (ax or plt.gca())
        #
        for eq in self.catalog:
            if eq['mag']<m0 or eq['event_date']<ms['event_date']: continue
            x = eq['lon'] - self.lon0
            y = eq['lat']
            if eq==ms:
                #
                ax.plot([x], [y], 'k*', zorder=7, ms=20, alpha=.8, transform=ax.projection)
                ax.plot([x], [y], 'r*', zorder=8, ms=18, label='mainshock', alpha=.8, transform=ax.projection)
            if eq['event_date']>eq['event_date']:
                ax.plot([x], [y], 'o', zorder=7, ms=20, alpha=.8, transform=ax.projection)
        #
        #return plt.gca()
        return ax
    #
