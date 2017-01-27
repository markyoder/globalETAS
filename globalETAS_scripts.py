'''
# globalETAS_scripts.py
# author: Mark R. Yoder, Ph.D.
#
# this code is free to use for non-commercial applications. it is research code, so it is probably buggy.
#
# some scripts for (global)ETAS stuff. mainly, ETAS product scripts like, 1) make an ETAS, 2) output contours, .xyz, .kml.
# later, integrate auto_etas() (find big earthquakes, run reports) and other automated bits.
'''

import globalETAS as gep
import pylab as plt
import multiprocessing as mpp
import os
import datetime as dtm
import pytz
#
n_contours = 25
fig_size   = (6.,6.)
color_map = 'spectral'

#
def etas_outputs(n_processes = None, output_path='etas_outputs/', kml_file='etas_kml.kml', png_file=None, xyz_file=None, fignum=0, color_map='jet', *args, **kwargs):
	#
	# standard set of output files for ETAS, including kml, .png, and .xyz files.
	#
	if not os.path.isdir(output_path): os.makedirs(output_path)
	#
	n_processes = (n_processes or mpp.cpu_count())
	# get a base string for output files.
	base_path = ((kml_file or png_file) or xyz_file)
	pth, fl = os.path.split(base_path)
	fl_root, ext = os.path.splitext(fl)
	#
	#
	kml_file = (kml_file or fl_root + '_kml.kml')
	png_file = (png_file or fl_root + '_png.png')
	xyz_file = (xyz_file or fl_root + '_xyz.xyz')
	#
	kml_file = os.path.join(output_path, kml_file)
	png_file = os.path.join(output_path, png_file)
	xyz_file = os.path.join(output_path, xyz_file)
	#
	#kml_file = (kml_file or os.path.join(pth, fl_root + '_kml.kml'))
	#png_file = (png_file or os.path.join(pth, fl_root + '_png.png'))
	#xyz_fyle = (xyz_file or os.path.join(pth, fl_root + '_xyz.xyz'))
	print('files: ', kml_file, png_file, xyz_file)
	#
	ETAS = gep.ETAS_mpp_handler_xyz(n_processes = n_processes, *args, **kwargs)
	#
	ETAS.make_etas_contour_map(n_contours=n_contours, fignum=fignum, fig_size=fig_size, contour_fig_file=png_file, contour_kml_file=kml_file, kml_contours_bottom=0., kml_contours_top=1.0, alpha=.6, alpha_kml=.6, refresh_etas=False, map_resolution='i', map_projection='cyl', map_cmap=color_map)
	#
	plt.figure(fignum)
	plt.title('ETAS: \n%s\n' % str(kwargs.get('t_now', dtm.datetime.now(pytz.timezone('UTC')))))
	plt.savefig(png_file)
	#
	ETAS.export_xyz(xyz_file)
	ETAS.export_kml(kml_file)
	#
	return ETAS
#
# we need a full diagnostic script here for the radial distributions. not sure if the local intensitis, particularly spatial distributions,
# are being calculated correctly... or what the consequence of that is. namely, we need to see spatial distributions...
#
# ... and it looks like they probably are; see one of the *.ipynb notebook worksheets.

"""
~~~~Defaults args~~~~
catalog=None, lats=[32., 36.], lons=[-117., -114.], mc=2.5, mc_etas=None, d_lon=.1, d_lat=.1, bin_lon0=0., bin_lat0=0., etas_range_factor=10.0, etas_range_padding=.25, etas_fit_factor=1.0, t_0=dtm.datetime(1990,1,1, tzinfo=tz_utc), t_now=dtm.datetime.now(tzutc), transform_type='equal_area', transform_ratio_max=2.5, cat_len=5.*365., calc_etas=True, n_contours=15, etas_cat_range=None, etas_xyz_range=None, p_cat=1.1, q_cat=1.5, p_etas=None
"""


etas_outputs(kml_file='etas_kml.kml', png_file='etas.png', xyz_file='etas_xyz.xyz', fignum=0, lats=[31.5, 42.2], lons=[-125., -114.], mc=2.5)

