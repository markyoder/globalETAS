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
	ETAS.make_etas_contour_map(n_contours=n_contours, fignum=fignum, fig_size=fig_size, contour_fig_file=png_file, contour_kml_file=kml_file, kml_contours_bottom=0., kml_contours_top=1.0, alpha=.6, alpha_kml=.5, refresh_etas=False, map_resolution='i', map_projection='cyl', map_cmap=color_map)
	#
	plt.figure(fignum)
	plt.title('ETAS: \n%s\n' % str(kwargs.get('t_now', dtm.datetime.now(pytz.timezone('UTC')))))
	plt.savefig(png_file)
	#
	ETAS.export_xyz(xyz_file)
	ETAS.export_kml(kml_file)
	#
