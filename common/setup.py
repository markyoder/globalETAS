#!/usr/bin/env python

from distutils.core import setup

setup(name='common',
      version='1.0',
      description='Common research utilities',
      author='Mark R. Yoder',
      author_email='mark.yoder@gmail.com',
      url='',
      py_modules = ['ANSStools', 'clusterutils', 'contours', 'eqcatalog', 'eqcataloglite', 'kmlparser', 'linefit', 'yodapy', 'rbIntervals', 'ygmapbits', 'pca_tools'],
      #packages=['../common'],
     )
