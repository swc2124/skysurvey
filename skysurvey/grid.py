"""
TODO
-------------------------------------------------------------------------------
 Author     : Sol W. Courtney
 Location   : Columbia U Dept. of Astronomy & Astrophysics NYC, Ny
 Email      : swc2124@Columbia.edu
 Date       : June 2017
 Title      : skysurvey/skysurvey/grid.py
-------------------------------------------------------------------------------
"""

from __future__ import division, absolute_import, print_function

import os

import ConfigParser
from numpy import load
from numpy import save
from numpy import ndarray

Config = ConfigParser.ConfigParser()
Config.read('../setup.cfg')

class Grid(ndarray):
	"""The output grid from spinbin binning halo data"""
	def __init__(self, 
			name=halo, 
			d_mpc=Config.getint('Distance', 'd_mpc'), 
			f_type=Config.get('Filter', 'filter_type'),
			f_prfx=Config.get('Filter', 'filter_prfx'),
			grid_dir=Config.get('PATH', 'grid_dir')):

		super(Grid, self).__init__()
		
		self.name = name
		self.d_mpc = d_mpc
		self.f_type = f_type
		self.f_prfx = f_prfx
		self.grid_dir = grid_dir

		self.file_handel = self.name + '_' + str(self.d_mpc) + 'Mpc_' + self.f_type + '_grid.npy'
		self.PATH = os.path.join(self.grid_dir, self.file_handel)

		if os.path.isfile(self.PATH):
			self.grid = load(self.PATH)

	