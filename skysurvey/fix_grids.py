from __future__ import division, absolute_import, print_function
import os
import sys
import numpy as np
from astropy.table import Table
from astropy.table import Column

import ConfigParser
from .new_config import SYS_CFG_FNAME
import skysurvey

sys_config_fh = os.path.join(os.path.dirname(
    os.path.realpath(skysurvey.__file__)), SYS_CFG_FNAME)
SysConfig = ConfigParser.ConfigParser()
SysConfig.read(sys_config_fh)
config_fh = SysConfig.get('skysurvey_global_settings', 'config_fh')
Config = ConfigParser.ConfigParser()
Config.read(config_fh)


def fix_rslice(grid, rslices=[9, 10]):
    center = grid.shape[0] / 2
    ratio = (1.0 * Config.getint('grid_options', 'size')) / \
        (2.0 * Config.getfloat('halo_data_settings', 'radial_cut'))
    print('fixing radial data slice')
    print('ratio: ', ratio)
    print('slices: ', rslices)
    print('center:', center)
    for r in rslices:
        for i in range(grid.shape[0]):
            for q in range(grid.shape[1]):
                value = np.sqrt(
                    (np.square(i - center) + np.square(q - center)))
                value /= ratio
                if value > 300.0:
                    value = 0.0
                grid[i, q, r] = value
    return grid


def fix_grids():
    grid_dir = Config.get('PATH', 'grid_dir')
    for i, grid_name in enumerate(os.listdir(grid_dir)):
        if not grid_name.endswith('.npy'):
            continue
        grid_fh = os.path.join(grid_dir, grid_name)
        print('fixing :', grid_fh)
        np.save(grid_fh, fix_rslice(np.load(grid_fh)))


def fix_grids2():
    slice_dict = {

        0: 'n_stars',
        1: 'apparent_mag',
        2: 'absolute_mag',
        3: 'mag_limit_min',
        4: 'mag_lim_med',
        5: 'mag_lim_max',
        6: 'sat_age',
        7: 'n_stars_total',
        8: 'sat_age_total'

    }
    table_path = Config.get('PATH', 'table_dir')
    table_bin_output_dir = os.path.join(table_path, 'table_bin_output')
    table_dir = os.path.join(table_path, 'fix_grids2_output')
    if not os.path.isdir(table_dir):
        os.mkdir(table_dir)

    grid_dir = Config.get('PATH', 'grid_dir')
    for grid_name in os.listdir(grid_dir):
        if not grid_name.endswith('.npy'):
            continue

        grid_fh = os.path.join(grid_dir, grid_name)
        grid = np.load(grid_fh)
        print('starting', grid_fh)
        data_table = Table()
        first_slice = True
        temp_x_idx_list = []
        temp_y_idx_list = []
        x_and_y = []
        zeros = []
        for slice_num in slice_dict.keys():

            print('slice:', slice_num, ':', slice_dict[slice_num])
            temp_value_list = []
            
            for y in range(grid.shape[0]):
                for x in range(grid.shape[1]):
                    temp_value_list.append(grid[y, x, slice_num])
                    if first_slice == True:
                        temp_x_idx_list.append(x)
                        temp_y_idx_list.append(y)
                        x_and_y.append((y, x))
                        zeros.append(0.)

            print('adding column data')
            data_table.add_column(Column(data=temp_value_list, name=slice_dict[slice_num], dtype=grid[
                                  :, :, slice_num].dtype, description=slice_dict[slice_num], meta={'slice': slice_num}))
            if first_slice == True:
                data_table.add_column(Column(data=temp_x_idx_list, name='X_idx', dtype=grid[
                                      :, :, slice_num].dtype, description='y indices'))
                data_table.add_column(Column(data=temp_y_idx_list, name='Y_idx', dtype=grid[
                                      :, :, slice_num].dtype, description='x indices'))
                first_slice = False

        table_name = grid_name.rstrip('.npy') + '_' + 'fix_grids2_table.hdf5'
        table_fh = os.path.join(table_dir, table_name)
        print('saving table:', table_fh)
        data_table.write(table_fh, format='hdf5', path='data',
                         overwrite=True, append=True)
    print('done\n')
