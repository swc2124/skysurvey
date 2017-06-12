
from __future__ import division, absolute_import, print_function

import os
import sys
import ConfigParser
import time

from astropy.table import Table
from astropy.table import Column

import numpy as np

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


def table_bin(radius_start=5, radius_end=250, step=1, percent=0.1):
    grid_dir = Config.get('PATH', 'grid_dir')
    table_dir = Config.get('PATH', 'table_dir')
    grid_files = os.listdir(grid_dir)
    for grid_fh in grid_files:
        if not grid_fh.endswith('.npy'):
            continue
        halo, d_mpc, f_type, ending = grid_fh.split(os.path.sep)[-1].split('_')
        d_mpc = float(d_mpc[:3])
        halo_number = int(halo[-2:])
        _grid_fh_ = os.path.join(grid_dir, grid_fh)
        print('laoding grid:', _grid_fh_)
        _grid_ = np.load(_grid_fh_)
        GRID = fix_rslice(_grid_)
        center = GRID.shape[0] / 2.0
        n_stars_total = len(np.load(os.path.join(
            Config.get('PATH', 'halo_dir'), halo, 'px.npy')))
        n_stars_in_grid = int(GRID[:, :, 0].sum())
        first_time = True
        for radius_Kpc in range(radius_start, radius_end, step):

            region = (radius_Kpc * percent)
            r_in = radius_Kpc - region
            r_out = radius_Kpc + region

            selected_boxes = np.nonzero(np.logical_and(
                GRID[:, :, 8] < r_out, GRID[:, :, 8] >= r_in))
            x, y = selected_boxes

            PHI = np.arctan2(y - center, x - center)
            BOXES = GRID[:, :, 0][selected_boxes]
            B_RADS = GRID[:, :, 9][selected_boxes]
            ALL_RADS = np.sqrt((np.square(x - center) + np.square(y - center)))
            AGE = GRID[:, :, 6][selected_boxes]

            mu = BOXES.mean()
            sigma_mu = mu
            XBOX = (BOXES - mu) / sigma_mu

            sigma_sqrt = np.sqrt(BOXES.mean())
            XBOX_sqrt = (BOXES - mu) / sigma_sqrt

            sigma_std = BOXES.std()
            XBOX_std = (BOXES - mu) / sigma_std

            n_boxes = len(BOXES)
            str_nboxes = str(n_boxes)

            if first_time == True:
                Boxes = Column(data=BOXES,
                               name='Boxes', dtype=np.float64, length=1,
                               description='Number of stars per grid box')
                Percent_Boxes = Column(data=BOXES / float(n_boxes),
                                       name='Percent_Boxes', dtype=np.float64, length=1,
                                       description='Boxes / n_boxes')
                LogBoxes = Column(data=np.log10(BOXES),
                                  name='LogBoxes', dtype=np.float64, length=1,
                                  description='Log(number of stars) per grid box')
                Xbox = Column(data=XBOX,
                              name='Xbox', dtype=np.float64, length=1,
                              description='(Boxes - mu) / sigma where sigma=Boxes.mean()=mu')
                Xbox_sqrt = Column(data=XBOX_sqrt,
                                   name='Xbox_sqrt', dtype=np.float64, length=1,
                                   description='(Boxes - mu) / sigma where sigma=sqrt(Boxes.mean())')
                Xbox_std = Column(data=XBOX_std,
                                  name='Xbox_std', dtype=np.float64, length=1,
                                  description='(Boxes - mu) / sigma where sigma=Boxes.std()')
                Phi = Column(data=PHI,
                             name='Phi', dtype=np.float64, length=1,
                             description='arctan(y - center / x - center)')
                Rads = Column(data=B_RADS,
                              name='Rads', dtype=np.float64, length=1,
                              description='Radius of box in Kpc')
                All_rads = Column(data=ALL_RADS,
                                  name='All_rads', dtype=np.float64, length=1,
                                  description='Radius of all boxes in Kpc')
                Age = Column(data=AGE,
                             name='t_accretion', dtype=np.float64, length=1,
                             description='Time since accretion Gyr')
                X_idx = Column(data=x - center,
                               name='X_idx', dtype=np.int64, length=1,
                               description='x indices on the grid')
                Y_idx = Column(data=y - center,
                               name='Y_idx', dtype=np.int64, length=1,
                               description='y indices on the grid')
                R_kpc = Column(data=np.ones_like(x) * radius_Kpc,
                               name='R_kpc', dtype=np.float64, length=1,
                               description='Selected radius in Kpc')
                Mu = Column(data=np.ones_like(x) * mu,
                            name='Mu', dtype=np.float64, length=1,
                            description='mu = Boxes.mean()')
                Halo = Column(data=np.ones_like(x) * halo_number,
                              name='Halo', dtype=np.int64, length=1,
                              description='Halo number ie. halo02 -> 2')
                Nstars_in_grid = Column(data=np.ones_like(x) * n_stars_in_grid,
                                        name='Nstars_in_grid', dtype=np.int64, length=1,
                                        description='number of stars present it the whole grid')
                Nstars_total = Column(data=np.ones_like(x) * n_stars_total,
                                      name='Nstars_total', dtype=np.int64, length=1,
                                      description='Total strs in halo')
                Nstars_r = Column(data=np.ones_like(x) * int(BOXES.sum()),
                                  name='Nstars_r', dtype=np.int64, length=1,
                                  description='Total strs in radial selection')
                N_boxes = Column(data=np.ones_like(x) * len(BOXES),
                                 name='N_boxes', dtype=np.int64, length=1,
                                 description='Total number of grid boxes at this R')
                table = Table()
                table.add_columns(
                    [Boxes, Percent_Boxes, LogBoxes, Xbox, Xbox_sqrt, Xbox_std,
                     Phi, Rads, All_rads, Age, X_idx, Y_idx, R_kpc, Mu, Halo, Nstars_in_grid,
                     Nstars_total, Nstars_r, N_boxes], indexes=None,
                    copy=False, rename_duplicate=False)

                #sys.stdout.write('table started\n')
                # sys.stdout.flush()
                first_time = False

            else:
                Boxes = Column(data=np.concatenate((table['Boxes'], BOXES)),
                               name='Boxes', dtype=np.float64, length=1,
                               description='Number of stars per grid box')
                Percent_Boxes = Column(data=np.concatenate((table['Percent_Boxes'], BOXES / float(n_boxes))),
                                       name='Percent_Boxes', dtype=np.float64, length=1,
                                       description='Boxes / n_boxes')
                LogBoxes = Column(data=np.concatenate((table['LogBoxes'], np.log10(BOXES))),
                                  name='LogBoxes', dtype=np.float64, length=1,
                                  description='Log(number of stars) per grid box')
                Xbox = Column(data=np.concatenate((table['Xbox'], XBOX)),
                              name='Xbox', dtype=np.float64, length=1,
                              description='(Boxes - mu) / sigma where sigma=Boxes.mean()=mu')
                Xbox_sqrt = Column(data=np.concatenate((table['Xbox_sqrt'], XBOX_sqrt)),
                                   name='Xbox_sqrt', dtype=np.float64, length=1,
                                   description='(Boxes - mu) / sigma where sigma=sqrt(Boxes.mean())')
                Xbox_std = Column(data=np.concatenate((table['Xbox_std'], XBOX_std)),
                                  name='Xbox_std', dtype=np.float64, length=1,
                                  description='(Boxes - mu) / sigma where sigma=Boxes.std()')
                Phi = Column(data=np.concatenate((table['Phi'], PHI)),
                             name='Phi', dtype=np.float64, length=1,
                             description='arctan(y - center / x - center)')
                Rads = Column(data=np.concatenate((table['Rads'], B_RADS)),
                              name='Rads', dtype=np.float64, length=1,
                              description='Radius of box in Kpc')
                All_rads = Column(data=np.concatenate((table['All_rads'], ALL_RADS)),
                                  name='All_rads', dtype=np.float64, length=1,
                                  description='Radius of all boxes in Kpc')
                Age = Column(data=np.concatenate((table['t_accretion'], AGE)),
                             name='t_accretion', dtype=np.float64, length=1,
                             description='Time since accretion Gyr')
                X_idx = Column(data=np.concatenate((table['X_idx'], (x - center))),
                               name='X_idx', dtype=np.int64, length=1,
                               description='x indices on the grid')
                Y_idx = Column(data=np.concatenate((table['Y_idx'], (y - center))),
                               name='Y_idx', dtype=np.int64, length=1,
                               description='y indices on the grid')
                R_kpc = Column(data=np.concatenate((table['R_kpc'], np.ones_like(x) * radius_Kpc)),
                               name='R_kpc', dtype=np.float, length=1,
                               description='Base radius in Kpc')
                Mu = Column(data=np.concatenate((table['Mu'], np.ones_like(x) * mu)),
                            name='Mu', dtype=np.float64, length=1,
                            description='mu = Boxes.mean()')
                Halo = Column(data=np.concatenate((table['Halo'], np.ones_like(x) * halo_number)),
                              name='Halo', dtype=np.int64, length=1,
                              description='Halo number ie. halo02 -> 2')
                Nstars_in_grid = Column(data=np.concatenate((table['Nstars_in_grid'], np.ones_like(x) * n_stars_in_grid)),
                                        name='Nstars_in_grid', dtype=np.int64, length=1,
                                        description='number of stars present it the whole grid')
                Nstars_total = Column(data=np.concatenate((table['Nstars_total'], np.ones_like(x) * n_stars_total)),
                                      name='Nstars_total', dtype=np.int64, length=1,
                                      description='Total strs in halo')
                Nstars_r = Column(data=np.concatenate((table['Nstars_r'], np.ones_like(x) * int(BOXES.sum()))),
                                  name='Nstars_r', dtype=np.int64, length=1,
                                  description='Total strs in radial selection')
                N_boxes = Column(data=np.concatenate((table['N_boxes'], np.ones_like(x) * len(BOXES))),
                                 name='N_boxes', dtype=np.int64, length=1,
                                 description='Total number of grid boxes at this R')

                temp_table = Table()
                temp_table.add_columns(
                    [Boxes, Percent_Boxes, LogBoxes, Xbox, Xbox_sqrt, Xbox_std,
                     Phi, Rads, All_rads, Age, X_idx, Y_idx, R_kpc, Mu, Halo, Nstars_in_grid,
                     Nstars_total, Nstars_r, N_boxes], indexes=None,
                    copy=False, rename_duplicate=False)

                table = temp_table
                del temp_table

            msg = '\r\r[' + halo + '] radius: ' + str(radius_Kpc) + '/' + str(radius_end) + ' Kpc - [ ' + str(
                float((1e2 * (radius_Kpc - radius_start)) / radius_end - radius_start)) + '% ] '
            sys.stdout.write(msg)
            sys.stdout.flush()

        table_fh = str(d_mpc) + 'Mpc_' + f_type + \
            '_' + halo + '_table_bin_table'
        target_dir = os.path.join(table_dir, 'table_bin_output')
        file_name = table_fh + '.hdf5'
        if not os.path.isdir(target_dir):
            os.mkdir(target_dir)
        fh = os.path.join(target_dir, file_name)
        table.write(fh, format='hdf5', path='data',
                    overwrite=True, append=True)
        sys.stdout.write(
            msg + '[' + halo.upper() + ' DONE] [' + time.ctime() + '] - saved to:' + fh + '\n')
        sys.stdout.flush()

        del table
