
from __future__ import division, absolute_import, print_function

import os
import sys
import ConfigParser
import time

from astropy.table import Table
from astropy.table import Column
from astropy.table import vstack

from astropy.units import kiloparsec
from astropy.units import radian
from astropy.units import degree
from astropy.units import gigayear
from astropy.units import ABmag
from astropy.units import mag

import numpy as np

from skysurvey.new_config import SYS_CFG_FNAME
import skysurvey

sys_config_fh = os.path.join(os.path.dirname(
    os.path.realpath(skysurvey.__file__)), SYS_CFG_FNAME)
SysConfig = ConfigParser.ConfigParser()
SysConfig.read(sys_config_fh)
config_fh = SysConfig.get('skysurvey_global_settings', 'config_fh')
Config = ConfigParser.ConfigParser()
Config.read(config_fh)

VERBOSE = Config.get('Global', 'verbose')


def window_size():

    from ctypes import windll, create_string_buffer

    # stdin handle is -10
    # stdout handle is -11
    # stderr handle is -12

    h = windll.kernel32.GetStdHandle(-12)
    csbi = create_string_buffer(22)
    res = windll.kernel32.GetConsoleScreenBufferInfo(h, csbi)

    if res:
        import struct
        (bufx, bufy,
         curx, cury,
         wattr,
         left, top, right, bottom,
         maxx, maxy) = struct.unpack("hhhhHhhhhhh", csbi.raw)
        sizex = right - left + 1
        sizey = bottom - top + 1
    else:
        # can't determine actual size - return default values
        sizex, sizey = 80, 25

    return sizex, sizey


def fix_rslice(grid, rslices=[4]):

    center = grid.shape[0] / 2
    ratio = (1.0 * Config.getint('grid_options', 'size')) / \
            (2.0 * Config.getfloat('halo_data_settings', 'radial_cut'))
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


def table_merge():

    from astropy.table import join as tjoin
    import ebf
    import sys

    grid_size = Config.get('grid_options', 'size')

    table_dir = Config.get('PATH', 'table_dir')

    print('loading satprop.ebf')
    satprop = ebf.read(os.path.join(
        Config.get('PATH', 'data_dir'),
        'satprop.ebf'))

    for halo in Config.get('Default_halos', 'halos').split(','):

        print(halo)

        bin_table_fh = os.path.join(
            table_dir,
            'table_bin_output',
            grid_size,
            '4.0Mpc_h158_' + halo + '_table_bin_table.hdf5')

        bin_table = Table.read(
            bin_table_fh,
            format='hdf5',
            path='data')

        print(bin_table.info)

        bin_table.keep_columns(
            ['x_int', 'y_int', 'Xbox',
             't_accretion', 'Rads', 'Phi', 'R_kpc'])

        print(bin_table.info(['attributes', 'stats']))
        bin_table.pprint(
            max_lines=100,
            max_width=window_size()[0],
            show_name=True,
            show_unit=True,
            show_dtype=True,
            align=None)

        stellar_table_fh = os.path.join(
            table_dir,
            'spinbin_output',
            grid_size,
            halo + '_4.0Mpc_h158_table.hdf5')

        stellar_table = Table.read(
            stellar_table_fh,
            format='hdf5',
            path='data')

        print(stellar_table.info)
        d_keys = ['x_int', 'y_int', 'satids', 'feh', 'age']
        stellar_table.keep_columns(d_keys)

        print('removing bound sats')
        satids = stellar_table.meta['satids']
        for satid in satids:
            if satprop['bsat'][satid]:
                stellar_table.remove_rows(
                    np.nonzero(stellar_table['satids'] == satid))
                print(' -> removing satid:', satid)
        print('done')

        print(stellar_table.info(['attributes', 'stats']))
        stellar_table.pprint(
            max_lines=100,
            max_width=window_size()[0],
            show_name=True,
            show_unit=True,
            show_dtype=True,
            align=None)

        print('joining tables for', halo, '...')
        new_table = tjoin(
            left=stellar_table,
            right=bin_table,
            join_type='right',
            keys=['x_int', 'y_int'])
        print('done')

        print(new_table.info(['attributes', 'stats']))
        new_table.pprint(
            max_lines=100,
            max_width=window_size()[0],
            show_name=True,
            show_unit=True,
            show_dtype=True,
            align=None)

        new_table_fh = os.path.join(
            table_dir,
            'merged_tables',
            halo + '.hdf5')

        print('writing new table to: ', new_table_fh)
        new_table.write(
            new_table_fh,
            format='hdf5',
            path='data',
            compression=True,
            overwrite=True,
            serialize_meta=True)
        print('done')


def table_bin(size=None, radius_start=1, radius_end=275,
              step=1, percent=0.025):

    if size == None:
        size = Config.get('grid_options', 'size')

    grid_dir = os.path.join(
        Config.get('PATH', 'grid_dir'),
        size)

    table_dir = os.path.join(
        Config.get('PATH', 'table_dir'),
        'table_bin_output', str(size))

    if not os.path.isdir(table_dir):
        os.mkdir(table_dir)

    grid_files = os.listdir(grid_dir)
    for grid_fh in grid_files:
        if not grid_fh.endswith('.npy'):
            continue

        fh_list = grid_fh.split(os.path.sep)[-1].split('_')
        halo, d_mpc, f_type, ending = fh_list
        d_mpc = float(d_mpc[:3])
        halo_number = int(halo[-2:])
        _grid_fh_ = os.path.join(grid_dir, grid_fh)

        _grid_ = np.load(_grid_fh_)
        GRID = fix_rslice(_grid_)
        center = GRID.shape[0] / 2.0

        n_stars_total = len(
            np.load(
                os.path.join(
                    Config.get('PATH', 'halo_dir'),
                    halo,
                    'px.npy')))

        n_stars_in_grid = GRID[:, :, 0].sum()
        prefix = str(d_mpc) + 'Mpc_' + f_type + '_' + halo
        table_fh = prefix + '_table_bin_table'
        target_dir = table_dir
        file_name = table_fh + '.hdf5'

        if not os.path.isdir(target_dir):
            os.mkdir(target_dir)

        fh = os.path.join(target_dir, file_name)

        table = Table()
        table.meta['_grid_fh_'] = _grid_fh_
        table.meta['grid_shape'] = GRID.shape
        table.meta['halo'] = halo
        table.meta['table_dir'] = table_dir
        table.meta['f_type'] = f_type
        table.meta['d_mpc'] = d_mpc
        table.meta['table_bin_output_fh'] = table_fh
        table.meta['full_path'] = fh
        table.meta['table_bin_creation_time'] = time.ctime()

        for radius_Kpc in range(radius_start, radius_end, step):
            region = (radius_Kpc * percent)
            r_in = radius_Kpc - region
            r_out = radius_Kpc + region

            selected_boxes = np.nonzero(
                np.logical_and(
                    GRID[:, :, 14] < r_out,
                    GRID[:, :, 14] >= r_in))

            x, y = selected_boxes
            PHI = np.arctan2(y - center, x - center)
            BOXES = GRID[:, :, 0][selected_boxes]
            B_RADS = GRID[:, :, 14][selected_boxes]
            AGE = GRID[:, :, 6][selected_boxes]
            mu = BOXES.mean()
            sigma_mu = mu
            XBOX = (BOXES - mu) / sigma_mu
            sigma_sqrt = np.sqrt(BOXES.mean())
            XBOX_sqrt = (BOXES - mu) / sigma_sqrt
            sigma_std = BOXES.std()
            XBOX_std = (BOXES - mu) / sigma_std
            n_boxes = len(BOXES)

            columns = [

                Column(
                    data=BOXES.astype(np.float32),
                    name='Boxes',
                    description='Number of stars per grid box'),

                Column(
                    data=np.divide(
                        BOXES,
                        np.float64(n_boxes)).astype(np.float32),
                    name='Percent_Boxes',
                    description='Boxes / n_boxes'),

                Column(
                    data=np.log10(BOXES).astype(np.float16),
                    name='LogBoxes',
                    description='Log(number of stars) per grid box'),

                Column(
                    data=XBOX.astype(np.float16),
                    name='Xbox',
                    description='(Boxes - mu) / mu'),

                Column(
                    data=XBOX_sqrt.astype(np.float16),
                    name='Xbox_sqrt',
                    description='(Boxes - mu) / sigma where sigma=sqrt(Boxes.mean())',
                    unit='float16'),

                Column(
                    data=XBOX_std.astype(np.float16),
                    name='Xbox_std',
                    description='(Boxes - mu) / Boxes.std()',
                    unit='float16'),

                Column(
                    data=PHI.astype(np.float16),
                    name='Phi',
                    description='arctan(y - center / x - center)',
                    unit='radian'),

                Column(
                    data=B_RADS.astype(np.uint16),
                    name='Rads',
                    description='Radius of box in Kpc',
                    unit='kiloparsec'),

                Column(
                    data=AGE.astype(np.float16),
                    name='t_accretion',
                    description='Time since accretion Gyr',
                    unit='gigayear'),

                Column(
                    data=(x - center).astype(np.int16),
                    name='X',
                    description='x indices on the grid - center',
                    unit='int16'),

                Column(
                    data=(y - center).astype(np.int16),
                    name='Y',
                    description='y indices on the grid - center',
                    unit='int16'),

                Column(
                    data=x.astype(np.uint16),
                    name='x_int',
                    description='x indices on the grid',
                    unit='uint16'),

                Column(
                    data=y.astype(np.uint16),
                    name='y_int',
                    description='y ints for merging table',
                    unit='uint16'),

                Column(
                    data=np.ones(
                        n_boxes,
                        dtype=np.uint16) * radius_Kpc,
                    name='R_kpc',
                    description='Selected radius in Kpc',
                    unit='kiloparsec'),

                Column(
                    data=np.ones(
                        n_boxes,
                        dtype=np.float32) * mu,
                    name='Mu',
                    description='mu = Boxes.mean()',
                    unit='float32'),

                Column(
                    data=np.ones(
                        n_boxes,
                        dtype=np.uint32) * BOXES.sum(),
                    name='Nstars_r',
                    description='Total strs in radial selection',
                    unit='uint32'),

                Column(
                    data=np.ones(
                        n_boxes,
                        dtype=np.float16) * np.float16(
                        np.log10(
                            BOXES.sum())),
                    name='Log10(Nstars)_r',
                    description='Total strs in radial selection',
                    unit='float16'),

                Column(
                    data=np.ones(
                        n_boxes,
                        dtype=np.uint32) * n_boxes,
                    name='N_boxes',
                    description='Total number of grid boxes at this R',
                    unit='uint32'),

                Column(
                    data=np.ones(
                        n_boxes,
                        dtype=np.float16) * np.float16(
                        np.log10(n_boxes)),
                    name='Log10(N_boxes_tot)',
                    description='Total number of grid boxes at this R',
                    unit='float16'),

                Column(
                    data=np.ones(
                        n_boxes,
                        dtype=np.uint32) * len(
                        np.nonzero(BOXES == 0.0)[0]),
                    name='N_boxes_empty',
                    unit='uint32'),

                Column(
                    data=np.ones(
                        n_boxes,
                        dtype=np.float16) * np.float16(
                        np.log10(
                            len(np.nonzero(BOXES == 0.0)[0]))),
                    name='Log10(N_boxes)_empty',
                    unit='float16'),

                Column(
                    data=np.ones(
                        n_boxes,
                        dtype=np.float16) * np.float16(
                        np.log10(
                            len(np.nonzero(BOXES >= 1.0)[0]))),
                    name='Log10(N_boxes)_full',
                    unit='float16'),

                Column(
                    data=np.ones(
                        n_boxes,
                        dtype=np.uint32) * len(
                        np.nonzero(BOXES >= 1.0)[0]),
                    name='N_boxes_full',
                    unit='uint32')]

            temp_table = Table()
            temp_table.add_columns(columns)
            table = vstack([table, temp_table])

            if VERBOSE:
                dr = radius_Kpc - radius_start
                dr_tot = radius_end - radius_start
                sfig = (1.0 * (dr + 1) / (1.0 * (dr_tot)))
                adj_fig = round(1e2 * sfig, 2)

                msg_00 = '\r\r' + '[' + halo + '] radius: '
                msg_01 = str(radius_Kpc + 1)
                msg_02 = '/' + str(radius_end)
                msg_03 = ' Kpc - [ ' + str(adj_fig) + '% ] '
                msg = msg_00 + msg_01 + msg_02 + msg_03

                sys.stdout.write(msg)
                sys.stdout.flush()

        sys.stdout.write('\n')
        sys.stdout.flush()

        log_n_stars_in_grid = str(np.log10(n_stars_in_grid))

        table.meta['radius_start'] = str(radius_start)
        table.meta['radius_end'] = str(radius_end)
        table.meta['log_n_stars_in_grid'] = log_n_stars_in_grid

        table.write(
            fh,
            format='hdf5',
            path='data',
            overwrite=True,
            serialize_meta=True)

        table.pprint(
            max_lines=100,
            max_width=window_size()[0],
            show_name=True,
            show_unit=True,
            show_dtype=True,
            align=None)

        sys.stdout.write('\n')
        sys.stdout.flush()

        if VERBOSE:

            m0 = '[' + halo.upper() + ' DONE] ['
            m1 = time.ctime() + '] - saved to:'
            m2 = fh + '\n'

            sys.stdout.write(msg + m0 + m1 + m2)
            sys.stdout.flush()


def table_binlite(size=None, radius_start=0, radius_end=220, step=1, percent=0.1):
    if size == None:
        size = Config.get('grid_options', 'size')
    grid_dir = os.path.join(Config.get('PATH', 'grid_dir'), str(size))
    table_dir = os.path.join(
        Config.get('PATH', 'table_dir'),
        'table_bin_output', str(size))
    if not os.path.isdir(table_dir):
        os.mkdir(table_dir)
    grid_files = os.listdir(grid_dir)
    for grid_fh in grid_files:
        if not grid_fh.endswith('.npy'):
            continue

        fh_list = grid_fh.split(os.path.sep)[-1].split('_')
        halo, d_mpc, f_type, ending = fh_list
        d_mpc = float(d_mpc[:3])
        halo_number = int(halo[-2:])

        _grid_fh_ = os.path.join(grid_dir, grid_fh)
        GRID = fix_rslice(np.load(_grid_fh_))

        center = GRID.shape[0] / 2.0
        n_stars_total = len(
            np.load(
                os.path.join(
                    Config.get('PATH', 'halo_dir'),
                    halo,
                    'px.npy')))
        n_stars_in_grid = GRID[:, :, 0].sum()

        table_fh = str(d_mpc) + 'Mpc_' + f_type + \
            '_' + halo + '_table_bin_table'
        target_dir = table_dir
        file_name = table_fh + '.hdf5'
        if not os.path.isdir(target_dir):
            os.mkdir(target_dir)
        fh = os.path.join(target_dir, file_name)
        table = Table()
        table.meta['_grid_fh_'] = _grid_fh_
        table.meta['grid_shape'] = GRID.shape
        table.meta['halo'] = halo
        table.meta['table_dir'] = table_dir
        table.meta['f_type'] = f_type
        table.meta['d_mpc'] = d_mpc
        table.meta['table_bin_output_fh'] = table_fh
        table.meta['full_path'] = fh
        table.meta['table_bin_creation_time'] = time.ctime()

        for radius_Kpc in range(radius_start, radius_end, step):
            region = (radius_Kpc * percent)
            r_in = radius_Kpc - region
            r_out = radius_Kpc + region

            selected_boxes = np.nonzero(
                np.logical_and(
                    GRID[:, :, 4] < r_out,
                    GRID[:, :, 4] >= r_in))

            #print(selected_boxes )
            x, y = selected_boxes

            PHI = np.arctan2(y - center, x - center)
            BOXES = GRID[:, :, 0][selected_boxes]
            B_RADS = GRID[:, :, 4][selected_boxes]
            AGE = GRID[:, :, 2][selected_boxes]
            mu = BOXES.mean()
            sigma_mu = mu
            XBOX = (BOXES - mu) / sigma_mu

            n_boxes = len(BOXES)

            columns = [
                Column(
                    data=BOXES.astype(np.float32),
                    name='Boxes',
                    description='Number of stars per grid box'),
                Column(
                    data=np.divide(
                        BOXES,
                        np.float64(n_boxes)).astype(np.float32),
                    name='Percent_Boxes',
                    description='Boxes / n_boxes'),
                Column(
                    data=np.log10(BOXES).astype(np.float16),
                    name='LogBoxes',
                    description='Log(number of stars) per grid box'),
                Column(
                    data=XBOX.astype(np.float16),
                    name='Xbox',
                    description='(Boxes - mu) / sigma where sigma=Boxes.mean()=mu'),
                Column(
                    data=np.log10(XBOX.astype(np.float16)),
                    name='LogXbox',
                    description='Log10((Boxes - mu) / sigma) where sigma=Boxes.mean()=mu'),
                Column(
                    data=PHI.astype(np.float16),
                    name='Phi',
                    description='arctan(y - center / x - center)',
                    unit='radian'),
                Column(
                    data=B_RADS.astype(np.uint16),
                    name='Rads',
                    description='Radius of box in Kpc',
                    unit='kiloparsec'),
                Column(
                    data=AGE.astype(np.float16),
                    name='t_accretion',
                    description='Time since accretion Gyr',
                    unit='gigayear'),
                Column(
                    data=x.astype(np.uint16),
                    name='x_int',
                    description='x indices on the grid',
                    unit='uint16'),
                Column(
                    data=y.astype(np.uint16),
                    name='y_int',
                    description='y ints for merging table',
                    unit='uint16'),
                Column(
                    data=np.ones(n_boxes).astype(np.uint16) * radius_Kpc,
                    name='R_kpc',
                    description='Selected radius in Kpc',
                    unit='kiloparsec')]

            temp_table = Table()
            temp_table.add_columns(columns)
            table = vstack([table, temp_table])

            if VERBOSE:
                msg = '\r\r[' + halo + '] radius: ' + str(radius_Kpc + 1) + '/' + str(
                    radius_end) + ' Kpc - [ ' + str(round(1e2 * (1.0 * (radius_Kpc - radius_start + 1) / (1.0 * (radius_end - radius_start))), 2)) + '% ] '
                sys.stdout.write(msg)
                sys.stdout.flush()
                #print(radius_Kpc, n_boxes, len(np.unique(table['R_kpc'])))
        sys.stdout.write('\n')
        sys.stdout.flush()

        table.meta['radius_start'] = str(radius_start)
        table.meta['radius_end'] = str(radius_end)
        table.meta['log_n_stars_in_grid'] = str(np.log10(n_stars_in_grid))
        table.write(fh, format='hdf5', path='data',
                    overwrite=True, serialize_meta=True)
        table.pprint(max_lines=100, max_width=window_size()[
            0], show_name=True, show_unit=True, show_dtype=True, align=None)
        sys.stdout.write('\n')
        sys.stdout.flush()
        if VERBOSE:
            sys.stdout.write(
                msg + '[' + halo.upper() + ' DONE] [' + time.ctime() + '] - saved to:' + fh + '\n')
            sys.stdout.flush()
