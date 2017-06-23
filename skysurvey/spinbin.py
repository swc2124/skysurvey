"""

TODO

-------------------------------------------------------------------------------
 Author     : Sol W. Courtney
 Location   : Columbia U Dept. of Astronomy & Astrophysics NYC, Ny
 Email      : swc2124@Columbia.edu
 Date       : Jan 2017
 Title      : skysurvey/skysurvey/spinbin.py
-------------------------------------------------------------------------------
"""

from __future__ import division, absolute_import, print_function

import os

import time

from numpy import sqrt
from numpy import square
from numpy import asarray
from numpy import int64
from numpy import uint16
from numpy import int16
from numpy import float64
from numpy import save
from numpy import int32
from numpy import float32
from numpy import float16
from numpy import unique

try:
    from .functions import load_ab_mags
    from .functions import load_positions
    from .functions import apparent_magnitude
    from .functions import calculate_abs_mag
    from .functions import load_satid
    from .new_config import SYS_CFG_FNAME

except ValueError, e:
    from skysurvey.functions import load_ab_mags
    from skysurvey.functions import load_positions
    from skysurvey.functions import apparent_magnitude
    from skysurvey.functions import calculate_abs_mag
    from skysurvey.functions import load_satid
    from skysurvey.new_config import SYS_CFG_FNAME

from c_functions import bin as bf
from c_functions import trippel_rotate
from c_functions import rotate
from c_functions import find_dlims
from c_functions import integerize

from astropy.table import Table
from astropy.table import Column

import ConfigParser

import skysurvey

sys_config_fh = os.path.join(os.path.dirname(
    os.path.realpath(skysurvey.__file__)), SYS_CFG_FNAME)
SysConfig = ConfigParser.ConfigParser()
SysConfig.read(sys_config_fh)
config_fh = SysConfig.get('skysurvey_global_settings', 'config_fh')
Config = ConfigParser.ConfigParser()
Config.read(config_fh)


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
        (bufx, bufy, curx, cury, wattr,
         left, top, right, bottom, maxx, maxy) = struct.unpack("hhhhHhhhhhh", csbi.raw)
        sizex = right - left + 1
        sizey = bottom - top + 1
    else:
        sizex, sizey = 80, 25  # can't determine actual size - return default values

    return sizex, sizey


def display(arr, to_print=False, log=False):
    if log:

        if not 'log10' in dir():
            from numpy import log10
        arr = log10(arr[arr != 0.0])
    msg_0 = 'min:' + str(round(arr.min(), 2))
    msg_1 = ' - mean:' + str(round(arr.mean(), 2))
    msg_2 = ' - max:' + str(round(arr.max(), 2))
    msg_3 = ' - std:' + str(round(arr.std(), 2))
    msg_4 = ' - len:' + str(len(arr.flatten()))
    msg = msg_0 + msg_1 + msg_2 + msg_3 + msg_4
    if to_print:
        print(msg)
    else:
        return str(msg)


def _spin(halo, _m_lims=None, _distance_mpc=None, _filter_type=None, table=True):
    '''[summary]

    [description]

    Arguments:
        halo {[type]} -- [description]

    Keyword Arguments:
        _m_lims {[type]} -- [description] (default: {apparent_mag_limits})
        _distance_mpc {[type]} -- [description] (default: {distance_Mpc})
        _filter_type {[type]} -- [description] (default: {filter_type})

    Returns:
        [type] -- [description]
    '''
    if _m_lims == None:
        _m_lims = asarray([float(lim) for name, lim in Config.items(
            'Default_magnitude_limits')], dtype=float64)
    if _distance_mpc == None:
        _distance_mpc = Config.getfloat('Distance', 'd_mpc')
    if _filter_type == None:
        _filter_type = Config.get('Filter', 'filter_type')

    print('prepairing arrays')
    ab_mag_arr = load_ab_mags(halo, f_type=_filter_type)
    abs_mag_limit = calculate_abs_mag(
        distance=_distance_mpc, _f_type=_filter_type)
    d_limits = find_dlims(ab_mag_arr, asarray([abs_mag_limit], dtype=float64))
    app_mags = apparent_magnitude(ab_mag_arr, _distance_mpc)
    r_px, r_py, r_pz = trippel_rotate(load_positions(halo, d_limits))
    integer_x_arr, intiger_y_arr = integerize(r_px, r_py)
    proj_rads = sqrt(square(px) + square(py))
    satids = load_satid(halo, d_limits)
    print('done')

    if table:

        print('starting table')
        d_table = Table()
        d_table.meta['spin_bin_creation_time'] = time.ctime()
        d_table.meta['grid_fh'] = filename
        d_table.meta['halo'] = halo
        d_table.meta['abs_mag_limit'] = str(round(abs_mag_limit, 2))
        d_table.meta['m_lims'] = m_lims
        d_table.meta['d_mpc'] = d_mpc
        d_table.meta['f_type'] = f_type
        d_table.meta['satids'] = unique(satids).tolist()
        d_table.meta['satid_start'] = str(satids.min())
        d_table.meta['satid_end'] = str(satids.max())
        d_table.meta['n_sats'] = len(unique(satids))

        print('table created')
        for k in d_table.meta.keys():
            print('  --> ', k, ':', d_table.meta[k])
        print('loading column data')
        d_table.add_columns([
            Column(data=ab_mag_arr[d_limits].astype(float16), name='ab_mag_arr', description=display(ab_mag_arr[d_limits]), unit='ABmag'),
            Column(data=app_mags[d_limits].astype(float16), name='app_mags', description=display(app_mags[d_limits]), unit='mag'),
            Column(data=r_px.astype(float16), name='r_px', description=display(r_px), unit='kiloparsec'),
            Column(data=r_py.astype(float16), name='r_py', description=display(r_py), unit='kiloparsec'),
            Column(data=r_pz.astype(float16), name='r_pz', description=display(r_pz), unit='kiloparsec'),
            Column(data=integer_x_arr.astype(uint16), name='x_int', unit='int16'),
            Column(data=intiger_y_arr.astype(uint16), name='y_int', unit='int16'),
            Column(data=satids.astype(uint16), name='satids', description=display(satids), unit='int16')
        ])
        print('column data loaded')
        for k in d_table.keys():
            print('  -->', k, ':', d_table[k].description)

        table_save_path = os.path.join(Config.get(
            'PATH', 'table_dir'), 'spinbin_output')
        print('save path :', table_save_path)
        if not os.path.isdir(table_save_path):
            print('making new dir:', table_save_path)
            os.mkdir(table_save_path)
            print('done')

        table_fh = os.path.join(
            table_save_path, halo + '_' + str(d_mpc) + 'Mpc_' + f_type + '_table.hdf5')
        print('table file handel:', table_fh)
        d_table.meta['spinbin_output_fh'] = table_fh
        print('writing table')
        d_table.write(table_fh, format='hdf5', path='data', compression=True,
                      overwrite=True, serialize_meta=True)
        d_table.pprint(max_lines=25, max_width=window_size()[
            0], show_name=True, show_unit=True, show_dtype=True, align=None)
        print('done')
        print(d_table.info())

    return bf(integer_x_arr, intiger_y_arr, ab_mag_arr, app_mags, proj_rads, _m_lims, satids)


def spinall(path=None, m_lims=None, d_mpc=None, f_type=None):
    '''[summary]

    [description]

    Keyword Arguments:
        path {[type]} -- [description] (default: {grid_dir})
        m_lims {[type]} -- [description] (default: {apparent_mag_limits})
        d_mpc {[type]} -- [description] (default: {distance_Mpc})
        f_type {[type]} -- [description] (default: {filter_type})
    '''
    if path == None:
        path = Config.get('PATH', 'grid_dir')
    if m_lims == None:
        m_lims = asarray([float(lim) for name, lim in Config.items(
            'Default_magnitude_limits')], dtype=float64)
    if d_mpc == None:
        d_mpc = Config.getfloat('Distance', 'd_mpc')
    if f_type == None:
        f_type = Config.get('Filter', 'filter_type')

    for halo in Config.get('Default_halos', 'halos').split(','):
        filename = os.path.join(
            path, halo + '_' + str(d_mpc) + 'Mpc_' + f_type + '_grid')
        print('--> [ STARTING ' + halo + ' ]\nspining and binning halo ' + halo)
        save(filename, _spin(halo, _m_lims=m_lims,
                             _distance_mpc=d_mpc, _filter_type=f_type))
        print(halo + ' saved : ' + filename +
              '\n---------------------------------\n')


def spinone(halo, path=None, m_lims=None, d_mpc=None, f_type=None):
    '''[summary]

    [description]

    Arguments:
        halo {[type]} -- [description]

    Keyword Arguments:
        path {[type]} -- [description] (default: {grid_dir})
        m_lims {[type]} -- [description] (default: {apparent_mag_limits})
        d_mpc {[type]} -- [description] (default: {distance_Mpc})
        f_type {[type]} -- [description] (default: {filter_type})
    '''
    if path == None:
        path = Config.get('PATH', 'grid_dir')
    if m_lims == None:
        m_lims = asarray([float(lim) for name, lim in Config.items(
            'Default_magnitude_limits')], dtype=float64)
    if d_mpc == None:
        d_mpc = Config.getfloat('Distance', 'd_mpc')
    if f_type == None:
        f_type = Config.get('Filter', 'filter_type')

    filename = os.path.join(
        path, halo + '_' + str(d_mpc) + 'Mpc_' + f_type + '_grid')
    print('spining and binning halo ' + halo)
    save(filename, _spin(halo, _m_lims=m_lims,
                         _distance_mpc=d_mpc, _filter_type=f_type))
    print(halo + ' saved : ' + filename +
          '\n---------------------------------\n')


def nospin_binall(path=None, m_lims=None, d_mpc=None, f_type=None, table=True):
    '''[summary]

    [description]

    Keyword Arguments:
        path {[type]} -- [description] (default: {grid_dir})
        m_lims {[type]} -- [description] (default: {apparent_mag_limits})
        d_mpc {[type]} -- [description] (default: {distance_Mpc})
        f_type {[type]} -- [description] (default: {filter_type})
    '''
    if path == None:
        path = os.path.join(Config.get('PATH', 'grid_dir'), Config.get('grid_options', 'size'))
        if not os.path.isdir(path):
            os.mkdir(path)
    if m_lims == None:
        m_lims = asarray([float(lim) for name, lim in Config.items(
            'Default_magnitude_limits')], dtype=float64)
    if d_mpc == None:
        d_mpc = Config.getfloat('Distance', 'd_mpc')
    if f_type == None:
        f_type = Config.get('Filter', 'filter_type')

    for halo in Config.get('Default_halos', 'halos').split(','):

        filename = os.path.join(
            path, halo + '_' + str(d_mpc) + 'Mpc_' + f_type + '_grid')
        print('filename:', filename)
        print('prepairing arrays')
        ab_mag_arr = load_ab_mags(halo)
        abs_mag_limit = calculate_abs_mag(distance=d_mpc, _f_type=f_type)
        d_limits = find_dlims(ab_mag_arr, asarray(
            [abs_mag_limit], dtype=float64))
        app_mags = apparent_magnitude(ab_mag_arr, d_mpc)
        px, py, pz = load_positions(halo, d_limits)
        integer_x_arr, intiger_y_arr = integerize(px, py)
        proj_rads = sqrt(square(px) + square(py))
        satids = load_satid(halo, d_limits)
        save(filename, bf(integer_x_arr, intiger_y_arr,
                          ab_mag_arr, app_mags, proj_rads, m_lims, satids))
        print(halo + ' saved : ' + filename +
              '\n---------------------------------\n')

        if table:

            print('starting table')
            d_table = Table()
            d_table.meta['spin_bin_creation_time'] = time.ctime()
            d_table.meta['grid_fh'] = filename
            d_table.meta['halo'] = halo
            d_table.meta['abs_mag_limit'] = str(round(abs_mag_limit, 2))
            d_table.meta['m_lims'] = m_lims
            d_table.meta['d_mpc'] = d_mpc
            d_table.meta['f_type'] = f_type
            d_table.meta['satids'] = unique(satids).tolist()
            d_table.meta['satid_start'] = str(satids.min())
            d_table.meta['satid_end'] = str(satids.max())
            d_table.meta['n_sats'] = len(unique(satids))

            print('table created')
            for k in d_table.meta.keys():
                print('  --> ', k, ':', d_table.meta[k])
            print('loading column data')
            d_table.add_columns([
                Column(data=ab_mag_arr[d_limits].astype(float16), name='ab_mag_arr', description=display(ab_mag_arr[d_limits]), unit='ABmag'),
                Column(data=app_mags[d_limits].astype(float16), name='app_mags', description=display(app_mags[d_limits]), unit='mag'),
                Column(data=px.astype(float16), name='px', description=display(px), unit='kiloparsec'),
                Column(data=py.astype(float16), name='py', description=display(py), unit='kiloparsec'),
                Column(data=pz.astype(float16), name='pz', description=display(pz), unit='kiloparsec'),
                Column(data=integer_x_arr.astype(uint16), name='x_int', unit='int16'),
                Column(data=intiger_y_arr.astype(uint16), name='y_int', unit='int16'),
                Column(data=proj_rads.astype(float16), name='r_proj', unit='kpc'),
                Column(data=satids.astype(uint16), name='satids', description=display(satids), unit='int16')
            ])
            print('column data loaded')
            for k in d_table.keys():
                print('  -->', k, ':', d_table[k].description)

            table_save_path = os.path.join(Config.get(
                'PATH', 'table_dir'), 'spinbin_output')
            print('save path :', table_save_path)
            if not os.path.isdir(table_save_path):
                print('making new dir:', table_save_path)
                os.mkdir(table_save_path)
                print('done')

            table_dir = os.path.join(table_save_path, str(Config.get('grid_options', 'size')))
            if not os.path.isdir(table_dir):
                os.mkdir(table_dir)
            table_fh = os.path.join(
                table_save_path, halo + '_' + str(d_mpc) + 'Mpc_' + f_type + '_table.hdf5')
            print('table file handel:', table_fh)
            d_table.meta['spinbin_output_fh'] = table_fh
            print('writing table')
            d_table.write(table_fh, format='hdf5', path='data', compression=True,
                          overwrite=True, serialize_meta=True)
            d_table.pprint(max_lines=25, max_width=window_size()[
                0], show_name=True, show_unit=True, show_dtype=True, align=None)
            print('done')
            print(d_table.info())
