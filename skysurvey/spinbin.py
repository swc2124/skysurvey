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
from numpy import float64
from numpy import save
from numpy import int32
from numpy import float64

from .functions import load_ab_mags
from .functions import load_positions
from .functions import apparent_magnitude
from .functions import calculate_abs_mag
from .functions import load_satid

from c_functions import bin as bf
from c_functions import trippel_rotate
from c_functions import rotate
from c_functions import find_dlims
from c_functions import integerize

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


def display(arr, to_print=False, log=False):
    if log:
        arr = np.log10(arr[arr != 0.0])
    msg_0 = '\nmin:' + str(round(arr.min(), 2))
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
        d_table.meta['creation_time'] = time.ctime()
        d_table.meta['halo'] = halo
        d_table.meta['abs_mag_limit'] = abs_mag_limit
        d_table.meta['m_lims'] = _m_lims
        d_table.meta['d_mpc'] = _distance_mpc
        d_table.meta['f_type'] = _filter_type

        d_table.add_columns([
            Column(data=ab_mag_arr[d_limits], name='ab_mag_arr', dtype=ab_mag_arr.dtype,
                   length=1, description=display(ab_mag_arr[d_limits])),
            Column(data=app_mags[d_limits], name='app_mags', dtype=app_mags.dtype,
                   length=1, description=display(app_mags[d_limits])),
            Column(data=r_px, name='px', dtype=r_px.dtype,
                   length=1, description=display(r_px)),
            Column(data=r_py, name='py', dtype=r_py.dtype,
                   length=1, description=display(r_py)),
            Column(data=r_pz, name='pz', dtype=r_pz.dtype,
                   length=1, description=display(r_pz)),
            Column(data=integer_x_arr, name='integer_x_arr', dtype=integer_x_arr.dtype,
                   length=1, description=display(integer_x_arr)),
            Column(data=intiger_y_arr, name='intiger_y_arr', dtype=intiger_y_arr.dtype,
                   length=1, description=display(intiger_y_arr)),
            Column(data=proj_rads, name='proj_rads', dtype=proj_rads.dtype,
                   length=1, description=display(proj_rads)),
            Column(data=satids, name='satids', dtype=satids.dtype,
                   length=1, description=display(satids))
        ])

        table_save_path = os.path.join(Config.get(
            'PATH', 'table_dir'), 'spinbin_output')
        if not os.path.isdir(table_save_path):
            os.mkdir(table_save_path)
        table_fh = os.path.join(
            table_save_path, halo + '_' + str(d_mpc) + 'Mpc_' + f_type + '_table.hdf5')

        d_table.meta['table_fh'] = table_fh
        d_table.write(d_table.meta['table_fh'], format='hdf5',
                      path='data', overwrite=True, append=True)
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
            d_table.meta['creation_time'] = time.ctime()
            d_table.meta['grid_fh'] = filename
            d_table.meta['halo'] = halo
            d_table.meta['abs_mag_limit'] = abs_mag_limit
            d_table.meta['m_lims'] = m_lims
            d_table.meta['d_mpc'] = d_mpc
            d_table.meta['f_type'] = f_type

            d_table.add_columns([
                Column(data=ab_mag_arr[d_limits], name='ab_mag_arr', dtype=ab_mag_arr.dtype,
                       length=1, description=display(ab_mag_arr[d_limits])),
                Column(data=app_mags[d_limits], name='app_mags', dtype=app_mags.dtype,
                       length=1, description=display(app_mags[d_limits])),
                Column(data=px, name='px', dtype=px.dtype,
                       length=1, description=display(px)),
                Column(data=py, name='py', dtype=py.dtype,
                       length=1, description=display(py)),
                Column(data=pz, name='pz', dtype=pz.dtype,
                       length=1, description=display(pz)),
                Column(data=integer_x_arr, name='integer_x_arr', dtype=integer_x_arr.dtype,
                       length=1, description=display(integer_x_arr)),
                Column(data=intiger_y_arr, name='intiger_y_arr', dtype=intiger_y_arr.dtype,
                       length=1, description=display(intiger_y_arr)),
                Column(data=proj_rads, name='proj_rads', dtype=proj_rads.dtype,
                       length=1, description=display(proj_rads)),
                Column(data=satids, name='satids', dtype=satids.dtype,
                       length=1, description=display(satids))
            ])

            table_save_path = os.path.join(Config.get(
                'PATH', 'table_dir'), 'spinbin_output')
            if not os.path.isdir(table_save_path):
                os.mkdir(table_save_path)
            table_fh = os.path.join(
                table_save_path, halo + '_' + str(d_mpc) + 'Mpc_' + f_type + '_table.hdf5')
            d_table.meta['table_fh'] = table_fh
            d_table.write(d_table.meta[
                          'table_fh'], format='hdf5', path='data', overwrite=True, append=True)
            print(d_table.info())
