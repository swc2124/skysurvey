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

from os.path import join as os_join
from os.path import abspath
from os.path import curdir
from warnings import warn

from numpy import sqrt
from numpy import square
from numpy import asarray
from numpy import int64
from numpy import float64
from numpy import save
from numpy import int32
from numpy import float64

from .options import halos
from .options import distance_Mpc
from .options import apparent_mag_limits
from .options import grid_dir
from .options import filter_type

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


def _spin(halo, _m_lims=apparent_mag_limits, _distance_mpc=distance_Mpc, _filter_type=filter_type):
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
    ab_mag_arr = load_ab_mags(
        halo, 
        f_type=_filter_type
        )
    d_limits = find_dlims(
        ab_mag_arr, 
        asarray(
            [
            calculate_abs_mag(
                distance=_distance_mpc, 
                _f_type=_filter_type
                )
            ], 
            dtype=float64
            )
        )
    r_px, r_py, r_pz = trippel_rotate(load_positions(halo, d_limits))
    integer_x_arr, intiger_y_arr = integerize(r_px, r_py)
    return bf(
        integer_x_arr, 
        intiger_y_arr, 
        ab_mag_arr, 
        apparent_magnitude(
            ab_mag_arr, 
            _distance_mpc
            ),
        sqrt(
            square(r_px) + square(r_py)
            ), 
        _m_lims, 
        load_satid(
            halo, 
            d_limits
            )
        )


def spinall(path=grid_dir, m_lims=apparent_mag_limits, d_mpc=distance_Mpc, f_type=filter_type):
    '''[summary]
    
    [description]
    
    Keyword Arguments:
        path {[type]} -- [description] (default: {grid_dir})
        m_lims {[type]} -- [description] (default: {apparent_mag_limits})
        d_mpc {[type]} -- [description] (default: {distance_Mpc})
        f_type {[type]} -- [description] (default: {filter_type})
    '''
    for halo in halos:
        filename = os_join(path, halo + '_' + str(d_mpc) + 'Mpc_' + f_type + '_grid')
        print('--> [ STARTING ' + halo + ' ]\nspining and binning halo ' + halo)
        grid = _spin(
            halo, 
            _m_lims=m_lims,
            _distance_mpc=d_mpc, 
            _filter_type=f_type
            )
        save(filename, grid)
        print(halo + ' saved : ' + filename + '\n---------------------------------\n')


def spinone(halo, path=grid_dir, m_lims=apparent_mag_limits, d_mpc=distance_Mpc, f_type=filter_type):
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
    filename = os_join(path, halo + '_' + str(d_mpc) + 'Mpc_' + f_type + '_grid')
    print('spining and binning halo ' + halo)
    grid = _spin(
        halo, 
        _m_lims=m_lims,
        _distance_mpc=d_mpc,  
        _filter_type=f_type
        )
    save(filename, grid)
    print(halo + ' saved : ' + filename + '\n---------------------------------\n')


def nospin_binall(path=grid_dir, m_lims=apparent_mag_limits, d_mpc=distance_Mpc, f_type=filter_type):
    '''[summary]
    
    [description]
    
    Keyword Arguments:
        path {[type]} -- [description] (default: {grid_dir})
        m_lims {[type]} -- [description] (default: {apparent_mag_limits})
        d_mpc {[type]} -- [description] (default: {distance_Mpc})
        f_type {[type]} -- [description] (default: {filter_type})
    '''
    for halo in halos:
        filename = os_join(path, halo + '_' + str(d_mpc) + 'Mpc_' + f_type + '_grid')
        ab_mag_arr = load_ab_mags(halo)
        d_limits = find_dlims(
            ab_mag_arr, 
            asarray(
                [
                calculate_abs_mag(
                    distance=d_mpc, 
                    _f_type=f_type
                    )
                ], 
                dtype=float64
                )
            )
        px, py, pz = load_positions(halo, d_limits)
        integer_x_arr, intiger_y_arr = integerize(px, py)
        save(
            filename, 
            bf(
                integer_x_arr, 
                intiger_y_arr, 
                ab_mag_arr, 
                apparent_magnitude(ab_mag_arr, d_mpc), 
                sqrt(square(px) + square(py)), 
                m_lims, 
                load_satid(halo, d_limits)
                )
            )
        print(halo + ' saved : ' + filename + '\n---------------------------------\n')