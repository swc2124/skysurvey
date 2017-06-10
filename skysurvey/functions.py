"""TODO"""
from __future__ import division, absolute_import, print_function

import os.path as ospth

import math
import numpy as np

from c_functions import find_dlims

__all__ = ['load_ab_mags', 'apparent_magnitude', 'load_positions', 'load_satid', 'load_allkeys',
           'calculate_abs_mag', 'calculate_nFov']

import ConfigParser
Config = ConfigParser.ConfigParser()
Config.read('../setup.cfg')

def calculate_app_mag(cutoff, t_fov=2000.0, t_exp=1000.0):
    '''
    Description : helper function for <calculate_abs_mag>.
            maybe will be expanded later

    Parameters
    ----------
    cutoff : float
            apparent magnitude limit / imaging depth
    t_fov : float
            time in seconds per FoV
    t_exp : float
            time in seconds per exposure

    Returns
    -------
    float
            returns the apparent magnitude limit given input parameters
    '''
    return np.float64(cutoff + 2.5 * np.log10(t_fov / t_exp))

def calculate_abs_mag(distance=1.0, _f_type=None):
    '''
    Description : used to calculate the intrinsic brightness a star must be in 
            magnitude in order to be seen at a given distance in Mpc

    Parameters
    ----------
    distance : float
            distance in Mpc to the star
    app_mag : float
            apparent magnitude limit of the equipment being used, aka imaging depth.

    Returns
    -------
    float
            absolute magnitude a star needs to be in order to be seen at the given 
                    distance with the given apparent magnitude limit
    '''
    #print('calculating absolute mag limit')
    if distance == None:
        distance = Config.getfloat('Distance', 'd_mpc')
    if _f_type == None:
        _f_type = Config.get('Filter', 'filter_type')
    app_mag = calculate_app_mag(filter_limits[_f_type])
    return np.float64(app_mag - (5.0 * np.log10(distance * 1e5)))


def load_ab_mags(halo_name, f_type=None, f_prfx=None):
    """
    Summary line.

    Extended description of function.

    Parameters
    ----------
    halo_name : str
        'halo02'

     halo_dir :  str
        'path/to/halo_name/'

    Returns
    -------
    np arr
        returns a magnitude numpy array 

    """
    if f_prfx == None:
        f_prfx = Config.get('Filter', 'filter_prfx')
    if f_type == None:
        f_type = Config.get('Filter', 'filter_type')
    _filter_ = f_prfx + f_type
    fh = ospth.join(Config.get('PATH', 'halo_dir') , halo_name, _filter_ + '.npy')
    #print('loading magnitude array: ' + fh)
    mags = np.load(fh)
    return mags.astype(np.float64)

def load_data_arr(halo_name, data_key, lim):
    """
    Summary line.

    Extended description of function.

    Parameters
    ----------
    halo_name : str
        'halo02'

     halo_dir :  str
        'path/to/halo_name/'

    Returns
    -------
    np arr
        returns a magnitude numpy array 

    """
    fh = ospth.join(Config.get('PATH', 'halo_dir') , halo_name, data_key + '.npy')
    #print('loading ' + data_key + ' array: ' + fh)
    mags = np.load(fh, mmap_mode='r')[lim]
    return mags.astype(np.float64)


def load_positions(halo, idx):
    '''
    Load indexed x,y,z position arrays for a gin halo 

    TODO <Extended description of function.>

    Parameters
    ----------
    halo : string
            Name of halo.

    idx : list
            List of desired index.

    Returns
    -------
    np.array,np.array,np.array
            Three numpy arrays representing x,y,z positions for a given halo.
    '''
    # use Config.get('PATH', 'halo_dir') from options to inform data directory
    # and concatenate it with the halo's name to form a PATH.
    path = ospth.join(Config.get('PATH', 'halo_dir') , halo)
    #print('loading position arrays from: ' + path)

    return np.asarray([
        np.load(ospth.join(path, 'px.npy'), mmap_mode='r')[idx],
        np.load(ospth.join(path, 'py.npy'), mmap_mode='r')[idx],
        np.load(ospth.join(path, 'pz.npy'), mmap_mode='r')[idx]],
        dtype=np.float64)

def load_satid(halo, idx):
    path = ospth.join(Config.get('PATH', 'halo_dir') , halo, 'satid.npy')
    print('loading satid from: ', path)
    return np.load(path, mmap_mode='r')[idx]
    
def load_allkeys(halo, idx):
    '''
    Load indexed x,y,z position arrays for a gin halo 

    TODO <Extended description of function.>

    Parameters
    ----------
    halo : string
            Name of halo.

    idx : list
            List of desired index.

    Returns
    -------
    np.array,np.array,np.array
            .
    '''
    # use Config.get('PATH', 'halo_dir') from options to inform data directory
    # and concatenate it with the halo's name to form a PATH.
    path = ospth.join(Config.get('PATH', 'halo_dir') , halo)
    #print('loading data arrays from: ' + path)
    return np.asarray([
        np.load(ospth.join(path, 'teff.npy'), mmap_mode='r')[idx],
        np.load(ospth.join(path, 'age.npy'), mmap_mode='r')[idx],
        np.load(ospth.join(path, 'feh.npy'), mmap_mode='r')[idx],
        np.load(ospth.join(path, 'alpha.npy'), mmap_mode='r')[idx],
        np.load(ospth.join(path, 'lum.npy'), mmap_mode='r')[idx],
        np.load(ospth.join(path, 'smass.npy'), mmap_mode='r')[idx],
        np.load(ospth.join(path, 'mact.npy'), mmap_mode='r')[idx],
        np.load(ospth.join(path, 'mtip.npy'), mmap_mode='r')[idx]],
        dtype=np.float64)


def apparent_magnitude(mag_abs, distance):
    '''
    Convert an array of absolute magnitudes into apparent magnitudes 
            for a given distance.

    TODO Extended description of function.

    Parameters
    ----------
    mag_abs : np.array
            array of absolute magnitudes.

    distance : float
            distance to target stars in units of Mpc.

    Returns
    -------
    np.array
            Returns the apparent magnitude for the given distance.
    '''
    #print('calculating apparent magnitudes')
    # make sure that the distance is greater than zero
    if distance < 0.0:
        distance = 0.1  # 100,000 Kpc
    # calculate the distance modulus.
    distance_modulus = 5. * np.log10(distance * 1e5)
    return mag_abs.astype(np.float64) + distance_modulus


def rotation_matrix(ax, th):
    '''
    Description : used to rotate x,y,z arrays. usually called by <rotate>.

    Parameters
    ----------
    ax : list
            an axis to rotate about.
    th : float
            the amount to rotate by in radians

    Returns
    -------
    list,list,list
            returns a rotation matrix
    '''
    #print('calculating rotation matrix')
    a = np.cos(np.asarray(th) / 2.0)
    b, c, d = -(np.asarray(ax) / (np.dot(ax, ax))**2) * \
        np.sin(np.asarray(th) / 2.0)
    aa, bb, cc, dd, bc = a * a, b * b, c * c, d * d, b * c
    ad, ac, ab, bd, cd = a * d, a * c, a * b, b * d, c * d
    px = [aa + bb - cc - dd, 2.0 * (bc + ad), 2.0 * (bd - ac)]
    py = [2.0 * (bc - ad), aa + cc - bb - dd, 2.0 * (cd + ab)]
    pz = [2.0 * (bd + ac), 2.0 * (cd - ab), aa + dd - bb - cc]
    return px, py, pz

def rotate(xyz, axis, theta):
    '''
    Consolidates data and rotation options into a single function.

    Parameters
    ----------
    xyz : list
            a list of numpy arrays [px,py,pz]
    axis : list
            an axis to rotate about [0,1,0]
    theta : float
            amount to rotate by in radians

    Returns
    -------
    array
            returns the rotated positions in the same form as they were input.
    '''
    #print('rotating x, y and z position arrays')
    return np.asarray(np.dot(rotation_matrix(axis, theta), xyz), dtype=np.float64)

def kpc_box_at_distance(distance, kpc, unit=3600.0):
    '''
    Input distance in Mpc
    Output square degrees if unit==3600.0
    ((206265.0*(_kpc_/1e2)/distance*1e1)/3600)^2
    '''
    return (((206265.0 * kpc) / (distance * 1e3)) / unit)**2

def calculate_nFov(distance, kpc):
    '''
    Input distance in Mpc
    Output number of FOV for a box covering _kpc_ radius at input distance 
    ((square deg of box) / WFIRST square degree FOV)
    '''
    WFIRST_FOV = 0.79 * 0.43
    n_fov = kpc_box_at_distance(distance, kpc) / WFIRST_FOV
    if n_fov < 1.0:
        return int(1)
    else:
        return math.ceil(n_fov)

