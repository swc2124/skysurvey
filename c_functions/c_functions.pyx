from __future__ import division, print_function

import sys
import os
from os.path import join

import cython
from cython.parallel import parallel
from cython.parallel import prange

import numpy as np
cimport numpy as np

import ebf

from libc.stdlib cimport rand
from libc.math cimport cos
from libc.math cimport sin
from libc.math cimport M_PI

import ConfigParser
from skysurvey.new_config import SYS_CFG_FNAME
import skysurvey

sys_config_fh = os.path.join(os.path.dirname(os.path.realpath(skysurvey.__file__)), SYS_CFG_FNAME)
SysConfig = ConfigParser.ConfigParser()
SysConfig.read(sys_config_fh)
config_fh = SysConfig.get('skysurvey_global_settings', 'config_fh')
Config = ConfigParser.ConfigParser()
Config.read(config_fh)

ebf_fh = os.path.join(Config.get('PATH', 'data_dir'), 'satprop.ebf')
SATPROP = ebf.read(ebf_fh)
#np.import_array()
if 'NUMBER_OF_PROCESSORS' in os.environ.keys():
    NUM_PROCESSORS = int(os.environ['NUMBER_OF_PROCESSORS'])
else:
    NUM_PROCESSORS = 1

@cython.boundscheck(False)
@cython.wraparound(False)
def rotation_matrix(np.ndarray[np.float64_t, ndim=1, mode='c'] ax, np.float64_t th):
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
    cdef np.float64_t a, b, c, d, aa, bb, cc, dd
    a = np.cos(th / 2.0)
    b, c, d = -(ax / (np.dot(ax, ax))**2) * np.sin(th / 2.0)
    aa, bb, cc, dd, bc = a * a, b * b, c * c, d * d, b * c
    ad, ac, ab, bd, cd = a * d, a * c, a * b, b * d, c * d
    return [aa + bb - cc - dd, 2.0 * (bc + ad), 2.0 * (bd - ac)], [2.0 * (bc - ad), aa + cc - bb - dd, 2.0 * (cd + ab)], [2.0 * (bd + ac), 2.0 * (cd - ab), aa + dd - bb - cc]


@cython.boundscheck(False)
@cython.wraparound(False)
def rotate(
    np.ndarray[np.float64_t, ndim=2, mode='c'] xyz,
    np.ndarray[np.float64_t, ndim=1, mode='c'] axis,
    np.float64_t theta):
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
    return np.asarray(
        np.dot(
            rotation_matrix(axis, theta),
            [xyz[0], xyz[1], xyz[2]]),
        dtype=np.float64)


@cython.boundscheck(False)
@cython.wraparound(False)
def trippel_rotate(
    np.ndarray[np.float64_t, ndim=2, mode='c'] xyz):
    '''
    '''
    #print('starting trippel rotatation')
    cdef np.float64_t theta_1 = np.float64(
        np.divide(
            (np.random.randint(360) * np.pi),
            180.0))
    cdef np.float64_t theta_2 = np.float64(
        np.divide(
            (np.random.randint(360) * np.pi),
            180.0))
    cdef np.float64_t theta_3 = np.float64(
        np.divide(
            (np.random.randint(360) * np.pi),
            180.0))
    return rotate(
        rotate(
            rotate(xyz,
                   np.asarray([0, 0, 1], dtype=np.float64), theta_1),
            np.asarray([0, 1, 0], dtype=np.float64), theta_2),
        np.asarray([1, 0, 0], dtype=np.float64), theta_3)

@cython.boundscheck(False)
@cython.wraparound(False)
def integerize(
    np.ndarray[np.float64_t, ndim=1] x,
    np.ndarray[np.float64_t, ndim=1] y):
    '''
    Convert np.array of float64 to int32 for binning into a grid.

    TODO Extended description of function.

    Parameters
    ----------
    x : np.array
            np.array of x coordinates.

    y : np.array
            np.array of x coordinates.

    center : float
            float representing the midpoint of the grid.

    scale : float
            float representing the scale factor for matching Kpc to the grid.
    Returns
    -------
    np.array,np.array
            Returns two equivalent length arrays in int32 form.
    '''
    center = Config.getint('grid_options', 'size') / 2.0
    scale = Config.getint('grid_options', 'size') / (x.max() + np.abs(x.min()))
    line = '-' * 85
    #print('converting px and py arrays to integers')

    print(line)
    print('[ integerize ]', '\n')

    print('[before ] px min, mean, max : ', x.min(), ',', x.mean(), ',', x.max())
    cdef np.ndarray[np.float64_t, ndim = 1, mode='c'] x1 = x.round(3)
    x1 *= scale
    x1 += center
    x2 = x1.astype(np.int32)
    print('[ after ] px min, mean, max : ', x2.min(), ',', x2.mean(), ',', x2.max())

    print('[before ] py min, mean, max : ', y.min(), ',', y.mean(), ',', y.max())
    cdef np.ndarray[np.float64_t, ndim = 1, mode='c'] y1 = y.round(3)
    y1 *= scale
    y1 += center
    y2 = y1.astype(np.int32)
    print('[ after ] py min, mean, max : ', y2.min(), ',', y2.mean(), ',', y2.max())
    print(line)

    return x2, y2

@cython.boundscheck(False)
@cython.wraparound(False)
def bin(
    np.ndarray[np.int32_t, ndim=1] px,
    np.ndarray[np.int32_t, ndim=1] py,
    np.ndarray[np.float64_t, ndim=1] ab_mags,
    np.ndarray[np.float64_t, ndim=1] ap_mags,
    np.ndarray[np.float64_t, ndim=1] r_proj,
    np.ndarray[np.float64_t, ndim=1] lims,
    np.ndarray[np.int32_t, ndim=1] satid):

    line = '-' * 85
    print(line)
    print('\n', '[ bin ]', '\n')
    cdef:

        # Itierable
        size_t i, ii
        #np.float64_t _one = 1.0

        # Magnitude limits (mlim_min, mlim_med, mlim_max).
        np.float64_t mlim_min = lims[0]
        np.float64_t mlim_med = lims[1]
        np.float64_t mlim_max = lims[2]

        # Apparent magnitude (apparent_mag).
        np.float64_t apparent_mag

        # The grid array (grid)

        np.ndarray[np.float64_t, ndim = 3, mode='c'] grid = np.zeros((Config.getint('grid_options', 'size'), Config.getint('grid_options', 'size'), Config.getint('grid_options', 'n_slices')), dtype=np.float64)

        # Satprop arrays. TODO package data file
        np.ndarray[np.float32_t, ndim=1] tsat = SATPROP['tsat']
        np.ndarray[np.int32_t, ndim=1] bsat = SATPROP['bsat']
        np.ndarray[np.int32_t, ndim=1] nsat = SATPROP['nsat']
        np.ndarray[np.int32_t, ndim=1] nsatc = SATPROP['nsatc']

        # Index arrays
        np.ndarray[np.int64_t, ndim=1] idx_1, idx_2
        #np.ndarray[np.int64_t, ndim=1] idx_1 = np.nonzero(bsat[satid + nsatc[2]] == 0)[0]
        #np.ndarray[np.int64_t, ndim=1] idx_2 = np.nonzero(bsat[satid + nsatc[2]] == 0)[1]

        # Number of stars (n_stars) to bin.
        np.int_t n_stars = len(px)

        # Number of threads
        np.int_t n_threads = int(px.shape[0]/(5 * 1e5))

        # X and Y grid boundaries (boundary_x, boundary_y)
        np.int_t boundary_y = grid.shape[0]
        np.int_t boundary_x = grid.shape[1]

        # Satellite number of star.
        np.int32_t sat_number

        # Age of sat (Gyr).
        np.float32_t sat_age

        # Bound or unbound status.
        np.int32_t sat_bound

        # Number of (Bound) or (unbound) stars.
        np.int32_t bound = 0
        np.int32_t unbound = 0

        # Number of stars outside (missed) of the superimposed grid.
        np.int32_t missed = 0

        np.int32_t bad = 0


    # Set (threads)
    if n_threads >= NUM_PROCESSORS:
        n_threads = NUM_PROCESSORS - 2
        if n_threads <= 0:
            n_threads = 1

    # Initial print statement.
    print('sample px and py before starting bin loop:')
    print('  -> ', np.random.choice(px, 10))
    print('  -> ', np.random.choice(py, 10))
    print('\n')
    print('threads : ', n_threads)
    print('mlim_min: ', mlim_min)
    print('mlim_med: ', mlim_med)
    print('mlim_max: ', mlim_max)

    with nogil, parallel(num_threads=n_threads):

        for i in prange(n_stars, schedule='dynamic'):

            # Check to see if star is in FOV of grid.
            if (px[i] >= boundary_x or px[i] <= 0 or py[i] >= boundary_y or py[i] <= 0):
                missed += 1

            # If star is in the grid.
            else:

                sat_number = satid[i]

                sat_age = tsat[sat_number]

                sat_bound = bsat[sat_number]

                apparent_mag = ap_mags[i]

                # If the star is unbound [0].
                if not sat_bound:
                    unbound += 1

                    # Bin stars.
                    grid[px[i], py[i], 0] += 1.0

                    # Magnitudes.
                    grid[px[i], py[i], 1] += ab_mags[i]

                    #grid[px[i], py[i], 2] += apparent_mag
                    #if apparent_mag < mlim_min:
                    #    grid[px[i], py[i], 3] += 1.0
                    #if apparent_mag < mlim_med:
                    #    grid[px[i], py[i], 4] += 1.0
                    #if apparent_mag < mlim_max:
                    #    grid[px[i], py[i], 5] += 1.0

                    # Accretion time (Gyr) of satellite.
                    grid[px[i], py[i], 2] += sat_age
                    # Satid.
                    grid[px[i], py[i], 3] += sat_number



                # If the star is bound [1].
                else:
                    bound += 1
                    # Bin stars.
                    #grid[px[i], py[i], 8] += 1.0
                    # Magnitudes.
                    #grid[px[i], py[i], 9] += ab_mags[i]
                    #grid[px[i], py[i], 10] += apparent_mag
                    # Accretion time (Gyr).
                    #grid[px[i], py[i], 11] += sat_age
                    # Satid.
                    #grid[px[i], py[i], 12] += sat_number
                    if not sat_bound:
                        bad += 1

                #grid[px[i], py[i], 13] += 1
                grid[px[i], py[i], 4] += r_proj[i]

    # Slices that need to be divided by the number of stars in each bin.
    idx_1, idx_2 = np.nonzero(grid[:, :, 0] > 0.)
    grid[idx_1, idx_2, 1] /= grid[idx_1, idx_2, 0]
    grid[idx_1, idx_2, 2] /= grid[idx_1, idx_2, 0]
    grid[idx_1, idx_2, 3] /= grid[idx_1, idx_2, 0]
    grid[idx_1, idx_2, 4] /= grid[idx_1, idx_2, 0]
    #grid[idx_1, idx_2, 14] /= grid[idx_1, idx_2, 0]

    #idx_1, idx_2 = np.nonzero(grid[:, :, 8] > 0.)
    #grid[idx_1, idx_2, 9] /= grid[idx_1, idx_2, 8]
    #grid[idx_1, idx_2, 10] /= grid[idx_1, idx_2, 8]
    #grid[idx_1, idx_2, 11] /= grid[idx_1, idx_2, 8]
    #grid[idx_1, idx_2, 12] /= grid[idx_1, idx_2, 8]

    #idx_1, idx_2 = np.nonzero(grid[:, :, 13] > 0.)
    #grid[idx_1, idx_2, 14] /= grid[idx_1, idx_2, 13]

    print('\n-------------------------------')
    print('    bound stars: ', bound)
    print('  unbound stars: ', unbound)
    print('percent unbound: ', round(1e2 * (float(unbound)/float(len(px))), 2), '%')
    print('   binned stars: ', (n_stars - missed))
    print('    total stars: ', n_stars)
    print('   missed stars: ', missed)
    print('-------------------------------\n')
    return grid

@cython.boundscheck(False)
@cython.wraparound(False)
def find_dlims(
    np.ndarray[np.float64_t, ndim = 1] mag_arr,
    double ab_mag_limit):
    """
    Find the indices of the stars in an array that are visible given a limit.


    Extended description of function.

    Parameters
    ----------
    mag_arr : np arr
        Numpy array of absolute magnitude.

    ab_mag_limit : float
        The minimum intrinsic brightness needed to be seen.

    Returns
    -------
    np arr
        Returns a Numpy array containing the indices of visible stars

    """
    line = '-' * 85
    print(line)
    print('\n[ find_dlims ]\nfinding indices of visible stars : ', ab_mag_limit, ' abs mag lim')
    print('mean of mag_arr:', mag_arr.mean())
    return np.nonzero(mag_arr < ab_mag_limit)[0]

@cython.boundscheck(False)
@cython.wraparound(False)
def box_lims(
    np.ndarray[np.float64_t, ndim = 1] px,
    np.ndarray[np.float64_t, ndim = 1] py,
    np.float64_t box_size,  np.float64_t box_step):
    cdef:
        size_t i
        np.int32_t n_stars = px.shape[0]
        np.float64_t outter_box = (box_size + box_step)
        np.float64_t neg_outter_box = -(outter_box)
        np.float64_t neg_box_size = -(box_size)
        np.ndarray[np.int64_t, ndim = 1, mode='c'] idx_arr = np.zeros((n_stars), dtype=np.int64)
        np.int_t n_threads = np.int(px.shape[0]/1e6)

    # Set (threads)
    if n_threads >= NUM_PROCESSORS:
        n_threads = NUM_PROCESSORS - 2
        if n_threads <= 0:
            n_threads = 1

        print('box_lims(): number of processors=' + str(NUM_PROCESSORS) + ' - number of threads=' + str(n_threads))
        with nogil, parallel(num_threads=n_threads):
            for i in prange(n_stars, schedule='dynamic'):
                if (px[i] > neg_outter_box and px[i] < outter_box):
                    if (py[i] >= box_size and py[i] < outter_box):
                        idx_arr[i] = 1
                    elif (py[i] <= neg_box_size and py[i] > neg_outter_box):
                        idx_arr[i] = 1
                if (py[i] > neg_outter_box and py[i] < outter_box):
                    if (px[i] >= box_size and px[i] < outter_box):
                        idx_arr[i] = 1
                    elif (px[i] <= neg_box_size and px[i] > neg_outter_box):
                        idx_arr[i] = 1
        return idx_arr

cdef extern from "math.h":
    double M_PI
@cython.boundscheck(False)
@cython.wraparound(False)
def rotate_positions(
    np.ndarray[np.float64_t, ndim = 2] xyz_arr,
    np.ndarray[np.float64_t, ndim = 2] out_arr):
    cdef:
        int _theta1 = rand() % 4000
        int _theta2 = rand() % 2000
        int _theta3 = rand() % 6000
        double theta1 = (_theta1 / 1e3) * M_PI # 0.01745
        double theta2 = (_theta2 / 1e3) * M_PI # 0.01745
        double theta3 = (_theta3 / 1e3) * M_PI # 0.01745
        size_t i, j, k
        double rot_matrix[3][3]
        int I = len(rot_matrix)
        int J = xyz_arr.shape[1]
        int K = xyz_arr.shape[0]
        np.float64_t s, a, aa, bb, cc, dd, bc, ad, ac, ab, bd, cd, b, c, d
    for a, b, c, d in [(cos(theta1), 0, 0, -sin(theta1)), (cos(theta2), 0, -sin(theta2), 0), (cos(theta3), -sin(theta3), 0, 0)]:
        aa, bb, cc, dd, bc = a * a, b * b, c * c, d * d, b * c
        ad, ac, ab, bd, cd = a * d, a * c, a * b, b * d, c * d
        rot_matrix[0][:] = [aa + bb - cc - dd, 2.0 * (bc + ad), 2.0 * (bd - ac)]
        rot_matrix[1][:] = [2.0 * (bc - ad), aa + cc - bb - dd, 2.0 * (cd + ab)]
        rot_matrix[2][:] = [2.0 * (bd + ac), 2.0 * (cd - ab), aa + dd - bb - cc]
        for i in range(I):
            for j in xrange(J):
                s = 0
                for k in range(K):
                    s += (xyz_arr[k, j] * rot_matrix[i][k])
                out_arr[i, j] = s

@cython.boundscheck(False)
@cython.wraparound(False)
def deposit_data(np.ndarray[np.float64_t, ndim = 2] in_arr,
    np.ndarray[np.float64_t, ndim = 1] out_arr,
    np.ndarray[np.int64_t, ndim = 1] bidx):
    cdef:
        size_t j, idx
        size_t n_data = 3
        np.float64_t n_stars = len(bidx)
    if n_stars <= 0:
        return 0
    out_arr[26] += n_stars
    for j in range(in_arr.shape[0]):
        idx = j * n_data
        out_arr[idx] += in_arr[j][bidx].min()
        out_arr[idx+1] += in_arr[j][bidx].mean()
        out_arr[idx+2] += in_arr[j][bidx].max()
    return n_stars

