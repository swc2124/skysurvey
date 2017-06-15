from __future__ import division, absolute_import, print_function

import os
import sys
import ConfigParser
from time import sleep as _sleep
from collections import OrderedDict

import skysurvey
from skysurvey.new_config import SYS_CFG_FNAME

from matplotlib import pyplot as plt

from astropy.table import Table
from astropy.table import Column

from astropy.units import kiloparsec
from astropy.units import radian
from astropy.units import degree
from astropy.units import gigayear
from astropy.units import ABmag
from astropy.units import mag

import numpy as np


# Table directory names and paths.
sys_config_fh = os.path.join(os.path.dirname(
    os.path.realpath(skysurvey.__file__)), SYS_CFG_FNAME)
SysConfig = ConfigParser.ConfigParser()
SysConfig.read(sys_config_fh)
config_fh = SysConfig.get('skysurvey_global_settings', 'config_fh')
Config = ConfigParser.ConfigParser()
Config.read(config_fh)
table_dir = Config.get('PATH', 'table_dir')
table_dirs_names = [name for name in os.listdir(table_dir)
                    if os.path.isdir(os.path.join(table_dir, name))
                    if not name.startswith('.')]


def new_feature(r0, r1, d0, d1, nbx, points):
    return {
        'r_in': r0,
        'r_out': r1,
        'degree0': d0,
        'degree1': d1,
        'nboxes': nbx,
        'points': set(points)}

def _plot(feat_dict, plot_number, halo):
    if feat_dict['nboxes'] < 10:
        return
    #print('plotting', halo, 'feature number:', plot_number)
    fig = plt.figure(figsize=(10, 10))
    title0 = 'group number:' + str(plot_number)
    title1 = '  r0: ' + str(feat_dict['r_in']) + \
        ' r1:' + str(feat_dict['r_out'])
    title2 = '  phi0: ' + \
        str(feat_dict['degree0']) + ' phi1:' + str(feat_dict['degree1'])
    title3 = '  nboxes: ' + str(feat_dict['nboxes'])
    title = title0 + title1 + title2 + title3
    fig.suptitle(title)
    rect = [0.1, 0.1, 0.8, 0.8]
    plr = fig.add_axes(rect, polar=True, alpha=.2,
                       facecolor=None, frame_on=False)
    plr.set_ylim([0, 300])
    ax = fig.add_axes(rect, alpha=.2, facecolor=None, frame_on=False)
    ax.axis
    ax.set_ylim([0, 500])
    ax.set_xlim([0, 500])
    ax.axis('off')
    ax.set_aspect('equal')
    points = feat_dict['points']
    x = [a for a, b in list(points)]
    y = [b for a, b in list(points)]
    ax.scatter(x, y)
    plot_fh = os.path.join(Config.get(
        'PATH', 'plot_dir'), 'groupfinder_testplots', halo + '_feature_' + str(plot_number))
    try:
        fig.savefig(plot_fh)
    except:
        import time
        time.sleep(1)
        fig.savefig(plot_fh)
    plt.close()


def _plot_halo(table):
    print('plotting', table.meta['halo'])
    fig = plt.figure(figsize=(10, 10))
    halo = table.meta['halo']
    fig.suptitle(halo)
    rect = [0.1, 0.1, 0.8, 0.8]
    plr = fig.add_axes(rect, polar=True, alpha=.2,
                       facecolor=None, frame_on=False)
    plr.set_ylim([0, 300])
    ax = fig.add_axes(rect, alpha=.2, facecolor=None, frame_on=False)
    ax.axis
    ax.set_ylim([0, 500])
    ax.set_xlim([0, 500])
    ax.axes.grid()
    ax.set_aspect('equal')
    ax.axis('off')
    ax.scatter(table['x_int'], table['y_int'])
    plot_fh = os.path.join(Config.get('PATH', 'plot_dir'),
                           'groupfinder_testplots', halo)
    try:
        fig.savefig(plot_fh)
    except:
        import time
        time.sleep(1)
        fig.savefig(plot_fh)
    plt.close()


def add_to_feature(fdict, r1, deg0, deg1, pts):
    fdict['r_out'] = r1

    if deg0 < fdict['degree0']:
        fdict['degree0'] = round(deg0, 2)

    if deg1 > fdict['degree1']:
        fdict['degree1'] = round(deg1, 2)

    fdict['points'].update(pts)
    fdict['nboxes'] = len(fdict['points'])
    return fdict

def merge_heavy(f_dict, halo):
    del_keys = []
    n_attempts = 0
    for f_id in f_dict.keys():
        feature = f_dict[f_id]
        if not 'points' in feature.keys():
            continue
        if not len(feature['points']):
            continue
        r_in = feature['r_in']
        r_out = feature['r_out']
        deg0 = feature['degree0']
        deg1 = feature['degree1']
        for _f_id in f_dict.keys():

            if _f_id == f_id:
                continue
            _feature = f_dict[_f_id]
            if not 'points' in _feature.keys():
                continue
            if not len(_feature['points']):
                continue
            _r_in = _feature['r_in']
            _r_out = _feature['r_out']
            _deg0 = _feature['degree0']
            _deg1 = _feature['degree1']
            center_deg = _deg0 + ((_deg1 - _deg0) / 2.0)
            merge = False
            if (deg0 - 20.) <= _deg0 <= (deg1 + 10.) or (deg0 - 10.) <= _deg1 <= (deg1 + 20.0):
                for x, y  in _feature['points']:
                    if (r_in - 15.0 <= np.sqrt(np.square(x) + np.square(y)) <= r_out + 15.0):
                        merge = True
                        break
            if merge:
                feature['points'].update(_feature['points'])
                if feature['r_in'] > _feature['r_in']:
                    feature['r_in'] = _feature['r_in']
                if feature['r_out'] < _feature['r_out']:
                    feature['r_out'] = _feature['r_out']
                if feature['degree0'] > _feature['degree0']:
                    feature['degree0'] = _feature['degree0']
                if feature['degree1'] < _feature['degree1']:
                    feature['degree1'] = _feature['degree1']
                _feature['points'] = set()

                plot_fh = os.path.join(Config.get(
                    'PATH', 'plot_dir'), 'groupfinder_testplots', halo + '_feature_' + str(_f_id) + '.png')
                if os.path.isfile(plot_fh):
                    os.remove(plot_fh)
                sys.stdout.write('\n')
                sys.stdout.flush()
                print('Heavy Merged feature', _f_id,  'into feature ', f_id)
                _plot(feature, f_id, halo)
                del_keys.append(_f_id)

    for key in del_keys:
        del f_dict[key]
    return f_dict


def merge(f_dict, halo):
    keys = f_dict.keys()
    del_keys = []
    for f_id in keys:
        feature = f_dict[f_id]
        if not len(feature['points']):
            continue
        for i, _f_id in enumerate(keys):
            if f_id == _f_id:
                continue
            merge = False
            _feature = f_dict[_f_id]
            if not len(_feature['points']):
                continue
            for point in _feature['points']:
                if point in feature['points']:
                    merge = True
                    break
            if merge:
                feature['points'].update(_feature['points'])
                if feature['r_in'] > _feature['r_in']:
                    feature['r_in'] = _feature['r_in']
                if feature['r_out'] < _feature['r_out']:
                    feature['r_out'] = _feature['r_out']
                if feature['degree0'] > _feature['degree0']:
                    feature['degree0'] = _feature['degree0']
                if feature['degree1'] < _feature['degree1']:
                    feature['degree1'] = _feature['degree1']
                _feature['points'] = set()
                plot_fh = os.path.join(Config.get(
                    'PATH', 'plot_dir'), 'groupfinder_testplots', halo + '_feature_' + str(_f_id) + '.png')
                if os.path.isfile(plot_fh):
                    os.remove(plot_fh)
                sys.stdout.write('\n')
                sys.stdout.flush()
                print('Merged feature', _f_id,  'into feature ', f_id)
                _plot(feature, f_id, halo)
                del_keys.append(_f_id)

    for key in del_keys:
        del f_dict[key]
    return f_dict

def resort_dict(f_dict):
    l0 = {}
    for key in f_dict.keys():
        l0[key] = f_dict[key]['nboxes']
    l0 = sorted(l0.items(), key=lambda x: x[1])
    l0.reverse()
    new_dict = {}
    for i, value in enumerate(l0):
        new_dict[i] = f_dict[value[0]]
    return new_dict

table_dirs = {}
for name in table_dirs_names:
    table_dirs[name] = os.path.join(table_dir, name)

percent = 0.1  # 10% of radius
annulus_degree_step = 1  # Within 1 degree
xbox_min_value = 10.0
filenames = os.listdir(table_dirs['merged_tables'])
filenames.reverse()

for name in filenames:
    table_fh = os.path.join(table_dirs['merged_tables'], name)
    table = Table.read(table_fh, path='data')
    table.remove_rows(np.nonzero(table['Xbox'] < xbox_min_value))
    table.add_column(Column(data=np.rad2deg(
        table['Phi']), name='Degree', description='degrees', dtype=np.float16))
    halo = table.meta['halo']
    _plot_halo(table)
    master_dict = {}
    feature_id = 0
    r_start = int(table['Rads'].min())
    r_stop =  int(table['Rads'].max())
    for i, annulus in enumerate(range(r_start, r_stop, 1)):
        region = annulus * 0.1
        r_in = annulus - region
        r_out = annulus + region
        deg_start = table['Degree'].min()
        deg_end = table['Degree'].max()
        for degree in np.linspace(deg_start, deg_end, 360):
            lims = np.nonzero(
                np.logical_and(
                    np.logical_and(
                        table['Degree'] >= degree - annulus_degree_step,
                        table['Degree'] < degree + annulus_degree_step),
                    np.logical_and(
                        table['Rads'] >= r_in,
                        table['Rads'] < r_out)))
            n_boxes = len(lims[0])
            if n_boxes:
                points = zip((table['x_int'][lims]).tolist(),
                             (table['y_int'][lims]).tolist())
                if not len(master_dict.keys()):
                    master_dict[feature_id] = new_feature(
                        r_in,
                        r_out,
                        round(degree, 2),
                        round(degree + annulus_degree_step, 2),
                        n_boxes,
                        points)
                    continue
                known_feature = False
                for known_feat_id in master_dict.keys():
                    feature = master_dict[known_feat_id]
                    new_points = set(points)
                    for point in new_points:
                        if point in feature['points']:
                            known_feature = True
                            break
                    if known_feature:
                        break
                if known_feature:
                    master_dict[known_feat_id] = add_to_feature(
                        master_dict[known_feat_id],
                        r_in,
                        degree,
                        degree + annulus_degree_step,
                        points)
                else:
                    feature_id += 1
                    master_dict[feature_id] = new_feature(
                        r_in,
                        r_out,
                        round(degree, 2),
                        round(degree + annulus_degree_step, 2),
                        n_boxes,
                        points)

            master_dict = merge(master_dict, halo)
            msg0 = '\r\r[ ' + halo + ' ] '
            msg1 = '[ ' + str(annulus) + ' Kpc ] '
            msg2 = '[ ' + str(feature_id) + ' ] '
            msg3 = '[ ' + str(n_boxes) + ' ] '
            msg4 = '[ ' + str(degree) + ' ] '
            msg = msg0 + msg1 + msg2 + msg3 + msg4
            sys.stdout.write(msg)
            sys.stdout.flush()

    if len(master_dict.keys()) > 20:
        master_dict = merge_heavy(master_dict, halo)
    
    master_dict = resort_dict(master_dict)
    for i in range(5):
        if i in master_dict.keys():
            _plot(master_dict[i], i, halo)

    logfile_fh = os.path.join(
        Config.get('PATH', 'plot_dir'),
        'groupfinder_testplots',
        halo + '_logfile')

    with open(logfile_fh, 'w') as logfile:
        logfile.write('\n' + '-' * 50 + '\n')
        for key in table.meta.keys():
            logfile.write('[table]' + key + ': ' + str(table.meta[key]) + '\n')
        logfile.write('\n' + '-' * 50 + '\n')
        logfile.write('-' * 50 + '\n')
        logfile.write('[Features]' + str(r_start) + '\n')
        logfile.write('[r_start]' + str(r_stop) +  ' Kpc\n')
        logfile.write('[r_stop]' + str() +  ' Kpc\n')
        logfile.write('[deg_start]' + str(deg_start) +  ' deg\n')
        logfile.write('[deg_stop]' + str(deg_end) +  ' deg\n')
        for f_id in master_dict.keys():
            feature = master_dict[f_id]
            logfile.write('\n' + '-' * 50 + '\n')
            logfile.write('[feature number: ' + str(f_id) + ' ]\n')
            for key in feature.keys():
                if key == 'points':
                    continue
                logfile.write('  -->' + key + ' : ' +
                              str(feature[key]) + '\n')
