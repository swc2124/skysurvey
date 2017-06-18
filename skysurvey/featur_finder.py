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

def _plot(feat_dict):
    if feat_dict['nboxes'] < 15:
        if os.path.isfile(feat_dict['plot_fh'] + '.png'):
            os.remove(feat_dict['plot_fh'] + '.png')
        return feat_dict
    if 'last_nboxes' in feat_dict.keys():
        if feat_dict['nboxes'] == feat_dict['last_nboxes']:
            return feat_dict
    else:
        feat_dict['last_nboxes'] = feat_dict['nboxes']
    plot_number = feat_dict['feature_id']
    sys.stdout.write('\nplotting ' + halo +
                     ' feature number: ' + str(plot_number))
    sys.stdout.flush()
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
    ax.set_ylim([0, 600])
    ax.set_xlim([0, 600])
    ax.axis('off')
    ax.set_aspect('equal')
    points = feat_dict['points']
    try:
        x = [a for a, b in list(points)]
        y = [b for a, b in list(points)]
        ax.scatter(x, y)
        plot_fh = feat_dict['plot_fh']
        try:
            fig.savefig(plot_fh)
        except:
            pass
        plt.close()
        return feat_dict
    except TypeError:
        plt.close()
        return feat_dict

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
    ax.set_ylim([0, 600])
    ax.set_xlim([0, 600])
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

def merge_heavy(f_dict):
    for f_id in f_dict.keys():
        feature = f_dict[f_id]
        if not 'points' in feature.keys():
            continue
        if not len(feature['points']):
            continue
        r_in = feature['r_in']
        r_out = feature['r_out']
        r_center = r_in + ((r_out - r_in) / 2.0)
        deg0 = feature['degree0'] - 0.5
        deg1 = feature['degree1'] + 0.5
        for _f_id in f_dict.keys():
            if _f_id <= f_id:
                continue
            _feature = f_dict[_f_id]
            if not 'points' in _feature.keys():
                continue
            if not len(_feature['points']):
                continue
            _r_in = _feature['r_in']
            _r_out = _feature['r_out']
            _r_center = _r_in + ((_r_out - _r_in) / 2.0)
            _deg0 = _feature['degree0']
            _deg1 = _feature['degree1']
            merge = False
            if (deg0 <= _deg0 <= deg1 or deg0 <= _deg1 <= deg1):
                if (r_in + 1.0 < _r_center < r_out - 1.0):
                    merge = True
            if merge:
                f_dict[f_id]['points'].update(_feature['points'])
                if f_dict[f_id]['r_in'] > _feature['r_in']:
                    f_dict[f_id]['r_in'] = _feature['r_in']
                if f_dict[f_id]['r_out'] < _feature['r_out']:
                    f_dict[f_id]['r_out'] = _feature['r_out']
                if f_dict[f_id]['degree0'] + 180.0 > _feature['degree0'] + 180.0:
                    f_dict[f_id]['degree0'] = _feature['degree0']
                if f_dict[f_id]['degree1'] + 180.0 < _feature['degree1'] + 180.0:
                    f_dict[f_id]['degree1'] = _feature['degree1']
                if f_dict[f_id]['max_xboxvalue'] < _feature['max_xboxvalue']:
                    f_dict[f_id]['max_xboxvalue'] = _feature['max_xboxvalue']
                if f_dict[f_id]['min_xboxvalue'] > _feature['min_xboxvalue']:
                    f_dict[f_id]['min_xboxvalue'] = _feature['min_xboxvalue']
                for satid in f_dict[f_id]['sats_book'].keys():
                    f_dict[f_id]['sats_book'][
                        satid] += f_dict[_f_id]['sats_book'][satid]
                f_dict[f_id]['nstars_total'] += _feature['nstars_total']
                f_dict[_f_id]['points'] = set()
                f_dict[_f_id]['nboxes'] = 0
                f_dict[_f_id]['nstars_total'] = 0
                sys.stdout.write(
                    '\n\n [<HEAVY> Merged feature ' + str(_f_id) + ' into feature ' + str(f_id) + '] \n')
                sys.stdout.flush()
        #f_dict[f_id] = _plot(f_dict[f_id])
    return f_dict

def merge(f_dict):
    keys = f_dict.keys()
    for f_id in keys:
        feature = f_dict[f_id]
        if not len(feature['points']):
            continue
        for i, _f_id in enumerate(keys):
            if _f_id == f_id:
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
                f_dict[f_id]['points'].update(_feature['points'])
                if f_dict[f_id]['r_in'] > _feature['r_in']:
                    f_dict[f_id]['r_in'] = _feature['r_in']
                if f_dict[f_id]['r_out'] < _feature['r_out']:
                    f_dict[f_id]['r_out'] = _feature['r_out']
                if f_dict[f_id]['degree0'] + 180.0 > _feature['degree0'] + 180.0:
                    f_dict[f_id]['degree0'] = _feature['degree0']
                if f_dict[f_id]['degree1'] + 180.0 < _feature['degree1'] + 180.0:
                    f_dict[f_id]['degree1'] = _feature['degree1']
                if f_dict[f_id]['max_xboxvalue'] < _feature['max_xboxvalue']:
                    f_dict[f_id]['max_xboxvalue'] = _feature['max_xboxvalue']
                if f_dict[f_id]['min_xboxvalue'] > _feature['min_xboxvalue']:
                    f_dict[f_id]['min_xboxvalue'] = _feature['min_xboxvalue']
                for satid in f_dict[f_id]['sats_book'].keys():
                    f_dict[f_id]['sats_book'][
                        satid] += f_dict[_f_id]['sats_book'][satid]
                f_dict[f_id]['nstars_total'] += _feature['nstars_total']
                f_dict[_f_id]['points'] = set()
                f_dict[_f_id]['nboxes'] = 0
                f_dict[_f_id]['nstars_total'] = 0
                sys.stdout.write(
                    '\n [Merged feature ' + str(_f_id) + ' into feature ' + str(f_id) + ']  \n')
                sys.stdout.flush()
        #f_dict[f_id] = _plot(f_dict[f_id])
    return f_dict

def new_feature(r0, r1, d0, d1, halo, plt_num, points, xboxes, sat_book):
    sys.stdout.write('\nstarting feature number: ' + str(plt_num))
    feature = {
        'r_in': r0,
        'r_out': r1,
        'degree0': d0,
        'degree1': d1,
        'nboxes': len(set(points)),
        'points': set(points),
        'nstars_total': len(points),
        'max_xboxvalue': max(xboxes),
        'min_xboxvalue': min(xboxes),
        'plot_fh': os.path.join(Config.get('PATH', 'plot_dir'), 'groupfinder_testplots', halo + '_feature_' + str(plt_num)),
        'halo': halo,
        'feature_id': plt_num,
        'sats_book': sat_book}
    return feature

def add_to_feature(fdict, r1, deg0, deg1, pts, xbx, sats):
    superset = False
    if fdict['points'].issuperset(pts):
        superset = True
    if fdict['r_out'] < r1:
        fdict['r_out'] = r1
    if deg0 + 180.0 < fdict['degree0'] + 180.0:
        fdict['degree0'] = round(deg0, 2)
    if deg1 + 180.0 > fdict['degree1'] + 180.0:
        fdict['degree1'] = round(deg1, 2)
    if fdict['max_xboxvalue'] < max(xbx):
        fdict['max_xboxvalue'] = max(xbx)
    if fdict['min_xboxvalue'] > min(xbx):
        fdict['min_xboxvalue'] = min(xbx)
    added = 0
    new_points = set(pts)
    unknown_points = new_points.difference(fdict['points'])
    for i, point in enumerate(pts):
        if point in unknown_points:
            fdict['nstars_total'] += 1
            added += 1
            fdict['sats_book'][sats[i]] += 1
    fdict['points'].update(pts)
    fdict['nboxes'] = len(fdict['points'])
    if not superset:
        sys.stdout.write('added ' + str(added) + ' out of ' + str(len(pts)) +
                         ' stars to feature ' + str(fdict['feature_id']) + '\n')
    if superset and added:
        sys.exit(1)
    sys.stdout.flush()
    return fdict

def count_satids(sats, sats_book):
    sats = sats.tolist()
    for satid in sats_book.keys():
        sats_book[satid] += sats.count(satid)
    return sats_book

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
    table = table[::10]
    table.remove_rows(np.nonzero(table['Xbox'] < xbox_min_value))
    halo = table.meta['halo']
    _plot_halo(table)
    master_dict = {}
    feature_id = 0
    r_start = int(table['Rads'].min())
    r_stop = int(table['Rads'].max())
    cycle = 0
    for i, annulus in enumerate(range(r_start, r_stop, 1)):
        cycle += 1
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
            points = zip((table['x_int'][lims]).tolist(),
                         (table['y_int'][lims]).tolist())
            satids = table['satids'][lims]
            xbox_values = table['Xbox'][lims]
            m5 = ' [nothing happend] '
            if n_boxes:
                if not feature_id in master_dict.keys():
                    sats_book = {}
                    for satid in table.meta['satids']:
                        sats_book[satid] = 0
                    feature_id += 1
                    # new_feature(r0, r1, d0, d1, halo, plt_num, points=[0],
                    master_dict[feature_id] = new_feature(
                        r_in,
                        r_out,
                        round(degree, 2),
                        round(degree + annulus_degree_step, 2),
                        halo,
                        feature_id,
                        points,
                        xbox_values,
                        count_satids(satids, sats_book))
                    m5 = '[started feature] '
                current_feature = False
                known_feature = False
                for known_feat_id in master_dict.keys():
                    feature = master_dict[known_feat_id]
                    for point in set(points):
                        if point in feature['points']:
                            m5 = '[known feature] '
                            known_feature = True
                            break
                    if known_feature:
                        break
                if not known_feature:
                    for unknown_point in set(points):
                        x0, y0 = unknown_point
                        known_points = master_dict[feature_id]['points']
                        for known_point in known_points:
                            x1, y1 = known_point
                            distance = np.sqrt(
                                np.square(x1 - x0) + np.square(y1 - y0))
                            if distance < 1.5:
                                m5 = '[current feature] '
                                current_feature = True
                                break
                        if current_feature:
                            break
                if current_feature:
                    master_dict[feature_id] = add_to_feature(
                        master_dict[feature_id],
                        r_out,
                        degree,
                        degree + annulus_degree_step,
                        points,
                        xbox_values,
                        satids)
                elif known_feature:
                    # add_to_feature(fdict, r1, deg0, deg1, pts, xbx)
                    master_dict[known_feat_id] = add_to_feature(
                        master_dict[known_feat_id],
                        r_out,
                        degree,
                        degree + annulus_degree_step,
                        points,
                        xbox_values,
                        satids)
                else:
                    feature_id += 1
                    # new_feature(r0, r1, d0, d1, halo, plt_num, points=[0],
                    # xboxes=np.array([0]))
                    sats_book = {}
                    for satid in table.meta['satids']:
                        sats_book[satid] = 0
                    master_dict[feature_id] = new_feature(
                        r_in,
                        r_out,
                        round(degree, 2),
                        round(degree + annulus_degree_step, 2),
                        halo,
                        feature_id,
                        points,
                        xbox_values,
                        count_satids(satids, sats_book))
                    m5 = '[else started feature] '
                    master_dict = merge(master_dict)
            msg0 = '\r\r[ ' + halo + ' ] '
            msg1 = '[ r:' + str(annulus) + ' Kpc ] '
            msg2 = '[ id: ' + str(feature_id) + ' ] '
            msg3 = '[ nbx:' + str(n_boxes) + ' ] '
            msg4 = '[ deg:' + str(round(degree, 2)) + ' ] '
            msg = msg0 + msg1 + msg2 + msg3 + msg4 + m5
            sys.stdout.write(msg)
            sys.stdout.flush()
        if m5 == ' [nothing happend] ':
            continue
        master_dict = merge(master_dict)
        if len(master_dict.keys()) > 100:
            master_dict = merge(master_dict)
            master_dict = merge(master_dict)
    for key in master_dict.keys():
        master_dict[key] = _plot(master_dict[key])
    logfile_fh = os.path.join(
        Config.get('PATH', 'plot_dir'),
        'groupfinder_testplots',
        halo + '_logfile')
    with open(logfile_fh, 'w') as logfile:
        logfile.write('\n' + '-' * 70 + '\n')
        for key in table.meta.keys():
            if key == 'satids':
                logfile.write(
                    ' | ----------------------------------------------------------- \n')
                logfile.write(' | [table] ' + key + '\n')
                lines = np.array_split(np.asarray(table.meta['satids']), 15)
                for line in lines:
                    l = line.tolist()
                    logfile.write(' | ' + str(l) + '\n')
                logfile.write(
                    ' | ----------------------------------------------------------- \n')
            else:
                logfile.write('[table] ' + key + ': ' +
                              str(table.meta[key]) + '\n')
        logfile.write('\n' + '-' * 70 + '\n')
        logfile.write('-' * 70 + '\n')
        logfile.write('[Features]            : ' +
                      str(len(master_dict.keys())) + '\n')
        logfile.write('[r_start]             : ' + str(r_start) + ' Kpc\n')
        logfile.write('[r_stop]              : ' + str(r_stop) + ' Kpc\n')
        logfile.write('[deg_start]           : ' + str(deg_start) + ' deg\n')
        logfile.write('[deg_stop]            : ' + str(deg_end) + ' deg\n')
        logfile.write('[xbox_min_value]      : ' + str(xbox_min_value) + '\n')
        logfile.write('[percent]             : ' + str(percent) + ' %\n')
        logfile.write('[annulus_degree_step] : ' +
                      str(annulus_degree_step) + ' degrees\n')
        reports = 0
        for f_id in master_dict.keys():
            if reports > 25 or master_dict[f_id]['nboxes'] < 10:
                continue
            reports += 1
            feature = master_dict[f_id]
            logfile.write('\n' + '-' * 70 + '\n')
            logfile.write('[feature number: ' + str(f_id) + ' ]\n')
            for key in feature.keys():
                if key == 'points':
                    continue
                if key == 'sats_book':
                    logfile.write('\n  ---------------------\n')
                    logfile.write('        SATIDS\n\n')
                    for satid in feature['sats_book'].keys():
                        if not feature['sats_book'][satid]:
                            continue
                        percentage = round(
                            (1e2 * feature['sats_book'][satid]) / float(feature['nstars_total']), 2)
                        if percentage > 75.0:
                            logfile.write('  ->[' + str(satid) + ']: ' + str(feature['sats_book'][
                                          satid]) + '         :' + str(percentage) + ' %' + '\n')
                        else:
                            logfile.write('  -> ' + str(satid) + ' : ' + str(feature['sats_book'][
                                          satid]) + '         :' + str(percentage) + ' %' + '\n')
                    logfile.write('\n  ---------------------\n')
                else:
                    logfile.write('  --> ' + key + ' : ' +
                                  str(feature[key]) + '\n')
