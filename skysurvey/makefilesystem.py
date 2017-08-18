"""[summary]

[description]
"""
from __future__ import division, absolute_import, print_function

import os
from os.path import join as join
import sys
from shutil import copyfile

import ebf
import numpy as np

import ConfigParser
from .new_config import SYS_CFG_FNAME
import skysurvey


__all__ = ['fsinit']


def fsinit():
    '''[summary]

    [description]

    Keyword Arguments:
        halo_ebf_dir {[type]} -- [description] (default: {Config.get('PATH','halo_ebf_dir')})
    '''
    sys_config_fh = os.path.join(os.path.dirname(
        os.path.realpath(skysurvey.__file__)), SYS_CFG_FNAME)
    SysConfig = ConfigParser.ConfigParser()
    SysConfig.read(sys_config_fh)
    config_fh = SysConfig.get('skysurvey_global_settings', 'config_fh')

    Config = ConfigParser.ConfigParser()
    Config.read(config_fh)

    for section in Config.sections():
        print(section)
        for item in Config.items(section):
            print(item)
    halo_ebf_dir = Config.get('PATH', 'halo_ebf_dir')
    print('\n\n')
    for key in Config.get('Sys_keys', 'keys').split(','):
        print(key + ':', os.environ[key])
    print('\n\n')
    halos = Config.get('Default_halos', 'halos').split(',')
    for name, _path in Config.items('PATH'):

        if not os.path.isdir(_path):
            print('Making PATH: [' + name + ']', _path)
            os.mkdir(_path)

            if name == 'table_dir':
                src_table_path = join(Config.get(
                    'PATH', 'src_data_path'), 'tables')
                for file in os.listdir(src_table_path):
                    src = join(src_table_path, file)
                    dst = join(_path, file)
                    print('copy:', src)
                    print('to  :', dst)
                    copyfile(src, dst)
                print('[' + name + '] done')

            if name == 'data_dir':
                src_table_path = join(Config.get(
                    'PATH', 'src_data_path'), 'data')
                for file in os.listdir(src_table_path):
                    src = join(src_table_path, file)
                    dst = join(_path, file)
                    print('copy:', src)
                    print('to  :', dst)
                    copyfile(src, dst)
                print('[' + name + '] done')

            if name == 'paper':
                src_paper_dir = join(Config.get(
                    'PATH', 'src_data_path'), 'paper')
                for src_dir in os.listdir(src_paper_dir):
                    new_paper_dir = join(_path, src_dir)
                    old_paper_dir = join(src_paper_dir, src_dir)
                    print('Making PATH: [' + name + ']', new_paper_dir)
                    os.mkdir(new_paper_dir)
                    for _file in os.listdir(old_paper_dir):
                        src = join(old_paper_dir, _file)
                        dst = join(new_paper_dir, _file)
                        print('copy:', src)
                        print('to  :', dst)
                        copyfile(src, dst)
                print('[' + name + '] done')
        else:
            print('PATH Exists: [' + name + ']', _path)

    for _halo in os.listdir(halo_ebf_dir):
        halo_fh = _halo.split('.')[0]
        if not halo_fh in halos:
            print('    ' + 'skipping ' + halo_fh)
            continue

        print('    ' + 'starting ' + halo_fh)
        new_path = join(Config.get('PATH', 'halo_dir'), halo_fh)
        if not os.path.isdir(new_path):
            print('    ' + 'Making PATH: [' + halo_fh + ']', new_path)
            os.mkdir(new_path)
        else:
            print('    ' + 'PATH Exists: [' + halo_fh + ']', new_path)

        data_keys = []
        for k in [value for item in Config.items('Data_keys') for value in item[1].split(',')]:
            if not os.path.exists(join(new_path, k + '.npy')):
                data_keys.append(k)
                print('    ' + 'Adding Key: [' + k + ']')
            else:
                print('    ' + 'Key Exists: [' + k + ']')
        if data_keys:
            print('    ' + 'loading and saving new keys')
            halo_array = ebf.read(join(halo_ebf_dir, _halo))
            for data_key in data_keys:
                key_fh = join(new_path, data_key)
                print('    --> ' + data_key + ': ' + key_fh)
                np.save(key_fh, halo_array[data_key][
                        ::Config.getint('Data_cut', 'cut')])
        else:
            print('    ' + 'nothing to do here')
        print('    ' + 'done' + '\n')

    print('\n[ process complete ]\n')
