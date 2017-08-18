"""
TODO
-------------------------------------------------------------------------------
 Author     : Sol W. Courtney
 Location   : Columbia U Dept. of Astronomy & Astrophysics NYC, Ny
 Email      : swc2124@Columbia.edu
 Date       : June 2017
 Title      : skysurvey/skysurvey/new_config.py
-------------------------------------------------------------------------------
"""
from __future__ import division, absolute_import, print_function

import os
import platform
import ConfigParser
import skysurvey
_my_halos = 'C:\Users\swc21\Halos - Copy'
SYS_CFG_FNAME = 'sys_conf.cfg'


def _select_config_file_location(prfx):
    selection = raw_input('interactive mode? yes/[NO]')
    if selection in ['yes', 'YES', 'y', 'Y']:
        print('prefix = ', prefix)
        dir_name = raw_input('enter new top leve dir name for skysurvey:')
    else:
        dir_name = 'skysurvey'
    return dir_name


def ask_halo_ebf_path():
    return raw_input('enter the location of the halo .ebf files:')


def new_cfg(halo_ebf_path=_my_halos, sys_config_filename=SYS_CFG_FNAME):
    '''Make a new configuration file <setup.cfg> in the top level skysurvey directory.

    This will erase and replace any existing cfg file

    Keyword Arguments:
        defaults {[dict]} -- [dic of values] (default: {_defaults})
        cfg_fh {[str]} -- [path to new cfg file] (default: {config_fh})
    '''
    config_filename = 'setup.cfg'
    config = ConfigParser.RawConfigParser()
    if 'HOME' in os.environ.keys():
        prefix = os.environ['HOME']
    if 'USERPROFILE' in os.environ.keys():
        prefix = os.environ['USERPROFILE']
    package_path = os.path.dirname(os.path.realpath(skysurvey.__file__))
    sys_config_fh = os.path.join(package_path, sys_config_filename)
    skysurvey_topdir = os.path.join(prefix, 'skysurvey')
    if not os.path.isdir(skysurvey_topdir):
        new_dir_name = _select_config_file_location(prefix)
        skysurvey_topdir = os.path.join(prefix, new_dir_name)
        print('making skysurvey top directory:')
        os.mkdir(skysurvey_topdir)

    config_fh = os.path.join(skysurvey_topdir, config_filename)
    print('new config file :', config_fh)

    if not os.path.isdir(halo_ebf_path):
        while not os.path.isdir(halo_ebf_path):
            halo_ebf_path = ask_halo_ebf_path()

    data_path = os.path.join(skysurvey_topdir)
    defaults = {
        'halo_data_settings': [
            ('radial_cut', 300.)
        ],
        'grid_options': [
            ('size', 600),
            ('n_slices', 15)
        ],
        'options': [
            ('zip_safe', False),
            ('include_package_data', True)
        ],
        'metadata': [
            ('name', 'skysurvey'),
            ('version', 'attr: 0.0.1'),
            ('license', 'MIT'),
            ('classifiers', 'Programming Language :: Python :: 2.7')
        ],
        'Global': [
            ('verbose', True)
        ],

        'Distance': [
            ('d_mpc', 4.0)
        ],

        'Filter': [
            ('filter_prfx', 'wfirst-hst_'),
            ('filter_type', 'h158')
        ],

        'Filter_limits': [
            ('z087', 27.15), ('y106', 27.13),
            ('j129', 27.14), ('h158', 27.12),
            ('f184', 26.15), ('w149', 27.67)
        ],

        'Default_magnitude_limits': [
            ('mag_lim_max', 26.7),
            ('mag_lim_mid', 23.5),
            ('mag_lim_low', 22.0)
        ],

        'Default_halos': [
            ('halos', 'halo02,halo05,halo07,halo08,halo09,halo10,halo12,halo14,halo15,halo17,halo20')
        ],

        'Platform': [
            ('platform', platform.platform(aliased=0, terse=0))
        ],

        'PATH': [
            ('prefix', prefix),
            ('skysurvey_topdir', skysurvey_topdir),
            ('src_data_path', os.path.join(package_path, 'data')),
            ('data_path', data_path),
            ('halo_ebf_dir', halo_ebf_path),
            ('data_dir', os.path.join(data_path, 'data')),
            ('plot_dir', os.path.join(data_path, 'plots')),
            ('text_dir', os.path.join(data_path, 'texts')),
            ('grid_dir', os.path.join(data_path, 'grids')),
            ('halo_dir', os.path.join(data_path, 'halos')),
            ('table_dir', os.path.join(data_path, 'tables')),
            ('paper', os.path.join(skysurvey_topdir, 'paper')),
            ('notebook_dir', os.path.join(skysurvey_topdir, 'notebooks'))

        ],

        'Sys_keys': [
            ('keys', 'OS,USERNAME,USERPROFILE,PROCESSOR_IDENTIFIER,PROCESSOR_ARCHITECTURE,NUMBER_OF_PROCESSORS')
        ],

        'Data_keys': [
            ('positions', 'px,py,pz'),
            ('mags', 'wfirst-hst_f110w,wfirst-hst_f160w,wfirst-hst_f184,wfirst-hst_f475w,wfirst-hst_f555w,wfirst-hst_f606w,wfirst-hst_f814w,wfirst-hst_h158,wfirst-hst_j129,wfirst-hst_w149,wfirst-hst_y106,wfirst-hst_z087'),
            ('data', 'teff,age,satid,mact,smass,feh,alpha,lum')
        ],

        'Data_cut': [
            ('cut', 1)
        ]

    }
    print('-' * 25 + '\nnew configuration file:\n' + '-' * 25)
    for key in defaults.keys():
        print('[' + key + ']')
        config.add_section(key)
        for option, value in defaults[key]:
            print(option, '=', value)
            config.set(key, option, value)
        print('')
    print('saving', config_filename)
    with open(config_fh, 'w') as configfile:
        config.write(configfile)
    sys_config = ConfigParser.RawConfigParser()
    sys_config.add_section('skysurvey_global_settings')
    sys_config.set('skysurvey_global_settings', 'package_path', package_path)
    sys_config.set('skysurvey_global_settings',
                   'skysurvey_topdir', skysurvey_topdir)
    sys_config.set('skysurvey_global_settings',
                   'sys_config_filename', sys_config_filename)
    sys_config.set('skysurvey_global_settings', 'sys_config_fh', sys_config_fh)
    sys_config.set('skysurvey_global_settings',
                   'config_filename', config_filename)
    sys_config.set('skysurvey_global_settings', 'config_fh', config_fh)
    with open(sys_config_fh, 'wb') as configfile:
        sys_config.write(configfile)
    print('done')
