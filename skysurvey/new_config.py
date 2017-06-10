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

from os.path import join
import os
import platform

try:
    prefix = os.environ['USERPROFILE']
except Exception, e:
    raise e

config_filename = 'setup.cfg'
main_path = join(prefix, 'GitHub', 'skysurvey')
data_path = join(main_path, 'data')
config_fh = join(main_path, config_filename)
src_data_path = join(main_path, 'skysurvey', 'data')
halo_ebf_path = 'C:\Users\swc21\Halos - Copy'#raw_input('--> path to halo ebf directory:')

_defaults = {
    'Global': [
        ('verbose', 0)
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
            ('halos','halo02,halo05,halo07,halo08,halo09,halo10,halo12,halo14,halo15,halo17,halo20')
    ],
    
    'Platform': [
        ('platform', platform.platform(aliased=0, terse=0))
    ],

    'PATH': [
        ('prefix', prefix),
        ('main_path', main_path),
        ('src_data_path', src_data_path),
        ('data_path', data_path),
        ('halo_ebf_dir', halo_ebf_path),
        ('plot_dir', join(data_path, 'plots')),
        ('notebook_dir', join(main_path, 'notebooks')),
        ('grid_dir', join(data_path, 'grids')),
        ('halo_dir', join(data_path, 'halos')),
        ('table_dir', join(data_path, 'tables')),
        ('paper', join(main_path, 'paper')),
        ('text_dir', join(data_path, 'texts'))
    ],

    'Sys_keys': [
        ('keys', 'OS,USERNAME,USERPROFILE,PROCESSOR_IDENTIFIER,PROCESSOR_ARCHITECTURE,NUMBER_OF_PROCESSORS')
    ],

    'Data_keys': [
        ('positions','px,py,pz'), 
        ('mags', 'wfirst-hst_f110w,wfirst-hst_f160w,wfirst-hst_f184,wfirst-hst_f475w,wfirst-hst_f555w,wfirst-hst_f606w,wfirst-hst_f814w,wfirst-hst_h158,wfirst-hst_j129,wfirst-hst_w149,wfirst-hst_y106,wfirst-hst_z087'),
        ('data', 'teff,age,satid,mact,smass')
    ],

    'Data_cut': [
        ('cut', 10)
    ]

    }

def new_cfg(config=ConfigParser.RawConfigParser(), defaults=_defaults, cfg_fh=config_fh):
    '''Make a new configuration file <setup.cfg> in the top level skysurvey directory.
    
    This will erase and replace any existing cfg file
    
    Keyword Arguments:
        defaults {[dict]} -- [dic of values] (default: {_defaults})
        cfg_fh {[str]} -- [path to new cfg file] (default: {config_fh})
    '''
    print('-'*25 + '\nnew configuration file:\n' + '-'*25)
    for key in _defaults.keys():
        print('['+ key + ']')
        config.add_section(key)
        for option, value in _defaults[key]:
            print(option,'=', value)
            config.set(key, option, value)
        print('')
    print('saving ', config_fh)
    with open(config_fh, 'w') as configfile:
        config.write(configfile)
    print('done')
if __name__ == '__main__':
    import ConfigParser
    import sys
    new_cfg()
    sys.exit(0)