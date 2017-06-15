"""
-------------------------------------------------------------------------------
 Author     : Sol W. Courtney
 Location   : Columbia U Dept. of Astronomy & Astrophysics NYC, Ny
 Email      : swc2124@Columbia.edu
 Date       : June 2017
 Title      : skysurvey/setup.py
-------------------------------------------------------------------------------
"""
import os


from setuptools.extension import Extension
from setuptools import setup
from setuptools import find_packages
import platform
from numpy import get_include
try:
    from Cython.Build import cythonize
    from Cython.Distutils import build_ext
except:
    print("You don't seem to have Cython installed. Please get a")
    print("copy from www.cython.org and install it")
    import sys
    sys.exit(1)

with open('./README.md') as f:
    Long_description = f.readlines()

ext_modules = Extension('c_functions',
                        [os.path.join('c_functions', 'c_functions.pyx')],
                        extra_compile_args=['-fopenmp'],
                        extra_link_args=['-fopenmp'])


def list_data_files():
    _include_list = []
    _data_dirs = ['data', 'paper', 'tables', 'text']
    _file_types = ['.pdf', '.bib', '.tex', '.ebf', '.hdf5', '.txt']
    _pkg_data_path = os.path.join(os.path.curdir, 'skysurvey', 'data')
    for directory in _data_dirs:
        for name in os.listdir(os.path.join(_pkg_data_path, directory)):
            lcl_pth = os.path.join(_pkg_data_path, directory, name)
            if os.path.isdir(lcl_pth):
                for name_ in os.listdir(lcl_pth):
                    lcl_lcl_pth = os.path.join(lcl_pth, name_)
                    if os.path.isfile(lcl_lcl_pth):
                        if os.path.splitext(lcl_lcl_pth)[-1].lower() in _file_types:
                            _include_list.append(
                                lcl_lcl_pth.replace('.\\skysurvey\\', ''))
            else:
                if os.path.splitext(lcl_pth)[-1].lower() in _file_types:
                    _include_list.append(lcl_pth.replace('.\\skysurvey\\', ''))
    return _include_list

setup(

    name='skysurvey',
    version='0.1.0',

    author='Sol W. Courtney',
    author_email='swc2124@Columbia.edu',

    maintainer='Sol W. Courtney',
    maintainer_email='swc2124@Columbia.edu',

    url='https://github.com/swc2124/skysurvey',

    description=(
        'Python tools for working with stellar data '
        'in the form of numpy arrays.  Intended to '
        'be used for defining WFIRST survey strategies'
    ),
    download_url='https://github.com/swc2124/skysurvey.git',
    license='MIT',
    include_package_data=True,
    include_dirs=[get_include()],

    packages=['skysurvey', 'c_functions'],
    package_dir={'skysurvey': 'skysurvey', 'c_functions': 'c_functions'},
    package_data={'skysurvey': list_data_files()},
    install_requires=['docutils>=0.3', 'numpy>=0.x',
                      'ebfpy>=0.x', 'astropy>=0.x', 'cython>=0.x'],
    cmdclass={'build_ext': build_ext},
    ext_modules=cythonize(ext_modules, include_path=[get_include()]),
    entry_points={
        'console_scripts': [
            'skysurvey-mkfs = skysurvey.makefilesystem:fsinit',
            'skysurvey-new-cfg = skysurvey.new_config:new_cfg',
            'skysurvey-nospin-binall = skysurvey.spinbin:nospin_binall',
            'skysurvey-table-bin = skysurvey.table_bin:table_bin',
            'skysurvey-table-merge = skysurvey.table_bin:table_merge',


        ]
    },


)

#-----------------------------------------------------------------------------#
