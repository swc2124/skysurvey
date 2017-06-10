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
from os.path import join

from setuptools import Extension
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

nthreads = int(int(os.environ['NUMBER_OF_PROCESSORS']) / 4.0)
if nthreads < 1.0:
    nthreads = 1

key_words = [
    'halo', 
    'astrophysics', 
    'nbody', 
    'sol', 
    'courtney'
    ]
packages = find_packages()
package_dir = {
    '': ''
    }
requirments = [
    'docutils>=0.3', 
    'numpy', 
    'ebfpy', 
    'astropy'
    ]

platform=platform.platform()

setup(

    name='skysurvey',
    version= '0.1.0',
    license= 'MIT',
    author= 'Sol W. Courtney',
    author_email= 'swc2124@Columbia.edu',
    url= 'https://github.com/swc2124/skysurvey.git',
    description= (
        'Python tools for working with stellar data '
        'in the form of numpy arrays.  Intended to '
        'be used for defining WFIRST survey strategies'
        ),
    keywords= key_words,
    packages= packages,
    install_requires= requirments,
    cmdclass= {'build_ext': build_ext},
    ext_modules= cythonize('c_functions.pyx'),
    include_dirs= [get_include()],
    scripts=[join('bin','skysurvey-mkfs'), join('bin','skysurvey-newcfg')],


    )

#-----------------------------------------------------------------------------#
