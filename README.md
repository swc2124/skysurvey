# skysurvey.
Python tools to make WFIRST survey strategies.  Library full of easy plotting functions.  Lots of automatic functions and features for working stellar data.

# basic definitions.

d_mpc = distance to target in Mpc. Default value: 4.0
f_type = WFIRST filter type

# PATH options

# Setting and resetting the config file


ebf.py is needed to handle .ebf files output from Galaxia [link](http://galaxia.sourceforge.net/Galaxia3pub.html)
# Skysurvey's file system.
'''
|-- data
    |-- halos
        |-- halo02
            |--  teff
            |--  age
            |--  px
            |--  etc.
        |-- halo05
            |--  teff
            |--  age
            |--  px
            |--  etc.
        |-- etc.
    |-- grids
        |-- fh_halo02_grid.npy
        |-- fh_halo03_grid.npy
        |-- etc.
    |-- plots
        |-- type_1
            |-- plot_1.png
            |-- plot_2.png
            |-- etc.
        |-- type_2
            |-- plot_1.png
            |-- plot_2.png
            |-- etc.
        |-- etc.
    |-- tables
        |-- table_1.hdf5
        |-- table_2.tex
        |-- etc
    |-- text
        |-- info.pdf
        |-- reference.pdf
    |-- etc.
|-- skysurvey
    |-- __init__.py
    |-- functions.py
    |-- grid.py
    |-- gridplot.py
    |-- options.py
    |-- spinbin.py
    |-- targetlist.py
    |-- etc.
|-- notebooks
    |-- Various IPython notebooks.
    |-- nb_1.ipynb
    |-- nb_2.ipynb
    |-- etc.
|-- paper
    |-- paper
        |-- paper.tex
        |-- paper.pdf
        |-- lib.bib
        |-- etc.
    |-- tables
        |-- tab_1.tex
        |-- tab_2.tex
        |-- etc.
    |-- etc.
'''
