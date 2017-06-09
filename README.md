## Welcome to Skysurvey

Skysurvey is my collection of Python tools for making mock survey strategies and outputting information pertaining to that survey.  Skysurvey has a good sized library full of easy plotting functions for making useful data comparisons and examples.  It also has lots of automatic functions and features for working with stellar data in the context of observation.  I manage and maintain this package alone so please feel free to email me at swc2124@columbia.edu with any comments or questions.

TODO...

### Installation & requirements.

    Version: Python 2.7 for now.

#### Input data.

You will need preexisting halo data from [Galaxia](http://galaxia.sourceforge.net/) in a single folder like so:
        
    ~/PATH/to/example_halo_folder/

        halo02.ebf
        halo05.ebf
        halo07.ebf
        halo09.ebf
        halo11.ebf
        halo12.ebf
        halo[...].ebf

TODO...
#### How to setup Galaxia and make BJ halos.
Downloads:
[1.) Galaxia Download](https://sourceforge.net/projects/galaxia/files).
[2.) Original Bullock and Johnston stellar halo site](http://www.astro.columbia.edu/~kvj/halos/).
[3.) Bullock and Johnston simulated N-body stellar halos (pluggable in Galaxia)](http://www.physics.usyd.edu.au/~sanjib/code/).

TODO...

#### Required Python packages.
            
    Package:
        - ebf
        - numpy
        - astropy
        - matplotlib
        - cython
        - setuptools

##### Package descriptions:

1. [ebf](http://pythonhosted.org/ebfpy/), [Galaxia](https://sourceforge.net/projects/galaxia/)

        needed to handle Galaxia's binary file output.

2. [Numpy](http://www.numpy.org/), [Memory mapping](https://docs.scipy.org/doc/numpy/reference/generated/numpy.memmap.html)

        is obviously needed.  
        Lot's of memory mapping going on (TODO... memmap options).

3. [Astropy](http://www.astropy.org/), [Data Tables](http://docs.astropy.org/en/stable/table/)

        is needed mostly for its data tables.

4. [Matplotlib](https://matplotlib.org/), [Pyplot](https://matplotlib.org/api/pyplot_api.html)

        for plotting, mostly Pyplot.

5. [Cython](http://cython.readthedocs.io/en/latest/)

        to speed things up.

6. [Setuptools](https://setuptools.readthedocs.io/en/latest/)
    
        for installation.

TODO...
#### Setting PATH options for Windows, Mac or Linux.

    Windows, Mac, Linux OS: TODO... (Table)

        External halo directory:

        User/
         └── HALO_DATA_DIR/
        
        Package-wide PATH values:

        User/
         └──main_path
             |
             └── data_path
                  |
                  ├── halo_dir/
                  ├── grid_dir/
                  ├── plot_dir/
                  ├── table_dir/
                  └── text_dir/

TODO...
#### Setting and resetting the configuration file.


TODO...
#### Installation steps.

1. install skysurvey

        pip install skysurvey -v
            or
        cd ~/
        git clone https://github.com/swc2124/skysurvey.git
        cd skysurvey
        pip install ./ -v

2. Overview of command line tools.
    
    make a new default configuration file:
    
        skysurvey -new_cfg --prefix=~/

    read a configuration file into skysurvey:

        skysurvey -configure ~/setup.cfg

    make the file system:

        skysurvey -make_fs

    make a new set of empty output grids for all halos:

        skysurvey -make_grids

    make a single new grid for one halo:

        skysurvey -make_grid halo_name

    make an example/test plot to see if everything is working:

        skysurvey -plot_example

    bin halo data into grids for all halos:

        skysurvey -binall

    TODO...

        TODO...

3. set or reset the config file setup.cfg
        
    The config file setup.cfg is where all the system defaults will be stored.
    To generate a new config file a user may run the following:

        skysurvey -new_cfg --prefix=/path/to/where/you/want/the/cfg/file

    Within the new congih file look for the following headers and select your OS by setting its value to 1.
    You may also set the path to the skysurvey file system if you wish.
    default path if value = False

        ~/setup.cfg

            [os]
            windows = 1
            mac = 0
            linux = 0

            [path]
            prefix = ~/your/new/path

            [...]


4. make file system:

        skysurvey -make_fs -v 

5. make halo output grids:

        skysurvey -make_grids -v 

6. bin halo data into grids:

        skysurvey -binall -v

            or

        skysurvey -spinnall -binall -v 

7. produce plots of binned data:

        skysurvey -plotall -v

8. spin and re-bin all halos and re-plot all grids:
    
        skysurvey -spinnall -binnall -plotall -v

TODO...
#### Universal/Default definitions/values used throughout the package.

    ~/setup.cfg

        [...]

        [Default distance]
        distance_Mpc = 4.0 

        [Default filter]
        filter_prfx = 'wfirst-hst_'
        filter_type = 'h158'

        [Other filters]
        filter_types = [
            'z087', 'y106', 'j129', 
            'h158', 'f184', 'w149'
            ]
        filter_limits = {
            'z087': 27.15,
            'y106': 27.13,
            'j129': 27.14,
            'h158': 27.12,
            'f184': 26.15,
            'w149': 27.67
            }

        [Default magnitude limits]
        mag_lim_max = 26.7
        mag_lim_mid = 23.5
        mag_lim_low = 22.0
        
        [Default halos]
        halos = [
            'halo02', 'halo05', 'halo07', 'halo08',
            'halo09', 'halo10', 'halo12', 'halo14',
            'halo15', 'halo17', 'halo20'
            ]

        [...]

The 11 [Bullock and Johnston](http://user.astro.columbia.edu/~kvj/halos/) simulation halos.

TODO...
### Getting started.

#### File naming system in skysurvey.

Skysurvey uses a single file handle naming system for all output files.  The basic order is such:

    any_file_name = '<d_mpc>Mpc_<f_type>_<halo*>_<file_type>_<ext>'

    where:
        d_mpc      = The current distance setting in Mpc ('4'). No decimal places, so '4.5' => '4_5'.
        f_type     = The current filter type ('h158').
        halo       = The name of the halo ('halo02'). This is not always needed, so 'halo## => 'allhalos' or 'NA'.
        file_type  = The tag from the function which created the file ('grid', 'plot_<type>', 'array_<type>', etc.).
        ext        = File extension ('.png', '.npy', '.pdf', etc.).

TODO...
#### Overview of skysurvey function usage.

TODO...
### Overview of the skysurvey file system.
    skysurvey/
        |
        ├── setup.py
        ├── setup.cfg
        ├── MANIFEST.in
        ├── README.md
        ├── LICENCE
        |
        ├─── data/
        |      |
        |      ├── halos/
        |      |     | 
        |      |     ├── halo02/
        |      |     |     |
        |      |     |     ├── teff.npy
        |      |     |     ├── age.npy
        |      |     |     ├── px.npy
        |      |     |     └── [...].npy
        |      |     |
        |      |     └── halo05/
        |      |           |
        |      |           ├── teff.npy
        |      |           ├── age.npy
        |      |           ├── px.npy
        |      |           └── [...].npy
        |      |
        |      ├── grids/
        |      |     |
        |      |     ├── halo_<fh>_grid.npy
        |      |     ├── halo_<fh>_grid.npy
        |      |     └── [...]_<fh>_grid.npy
        |      |
        |      ├── plots/
        |      |    |
        |      |    ├── plot_type_01/
        |      |    |      |
        |      |    |      ├── plot_01.png
        |      |    |      └── [...].png
        |      |    |
        |      |    └── plot_type_02/
        |      |           |
        |      |           ├── plot_01.png
        |      |           └── [...].png
        |      |
        |      ├── tables/
        |      |     |
        |      |     ├── target_table
        |      |     └── neargalcat.txt
        |      |
        |      └─── text/
        |            |
        |            ├── reference.pdf
        |            └── WFIRST.pdf
        |
        ├── skysurvey/
        |       |
        |       ├── __init__.py
        |       ├── functions.py
        |       ├── grid.py
        |       ├── spinbin.py
        |       ├── gridplot.py
        |       ├── targetlist.py
        |       └── [...].py
        |
        ├── paper/
        |     |
        |     ├── paper_text/
        |     |       |
        |     |       ├── paper.tex
        |     |       └── lib.bib
        |     |       
        |     |
        |     └── paper_tables/
        |             |
        |             └── tab_01.tex
        |
        └── notebooks/
                |
                └── [...].ipynb