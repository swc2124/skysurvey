"""TODO"""
from __future__ import division, absolute_import, print_function

from warnings import warn
from os import listdir
from os.path import join as join

from numpy import float64
from numpy import array
from numpy import load
from numpy import save
from numpy import arctan
from numpy import square
from numpy import pi
from numpy import log10
from numpy import linspace
from numpy import logical_and


import matplotlib
matplotlib.use('AGG')
from matplotlib import pyplot as plt
'''
from matplotlib.pyplot import figure
from matplotlib.pyplot import close
from matplotlib.pyplot import grid
from matplotlib.pyplot import legend
from matplotlib.pyplot import cm
from matplotlib.pyplot import pcolormesh
from matplotlib.pyplot import colorbar
'''

from .options import HALO_DATA_DIR
from .options import savedplot_path
from .options import data_dir
from .options import halos
from .options import distance_Mpc
from .options import grid_dir
from .options import filter_type
from .options import mag_lim_max
from .options import mag_lim_mid
from .options import mag_lim_low

from .spinbin import spinone


def sqrdeg(d_Mpc=distance_Mpc):
    '''
    conversion mod to deg^2
    '''
    # edge length of box _Grid - 6 Kpc
    d = 6 * 1e3
    # distance in Mpc - 0.75 ~ 11.0
    D = d_Mpc * 1e6
    # small angle formula
    rad = 2 * (arctan((d / (2 * D))))
    # convert radians to degrees
    deg = rad * (180 / pi)
    # square for area
    return square(deg)

def sqr_arcmin(d_Mpc=distance_Mpc):
    '''
    conversion mod to arcmin^2
    '''
    # edge length of box _Grid - 6 Kpc
    d = 6.0 * 1e3
    # distance in Mpc - 0.75 ~ 11.0
    
    D = d_Mpc * 1e6

    return square(3437.75 * (2 * (arctan((d / (2 * D))))))

def _Grid_list(_Grid_path=grid_dir, _d_Mpc=distance_Mpc, _f_type=filter_type):
    '''
    Make a list of output files writing to be plotted.  These files are
        the output of spinbin.spinall ro spin.
    '''
    print('reading arrays from : ' + _Grid_path)
    print('[' + str(_d_Mpc) + ']')
    _Grids = [i for i in listdir(_Grid_path) if i.endswith(
        '_' + str(_d_Mpc) + 'Mpc_' + _f_type + '_grid.npy')]
    if len(_Grids) > 0:
        _Grids.sort()
        for _Grid in _Grids:
            print(' -> ' + _Grid)
        return _Grids
    else:
        warn((
            'there are no arrays to be plotted!'
            'Have you run < spinbin.spinall() > yet?'))
        return []


def nstars_per_sqdeg(d_Mpc=distance_Mpc, plot_path=savedplot_path, f_type=filter_type):
    """
    Render and save the main plot pannal.

    Extended description of function.

    Parameters
    ----------
    d_Mpc : float
        The distance to the halo.

    plot_path : str
        /path/to/plotdir

    Returns
    -------
    plot
        Description of return value

    """
    # figure dimensions
    fig_x = 31
    fig_y = 23
    plotname = 'log(n_stars(m_app<26.87))/arcmin^2 vs. radial separation Kpc'
    x0 = 0
    x1 = 301

    xticks = [str(round(i)) for i in linspace(0, 300, 9)]
    x_label = 'Radial Distance (Kpc)'
    y_label = 'Log(N_stars)/arcmin^2'
    mod = sqr_arcmin(d_Mpc=d_Mpc)
    fig = plt.figure(figsize=(fig_x, fig_y))
    fig.suptitle(plotname + '\n' + str(d_Mpc) +
                 ' Mpc  ' + f_type, fontsize=40)
    _Grids = _Grid_list(_d_Mpc=d_Mpc, _f_type=f_type)
    if len(_Grids) == 0:
        warn('There are no arrays to be plotted')
        return
    min_strs = []
    max_strs = []
    print('finding min and max for y-axis')
    for fh in _Grids:
        filename = join(grid_dir, fh)
        arr = load(filename)
        strs = log10(arr[:, :, 3:6][arr[:, :, 0] > 0.].flatten() / mod)
        min_strs.append(strs[strs != -float('Inf')].min())
        max_strs.append(strs.max())
    y0 = min(min_strs)
    y1 = max(max_strs)
    print('y-axis min: ', round(y0, 2))
    print('y-axis max: ', round(y1, 2))
    print('plotting')
    for i, _Grid in enumerate(_Grids):
        print(' -> ' + _Grid)
        filename = join(grid_dir, _Grid)
        array = load(filename)
        ax = fig.add_subplot(3, 4, i + 1)
        header = _Grid.split('_')
        title = header[0]
        ax.set_title(title, fontsize=30)
        #mlim_min = log10(array[:, :, 3][array[:, :, 0] > 0.].flatten() / mod)
        #mlim_med = log10(array[:, :, 4][array[:, :, 0] > 0.].flatten() / mod)
        mlim_max = log10(array[:, :, 5][array[:, :, 0] > 0.].flatten() / mod)
        rads = array[:, :, 6][array[:, :, 0] > 0.].flatten()
        ax.scatter(rads, mlim_max, s=5, alpha=.75, color='k',
                   marker='o')#, label='(10^5) blue: m<' + str(mag_lim_max))
        '''
        ax.axhline(y=1.44, xmin=0, xmax=301, c="blue",
                   linewidth=7, alpha=.2, zorder=0)
        ax.scatter(rads, mlim_med, s=9, alpha=.5, color='red',
                   marker='o', label='(10^4.6) red: m<' + str(mag_lim_mid))
        ax.axhline(y=1.04, xmin=0, xmax=301, c="red",
                   linewidth=7, alpha=.2, zorder=0)
        ax.scatter(rads, mlim_min, s=9, alpha=.5, color='black',
                   marker='o', label='(10^4) black: m<' + str(mag_lim_low))
        ax.axhline(y=.44, xmin=0, xmax=301, c="black",
                   linewidth=7, alpha=.2, zorder=0)
        #plt.legend(fontsize='large', numpoints=1)
        '''
        plt.grid()
        ax.set_xlim([x0, x1])
        # ax.set_xticklabels(xticks)
        ax.set_xlabel(x_label, fontsize=25)
        ax.set_ylim([y0, y1])
        ax.set_ylabel(y_label, fontsize=25)

    #ax12 = fig.add_subplot(3, 4, 12)
    #ax12.axis('off')
    #ax12.legend(['(10^5) blue: m<' + str(mag_lim_max), '(10^4.6) red: m<' + str(mag_lim_mid), '(10^4) black: m<' + str(mag_lim_low)],fontsize=30, shadow=True)
    full_plot_filename = savedplot_path + 'sqdeg_' + \
        str(d_Mpc).split('.')[0] + 'Mpc_' + f_type
    fig.savefig(full_plot_filename)
    plt.close()


def nstars(d_Mpc=distance_Mpc, plot_path=savedplot_path, f_type=filter_type, lims=[100.0, 300.0]):
    fig_x = 31
    fig_y = 23
    plotname = 'Log number of stars per square arcmin\n'
    x0 = 0
    x1 = 100
    y0 = 0
    y1 = 100
    x_label = 'Kpc'
    y_label = 'Kpc'
    ticks = linspace(-300, 300, 6)
    x_ticklabels = [str(round(i)) for i in ticks]
    y_ticklabels = [str(round(i)) for i in ticks]
    x_ticklabels[0] = ''
    r0, r1 = lims
    mod = sqr_arcmin(d_Mpc=d_Mpc)
    fig = plt.figure(figsize=(fig_x, fig_y))
    fig.suptitle(plotname + '   ' + str(d_Mpc) +
                 ' Mpc  ' + f_type + ' m_app < 26.78', fontsize=50)
    _Grids = _Grid_list(_d_Mpc=d_Mpc, _f_type=f_type)
    if len(_Grids) == 0:
        warn('There are no arrays to be plotted')
        return
    print('plotting')
    for i, _Grid in enumerate(_Grids):
        filename = join(grid_dir, _Grid)
        print(' -> ' + filename)
        array = load(filename)
        ax = fig.add_subplot(3, 4, i + 1)

        header = _Grid.split('_')
        title = header[0]
        ax.set_title(title, fontsize=30)
        hm = plt.pcolormesh(
            log10(array[:, :, 5] / mod), cmap=plt.cm.plasma, vmin=-1, vmax=3.5)
        cp = ax.contour(array[:, :, 6], [r0, r1], colors='white',
                        linewidths=3, alpha=.6, linestyles='dashed')
        ax.clabel(cp, [r0, r1], inline=1, fmt='%s Kpc',
                  fontsize=15, color='white', linewidth=7, alpha=1)
        # hm.set_linewidth(0.01)
        if i == 10:
            ax1 = fig.add_subplot(3, 4, i + 2)
            cb = plt.colorbar(hm, ax1)
            ax1.axes.grid(color='white', alpha=.8,
                          linewidth=5, linestyle='dashed')
            cb.ax.tick_params(labelsize=28)
            ax1.yaxis.set_label_position("right")
            ax1.set_ylabel("Log Nstars/arcmin^2",fontsize=40)

        ax.axes.grid(color='white', alpha=.4, linewidth=2, linestyle='dashed')
        ax.axes.set_yticklabels(y_ticklabels, size=20, rotation=-30)
        ax.axes.set_xticklabels(x_ticklabels, size=20, rotation=-30)
        ax.axes.set_xlim([x0, x1])
        ax.axes.set_ylim([y0, y1])
        if i in [0, 4, 8]:
            ax.axes.set_ylabel(y_label, fontsize=40)
        if i in [8, 9, 10]:
            ax.axes.set_xlabel(x_label, fontsize=40)

    full_plot_filename = savedplot_path + 'nstars_' + \
        str(d_Mpc).split('.')[0] + 'Mpc_' + f_type + '.png'
    fig.savefig(full_plot_filename, dpi=fig.dpi)
    plt.close()


def mixplot(plot_halos=['halo12', 'halo15', 'halo20'],
            d_Mpc=distance_Mpc, plot_path=savedplot_path,
            f_type=filter_type, radius=None):
    plt_halos = []
    print('loading filenames')
    for i in _Grid_list(_d_Mpc=d_Mpc, _f_type=f_type):
        breakdown = i.split('_')
        if breakdown[0] in plot_halos:
            plt_halos.append(i)
            print(' -> [ selected ] ', i)
    fig_x = 42
    fig_y = 40
    plotname = 'Log Number of stars per square arcmin'
    fig = plt.figure(figsize=(fig_x, fig_y))
    fig.suptitle(plotname + '\n' + str(d_Mpc) +
                 ' Mpc  ' + f_type + ' m_app < 26.78', fontsize=80)
    mod = sqr_arcmin(d_Mpc=d_Mpc)
    min_strs = []
    max_strs = []
    print('finding min and max for y-axis')
    for fh in plt_halos:
        filename = join(grid_dir, fh)
        arr = load(filename)
        strs = log10(arr[:, :, 3:6][arr[:, :, 0] > 0.].flatten() / mod)
        min_strs.append(strs[strs != -float('Inf')].min())
        max_strs.append(strs.max())
    _y0 = min(min_strs)
    _y1 = max(max_strs)
    print('y-axis min: ', round(_y0, 2))
    print('y-axis max: ', round(_y1, 2))
    for p_num, _file in enumerate(plt_halos):
        if radius:
            r0, r1 = radius[p_num]
        else:
            r0, r1 = 0, 1
        header = _file.split('_')
        title = header[0]
        filename = join(grid_dir, _file)
        print(' -> ' + filename)
        array = load(filename)
        for i in range(3):
            
            plot_number = (p_num * len(plt_halos)) + (i + 1)
            
            print('plotting ', plot_number, title)

            ax = fig.add_subplot(3, 3, plot_number)

            if i in [0, 3, 6]:
                plotname = 'Log(N_stars)/arcmin^2'
                x0 = 0
                x1 = 301
                y0 = _y0
                y1 = _y1
                x_ticklabels = [str(int(i)) for i in linspace(0, 300, 7)]
                #y_ticklabels = [str(i) for i in [-1,0,1,2,3]]
                x_label = 'Radial Distance (Kpc)'
                y_label = title + '\nlog(n_stars)_vs_arcmin^2'
                '''
                mlim_min = log10(
                    array[:, :, 3][array[:, :, 0] > 0.].flatten() / mod)
                mlim_med = log10(
                    array[:, :, 4][array[:, :, 0] > 0.].flatten() / mod)
                '''
                mlim_max = log10(
                    array[:, :, 5][array[:, :, 0] > 0.].flatten() / mod)

                rads = array[:, :, 6][array[:, :, 0] > 0.].flatten()
                idx = logical_and(rads > r0, rads < r1)
                
                ax.scatter(rads[idx], mlim_max[idx], s=80, alpha=.6, color='orange',
                           marker='o')
                '''
                ax.scatter(rads[idx], mlim_med[idx], s=80, alpha=.6, color='orange',
                           marker='o')
                ax.scatter(rads[idx], mlim_min[idx], s=80, alpha=.6, color='orange',
                           marker='o')
                '''
                ax.scatter(rads, mlim_max, s=15, alpha=.5, color='k',
                           marker='o')#, label='(10^5) blue: m<' + str(mag_lim_max))
                '''
                ax.axhline(y=1.44, xmin=0, xmax=301, c="blue",
                           linewidth=7, alpha=.2, zorder=0)
                ax.scatter(rads, mlim_med, s=9, alpha=.5, color='red',
                           marker='o', label='(10^4.6) red: m<' + str(mag_lim_mid))
                ax.axhline(y=1.04, xmin=0, xmax=301, c="red",
                           linewidth=7, alpha=.2, zorder=0)
                ax.scatter(rads, mlim_min, s=9, alpha=.5, color='black',
                           marker='o', label='(10^4) black: m<' + str(mag_lim_low))
                ax.axhline(y=.44, xmin=0, xmax=301, c="black",
                           linewidth=7, alpha=.2, zorder=0)
                '''
                ax.axes.legend(fontsize='xx-large', numpoints=1, shadow=True)
                ax.axes.grid(color='black', alpha=.5,
                                 linewidth=3, linestyle='dashed')
                y_ticklabels = [str(i) for i in range(-1,4)]
                ax.axes.set_yticks([-1.0,0.0,1.0,2.0,3.0])
                #ax.axes.tick_params(axis=)

            else:

                hm = ax.pcolormesh(
                    log10(array[:, :, 0] / mod), cmap=plt.cm.plasma, vmin=-1, vmax=3.5)
                cp = ax.contour(array[:, :, 6], [r0, r1], colors='white',
                                linewidths=7, alpha=.5, linestyles='dashed')
                ax.clabel(cp, [r0, r1], inline=1, fmt='%s Kpc',
                          fontsize=25, color='white', linewidth=7, alpha=.9)
                if i in [1, 4, 7]:
                    plotname = 'n_stars'
                    x0 = 0
                    x1 = 100
                    y0 = 0
                    y1 = 100
                    x_label = 'Kpc'
                    y_label = 'Kpc'
                    x_ticks = linspace(-300, 300, 6)
                    x_ticklabels = [str(int(i)) for i in x_ticks]
                    y_ticks = linspace(-300, 300, 6)
                    y_ticklabels = [str(int(i)) for i in x_ticks]
                    x_ticklabels[0] = ''
                    ax.axes.grid(color='white', alpha=.2,
                                 linewidth=3, linestyle='dashed')

                else:
                    plotname = 'n_stars_zoomed'
                    x0 = 20
                    x1 = 80
                    y0 = 20
                    y1 = 80
                    x_label = 'Kpc'
                    y_label = 'Kpc'
                    x_ticks = linspace(-180, 180, 7)
                    x_ticklabels = [str(int(i)) for i in x_ticks]
                    y_ticks = linspace(-180, 180, 7)
                    y_ticklabels = [str(int(i)) for i in y_ticks]
                    x_ticklabels[0] = ''

                    cax = plt.axes([0.92, 0.15, 0.05, 0.7])
                    cbar = plt.colorbar(hm, cax=cax)
                    cbar.ax.tick_params(labelsize=40)
                    cbar.ax.grid(color='white', alpha=.8,
                                 linewidth=6, linestyle='dashed')

            ax.axes.set_yticklabels(y_ticklabels, size=40, rotation=-30)

            ax.set_title(title, fontsize=30)
            
            ax.axes.set_xticklabels(x_ticklabels, size=30)  # , rotation=-30)
            ax.axes.set_title(title, fontsize=50)
            ax.axes.set_xlim([x0, x1])
            ax.axes.set_xlabel(x_label, fontsize=30)
            ax.axes.set_ylim([y0, y1])
            ax.axes.set_ylabel(y_label, fontsize=30)
            if i in [0, 3, 6]:
                ax.set_ylabel(y_label, fontsize=50)
            print('done')
    full_plot_filename = savedplot_path + 'mixplt_' + \
        str(d_Mpc).split('.')[0] + 'Mpc_' + f_type + '.png'
    fig.savefig(full_plot_filename, dpi=fig.dpi)
    plt.close()

def make_plots(distances=[2.0,5.0], f_types=['h158'],
    plot_path=savedplot_path,
    plot_halos=['halo08', 'halo15', 'halo20'],
    radius=None):

    for f_typ in f_types:

        for dist in distances:

            nstars_per_sqdeg(d_Mpc=dist, plot_path=plot_path, f_type=f_typ)
            nstars(d_Mpc=dist, plot_path=plot_path, f_type=f_typ)
            mixplot(plot_halos=plot_halos,
                d_Mpc=dist, plot_path=plot_path,
                f_type=f_typ, radius=radius)

def run_limits(target_n=85.0, radius=75.0, n_attemts=100, distances=[2.0,5.0], f_typ='h158'):

    for distance in distances:
        mod = sqr_arcmin(d_Mpc=distance)
    
        for halo in _Grid_list(_d_Mpc=distance, _f_type=f_typ):
            attemts = 0
            best_arr = load(grid_dir+halo)
            best_stars = 0.0

            for ii in range(n_attemts):
                stars = 0
                boxes = 0
                halo_name = halo.split('_')[0]
                arr = load(grid_dir+halo)

                for i,rad in enumerate(arr[:,:,6].flatten()):
                    if rad >= radius:
                        if arr[:,:,5].flatten()[i]/mod>=target_n:
                            stars+=arr[:,:,5].flatten()[i]
                            boxes+=1.0
                curent_n = round((stars/mod)/(boxes),2)

                if curent_n>best_stars:
                    print('\nhalo: ', halo_name, '| stars: ', round(stars/mod,2), '| boxes: ', boxes, '| best n: ', best_stars, '| current n: ', curent_n, '| attempt [ ' + str(ii) + '/' + str(n_attemts)+ ' ]', '\n')
                    best_stars = curent_n
                    best_arr = arr
                    save(grid_dir+halo, best_arr)

                spinone(halo_name,m_lims=array([27.12,27.87,26.87],dtype=float64),d_mpc=distance,f_type=f_typ)
                
            save(grid_dir+halo, best_arr)
            nstars_per_sqdeg(d_Mpc=distance, f_type=f_typ)