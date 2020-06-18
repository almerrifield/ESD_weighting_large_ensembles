#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""""

(c) 2019 under a MIT License (https://mit-license.org)

To run:

python ESD_fig3_med_box_DA_fix.py -p /home/meranna/python/esd_weighting_large_ensembles/ CMIP5_med_9P_DA_fix.nc CMIP5-ALL_med_9P_DA_fix.nc CMIP5_med_9P_DA_fix.nc CMIP5-ALL_med_9P_DA_fix.nc CMIP5_med_9P_DA_fix.nc CMIP5-ALL_med_9P_DA_fix.nc CMIP5_med_9P_DA_fix.nc CMIP5-ALL_med_9P_DA_fix.nc CMIP5_med_9P_DA_fix.nc CMIP5-ALL_med_9P_DA_fix.nc -u -s box_MED_9P_DA_fix.png

"""


import os
import argparse
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import seaborn as sns

from utils_python.xarray import area_weighted_mean
from boxplot import boxplot

SAVEPATH = '../../'


def read_input():
    """Read the given configuration from the config file"""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        dest='filenames', nargs='+', type=str,
        help='')
    parser.add_argument(
        '--path', '-p', dest='path', type=str, default='',
        help='')
    parser.add_argument(
        '--no-unweighted', '-u', dest='unweighted', action='store_false',
        help='')
    parser.add_argument(
        '--title', '-t', dest='title', type=str, default=None,
        help='')
    parser.add_argument(
        '--savename', '-s', dest='savename', type=str, default=None,
        help='')
    args = parser.parse_args()
    return args


def read_data(filename, path):
    ds = xr.open_dataset(os.path.join(path, filename))
    ds = area_weighted_mean(ds)

    # for CMIP5-ALL_med_9P_DA_fix.nc
    ds['weights_ic_all'] = []
    ds['weights_ic_all'] = ds['weights_ic_all']/ds['weights_ic_all'].sum()

    ds['weights_mc_all'] = []
    ds['weights_mc_all'] = ds['weights_mc_all']/ds['weights_mc_all'].sum()

    ds['weights_q_norm_all'] = ds['weights_q']/ds['weights_q'].sum()
    ds['weights_none_all'] = np.ones(288)/288

    ds['weights_ic'] = []
    ds['weights_ic'] = ds['weights_ic']/ds['weights_ic'].sum()

    ds['weights_mc'] = []
    ds['weights_mc'] = ds['weights_mc']/ds['weights_mc'].sum()

    ds['weights_q_norm'] = ds['weights_q']/ds['weights_q'].sum()
    ds['weights_none'] = np.ones(88)/88


    return ds


def main():
    args = read_input()

    fig, ax = plt.subplots(figsize=(7, 5))
    fig.subplots_adjust(left=0.2, right=0.8, bottom=.22, top=.91)

    xticks = []
    xticklabels = ['equal','equal','performance','performance', '1/N, ic members','1/N, ic members', '1/N, models','1/N, models', 'RMSE','RMSE']
    colors = ['blue','0.65','blue','0.65','blue','0.65','blue','0.65','blue','0.65']
    weights = ['weights_none','weights_none_all','weights_q_norm','weights_q_norm_all','weights_ic','weights_ic_all','weights_mc','weights_mc_all','weights','weights']

    for xx, filename in enumerate(args.filenames):
        if filename == '':
            continue

        ds = read_data(filename, args.path)

        varn = ds.attrs['target']
        xticks.append(xx)
        xticklabels.append(ds.attrs['config'])

        if args.unweighted:
            h1 = boxplot(
                ax, xx,
                mean=ds[varn],
                box=ds[varn],
                whis=ds[varn],  # (ds[varn].min(), ds[varn].max()),
                width=.8,
                color=colors[xx], # sns.xkcd_rgb['greyish'],
                alpha=.3,
            )

        h2 = boxplot(
            ax, xx,
            mean=ds[varn],
            box=ds[varn],
            whis=ds[varn],
            weights=ds[weights[xx]],
            width=.6,
            color=colors[xx],
            alpha=1,
        )

    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels, rotation=30, ha='center')

    try:
        unit = f' ({ds[varn].attrs["units"]})'
    except KeyError:
        unit = ''
    ax.set_ylabel('SAT $\Delta$ ($\degree$C)')
    ax.set_ylim([2, 9])
    #plt.gca().set_aspect('equal', adjustable='box')

    if args.title is not None:
        plt.title(args.title)

    if args.savename is None:
        plt.show()
    else:
        plt.savefig(os.path.join(SAVEPATH, args.savename), dpi=300)


if __name__ == '__main__':
    main()
