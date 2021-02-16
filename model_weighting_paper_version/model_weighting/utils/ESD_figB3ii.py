#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Time-stamp: <2019-09-12 09:55:13 lukbrunn>

(c) 2019 under a MIT License (https://mit-license.org)

Authors:
- Lukas Brunner || lukas.brunner@env.ethz.ch

Abstract:

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



    return ds


def main():
    args = read_input()

    fig, ax = plt.subplots(figsize=(7, 5))
    fig.subplots_adjust(left=0.2, right=0.8, bottom=.22, top=.91)
    xticks = []
    xticklabels = ['0.05', '0.1', '0.2', '0.3','0.4','0.5','0.6','0.7','0.8']
    colors = ['0.65','0.65','0.65','0.65','0.65','0.65','0.65','0.65','0.65']


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
            weights=ds['weights'],
            width=.6,
            color=colors[xx], #sns.xkcd_rgb['greyish'],
            alpha=1,
            )



    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels, rotation=0, ha='center')

    try:
        unit = f' ({ds[varn].attrs["units"]})'
    except KeyError:
        unit = ''
    ax.set_ylabel('SAT $\Delta$ ($\degree$C)')
    ax.set_xlabel('Independence shape parameter, $\sigma_s$')   #f'{varn}{unit}'
    ax.set_ylim([3, 9])
    plt.gca().set_aspect('equal', adjustable='box')


    if args.title is not None:
        plt.title(args.title)

    if args.savename is None:
        plt.show()
    else:
        plt.savefig(os.path.join(SAVEPATH, args.savename), dpi=300)


if __name__ == '__main__':
    main()
