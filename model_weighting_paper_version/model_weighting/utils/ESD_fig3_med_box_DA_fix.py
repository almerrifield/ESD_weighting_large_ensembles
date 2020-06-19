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
    ds['weights_ic_all'] = [0.0173645245053380, 0.0203569712595225, 0.000386623999800020, 0.00175385476816864, 0.00175385476816864, 0.00175385476816864, 0.00175385476816864, 0.00175385476816864, 0.00175385476816864, 0.0149877194350907, 0.00534133867487110, 0.00534133867487110, 0.00534133867487110, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.000150857296767254, 0.0111500739284606, 0.0173243017086583, 0.0254478697381020, 0.00660444157632810, 0.00660444157632810, 0.00660444157632810, 0.00660444157632810, 0.00660444157632810, 0.000789256962921410, 0.000789256962921410, 0.000789256962921410, 0.000789256962921410, 0.000789256962921410, 0.000789256962921410, 0.000789256962921410, 0.000789256962921410, 0.000789256962921410, 0.000789256962921410, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 7.81218954017782e-05, 0.000657141865269707, 0.000657141865269707, 0.000657141865269707, 0.000657141865269707, 0.000657141865269707, 0.00975293795439704, 0.00975293795439704, 0.00975293795439704, 0.00975293795439704, 0.00975293795439704, 0.0205828413526940, 0.000777714869769811, 0.000777714869769811, 0.000777714869769811, 0.00936654632862810, 0.00150766657362369, 0.000257303183878334, 0.0101293911099127, 0.00634237391463581, 0.0136328786842070, 0.0108653021490644, 0.00634237391463581, 0.0108653021490644, 0.0149985193452363, 0.00590432468557565, 0.0149412285347167, 0.00592697102067078, 0.00590432468557565, 0.00592697102067078, 0.0163786490319080, 0.0266036508275658, 0.00591207366188384, 0.00591207366188384, 0.00591207366188384, 0.00591207366188384, 0.00113746952131555, 0.00113746952131555, 0.00113746952131555, 0.00113746952131555, 0.0130571399126652, 0.00237580290248032, 0.00438123681715656, 0.00592322682981926, 0.00239185126410343, 0.00239185126410343, 0.00239185126410343, 0.00646243993086172, 0.00646243993086172, 0.00646243993086172, 0.0189238866516331, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 7.69126865451268e-05, 0.00807208094414755, 0.00748992176859595, 0.0293103380988083, 0.0200810011517527, 0.00433793743220949, 0.0119406692125648, 0.0108898285395486]
    ds['weights_ic_all'] = ds['weights_ic_all']/ds['weights_ic_all'].sum()

    ds['weights_mc_all'] = [0.00943037394121511, 0.00943037394121511, 0.000386623999800020, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.000139817161411427, 0.00599136059724678, 0.00599136059724678, 0.00599136059724678, 0.00660444157632810, 0.00660444157632810, 0.00660444157632810, 0.00660444157632810, 0.00660444157632810, 0.000789256962921410, 0.000789256962921410, 0.000789256962921410, 0.000789256962921410, 0.000789256962921410, 0.000789256962921410, 0.000789256962921410, 0.000789256962921410, 0.000789256962921410, 0.000789256962921410, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 6.99944744251862e-05, 0.00975293795439704, 0.00975293795439704, 0.00975293795439704, 0.00975293795439704, 0.00975293795439704, 0.0205828413526940, 0.000777714869769811, 0.000777714869769811, 0.000777714869769811, 0.00123683512068112, 0.00123683512068112, 0.00123683512068112, 0.00117956878301291, 0.00117956878301291, 0.00117956878301291, 0.00117956878301291, 0.00117956878301291, 0.00117956878301291, 0.00117956878301291, 0.00117956878301291, 0.00117956878301291, 0.00117956878301291, 0.00117956878301291, 0.00117956878301291, 0.00382154106804487, 0.00382154106804487, 0.00382154106804487, 0.00382154106804487, 0.00382154106804487, 0.00382154106804487, 0.000934234865449844, 0.000934234865449844, 0.000934234865449844, 0.000934234865449844, 0.000934234865449844, 0.000934234865449844, 0.00127324500095627, 0.00127324500095627, 0.00127324500095627, 0.00127324500095627, 0.00127324500095627, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 7.82371219934039e-05, 0.00389050067818587, 0.00389050067818587, 0.0123478348126402, 0.0123478348126402, 0.00406965166119357, 0.00406965166119357, 0.0108898285395486]
    ds['weights_mc_all'] = ds['weights_mc_all']/ds['weights_mc_all'].sum()

    ds['weights_q_norm_all'] = ds['weights_q']/ds['weights_q'].sum()
    ds['weights_none_all'] = np.ones(288)/288

    ds['weights_ic'] = [0.0174718770203454, 0.0205613863535534, 0.000412558961838729, 0.00181139252421693, 0.00181139252421693, 0.00181139252421693, 0.00181139252421693, 0.00181139252421693, 0.00181139252421693, 0.0152283273055325, 0.00550018051180113, 0.00550018051180113, 0.00550018051180113, 0.0112088503665843, 0.0176843558125814, 0.0260773435951545, 0.00665451624081320, 0.00665451624081320, 0.00665451624081320, 0.00665451624081320, 0.00665451624081320, 0.000801662917011347, 0.000801662917011347, 0.000801662917011347, 0.000801662917011347, 0.000801662917011347, 0.000801662917011347, 0.000801662917011347, 0.000801662917011347, 0.000801662917011347, 0.000801662917011347, 0.000656719634757235, 0.000656719634757235, 0.000656719634757235, 0.000656719634757235, 0.000656719634757235, 0.00899614071852682, 0.00899614071852682, 0.00899614071852682, 0.00899614071852682, 0.00899614071852682, 0.00899614071852682, 0.000781379408867015, 0.000781379408867015, 0.000781379408867015, 0.00972743527533972, 0.00163266888094045, 0.000301337471280230, 0.0105003069513024, 0.00651606493379239, 0.0139223845374914, 0.0111001101725287, 0.00651606493379239, 0.0111001101725287, 0.0152014391759463, 0.00607054179337712, 0.0154118121987557, 0.00616597283229141, 0.00607054179337712, 0.00616597283229141, 0.0164066070898815, 0.0265744029546302, 0.00595595344081403, 0.00595595344081403, 0.00595595344081403, 0.00595595344081403, 0.00112909631444207, 0.00112909631444207, 0.00112909631444207, 0.00112909631444207, 0.0130100732010882, 0.00237656896992004, 0.00442777010500610, 0.00598900416553578, 0.00239570964391535, 0.00239570964391535, 0.00239570964391535, 0.00662762048171131, 0.00662762048171131, 0.00662762048171131, 0.0190513881039158, 0.00808316367371425, 0.00755266722505562, 0.0294790423288641, 0.0205224435384989, 0.00453997994364382, 0.0125080241634372, 0.0109314202361373]
    ds['weights_ic'] = ds['weights_ic']/ds['weights_ic'].sum()

    ds['weights_mc'] = [0.00950831584347471, 0.00950831584347471, 0.000412558961838729, 0.00129940082783552, 0.00129940082783552, 0.00129940082783552, 0.00129940082783552, 0.00129940082783552, 0.00129940082783552, 0.00129940082783552, 0.00129940082783552, 0.00129940082783552, 0.00129940082783552, 0.00610783886381335, 0.00610783886381335, 0.00610783886381335, 0.00665451624081320, 0.00665451624081320, 0.00665451624081320, 0.00665451624081320, 0.00665451624081320, 0.000801662917011347, 0.000801662917011347, 0.000801662917011347, 0.000801662917011347, 0.000801662917011347, 0.000801662917011347, 0.000801662917011347, 0.000801662917011347, 0.000801662917011347, 0.000801662917011347, 0.000656719634757235, 0.000656719634757235, 0.000656719634757235, 0.000656719634757235, 0.000656719634757235, 0.00899614071852682, 0.00899614071852682, 0.00899614071852682, 0.00899614071852682, 0.00899614071852682, 0.00899614071852682, 0.000781379408867015, 0.000781379408867015, 0.000781379408867015, 0.00129571573639560, 0.00129571573639560, 0.00129571573639560, 0.00121143542910732, 0.00121143542910732, 0.00121143542910732, 0.00121143542910732, 0.00121143542910732, 0.00121143542910732, 0.00121143542910732, 0.00121143542910732, 0.00121143542910732, 0.00121143542910732, 0.00121143542910732, 0.00121143542910732, 0.00384100736382045, 0.00384100736382045, 0.00384100736382045, 0.00384100736382045, 0.00384100736382045, 0.00384100736382045, 0.000929227311168926, 0.000929227311168926, 0.000929227311168926, 0.000929227311168926, 0.000929227311168926, 0.000929227311168926, 0.00127912644263120, 0.00127912644263120, 0.00127912644263120, 0.00127912644263120, 0.00127912644263120, 0.00491874827745735, 0.00491874827745735, 0.00491874827745735, 0.00491874827745735, 0.00390895772469247, 0.00390895772469247, 0.0125003714668407, 0.0125003714668407, 0.00426200102677025, 0.00426200102677025, 0.0109314202361373]
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
