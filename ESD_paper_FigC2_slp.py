#################################
# packages
#################################

import xarray as xr
import cartopy.crs as ccrs
import numpy as np
import matplotlib.pyplot as plt
import regionmask
import cdo

#################################
# load things in
#################################

dirT = '/net/h2o/climphys/meranna/files_for_CC_Dyn_Adj/cesm122_members/'
dsP = xr.open_dataset(dirT + 'psl_mon_CESM12-LE_rcp85_00-49_g025.nc',use_cftime = True)

dirTm = '/net/h2o/climphys/meranna/files_for_CC_Dyn_Adj/mpi_ge_members/'
dsPm = xr.open_dataset(dirTm + 'psl_mon_MPI-ESM_rcp85_001-100_g025.nc',use_cftime = True)

dirTn = '/net/h2o/climphys/meranna/files_for_CC_Dyn_Adj/canesm_members/'
dsPn = xr.open_dataset(dirTn + 'psl_mon_CanESM2-LE_rcp85_01-50_g025.nc',use_cftime = True)

dirTc = '/net/h2o/climphys/meranna/files_for_CC_Dyn_Adj/cmip5/temp/'
dsPc = xr.open_dataset(dirTc + 'psl_mon_cmip5_rcp85_g025.nc',use_cftime = True)

dirTobs = '/net/h2o/climphys/meranna/files_for_CC_Dyn_Adj/obs/'
dsPobs = xr.open_dataset(dirTobs + 'psl_mon_ERA-20C_g025.nc',use_cftime = True)
dsPobs1 = xr.open_dataset(dirTobs + 'psl_mon_BEST_g025.nc',use_cftime = True)


#################################
# convert Pa to hPa
#################################

def Pa_to_hPa(ds):
    if ds.psl.units == 'Pa':
        ds['psl_hPa'] = ds.psl/100
    else:
        ds['psl_hPa'] = ds.psl
    return ds

dsP = Pa_to_hPa(dsP)
dsPm = Pa_to_hPa(dsPm)
dsPn = Pa_to_hPa(dsPn)
dsPc = Pa_to_hPa(dsPc)
dsPobs = Pa_to_hPa(dsPobs)
dsPobs1 = Pa_to_hPa(dsPobs1)

#################################
# Use SREX regions or a Custom Mask for EU region
#################################

def mask_region(ds,region,wraplon):
    if region == 'NEU':
        mask = regionmask.defined_regions.srex.mask(ds, wrap_lon=wraplon)
        ds_r = ds.where(mask == 11)
    if region == 'MED':
        mask = regionmask.defined_regions.srex.mask(ds, wrap_lon=wraplon)
        ds_r = ds.where(mask == 13)
    if region == 'DA':
        da_msk = [[-60., 90.], [-60., 25.], [100., 25.], [100., 90.]]
        names = ['Adj. Region']
        abbrevs = ['DA']
        mask = regionmask.Regions([da_msk],names=names, abbrevs=abbrevs, name='DA_Region')
        mask = mask.mask(ds, wrap_lon=wraplon) == 0
        ds_r = ds.where(mask)
    return ds_r

dsP_DA = mask_region(dsP,'DA',True)
dsPm_DA = mask_region(dsPm,'DA',True)
dsPn_DA = mask_region(dsPn,'DA',True)
dsPc_DA = mask_region(dsPc,'DA',True)
dsPobs_DA = mask_region(dsPobs,'DA',True)
dsPobs1_DA = mask_region(dsPobs1,'DA',True)

#################################
# Seasonal-Averages
# Author:  Mario S. Könz <mskoenz@gmx.net>
# to distingush - JJA: time --> year, DJF: time remains
#################################

from xarray_season_select import select_season_and_groupby_season_year

dsP_DA_jja = select_season_and_groupby_season_year(dsP_DA, "JJA").mean('time')
dsPm_DA_jja = select_season_and_groupby_season_year(dsPm_DA, "JJA").mean('time')
dsPn_DA_jja = select_season_and_groupby_season_year(dsPn_DA, "JJA").mean('time')
dsPc_DA_jja = select_season_and_groupby_season_year(dsPc_DA, "JJA").mean('time')
dsPobs_DA_jja = select_season_and_groupby_season_year(dsPobs_DA, "JJA").mean('time')
dsPobs1_DA_jja = select_season_and_groupby_season_year(dsPobs1_DA, "JJA").mean('time')

#################################
# Plot 20-Year CLIM (update needed here)
#################################

def cos_lat_weighted_mean(ds):
  weights = np.cos(np.deg2rad(ds.lat))
  weights.name = "weights"
  ds_weighted = ds.weighted(weights)
  weighted_mean = ds_weighted.mean(('lon', 'lat'))
  return weighted_mean

import utils_python.xarray_lew as ut
from scipy.stats import linregress

def plot_scatter_CLIM(ds_dict, varn, year_start, year_end, fut_year_start, fut_year_end):
    color_dict = dict(
        CESM='r',
        MPI='g',
        CanESM='y',
        ERA20C='k',
        BEST='k',
        CMIP5='b')
    for key in ds_dict.keys():
        if key in ['ERA20C']:
            plt.axvline(
                cos_lat_weighted_mean(
                    ds_dict[key][varn].sel(year=slice(year_start, year_end)).mean('year')).data,
                color=color_dict[key],
                linestyle='solid',
                linewidth=1,
                label=key)
        if key in ['BEST']:
            plt.axvline(
                cos_lat_weighted_mean(
                    ds_dict[key][varn].sel(year=slice(year_start, year_end)).mean('year')).data,
                color=color_dict[key],
                linestyle='dashed',
                linewidth=1,
                label=key)
        if key in ['CESM','MPI','CanESM','CMIP5']:
            plt.scatter(
                cos_lat_weighted_mean(
                    ds_dict[key][varn].sel(year=slice(year_start, year_end)).mean('year')).data,
                cos_lat_weighted_mean(
                    ds_dict[key][varn].sel(year=slice(fut_year_start, fut_year_end)).mean('year')).data,
                s=5,
                color=color_dict[key],
                label=key)


def plot_fit_CLIM(ds_dict, varn, year_start, year_end, fut_year_start, fut_year_end, keys=None, **kwargs):
    if keys is None:
        keys = ['CESM', 'MPI', 'CanESM', 'CMIP5']
    elif isinstance(keys, str):
        keys = [keys]

    xx, yy = [], []
    for key in keys:
        xx += list(cos_lat_weighted_mean(ds_dict[key][varn].sel(year=slice(year_start, year_end)).mean('year')).data)
        yy += list(cos_lat_weighted_mean(ds_dict[key][varn].sel(year=slice(fut_year_start, fut_year_end)).mean('year')).data)

    a, b, r_value, p_value, std_err = linregress(xx, yy)
    f = lambda x: a*x + b
    x = np.array([np.min(xx), np.max(xx)])
    plt.plot(x, f(x), **kwargs)
    return r_value


ds_psl_JJA = dict(
        CESM=dsP_DA_jja,
        MPI=dsPm_DA_jja,
        CanESM=dsPn_DA_jja,
        ERA20C=dsPobs_DA_jja,
        BEST=dsPobs1_DA_jja,
        CMIP5=dsPc_DA_jja)


plot_scatter_CLIM(ds_psl_JJA, 'psl_hPa', 1950, 1969, 2080, 2099)
r_value_cmip = plot_fit_CLIM(ds_psl_JJA, 'psl_hPa', 1950, 1969, 2080, 2099,'CMIP5', color='b')
r_value_all = plot_fit_CLIM(ds_psl_JJA, 'psl_hPa', 1950, 1969, 2080, 2099,color='k')
plt.xlabel('JJA European Sector Mean SLP (hPa), 1950-1969')
plt.ylabel('JJA European Sector Mean SLP (hPa), 2080-2099')
plt.xlim(1009, 1019)
plt.ylim(1008, 1018)
plt.gca().set_aspect('equal', adjustable='box')
plt.figtext(.49, .9, f'R$^2$={r_value_cmip**2:.2f}', fontsize='large', color='b', ha='right')
plt.figtext(.51, .9, f'R$^2$={r_value_all**2:.2f}', fontsize='large', color='k', ha='left')
plt.savefig('ESD_figC2_slp_a.png', dpi=300)
plt.close()

plot_scatter_CLIM(ds_psl_JJA, 'psl_hPa', 1990, 2009, 2080, 2099)
r_value_cmip = plot_fit_CLIM(ds_psl_JJA, 'psl_hPa', 1990, 2009, 2080, 2099, 'CMIP5', color='b')
r_value_all = plot_fit_CLIM(ds_psl_JJA, 'psl_hPa', 1990, 2009, 2080, 2099, color='k')
plt.xlabel('JJA European Sector Mean SLP (hPa), 1990-2009')
plt.ylabel('JJA European Sector Mean SLP (hPa), 2080-2099')
plt.xlim(1009, 1019)
plt.ylim(1008, 1018)
plt.gca().set_aspect('equal', adjustable='box')
plt.figtext(.49, .9, f'R$^2$={r_value_cmip**2:.2f}', fontsize='large', color='b', ha='right')
plt.figtext(.51, .9, f'R$^2$={r_value_all**2:.2f}', fontsize='large', color='k', ha='left')
plt.savefig('ESD_figC2_slp_b.png', dpi=300)
plt.close()


#################################
# Plot 20-Year STD (update needed here)
#################################

def cos_lat_weighted_mean(ds):
  weights = np.cos(np.deg2rad(ds.lat))
  weights.name = "weights"
  ds_weighted = ds.weighted(weights)
  weighted_mean = ds_weighted.mean(('lon', 'lat'))
  return weighted_mean

import utils_python.xarray as ut
from scipy.stats import linregress

def plot_scatter_STD(ds_dict, varn, year_start, year_end, fut_year_start, fut_year_end):
    color_dict = dict(
        CESM='r',
        MPI='g',
        CanESM='y',
        ERA20C='k',
        BEST='k',
        CMIP5='b')
    for key in ds_dict.keys():
        if key in ['ERA20C']:
            plt.axvline(
                cos_lat_weighted_mean(ut.detrend(
                    ds_dict[key][varn].sel(year=slice(year_start, year_end)),time_name='year').std('year')).data,
                color=color_dict[key],
                linestyle='solid',
                linewidth=1,
                label=key)
        if key in ['BEST']:
            plt.axvline(
                cos_lat_weighted_mean(ut.detrend(
                    ds_dict[key][varn].sel(year=slice(year_start, year_end)),time_name='year').std('year')).data,
                color=color_dict[key],
                linestyle='dashed',
                linewidth=1,
                label=key)
        if key in ['CESM','MPI','CanESM','CMIP5']:
            plt.scatter(
                cos_lat_weighted_mean(ut.detrend(
                    ds_dict[key][varn].sel(year=slice(year_start, year_end)),time_name='year').std('year')).data,
                cos_lat_weighted_mean(ut.detrend(
                    ds_dict[key][varn].sel(year=slice(fut_year_start, fut_year_end)),time_name='year').std('year')).data,
                s=5,
                color=color_dict[key],
                label=key)


def plot_fit_STD(ds_dict, varn, year_start, year_end, fut_year_start, fut_year_end, keys=None, **kwargs):
    if keys is None:
        keys = ['CESM', 'MPI', 'CanESM', 'CMIP5']
    elif isinstance(keys, str):
        keys = [keys]

    xx, yy = [], []
    for key in keys:
        xx += list(cos_lat_weighted_mean(ut.detrend(ds_dict[key][varn].sel(year=slice(year_start, year_end)),time_name='year').std('year')).data)
        yy += list(cos_lat_weighted_mean(ut.detrend(ds_dict[key][varn].sel(year=slice(fut_year_start, fut_year_end)),time_name='year').std('year')).data)

    a, b, r_value, p_value, std_err = linregress(xx, yy)
    f = lambda x: a*x + b
    x = np.array([np.min(xx), np.max(xx)])
    plt.plot(x, f(x), **kwargs)
    return r_value

ds_psl_JJA = dict(
        CESM=dsP_DA_jja,
        MPI=dsPm_DA_jja,
        CanESM=dsPn_DA_jja,
        ERA20C=dsPobs_DA_jja,
        BEST=dsPobs1_DA_jja,
        CMIP5=dsPc_DA_jja)


plot_scatter_STD(ds_psl_JJA, 'psl_hPa', 1950, 1969, 2080, 2099)
r_value_cmip = plot_fit_STD(ds_psl_JJA, 'psl_hPa', 1950, 1969, 2080, 2099,'CMIP5', color='b')
r_value_all = plot_fit_STD(ds_psl_JJA, 'psl_hPa', 1950, 1969, 2080, 2099,color='k')
plt.xlabel('JJA European Sector SLP Std. Dev. (hPa), 1950-1969')
plt.ylabel('JJA European Sector SLP Std. Dev. (hPa), 2080-2099')
plt.xlim(0.8, 1.8)
plt.ylim(0.8, 1.8)
plt.gca().set_aspect('equal', adjustable='box')
plt.figtext(.49, .9, f'R$^2$={r_value_cmip**2:.2f}', fontsize='large', color='b', ha='right')
plt.figtext(.51, .9, f'R$^2$={r_value_all**2:.2f}', fontsize='large', color='k', ha='left')
plt.savefig('ESD_figC2_slp_c.png', dpi=300)
plt.close()

plot_scatter_STD(ds_psl_JJA, 'psl_hPa', 1990, 2009, 2080, 2099)
r_value_cmip = plot_fit_STD(ds_psl_JJA, 'psl_hPa', 1990, 2009, 2080, 2099, 'CMIP5', color='b')
r_value_all = plot_fit_STD(ds_psl_JJA, 'psl_hPa', 1990, 2009, 2080, 2099, color='k')
plt.xlabel('JJA European Sector SLP Std. Dev. (hPa), 1990-2009')
plt.ylabel('JJA European Sector SLP Std. Dev. (hPa), 2080-2099')
plt.xlim(0.8, 1.8)
plt.ylim(0.8, 1.8)
plt.gca().set_aspect('equal', adjustable='box')
plt.figtext(.49, .9, f'R$^2$={r_value_cmip**2:.2f}', fontsize='large', color='b', ha='right')
plt.figtext(.51, .9, f'R$^2$={r_value_all**2:.2f}', fontsize='large', color='k', ha='left')
plt.savefig('ESD_figC2_slp_d.png', dpi=300)
plt.close()
