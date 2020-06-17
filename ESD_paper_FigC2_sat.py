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
dsT = xr.open_dataset(dirT + 'tas_mon_CESM12-LE_rcp85_00-49_g025.nc',use_cftime = True)

dirTm = '/net/h2o/climphys/meranna/files_for_CC_Dyn_Adj/mpi_ge_members/'
dsTm = xr.open_dataset(dirTm + 'tas_mon_MPI-ESM_rcp85_001-100_g025.nc',use_cftime = True)

dirTn = '/net/h2o/climphys/meranna/files_for_CC_Dyn_Adj/canesm_members/'
dsTn = xr.open_dataset(dirTn + 'tas_mon_CanESM2-LE_rcp85_01-50_g025.nc',use_cftime = True)

dirTc = '/net/h2o/climphys/meranna/files_for_CC_Dyn_Adj/cmip5/temp/'
dsTc = xr.open_dataset(dirTc + 'tas_mon_cmip5_rcp85_g025.nc',use_cftime = True)

dirTobs = '/net/h2o/climphys/meranna/files_for_CC_Dyn_Adj/obs/'
dsTobs = xr.open_dataset(dirTobs + 'tas_mon_ERA-20C_g025.nc',use_cftime = True)
dsTobs1 = xr.open_dataset(dirTobs + 'tas_mon_BEST_g025.nc',use_cftime = True)

#################################
# Convert BEST from anomalies
#################################

clim = [dsTobs1.clim]*166
clim_obs1 = xr.concat(clim, "month_number")
clim_obs1 = clim_obs1.rename({'month_number': 'time'})

dsTobs1['tas_raw'] = dsTobs1.tas + clim_obs1

#################################
# convert K to C
#################################

def K_to_C(ds):
    if ds.tas.units == 'K':
        ds['tas_C'] = ds.tas - 273.15
    else:
        ds['tas_C'] = ds.tas_raw # BEST case
    return ds

dsT = K_to_C(dsT)
dsTm = K_to_C(dsTm)
dsTn = K_to_C(dsTn)
dsTc = K_to_C(dsTc)
dsTobs = K_to_C(dsTobs)
dsTobs1 = K_to_C(dsTobs1)

#################################
# mask ocean
#################################

def mask_ocean(ds,wraplon):
    land = regionmask.defined_regions.natural_earth.land_110
    mask_lnd = land.mask(ds,wrap_lon=wraplon)
    ds = ds.where(mask_lnd==0)
    return ds

dsT = mask_ocean(dsT,True)
dsTm = mask_ocean(dsTm,True)
dsTn = mask_ocean(dsTn,True)
dsTc = mask_ocean(dsTc,True)
dsTobs = mask_ocean(dsTobs,True)
dsTobs1 = mask_ocean(dsTobs1,True)

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


dsT_MED = mask_region(dsT,'MED',True)
dsTm_MED = mask_region(dsTm,'MED',True)
dsTn_MED = mask_region(dsTn,'MED',True)
dsTc_MED = mask_region(dsTc,'MED',True)
dsTobs_MED = mask_region(dsTobs,'MED',True)
dsTobs1_MED = mask_region(dsTobs1,'MED',True)


#################################
# Seasonal-Averages
# Author:  Mario S. Könz <mskoenz@gmx.net>
# to distingush - JJA: time --> year, DJF: time remains
#################################

from xarray_season_select import select_season_and_groupby_season_year

dsT_MED_jja = select_season_and_groupby_season_year(dsT_MED, "JJA").mean('time')
dsTm_MED_jja = select_season_and_groupby_season_year(dsTm_MED, "JJA").mean('time')
dsTn_MED_jja = select_season_and_groupby_season_year(dsTn_MED, "JJA").mean('time')
dsTc_MED_jja = select_season_and_groupby_season_year(dsTc_MED, "JJA").mean('time')
dsTobs_MED_jja = select_season_and_groupby_season_year(dsTobs_MED, "JJA").mean('time')
dsTobs1_MED_jja = select_season_and_groupby_season_year(dsTobs1_MED, "JJA").mean('time')

#################################
# Plot 20-Year CLIM (update needed here)
#################################

def cos_lat_weighted_mean(ds):
  weights = np.cos(np.deg2rad(ds.lat))
  weights.name = "weights"
  ds_weighted = ds.weighted(weights)
  weighted_mean = ds_weighted.mean(('lon', 'lat'))
  return weighted_mean

import linear_regression_functions as ut
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


ds_tas_JJA = dict(
        CESM=dsT_MED_jja,
        MPI=dsTm_MED_jja,
        CanESM=dsTn_MED_jja,
        ERA20C=dsTobs_MED_jja,
        BEST=dsTobs1_MED_jja,
        CMIP5=dsTc_MED_jja)


plot_scatter_CLIM(ds_tas_JJA, 'tas_C', 1950, 1969, 2080, 2099)
r_value_cmip = plot_fit_CLIM(ds_tas_JJA, 'tas_C', 1950, 1969, 2080, 2099,'CMIP5', color='b')
r_value_all = plot_fit_CLIM(ds_tas_JJA, 'tas_C', 1950, 1969, 2080, 2099,color='k')
plt.xlabel('JJA MED Mean SAT ($\degree$C), 1950-1969')
plt.ylabel('JJA MED Mean SAT ($\degree$C), 2080-2099')
plt.xlim(21, 30)
plt.ylim(26, 35)
plt.gca().set_aspect('equal', adjustable='box')
plt.figtext(.49, .9, f'R$^2$={r_value_cmip**2:.2f}', fontsize='large', color='b', ha='right')
plt.figtext(.51, .9, f'R$^2$={r_value_all**2:.2f}', fontsize='large', color='k', ha='left')
plt.savefig('ESD_figC2_sat_a.png', dpi=300)
plt.close()

plot_scatter_CLIM(ds_tas_JJA, 'tas_C', 1990, 2009, 2080, 2099)
r_value_cmip = plot_fit_CLIM(ds_tas_JJA, 'tas_C', 1990, 2009, 2080, 2099, 'CMIP5', color='b')
r_value_all = plot_fit_CLIM(ds_tas_JJA, 'tas_C', 1990, 2009, 2080, 2099, color='k')
plt.xlabel('JJA MED Mean SAT ($\degree$C), 1990-2009')
plt.ylabel('JJA MED Mean SAT ($\degree$C), 2080-2099')
plt.xlim(21, 30)
plt.ylim(26, 35)
plt.gca().set_aspect('equal', adjustable='box')
plt.figtext(.49, .9, f'R$^2$={r_value_cmip**2:.2f}', fontsize='large', color='b', ha='right')
plt.figtext(.51, .9, f'R$^2$={r_value_all**2:.2f}', fontsize='large', color='k', ha='left')
plt.savefig('ESD_figC2_sat_b.png', dpi=300)
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

import linear_regression_functions as ut
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
                cos_lat_weighted_mean(ut.linear_detrend(
                    ds_dict[key][varn].sel(year=slice(year_start, year_end)),time_name='year').std('year')).data,
                color=color_dict[key],
                linestyle='solid',
                linewidth=1,
                label=key)
        if key in ['BEST']:
            plt.axvline(
                cos_lat_weighted_mean(ut.linear_detrend(
                    ds_dict[key][varn].sel(year=slice(year_start, year_end)),time_name='year').std('year')).data,
                color=color_dict[key],
                linestyle='dashed',
                linewidth=1,
                label=key)
        if key in ['CESM','MPI','CanESM','CMIP5']:
            plt.scatter(
                cos_lat_weighted_mean(ut.linear_detrend(
                    ds_dict[key][varn].sel(year=slice(year_start, year_end)),time_name='year').std('year')).data,
                cos_lat_weighted_mean(ut.linear_detrend(
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
        xx += list(cos_lat_weighted_mean(ut.linear_detrend(ds_dict[key][varn].sel(year=slice(year_start, year_end)),time_name='year').std('year')).data)
        yy += list(cos_lat_weighted_mean(ut.linear_detrend(ds_dict[key][varn].sel(year=slice(fut_year_start, fut_year_end)),time_name='year').std('year')).data)

    a, b, r_value, p_value, std_err = linregress(xx, yy)
    f = lambda x: a*x + b
    x = np.array([np.min(xx), np.max(xx)])
    plt.plot(x, f(x), **kwargs)
    return r_value

ds_tas_JJA = dict(
        CESM=dsT_MED_jja,
        MPI=dsTm_MED_jja,
        CanESM=dsTn_MED_jja,
        ERA20C=dsTobs_MED_jja,
        BEST=dsTobs1_MED_jja,
        CMIP5=dsTc_MED_jja)


plot_scatter_STD(ds_tas_JJA, 'tas_C', 1950, 1969, 2080, 2099)
r_value_cmip = plot_fit_STD(ds_tas_JJA, 'tas_C', 1950, 1969, 2080, 2099,'CMIP5', color='b')
r_value_all = plot_fit_STD(ds_tas_JJA, 'tas_C', 1950, 1969, 2080, 2099,color='k')
plt.xlabel('JJA MED SAT Std. Dev. (˚C), 1950-1969')
plt.ylabel('JJA MED SAT Std. Dev. (˚C), 2080-2099')
plt.xlim(0.4, 1.6)
plt.ylim(0.4, 1.6)
plt.gca().set_aspect('equal', adjustable='box')
plt.figtext(.49, .9, f'R$^2$={r_value_cmip**2:.2f}', fontsize='large', color='b', ha='right')
plt.figtext(.51, .9, f'R$^2$={r_value_all**2:.2f}', fontsize='large', color='k', ha='left')
plt.savefig('ESD_figC2_sat_c.png', dpi=300)
plt.close()

plot_scatter_STD(ds_tas_JJA, 'tas_C', 1990, 2009, 2080, 2099)
r_value_cmip = plot_fit_STD(ds_tas_JJA, 'tas_C', 1990, 2009, 2080, 2099, 'CMIP5', color='b')
r_value_all = plot_fit_STD(ds_tas_JJA, 'tas_C', 1990, 2009, 2080, 2099, color='k')
plt.xlabel('JJA MED SAT Std. Dev. (˚C), 1990-2009')
plt.ylabel('JJA MED SAT Std. Dev. (˚C), 2080-2099')
plt.xlim(0.4, 1.6)
plt.ylim(0.4, 1.6)
plt.gca().set_aspect('equal', adjustable='box')
plt.figtext(.49, .9, f'R$^2$={r_value_cmip**2:.2f}', fontsize='large', color='b', ha='right')
plt.figtext(.51, .9, f'R$^2$={r_value_all**2:.2f}', fontsize='large', color='k', ha='left')
plt.savefig('ESD_figC2_sat_d.png', dpi=300)
plt.close()
