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
dsTadj = xr.open_dataset(dirT + 'tasAdj_mon_CESM12-LE_rcp85_00-49_g025_anom_cat.nc',use_cftime = True)

dirTm = '/net/h2o/climphys/meranna/files_for_CC_Dyn_Adj/mpi_ge_members/'
dsTm = xr.open_dataset(dirTm + 'tas_mon_MPI-ESM_rcp85_001-100_g025.nc',use_cftime = True)
dsTadjm = xr.open_dataset(dirTm + 'tasAdj_mon_MPI-ESM_rcp85_001-100_g025_cat.nc',use_cftime = True)

dirTn = '/net/h2o/climphys/meranna/files_for_CC_Dyn_Adj/canesm_members/'
dsTn = xr.open_dataset(dirTn + 'tas_mon_CanESM2-LE_rcp85_01-50_g025.nc',use_cftime = True)
dsTadjn = xr.open_dataset(dirTn + 'tasAdj_mon_CanESM2-LE_rcp85_01-50_g025_cat.nc',use_cftime = True)

dirTc = '/net/h2o/climphys/meranna/files_for_CC_Dyn_Adj/cmip5/temp/'
dsTc = xr.open_dataset(dirTc + 'tas_mon_cmip5_rcp85_g025.nc',use_cftime = True)
dsTadjc = xr.open_dataset(dirTc + 'tasAdj_mon_cmip5_rcp85_g025_cat.nc',use_cftime = True)

dirTobs = '/net/h2o/climphys/meranna/files_for_CC_Dyn_Adj/obs/'
dsTobs = xr.open_dataset(dirTobs + 'tas_mon_ERA-20C_g025.nc',use_cftime = True)
dsTobs_adj = xr.open_dataset(dirTobs + 'tasAdj_mon_ERA-20C_g025.nc',use_cftime = True)
dsTobs1 = xr.open_dataset(dirTobs + 'tas_mon_BEST_g025.nc',use_cftime = True)
dsTobs1_adj = xr.open_dataset(dirTobs + 'tasAdj_mon_BEST_g025.nc',use_cftime = True)

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

dsTadj = mask_ocean(dsTadj,False)
dsTadjm = mask_ocean(dsTadjm,False)
dsTadjn = mask_ocean(dsTadjn,False)
dsTadjc = mask_ocean(dsTadjc,False)
dsTobs_adj = mask_ocean(dsTobs_adj,False)
dsTobs1_adj = mask_ocean(dsTobs1_adj,False)

# cmap = plt.get_cmap('YlOrRd')
# ax = plt.axes(projection=ccrs.PlateCarree())
# dsTadj.isel(time=0,member=-1).tasAdj.plot.pcolormesh(ax=ax, transform=ccrs.PlateCarree(), cmap=cmap,
#                                             cbar_kwargs={'orientation': 'horizontal',
#                                             'label': 'SAT (C)',
#                                             'pad': .1})
# ax.coastlines();
# plt.title('Land Mask Check')
# plt.savefig('land_mask_test_obs.png', dpi=300)
# plt.close()

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
    if region == 'EUR':
        da_msk = [[-10., 76.25], [-10., 30.], [39., 30.], [39., 76.25]]
        # numbers = [0]
        names = ['Europe']
        abbrevs = ['EUR']
        da_mask = regionmask.Regions([da_msk],names=names, abbrevs=abbrevs, name='Europe')
        mask2 = da_mask.mask(ds,wrap_lon=wraplon)
        ds_r = ds.where(mask2==0)
    return ds_r

dsT_EU = mask_region(dsT,'EUR',True)
dsTm_EU = mask_region(dsTm,'EUR',True)
dsTn_EU = mask_region(dsTn,'EUR',True)
dsTc_EU = mask_region(dsTc,'EUR',True)
dsTobs_EU = mask_region(dsTobs,'EUR',True)
dsTobs1_EU = mask_region(dsTobs1,'EUR',True)

dsTadj_EU = mask_region(dsTadj,'EUR',False)
dsTadjm_EU = mask_region(dsTadjm,'EUR',False)
dsTadjn_EU = mask_region(dsTadjn,'EUR',False)
dsTadjc_EU = mask_region(dsTadjc,'EUR',False)
dsTobs_adj_EU = mask_region(dsTobs_adj,'EUR',False)
dsTobs1_adj_EU = mask_region(dsTobs1_adj,'EUR',False)

# cmap = plt.get_cmap('YlOrRd')
# ax = plt.axes(projection=ccrs.PlateCarree())
# dsTadjc_EU.isel(time=0,member=-1).tasAdj.plot.pcolormesh(ax=ax, transform=ccrs.PlateCarree(), cmap=cmap,
#                                             cbar_kwargs={'orientation': 'horizontal',
#                                             'label': 'SAT (C)',
#                                             'pad': .1})
# ax.coastlines();
# plt.title('Land Mask Check')
# plt.savefig('land_mask_test_obs.png', dpi=300)
# plt.close()

#################################
# Seasonal-Averages
# Author:  Mario S. Könz <mskoenz@gmx.net>
# to distingush - JJA: time --> year, DJF: time remains
#################################

from xarray_season_select import select_season_and_groupby_season_year

dsT_EU_djf = select_season_and_groupby_season_year(dsT_EU, "DJF").mean("time")
dsT_EU_jja = select_season_and_groupby_season_year(dsT_EU, "JJA").mean("time")

dsTadj_EU_djf = select_season_and_groupby_season_year(dsTadj_EU, "DJF").mean("time")
dsTadj_EU_jja = select_season_and_groupby_season_year(dsTadj_EU, "JJA").mean("time")

dsTm_EU_djf = select_season_and_groupby_season_year(dsTm_EU, "DJF").mean("time")
dsTm_EU_jja = select_season_and_groupby_season_year(dsTm_EU, "JJA").mean("time")

dsTadjm_EU_djf = select_season_and_groupby_season_year(dsTadjm_EU, "DJF").mean("time")
dsTadjm_EU_jja = select_season_and_groupby_season_year(dsTadjm_EU, "JJA").mean("time")

dsTn_EU_djf = select_season_and_groupby_season_year(dsTn_EU, "DJF").mean("time")
dsTn_EU_jja = select_season_and_groupby_season_year(dsTn_EU, "JJA").mean("time")

dsTadjn_EU_djf = select_season_and_groupby_season_year(dsTadjn_EU, "DJF").mean("time")
dsTadjn_EU_jja = select_season_and_groupby_season_year(dsTadjn_EU, "JJA").mean("time")

dsTc_EU_djf = select_season_and_groupby_season_year(dsTc_EU, "DJF").mean("time")
dsTc_EU_jja = select_season_and_groupby_season_year(dsTc_EU, "JJA").mean("time")

dsTadjc_EU_djf = select_season_and_groupby_season_year(dsTadjc_EU, "DJF").mean("time")
dsTadjc_EU_jja = select_season_and_groupby_season_year(dsTadjc_EU, "JJA").mean("time")

dsTobs_EU_djf = select_season_and_groupby_season_year(dsTobs_EU, "DJF").mean("time")
dsTobs_EU_jja = select_season_and_groupby_season_year(dsTobs_EU, "JJA").mean("time")

dsTobs_adj_EU_djf = select_season_and_groupby_season_year(dsTobs_adj_EU, "DJF").mean("time")
dsTobs_adj_EU_jja = select_season_and_groupby_season_year(dsTobs_adj_EU, "JJA").mean("time")

dsTobs1_EU_djf = select_season_and_groupby_season_year(dsTobs1_EU, "DJF").mean("time")
dsTobs1_EU_jja = select_season_and_groupby_season_year(dsTobs1_EU, "JJA").mean("time")

dsTobs1_adj_EU_djf = select_season_and_groupby_season_year(dsTobs1_adj_EU, "DJF").mean("time")
dsTobs1_adj_EU_jja = select_season_and_groupby_season_year(dsTobs1_adj_EU, "JJA").mean("time")

#################################
# Plot 50-Year Trends (!!update needed here!!)
#################################

def cos_lat_weighted_mean(ds):
  weights = np.cos(np.deg2rad(ds.lat))
  weights.name = "weights"
  ds_weighted = ds.weighted(weights)
  weighted_mean = ds_weighted.mean(('lon', 'lat'))
  return weighted_mean

import linear_regression_functions as ut
from scipy.stats import linregress


def plot_scatter(ds_dict, varn):
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
                    ut.linear_trend(ds_dict[key][varn].sel(year=slice(1960, 2009)),axis='year')).data * 50,
                color=color_dict[key],
                linestyle='solid',
                linewidth=1,
                label=key)
        if key in ['BEST']:
            plt.axvline(
                cos_lat_weighted_mean(
                    ut.linear_trend(ds_dict[key][varn].sel(year=slice(1960, 2009)),axis='year')).data * 50,
                color=color_dict[key],
                linestyle='dashed',
                linewidth=1,
                label=key)
        if key in ['CESM','MPI','CanESM','CMIP5']:
            plt.scatter(
                cos_lat_weighted_mean(ut.linear_trend(ds_dict[key][varn].sel(year=slice(1960, 2009)),axis='year')).data * 50,
                cos_lat_weighted_mean(ut.linear_trend(ds_dict[key][varn].sel(year=slice(2050, 2099)),axis='year')).data * 50,
                s=10,
                color=color_dict[key],
                label=key,
            )


def plot_fit(ds_dict, varn, keys=None, **kwargs):
    if keys is None:
        keys = ['CESM', 'MPI', 'CanESM', 'CMIP5']
    elif isinstance(keys, str):
        keys = [keys]

    xx, yy = [], []
    for key in keys:
        xx += list(cos_lat_weighted_mean(ut.linear_trend(ds_dict[key][varn].sel(year=slice(1960, 2009)),axis='year')).data * 50)
        yy += list(cos_lat_weighted_mean(ut.linear_trend(ds_dict[key][varn].sel(year=slice(2050, 2099)),axis='year')).data * 50)

    a, b, r_value, p_value, std_err = linregress(xx, yy)
    f = lambda x: a*x + b
    x = np.array([np.min(xx), np.max(xx)])
    plt.plot(x, f(x), **kwargs)
    return r_value



ds_tas_DJF = dict(
        CESM=dsT_EU_djf,
        MPI=dsTm_EU_djf,
        CanESM=dsTn_EU_djf,
        ERA20C=dsTobs_EU_djf,
        BEST=dsTobs1_EU_djf,
        CMIP5=dsTc_EU_djf)

ds_tasAdj_DJF = dict(
        CESM=dsTadj_EU_djf,
        MPI=dsTadjm_EU_djf,
        CanESM=dsTadjn_EU_djf,
        ERA20C=dsTobs_adj_EU_djf,
        BEST=dsTobs1_adj_EU_djf,
        CMIP5=dsTadjc_EU_djf)


plot_scatter(ds_tas_DJF, 'tas_C')
r_value_cmip = plot_fit(ds_tas_DJF, 'tas_C', 'CMIP5', color='b')
r_value_all = plot_fit(ds_tas_DJF, 'tas_C', color='k')
plt.xlabel('DJF EUR SAT Trend (˚C/50 years), 1960-2009')
plt.ylabel('DJF EUR SAT Trend (˚C/50 years), 2050-2099')
plt.xlim(-1, 4)
plt.ylim(-1, 6)
plt.axvline(0, color='k', linestyle='dotted', linewidth=0.5)
plt.axhline(0, color='k', linestyle='dotted', linewidth=0.5)
plt.gca().set_aspect('equal', adjustable='box')
plt.figtext(.49, .9, f'R$^2$={r_value_cmip**2:.2f}', fontsize='large', color='b', ha='right')
plt.figtext(.51, .9, f'R$^2$={r_value_all**2:.2f}', fontsize='large', color='k', ha='left')
plt.savefig('ESD_fig2a.png', dpi=300)
plt.close()

plot_scatter(ds_tasAdj_DJF, 'tasAdj')
r_value_cmip = plot_fit(ds_tasAdj_DJF, 'tasAdj','CMIP5', color='b')
r_value_all = plot_fit(ds_tasAdj_DJF, 'tasAdj',color='k')
plt.xlabel('DJF EUR ERT SAT Trend (˚C/50 years), 1960-2009')
plt.ylabel('DJF EUR ERT SAT Trend (˚C/50 years), 2050-2099')
plt.xlim(-1, 4)
plt.ylim(-1, 6)
plt.axvline(0, color='k', linestyle='dotted', linewidth=0.5)
plt.axhline(0, color='k', linestyle='dotted', linewidth=0.5)
plt.gca().set_aspect('equal', adjustable='box')
plt.figtext(.49, .9, f'R$^2$={r_value_cmip**2:.2f}', fontsize='large', color='b', ha='right')
plt.figtext(.51, .9, f'R$^2$={r_value_all**2:.2f}', fontsize='large', color='k', ha='left')
plt.savefig('ESD_fig2b.png', dpi=300)
plt.close()



ds_tas_JJA = dict(
        CESM=dsT_EU_jja,
        MPI=dsTm_EU_jja,
        CanESM=dsTn_EU_jja,
        ERA20C=dsTobs_EU_jja,
        BEST=dsTobs1_EU_jja,
        CMIP5=dsTc_EU_jja)

ds_tasAdj_JJA = dict(
        CESM=dsTadj_EU_jja,
        MPI=dsTadjm_EU_jja,
        CanESM=dsTadjn_EU_jja,
        ERA20C=dsTobs_adj_EU_jja,
        BEST=dsTobs1_adj_EU_jja,
        CMIP5=dsTadjc_EU_jja)


plot_scatter(ds_tas_JJA, 'tas_C')
r_value_cmip = plot_fit(ds_tas_JJA, 'tas_C', 'CMIP5', color='b')
r_value_all = plot_fit(ds_tas_JJA, 'tas_C', color='k')
plt.xlabel('JJA EUR SAT Trend (˚C/50 years), 1960-2009')
plt.ylabel('JJA EUR SAT Trend (˚C/50 years), 2050-2099')
plt.xlim(-1, 4)
plt.ylim(-1, 6)
plt.axvline(0, color='k', linestyle='dotted', linewidth=0.5)
plt.axhline(0, color='k', linestyle='dotted', linewidth=0.5)
plt.gca().set_aspect('equal', adjustable='box')
plt.figtext(.49, .9, f'R$^2$={r_value_cmip**2:.2f}', fontsize='large', color='b', ha='right')
plt.figtext(.51, .9, f'R$^2$={r_value_all**2:.2f}', fontsize='large', color='k', ha='left')
plt.savefig('ESD_fig2c.png', dpi=300)
plt.close()

plot_scatter(ds_tasAdj_JJA, 'tasAdj')
r_value_cmip = plot_fit(ds_tasAdj_JJA, 'tasAdj','CMIP5', color='b')
r_value_all = plot_fit(ds_tasAdj_JJA, 'tasAdj',color='k')
plt.xlabel('JJA EUR ERT SAT Trend (˚C/50 years), 1960-2009')
plt.ylabel('JJA EUR ERT SAT Trend (˚C/50 years), 2050-2099')
plt.xlim(-1, 4)
plt.ylim(-1, 6)
plt.axvline(0, color='k', linestyle='dotted', linewidth=0.5)
plt.axhline(0, color='k', linestyle='dotted', linewidth=0.5)
plt.gca().set_aspect('equal', adjustable='box')
plt.figtext(.49, .9, f'R$^2$={r_value_cmip**2:.2f}', fontsize='large', color='b', ha='right')
plt.figtext(.51, .9, f'R$^2$={r_value_all**2:.2f}', fontsize='large', color='k', ha='left')
plt.savefig('ESD_fig2d.png', dpi=300)
plt.close()
