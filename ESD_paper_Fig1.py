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

def mask_ocean(ds):
    land = regionmask.defined_regions.natural_earth.land_110
    mask_lnd = land.mask(ds,wrap_lon=True)
    ds = ds.where(mask_lnd==0)
    return ds

dsT = mask_ocean(dsT)
dsTm = mask_ocean(dsTm)
dsTn = mask_ocean(dsTn)
dsTc = mask_ocean(dsTc)
dsTobs = mask_ocean(dsTobs)
dsTobs1 = mask_ocean(dsTobs1)

#################################
# check mask
#################################

# cmap = plt.get_cmap('YlOrRd')
# ax = plt.axes(projection=ccrs.PlateCarree())
# dsTobs.isel(time=-4).tas_C.plot.pcolormesh(ax=ax, transform=ccrs.PlateCarree(), cmap=cmap,
#                                             cbar_kwargs={'orientation': 'horizontal',
#                                             'label': 'SAT (K)',
#                                             'pad': .1})
# ax.coastlines();
# plt.title('Land Mask Check')
# plt.savefig('land_mask_test_obs.png', dpi=300)
# plt.close()

# cmap = plt.get_cmap('YlOrRd')
# ax = plt.axes(projection=ccrs.PlateCarree())
# dsTn.isel(time=-1,member=0).tas_C.plot.pcolormesh(ax=ax, transform=ccrs.PlateCarree(), cmap=cmap,
#                                             cbar_kwargs={'orientation': 'horizontal',
#                                             'label': 'SAT (C)',
#                                             'pad': .1})
# ax.coastlines();
# plt.title('Land Mask Check')
# plt.savefig('land_mask_test.png', dpi=300)
# plt.close()


#################################
# Use SREX regions or a Custom Mask for EU region
#################################

def mask_region(ds,region):
    if region == 'NEU':
        mask = regionmask.defined_regions.srex.mask(ds, wrap_lon=True)
        ds_r = ds.where(mask == 11)
    if region == 'MED':
        mask = regionmask.defined_regions.srex.mask(ds, wrap_lon=True)
        ds_r = ds.where(mask == 13)
    return ds_r

dsT_NEU = mask_region(dsT,'NEU')
dsT_MED = mask_region(dsT,'MED')

dsTm_NEU = mask_region(dsTm,'NEU')
dsTm_MED = mask_region(dsTm,'MED')

dsTn_NEU = mask_region(dsTn,'NEU')
dsTn_MED = mask_region(dsTn,'MED')

dsTc_NEU = mask_region(dsTc,'NEU')
dsTc_MED = mask_region(dsTc,'MED')

dsTobs_NEU = mask_region(dsTobs,'NEU')
dsTobs_MED = mask_region(dsTobs,'MED')

dsTobs1_NEU = mask_region(dsTobs1,'NEU')
dsTobs1_MED = mask_region(dsTobs1,'MED')

#################################
# Seasonal-Averages
# Author:  Mario S. Könz <mskoenz@gmx.net>
# to distingush: JJA time --> year, DJF time remains
#################################

from xarray_season_select import select_season_and_groupby_season_year
# importlib.reload(xarray_season_select)

dsT_NEU_djf = select_season_and_groupby_season_year(dsT_NEU, "DJF").mean('time')
dsT_MED_jja = select_season_and_groupby_season_year(dsT_MED, "JJA").mean('time')

dsTm_NEU_djf = select_season_and_groupby_season_year(dsTm_NEU, "DJF").mean('time')
dsTm_MED_jja = select_season_and_groupby_season_year(dsTm_MED, "JJA").mean('time')

dsTn_NEU_djf = select_season_and_groupby_season_year(dsTn_NEU, "DJF").mean('time')
dsTn_MED_jja = select_season_and_groupby_season_year(dsTn_MED, "JJA").mean('time')

dsTc_NEU_djf = select_season_and_groupby_season_year(dsTc_NEU, "DJF").mean('time')
dsTc_MED_jja = select_season_and_groupby_season_year(dsTc_MED, "JJA").mean('time')

dsTobs_NEU_djf = select_season_and_groupby_season_year(dsTobs_NEU, "DJF").mean('time')
dsTobs_MED_jja = select_season_and_groupby_season_year(dsTobs_MED, "JJA").mean('time')

dsTobs1_NEU_djf = select_season_and_groupby_season_year(dsTobs1_NEU, "DJF").mean('time')
dsTobs1_MED_jja = select_season_and_groupby_season_year(dsTobs1_MED, "JJA").mean('time')

# cmap = plt.get_cmap('YlOrRd')
# ax = plt.axes(projection=ccrs.PlateCarree())
# dsTobs_MED_jja.isel(year=0).tas_C.plot.pcolormesh(ax=ax, transform=ccrs.PlateCarree(), cmap=cmap,
#                                             cbar_kwargs={'orientation': 'horizontal',
#                                             'label': 'SAT (C)',
#                                             'pad': .1})
# ax.coastlines();
# plt.title('Land Mask Check')
# plt.savefig('land_mask_test_obs.png', dpi=300)
# plt.close()

#################################
# Area-Averages
#################################

# importlib.reload(cos_lat_weighted_mean)

def cos_lat_weighted_mean(ds):
  weights = np.cos(np.deg2rad(ds.lat))
  weights.name = "weights"
  ds_weighted = ds.weighted(weights)
  weighted_mean = ds_weighted.mean(('lon', 'lat'))
  return weighted_mean

dsT_NEU_djf_ts = cos_lat_weighted_mean(dsT_NEU_djf.sel(year=slice('1950','2099')))
dsTm_NEU_djf_ts = cos_lat_weighted_mean(dsTm_NEU_djf.sel(year=slice('1950','2099')))
dsTn_NEU_djf_ts = cos_lat_weighted_mean(dsTn_NEU_djf.sel(year=slice('1950','2099')))
dsTc_NEU_djf_ts = cos_lat_weighted_mean(dsTc_NEU_djf.sel(year=slice('1950','2099')))
dsTobs_NEU_djf_ts = cos_lat_weighted_mean(dsTobs_NEU_djf.sel(year=slice('1950','2009')))
dsTobs1_NEU_djf_ts = cos_lat_weighted_mean(dsTobs1_NEU_djf.sel(year=slice('1950','2009')))

dsT_MED_jja_ts = cos_lat_weighted_mean(dsT_MED_jja.sel(year=slice('1950','2099')))
dsTm_MED_jja_ts = cos_lat_weighted_mean(dsTm_MED_jja.sel(year=slice('1950','2099')))
dsTn_MED_jja_ts = cos_lat_weighted_mean(dsTn_MED_jja.sel(year=slice('1950','2099')))
dsTc_MED_jja_ts = cos_lat_weighted_mean(dsTc_MED_jja.sel(year=slice('1950','2099')))

dsTobs_MED_jja_ts = cos_lat_weighted_mean(dsTobs_MED_jja.sel(year=slice('1950','2009')))
dsTobs1_MED_jja_ts = cos_lat_weighted_mean(dsTobs1_MED_jja.sel(year=slice('1950','2009')))

#################################
# Figure 1: panel a (legend created externally)
#################################

dsTobs_NEU_djf_cat = xr.concat([dsTobs_NEU_djf_ts, dsTobs1_NEU_djf_ts.drop(['clim','tas_raw'])], dim='member')

xc = dsTc_NEU_djf_ts.year
# yupc = dsTc_NEU_djf_ts.tas_C.mean('member') +  dsTc_NEU_djf_ts.tas_C.std('member')
# ydnc = dsTc_NEU_djf_ts.tas_C.mean('member') -  dsTc_NEU_djf_ts.tas_C.std('member')
yupc = np.quantile(dsTc_NEU_djf_ts.tas_C,0.95,1)
ydnc = np.quantile(dsTc_NEU_djf_ts.tas_C,0.05,1)

x = dsT_NEU_djf_ts.year
# yup = dsT_NEU_djf_ts.tas_C.mean('member') +  dsT_NEU_djf_ts.tas_C.std('member')
# ydn = dsT_NEU_djf_ts.tas_C.mean('member') -  dsT_NEU_djf_ts.tas_C.std('member')
yup = np.quantile(dsT_NEU_djf_ts.tas_C,0.95,1)
ydn = np.quantile(dsT_NEU_djf_ts.tas_C,0.05,1)

xm = dsTm_NEU_djf_ts.year
# yupm = dsTm_NEU_djf_ts.tas_C.mean('member') +  dsTm_NEU_djf_ts.tas_C.std('member')
# ydnm = dsTm_NEU_djf_ts.tas_C.mean('member') -  dsTm_NEU_djf_ts.tas_C.std('member')
yupm = np.quantile(dsTm_NEU_djf_ts.tas_C,0.95,1)
ydnm = np.quantile(dsTm_NEU_djf_ts.tas_C,0.05,1)

xn = dsTn_NEU_djf_ts.year
# yupn = dsTn_NEU_djf_ts.tas_C.mean('member') +  dsTn_NEU_djf_ts.tas_C.std('member')
# ydnn = dsTn_NEU_djf_ts.tas_C.mean('member') -  dsTn_NEU_djf_ts.tas_C.std('member')
yupn = np.quantile(dsTn_NEU_djf_ts.tas_C,0.95,1)
ydnn = np.quantile(dsTn_NEU_djf_ts.tas_C,0.05,1)

plt.fill_between(xc,ydnc,yupc, alpha=0.1, color="b",label="CMIP5 (N=88)")
plt.plot(dsTc_NEU_djf_ts.year,dsTc_NEU_djf_ts.tas_C.mean('member'),linewidth=2, color="w")
plt.plot(dsTc_NEU_djf_ts.year,dsTc_NEU_djf_ts.tas_C.mean('member'),alpha=0.8, color="b")
plt.fill_between(x,ydn,yup, alpha=0.1, color="r",label="CESM1.2.2 (N=50)")
plt.plot(dsT_NEU_djf_ts.year,dsT_NEU_djf_ts.tas_C.mean('member'),linewidth=2, color="w")
plt.plot(dsT_NEU_djf_ts.year,dsT_NEU_djf_ts.tas_C.mean('member'),alpha=0.8, color="r")
plt.fill_between(xm,ydnm,yupm, alpha=0.2, color="g",label="MPI-GE (N=100)")
plt.plot(dsTm_NEU_djf_ts.year,dsTm_NEU_djf_ts.tas_C.mean('member'),linewidth=2, color="w")
plt.plot(dsTm_NEU_djf_ts.year,dsTm_NEU_djf_ts.tas_C.mean('member'),alpha=0.8, color="g")
plt.fill_between(xn,ydnn,yupn, alpha=0.3, color="y",label="CanESM2 (N=50)")
plt.plot(dsTn_NEU_djf_ts.year,dsTn_NEU_djf_ts.tas_C.mean('member'),linewidth=2, color="w")
plt.plot(dsTn_NEU_djf_ts.year,dsTn_NEU_djf_ts.tas_C.mean('member'),alpha=0.8, color="y")
plt.plot(dsTobs_NEU_djf_ts.year,dsTobs_NEU_djf_ts.tas_C,"k",alpha=0.4,linewidth=1,linestyle="solid",label="OBS:ERA-20C")
plt.plot(dsTobs1_NEU_djf_ts.year,dsTobs1_NEU_djf_ts.tas_C,"k",alpha=0.4,linewidth=1,linestyle="dashed",label="OBS:BEST")
plt.plot(dsTobs_NEU_djf_cat.year,dsTobs_NEU_djf_cat.tas_C.mean('member'),linewidth=2, color="w")
plt.plot(dsTobs_NEU_djf_cat.year,dsTobs_NEU_djf_cat.tas_C.mean('member'),alpha=1, color="k",label="MEAN OBS")
plt.xlabel('Year')
plt.xlim([1950,2099])
plt.ylabel('DJF NEU SAT (˚C)')
plt.ylim([-14,3])
plt.savefig('ESD_fig1a', dpi=300)
plt.close()

#################################
# Figure 1: panel b (legend created externally)
#################################
dsTobs_MED_jja_cat = xr.concat([dsTobs_MED_jja_ts, dsTobs1_MED_jja_ts.drop(['clim','tas_raw'])], dim='member')

xc = dsTc_MED_jja_ts.year
# yupc = dsTc_MED_jja_ts.tas_C.mean('member') + dsTc_MED_jja_ts.tas_C.std('member')
# ydnc = dsTc_MED_jja_ts.tas_C.mean('member') - dsTc_MED_jja_ts.tas_C.std('member')
yupc = np.quantile(dsTc_MED_jja_ts.tas_C,0.95,1)
ydnc = np.quantile(dsTc_MED_jja_ts.tas_C,0.05,1)

x = dsT_MED_jja_ts.year
# yup = dsT_MED_jja_ts.tas_C.mean('member') + dsT_MED_jja_ts.tas_C.std('member')
# ydn = dsT_MED_jja_ts.tas_C.mean('member') - dsT_MED_jja_ts.tas_C.std('member')
yup = np.quantile(dsT_MED_jja_ts.tas_C,0.95,1)
ydn = np.quantile(dsT_MED_jja_ts.tas_C,0.05,1)

xm = dsTm_MED_jja_ts.year
# yupm = dsTm_MED_jja_ts.tas_C.mean('member') + dsTm_MED_jja_ts.tas_C.std('member')
# ydnm = dsTm_MED_jja_ts.tas_C.mean('member') - dsTm_MED_jja_ts.tas_C.std('member')
yupm = np.quantile(dsTm_MED_jja_ts.tas_C,0.95,1)
ydnm = np.quantile(dsTm_MED_jja_ts.tas_C,0.05,1)

xn = dsTn_MED_jja_ts.year
# yupn = dsTn_MED_jja_ts.tas_C.mean('member') + dsTn_MED_jja_ts.tas_C.std('member')
# ydnn = dsTn_MED_jja_ts.tas_C.mean('member') - dsTn_MED_jja_ts.tas_C.std('member')
yupn = np.quantile(dsTn_MED_jja_ts.tas_C,0.95,1)
ydnn = np.quantile(dsTn_MED_jja_ts.tas_C,0.05,1)

xo = dsTobs_MED_jja_cat.year

plt.fill_between(xc,ydnc,yupc, alpha=0.1, color="b",label="CMIP5 (N=88)")
plt.plot(dsTc_MED_jja_ts.year,dsTc_MED_jja_ts.tas_C.mean('member'),linewidth=2, color="w")
plt.plot(dsTc_MED_jja_ts.year,dsTc_MED_jja_ts.tas_C.mean('member'),alpha=0.8, color="b")
plt.fill_between(x,ydn,yup, alpha=0.1, color="r",label="CESM1.2.2-LE (N=50)")
plt.plot(dsT_MED_jja_ts.year,dsT_MED_jja_ts.tas_C.mean('member'),linewidth=2, color="w")
plt.plot(dsT_MED_jja_ts.year,dsT_MED_jja_ts.tas_C.mean('member'),alpha=0.8, color="r")
plt.fill_between(xn,ydnn,yupn, alpha=0.3, color="y",label="CanESM2-LE (N=50)")
plt.plot(dsTn_MED_jja_ts.year,dsTn_MED_jja_ts.tas_C.mean('member'),linewidth=2, color="w")
plt.plot(dsTn_MED_jja_ts.year,dsTn_MED_jja_ts.tas_C.mean('member'),alpha=0.8, color="y")
plt.fill_between(xm,ydnm,yupm, alpha=0.2, color="g",label="MPI-GE (N=100)")
plt.plot(dsTm_MED_jja_ts.year,dsTm_MED_jja_ts.tas_C.mean('member'),linewidth=2, color="w")
plt.plot(dsTm_MED_jja_ts.year,dsTm_MED_jja_ts.tas_C.mean('member'),alpha=0.8, color="g")
plt.plot(dsTobs_MED_jja_ts.year,dsTobs_MED_jja_ts.tas_C,"k",alpha=0.4,linewidth=1,label="OBS:ERA-20C")
plt.plot(dsTobs1_MED_jja_ts.year,dsTobs1_MED_jja_ts.tas_C,"k",alpha=0.4,linewidth=1,linestyle="dashed",label="OBS:BEST")
plt.plot(dsTobs_MED_jja_cat.year,dsTobs_MED_jja_cat.tas_C.mean('member'),linewidth=2, color="w")
plt.plot(dsTobs_MED_jja_cat.year,dsTobs_MED_jja_cat.tas_C.mean('member'),alpha=1, color="k",label="MEAN OBS")
plt.xlabel('Year')
plt.xlim([1950,2099])
plt.ylim([21,35])
plt.ylabel('JJA MED SAT (˚C)')
plt.savefig('ESD_fig1b', dpi=300)
plt.close()
