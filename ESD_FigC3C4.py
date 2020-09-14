import xarray as xr
import cartopy.crs as ccrs
import numpy as np
import matplotlib.pyplot as plt
import regionmask
import cdo

#################################
# load things in
#################################

dirTn = '/net/h2o/climphys/meranna/files_for_CC_Dyn_Adj/canesm_members/'
dsPn = xr.open_dataset(dirTn + 'psl_mon_CanESM2-LE_rcp85_01-50_g025.nc',use_cftime = True)

dirTc = '/net/h2o/climphys/meranna/files_for_CC_Dyn_Adj/cmip5/temp/'
dsPc = xr.open_dataset(dirTc + 'psl_mon_cmip5_rcp85_g025.nc',use_cftime = True)

dirTobs = '/net/h2o/climphys/meranna/files_for_CC_Dyn_Adj/obs/'
dsPobs = xr.open_dataset(dirTobs + 'psl_mon_ERA-20C_g025.nc',use_cftime = True)

#################################
# convert Pa to hPa
#################################

def Pa_to_hPa(ds):
    if ds.psl.units == 'Pa':
        ds['psl_hPa'] = ds.psl/100
    else:
        ds['psl_hPa'] = ds.psl
    return ds

dsPn = Pa_to_hPa(dsPn)
dsPc = Pa_to_hPa(dsPc)
dsPobs = Pa_to_hPa(dsPobs)

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

dsPn_DA = mask_region(dsPn,'DA',True)
dsPc_DA = mask_region(dsPc,'DA',True)
dsPobs_DA = mask_region(dsPobs,'DA',True)

#################################
# Seasonal-Averages
# Author:  Mario S. Könz <mskoenz@gmx.net>
# to distingush - JJA: time --> year, DJF: time remains
#################################

from xarray_season_select import select_season_and_groupby_season_year

dsPn_DA_djf = select_season_and_groupby_season_year(dsPn_DA, "DJF").mean('time')
dsPc_DA_djf = select_season_and_groupby_season_year(dsPc_DA, "DJF").mean('time')
dsPobs_DA_djf = select_season_and_groupby_season_year(dsPobs_DA, "DJF").mean('time')

dsPn_DA_jja = select_season_and_groupby_season_year(dsPn_DA, "JJA").mean('time')
dsPc_DA_jja = select_season_and_groupby_season_year(dsPc_DA, "JJA").mean('time')
dsPobs_DA_jja = select_season_and_groupby_season_year(dsPobs_DA, "JJA").mean('time')


#################################
# Plot 20-Year CLIM
#################################

def cos_lat_weighted_mean(ds):
  weights = np.cos(np.deg2rad(ds.lat))
  weights.name = "weights"
  ds_weighted = ds.weighted(weights)
  weighted_mean = ds_weighted.mean(('lon', 'lat'))
  return weighted_mean

dsPn_DA_djf_clm1 = dsPn_DA_djf.sel(year=slice('1950','1969')).mean('year')
dsPn_DA_djf_clm2 = dsPn_DA_djf.sel(year=slice('1990','2009')).mean('year')
dsPn_DA_djf_clm_lt = dsPn_DA_djf.sel(year=slice('2080','2099')).mean('year')

dsPc_DA_djf_clm1 = dsPc_DA_djf.sel(year=slice('1950','1969')).mean('year')
dsPc_DA_djf_clm2 = dsPc_DA_djf.sel(year=slice('1990','2009')).mean('year')
dsPc_DA_djf_clm_lt = dsPc_DA_djf.sel(year=slice('2080','2099')).mean('year')

ind_djf_min = cos_lat_weighted_mean(dsPc_DA_djf_clm2).argmin()
ind_djf_max = cos_lat_weighted_mean(dsPc_DA_djf_clm2).argmax()

dsPobs_DA_djf_clm1 = dsPobs_DA_djf.sel(year=slice('1950','1969')).mean('year')
dsPobs_DA_djf_clm2 = dsPobs_DA_djf.sel(year=slice('1990','2009')).mean('year')

cmap = plt.get_cmap('PRGn')
proj=ccrs.LambertConformal(central_longitude=15)
ax = plt.subplot(331, projection=proj)
dsPobs_DA_djf_clm1.psl_hPa.plot.pcolormesh(ax=ax,transform=ccrs.PlateCarree(),vmin=990,vmax=1040,cmap=cmap,add_colorbar=False)
ax.coastlines();
ax.set_extent([-60, 100, 20, 90], crs=ccrs.PlateCarree());

ax = plt.subplot(334, projection=proj)
dsPobs_DA_djf_clm2.psl_hPa.plot.pcolormesh(ax=ax,transform=ccrs.PlateCarree(),vmin=990,vmax=1040,cmap=cmap,add_colorbar=False)
ax.coastlines();
ax.set_extent([-60, 100, 20, 90], crs=ccrs.PlateCarree());

ax = plt.subplot(332, projection=proj)
dsPc_DA_djf_clm1.isel(member=ind_djf_min.psl.data).psl_hPa.plot.pcolormesh(ax=ax,transform=ccrs.PlateCarree(),vmin=990,vmax=1040,cmap=cmap,add_colorbar=False)
ax.coastlines();
ax.set_extent([-60, 100, 20, 90], crs=ccrs.PlateCarree());

ax = plt.subplot(335, projection=proj)
dsPc_DA_djf_clm2.isel(member=ind_djf_min.psl.data).psl_hPa.plot.pcolormesh(ax=ax,transform=ccrs.PlateCarree(),vmin=990,vmax=1040,cmap=cmap,add_colorbar=False)
ax.coastlines();
ax.set_extent([-60, 100, 20, 90], crs=ccrs.PlateCarree());

ax = plt.subplot(338, projection=proj)
dsPc_DA_djf_clm_lt.isel(member=ind_djf_min.psl.data).psl_hPa.plot.pcolormesh(ax=ax,transform=ccrs.PlateCarree(),vmin=990,vmax=1040,cmap=cmap,add_colorbar=False)
ax.coastlines();
ax.set_extent([-60, 100, 20, 90], crs=ccrs.PlateCarree());

ax = plt.subplot(333, projection=proj)
dsPc_DA_djf_clm1.isel(member=ind_djf_max.psl.data).psl_hPa.plot.pcolormesh(ax=ax,transform=ccrs.PlateCarree(),vmin=990,vmax=1040,cmap=cmap,add_colorbar=False)
ax.coastlines();
ax.set_extent([-60, 100, 20, 90], crs=ccrs.PlateCarree());

ax = plt.subplot(336, projection=proj)
dsPc_DA_djf_clm2.isel(member=ind_djf_max.psl.data).psl_hPa.plot.pcolormesh(ax=ax,transform=ccrs.PlateCarree(),vmin=990,vmax=1040,cmap=cmap,add_colorbar=False)
ax.coastlines();
ax.set_extent([-60, 100, 20, 90], crs=ccrs.PlateCarree());

ax = plt.subplot(339, projection=proj)
dsPc_DA_djf_clm_lt.isel(member=ind_djf_max.psl.data).psl_hPa.plot.pcolormesh(ax=ax,transform=ccrs.PlateCarree(),vmin=990,vmax=1040,cmap=cmap,add_colorbar=False)
ax.coastlines();
ax.set_extent([-60, 100, 20, 90], crs=ccrs.PlateCarree());

plt.savefig('clim_djf.png', dpi=300)
plt.close()

dsPn_DA_jja_clm1 = dsPn_DA_jja.sel(year=slice('1950','1969')).mean('year')
dsPn_DA_jja_clm2 = dsPn_DA_jja.sel(year=slice('1990','2009')).mean('year')
dsPn_DA_jja_clm_lt = dsPn_DA_jja.sel(year=slice('2080','2099')).mean('year')

ind_jja_min = cos_lat_weighted_mean(dsPn_DA_jja_clm2).argmin()

dsPc_DA_jja_clm1 = dsPc_DA_jja.sel(year=slice('1950','1969')).mean('year')
dsPc_DA_jja_clm2 = dsPc_DA_jja.sel(year=slice('1990','2009')).mean('year')
dsPc_DA_jja_clm_lt = dsPc_DA_jja.sel(year=slice('2080','2099')).mean('year')

ind_jja_max = cos_lat_weighted_mean(dsPc_DA_jja_clm2).argmax()

dsPobs_DA_jja_clm1 = dsPobs_DA_jja.sel(year=slice('1950','1969')).mean('year')
dsPobs_DA_jja_clm2 = dsPobs_DA_jja.sel(year=slice('1990','2009')).mean('year')

cmap = plt.get_cmap('PRGn')
proj=ccrs.LambertConformal(central_longitude=15)
ax = plt.subplot(331, projection=proj)
dsPobs_DA_jja_clm1.psl_hPa.plot.pcolormesh(ax=ax,transform=ccrs.PlateCarree(),vmin=990,vmax=1040,cmap=cmap,add_colorbar=False)
ax.coastlines();
ax.set_extent([-60, 100, 20, 90], crs=ccrs.PlateCarree());

ax = plt.subplot(334, projection=proj)
dsPobs_DA_jja_clm2.psl_hPa.plot.pcolormesh(ax=ax,transform=ccrs.PlateCarree(),vmin=990,vmax=1040,cmap=cmap,add_colorbar=False)
ax.coastlines();
ax.set_extent([-60, 100, 20, 90], crs=ccrs.PlateCarree());

ax = plt.subplot(332, projection=proj)
dsPn_DA_jja_clm1.isel(member=ind_jja_min.psl.data).psl_hPa.plot.pcolormesh(ax=ax,transform=ccrs.PlateCarree(),vmin=990,vmax=1040,cmap=cmap,add_colorbar=False)
ax.coastlines();
ax.set_extent([-60, 100, 20, 90], crs=ccrs.PlateCarree());

ax = plt.subplot(335, projection=proj)
dsPn_DA_jja_clm2.isel(member=ind_jja_min.psl.data).psl_hPa.plot.pcolormesh(ax=ax,transform=ccrs.PlateCarree(),vmin=990,vmax=1040,cmap=cmap,add_colorbar=False)
ax.coastlines();
ax.set_extent([-60, 100, 20, 90], crs=ccrs.PlateCarree());

ax = plt.subplot(338, projection=proj)
dsPn_DA_jja_clm_lt.isel(member=ind_jja_min.psl.data).psl_hPa.plot.pcolormesh(ax=ax,transform=ccrs.PlateCarree(),vmin=990,vmax=1040,cmap=cmap,add_colorbar=False)
ax.coastlines();
ax.set_extent([-60, 100, 20, 90], crs=ccrs.PlateCarree());

ax = plt.subplot(333, projection=proj)
dsPc_DA_jja_clm1.isel(member=ind_jja_max.psl.data).psl_hPa.plot.pcolormesh(ax=ax,transform=ccrs.PlateCarree(),vmin=990,vmax=1040,cmap=cmap,add_colorbar=False)
ax.coastlines();
ax.set_extent([-60, 100, 20, 90], crs=ccrs.PlateCarree());

ax = plt.subplot(336, projection=proj)
dsPc_DA_jja_clm2.isel(member=ind_jja_max.psl.data).psl_hPa.plot.pcolormesh(ax=ax,transform=ccrs.PlateCarree(),vmin=990,vmax=1040,cmap=cmap,add_colorbar=False)
ax.coastlines();
ax.set_extent([-60, 100, 20, 90], crs=ccrs.PlateCarree());

ax = plt.subplot(339, projection=proj)
dsPc_DA_jja_clm_lt.isel(member=ind_jja_max.psl.data).psl_hPa.plot.pcolormesh(ax=ax,transform=ccrs.PlateCarree(),vmin=990,vmax=1040,cmap=cmap,add_colorbar=False)
ax.coastlines();
ax.set_extent([-60, 100, 20, 90], crs=ccrs.PlateCarree());

plt.savefig('clim_jja.png', dpi=300)
plt.close()
