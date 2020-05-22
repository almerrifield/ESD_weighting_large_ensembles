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
dsP = xr.open_dataset(dirT + 'psl_mon_CESM12-LE_rcp85_00-49_g025.nc',use_cftime = True)

dirTm = '/net/h2o/climphys/meranna/files_for_CC_Dyn_Adj/mpi_ge_members/'
dsTm = xr.open_dataset(dirTm + 'tas_mon_MPI-ESM_rcp85_001-100_g025.nc',use_cftime = True)
dsPm = xr.open_dataset(dirTm + 'psl_mon_MPI-ESM_rcp85_001-100_g025.nc',use_cftime = True)

dirTn = '/net/h2o/climphys/meranna/files_for_CC_Dyn_Adj/canesm_members/'
dsTn = xr.open_dataset(dirTn + 'tas_mon_CanESM2-LE_rcp85_01-50_g025.nc',use_cftime = True)
dsPn = xr.open_dataset(dirTn + 'psl_mon_CanESM2-LE_rcp85_01-50_g025.nc',use_cftime = True)

dirTc = '/net/h2o/climphys/meranna/files_for_CC_Dyn_Adj/cmip5/temp/'
dsTc = xr.open_dataset(dirTc + 'tas_mon_cmip5_rcp85_g025.nc',use_cftime = True)
dsPc = xr.open_dataset(dirTc + 'psl_mon_cmip5_rcp85_g025.nc',use_cftime = True)

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
    if region == 'NH':
        da_msk = [[-180., 90.], [-180., 0.], [180., 0.], [180., 90.]]
        names = ['Adj. Region']
        abbrevs = ['DA']
        mask = regionmask.Regions([da_msk],names=names, abbrevs=abbrevs, name='NH_Region')
        mask = mask.mask(ds, wrap_lon=wraplon) == 0
        ds_r = ds.where(mask)
    return ds_r

dsP_NH = mask_region(dsP,'NH',True)
dsPm_NH = mask_region(dsPm,'NH',True)
dsPn_NH = mask_region(dsPn,'NH',True)
dsPc_NH = mask_region(dsPc,'NH',True)

# cmap = plt.get_cmap('YlOrRd')
# ax = plt.axes(projection=ccrs.PlateCarree())
# dsP_NH.isel(time=-1,member=0).psl_hPa.plot.pcolormesh(ax=ax, transform=ccrs.PlateCarree(), cmap=cmap,
#                                             cbar_kwargs={'orientation': 'horizontal',
#                                             'label': 'SLP (hPa)',
#                                             'pad': .1})
# ax.coastlines();
# plt.title('Land Mask Check')
# plt.savefig('land_mask_test.png', dpi=300)
# plt.close()

#################################
# Time-Aggregation: Annual
#################################

def annual_average(ds):
    ds_out = ds.groupby('time.year').mean('time')
    return ds_out

dsT_ANN = annual_average(dsT)
dsTm_ANN = annual_average(dsTm)
dsTn_ANN = annual_average(dsTn)
dsTc_ANN = annual_average(dsTc)

dsP_NH_ANN = annual_average(dsP_NH)
dsPm_NH_ANN = annual_average(dsPm_NH)
dsPn_NH_ANN = annual_average(dsPn_NH)
dsPc_NH_ANN = annual_average(dsPc_NH)

#################################
# Plotting
#################################

def cos_lat_weighted_mean(ds):
  weights = np.cos(np.deg2rad(ds.lat))
  weights.name = "weights"
  ds_weighted = ds.weighted(weights)
  weighted_mean = ds_weighted.mean(('lon', 'lat'))
  return weighted_mean

def plot_scatter(ds_dict1,varn1,ds_dict2,varn2):
    color_dict = dict(
        CESM='r',
        MPI='g',
        CanESM='y',
        CMIP5='b')

    for key in ds_dict1.keys():
            plt.scatter(
                cos_lat_weighted_mean(
                    ds_dict1[key][varn1].sel(year=slice(1950, 2009)).mean('year')).data,
                cos_lat_weighted_mean(
                    ds_dict2[key][varn2].sel(year=slice(1950, 2009)).mean('year')).data,
                s=5,
                color=color_dict[key],
                label=key)

ds_tas = dict(
        CESM=dsT_ANN,
        MPI=dsTm_ANN,
        CanESM=dsTn_ANN,
        CMIP5=dsTc_ANN)

ds_psl = dict(
        CESM=dsP_NH_ANN,
        MPI=dsPm_NH_ANN,
        CanESM=dsPn_NH_ANN,
        CMIP5=dsPc_NH_ANN)

plot_scatter(ds_psl, 'psl_hPa',ds_tas,'tas_C')
plt.xlabel('Annual Mean Northern Hemisphere SLP (hPa), 1950-2009')
plt.ylabel('Annual Mean Global Land SAT ($\degree$C), 1950-2009')
plt.savefig('ESD_fig5.png', dpi=300)
plt.close()
