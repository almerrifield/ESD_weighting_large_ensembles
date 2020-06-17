#################################
# Useful packages
#################################

import xarray as xr
import cartopy.crs as ccrs
import numpy as np
import matplotlib.pyplot as plt
import regionmask
import cdo

import datetime
from scipy import signal, stats


def linear_detrend(ds,axis='time'):
	""" remove linear trend from data """
	def detrend(data):
		""" function to apply """
		if np.any(np.isnan(data)):
			return data * np.nan
		assert len(data.shape) == 1 # insures one specified dimension
		return signal.detrend(data)

	return xr.apply_ufunc(
		detrend,ds,
		input_core_dims=[[axis]],
		output_core_dims=[[axis]],
		vectorize=True,keep_attrs=True)


def linear_trend(ds, axis='time'):
	"""Calculate the linear least-squares trend"""
	def trend_data(data):
		if np.any(np.isnan(data)):
			return np.nan
		return stats.linregress(np.arange(len(data)), data)[0]

	return xr.apply_ufunc(
		trend_data, ds,
		input_core_dims=[[axis]],
		vectorize=True,
		keep_attrs=False)
