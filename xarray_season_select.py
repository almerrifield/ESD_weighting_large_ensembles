#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author:  Mario S. Könz <mskoenz@gmx.net>
# Adapted by Anna Merrifield

import datetime
import numpy as np
import xarray as xr


def groupby_season_year_wrong(season_select):
    # wrong for DJF, ok for all other seasons
    return season_select.groupby("time.year")

def groupby_season_year(season_select, season):
    ###############################
    # Pre-condition: dataset with a season in it; season = DJF, JJA, MAM, or SON
    # Post-condition: groupby season/seasons present and year, with DJF in the JF year
    ###############################
    if season not in ['DJF', 'MAM', 'JJA',
                      'SON']:  # checks that season is valid
        raise NotImplementedError(season)
    times = season_select.coords['time'].data
    if season != 'DJF':
        return season_select.groupby("time.year")
    # If it's DJF...
    def get_season_label(t):
        if t.month == 12:
            return t.year + 1
        return t.year
    year_label = [get_season_label(t) for t in times]
    xr_year_label = xr.DataArray(year_label, dims=["time"], name="year")
    return season_select.groupby(xr_year_label)


def select_season_and_groupby_season_year(ds, season):
    season_select = ds.isel(time=ds['time.season'] == season)
    return groupby_season_year(season_select, season)


def main():
    times = xr.cftime_range(start="1999-01-01", end="2005-12-01", freq="MS")
    data = np.array(list(range(len(times))))
    ds = xr.DataArray(data, dims=['time'], coords=[times])
    print(ds)  #linear increment over months
    res = select_season_and_groupby_season_year(ds, "DJF")
    for x in res:
        print(x)


if __name__ == "__main__":
    main()
