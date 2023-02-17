# make average salinity files for the whole dang season

import xarray as xr
import datetime as dt
import numpy as np

path_files = "/ocean/mdunphy/CIOPSW-BC12/"
path_save = "/ocean/rbeutel/CIOPSBC12_hourly/seasonal_salinity/"

# set running dates:
startdayW16 = [dt.datetime(2016,10,2)+dt.timedelta(days=i) for i in range(int(34*7))]
folders = [dt.datetime(2016,10,2)+dt.timedelta(days=7*(i+1)) for i in range(int(34))]
folders = np.repeat(folders,7)
salW16 = xr.open_mfdataset([str(path_files)+'{:%Y%m%d}00/BC12_1d_grid_T_{:%Y%m%d}_{:%Y%m%d}.nc'.format(folders[i],startdayW16[i],startdayW16[i]) for i in range(len(startdayW16))]).so.mean(dim="time_counter")

startdayS17 = [dt.datetime(2017,5,7)+dt.timedelta(days=i) for i in range(int(21*7))]
folders = [dt.datetime(2017,5,7)+dt.timedelta(days=7*(i+1)) for i in range(int(21))]
folders = np.repeat(folders,7)
salS17 = xr.open_mfdataset([str(path_files)+'{:%Y%m%d}00/BC12_1d_grid_T_{:%Y%m%d}_{:%Y%m%d}.nc'.format(folders[i],startdayS17[i],startdayS17[i]) for i in range(len(startdayS17))]).so.mean(dim="time_counter")

startdayW17 = [dt.datetime(2017,10,1)+dt.timedelta(days=i) for i in range(int(26*7))]
folders = [dt.datetime(2017,10,1)+dt.timedelta(days=7*(i+1)) for i in range(int(26))]
folders = np.repeat(folders,7)
salW17 = xr.open_mfdataset([str(path_files)+'{:%Y%m%d}00/BC12_1d_grid_T_{:%Y%m%d}_{:%Y%m%d}.nc'.format(folders[i],startdayW17[i],startdayW17[i]) for i in range(len(startdayW17))]).so.mean(dim="time_counter")

# And save!
    
encoding={
      "so": {"zlib": True, "complevel": 4, "_FillValue": 0}
}
    
salW16.to_netcdf(str(path_save)+'S_avg_{:%Y%m%d}_{:%Y%m%d}.nc'.format(startdayW16[0],startdayW16[-1]), encoding=encoding)
salW17.to_netcdf(str(path_save)+'S_avg_{:%Y%m%d}_{:%Y%m%d}.nc'.format(startdayW17[0],startdayW17[-1]), encoding=encoding)
salS17.to_netcdf(str(path_save)+'S_avg_{:%Y%m%d}_{:%Y%m%d}.nc'.format(startdayS17[0],startdayS17[-1]), encoding=encoding)