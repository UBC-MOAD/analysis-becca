# make average salinity files for the whole dang season

import xarray as xr
import datetime as dt

path = "/ocean/rbeutel/data/"

datesW16 = [dt.datetime(2016,10,1)+dt.timedelta(days=i) for i in range(30*6+2)]
salW16 = xr.open_mfdataset([str(path)+'{:%Y%m}/S_new_{:%Y%m%d}.nc'.format(d,d) for d in datesW16]).mean(dim="time_counter")

datesW17 = [dt.datetime(2017,10,1)+dt.timedelta(days=i) for i in range(30*6+2)]
salW17 = xr.open_mfdataset([str(path)+'{:%Y%m}/S_new_{:%Y%m%d}.nc'.format(d,d) for d in datesW17]).mean(dim="time_counter")

datesS17 = [dt.datetime(2017,4,1)+dt.timedelta(days=i) for i in range(30*6+3)]
salS17 = xr.open_mfdataset([str(path)+'{:%Y%m}/S_new_{:%Y%m%d}.nc'.format(d,d) for d in datesS17]).mean(dim="time_counter")

# And save!
    
encoding={
      "vosaline": {"zlib": True, "complevel": 4, "_FillValue": 0}
}
    
salW16.to_netcdf(str(path)+'S_avg_{:%Y%m%d}_{:%Y%m%d}.nc'.format(datesW16[0],datesW16[-1]), encoding=encoding)
salW17.to_netcdf(str(path)+'S_avg_{:%Y%m%d}_{:%Y%m%d}.nc'.format(datesW17[0],datesW17[-1]), encoding=encoding)
salS17.to_netcdf(str(path)+'S_avg_{:%Y%m%d}_{:%Y%m%d}.nc'.format(datesS17[0],datesS17[-1]), encoding=encoding)