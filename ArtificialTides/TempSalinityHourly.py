# # Making 3D Temp and Salinity hourly in BC12 output

import xarray as xr
from pathlib import Path
import numpy as np
import datetime as dt
import gsw

# set runing dates:
n = 18
startday = [dt.datetime(2016,5,1)+dt.timedelta(days=i) for i in range(int(n*7))]
folders = [dt.datetime(2016,5,1)+dt.timedelta(days=7*(i+1)) for i in range(int(n))]
folders = np.repeat(folders,7)
print(startday[0])
print(startday[-1])

for i in range(len(startday)):
# all within a for loop so you dont have to restart the code every day for 4 years

    # dates for each run
    date_list = [startday[i],startday[i+1]]
#     print(date_list)
    folderday = [folders[i], folders[i+1]]
#     print(folderday)

    # In[3]:
    #load U
    path = Path("/ocean/mdunphy/CIOPSW-BC12/")

    drop_vars = (
        "deptht_bounds","time_counter_bounds","time_instant_bounds",
    )

    files = [sorted(path.glob("{:%Y%m%d}00/BC12_1d_grid_T_{:%Y%m%d}_{:%Y%m%d}.nc".format(folderday[i], date_list[i], date_list[i]))) for i in range(len(date_list))]
#     print(files)

    mydata = xr.open_mfdataset(files, drop_variables=drop_vars)
    sal = mydata['so']
    potT = mydata['thetao']
    
    # replace all land values with nan so that math isn't done on them
    sal = sal.where(sal != 0)
    potT = potT.where(potT != 0)

    #convert potential temp to in-situ temp using gws toolbox
    CT = gsw.CT_from_pt(sal,potT)
    T = gsw.t_from_CT(sal,CT,potT.deptht)

    # interpolate + resample uc_d to get it in an hourly format
    sal_interp = sal.resample(time_counter="1H", loffset=dt.timedelta(hours=1)).interpolate("linear")
    T_interp = T.resample(time_counter="1H", loffset=dt.timedelta(hours=1)).interpolate("linear")

    # trim the extra hour
    sal_new = sal_interp.isel(time_counter = np.arange(0,24,1))
    T_new = T_interp.isel(time_counter = np.arange(0,24,1))
    
    #naming format same as salishseacast
    sal_new = sal_new.rename('vosaline')
    T_new = T_new.rename('votemper')
    
    # order the variables the way you need em
    T_new = T_new.transpose('time_counter','deptht','y','x')
    sal_new = sal_new.transpose('time_counter','deptht','y','x')

    # And save!
    path = '/ocean/rbeutel/data/'
    
    encoding={
          "vosaline": {"zlib": True, "complevel": 4, "_FillValue": 0}
    }
    
    sal_new.to_netcdf(str(path)+'{:%Y%m}/S_new_{:%Y%m%d}.nc'.format(date_list[1],date_list[1]), encoding=encoding)
    print(str(path)+'{:%Y%m}/S_new_{:%Y%m%d}.nc'.format(date_list[1],date_list[1]))
    
    encoding={
          "votemper": {"zlib": True, "complevel": 4, "_FillValue": 0}
    }
    
    T_new.to_netcdf(str(path)+'{:%Y%m}/T_new_{:%Y%m%d}.nc'.format(date_list[1],date_list[1]), encoding=encoding)
    print(str(path)+'{:%Y%m}/T_new_{:%Y%m%d}.nc'.format(date_list[1],date_list[1]))
