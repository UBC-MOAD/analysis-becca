# # Making 3D Temp and Salinity hourly in BC12 output

import xarray as xr
from pathlib import Path
import numpy as np
import datetime as dt
import gsw

# set runing dates:
startday = [dt.datetime(2015,11,22)+dt.timedelta(days=i) for i in range(int(1498))]
# print(len(startday))
folders = [dt.datetime(2015,11,29)+dt.timedelta(days=7*i) for i in range(int(214))]
folders = np.repeat(folders,7)
# print(len(folders))

for i in range(866,len(startday)-1):
# all within a for loop so you dont have to restart the code every day for 4 years

    # dates for each run
    date_list = [startday[i],startday[i+1]]
    folderday = [folders[i], folders[i+1]]

    # In[3]:
    #load U
    path = Path("/ocean/mdunphy/CIOPSW-BC12/")

    drop_vars = (
        "deptht_bounds","time_counter_bounds","time_instant_bounds",
    )

    files = [sorted(path.glob("{:%Y%m%d}00/BC12_1d_grid_T_{:%Y%m%d}_{:%Y%m%d}.nc".format(folderday[i], date_list[i], date_list[i]))) for i in range(len(date_list))]

    mydata = xr.open_mfdataset(files, drop_variables=drop_vars)
    sal = mydata['so']
    potT = mydata['thetao']

    #convert potential temp to in-situ temp using gws toolbox
    CT = gsw.CT_from_pt(sal,potT)
    T = gsw.t_from_CT(sal,CT,potT.deptht)
    T = T.where(T>1,0) #the conversions mess with all the land values, convert these back to 0.. the region doesnt go lower than 1 so this is fine

    # interpolate + resample uc_d to get it in an hourly format
    sal_interp = sal.resample(time_counter="1H", loffset=dt.timedelta(hours=-23)).interpolate("linear")
    T_interp = T.resample(time_counter="1H", loffset=dt.timedelta(hours=-23)).interpolate("linear")

    # trim the extra hour
    sal_new = sal_interp.isel(time_counter = np.arange(0,24,1))
    T_new = T_interp.isel(time_counter = np.arange(0,24,1))
    
    #naming format same as salishseacast
    sal_new = sal_new.rename('vosaline')
    T_new = T_new.rename('votemper')

    # And save!
    path = '/ocean/rbeutel/data/'
    sal_new.to_netcdf(str(path)+'{:%Y%m}/S_new_{:%Y%m%d}.nc'.format(date_list[0],date_list[0]))
    print(str(path)+'{:%Y%m}/S_new_{:%Y%m%d}.nc'.format(date_list[0],date_list[0]))
    T_new.to_netcdf(str(path)+'{:%Y%m}/T_new_{:%Y%m%d}.nc'.format(date_list[0],date_list[0]))
    print(str(path)+'{:%Y%m}/S_new_{:%Y%m%d}.nc'.format(date_list[0],date_list[0]))