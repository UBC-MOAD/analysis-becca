# # Adding Back the Tides to CIOPS BC12 output
# An effort to make the daily files more accurate as they are currently lacking the tidal pumping that is so important to the flow of the Salish Sea

import xarray as xr
from pathlib import Path
import numpy as np
import datetime as dt

startday = [dt.datetime(2016,11,22)+dt.timedelta(days=i) for i in range(int(1498))]
folders = [dt.datetime(2015,11,29)+dt.timedelta(days=7*i) for i in range(int(214))]
folders = np.repeat(folders,7)

for i in range(131,len(startday)-1):
# all within a for loop so you dont have to restart the code every day for 4 years

    # dates for each run
    date_list = [startday[i],startday[i+1]]
    folderday = [folders[i], folders[i+1]]

    # In[3]:
    #load U
    path = Path("/ocean/mdunphy/CIOPSW-BC12/")

    drop_vars = (
        "sbu", "tauuo", "time_counter_bounds", "time_instant_bounds", "uos",
    )

    files = [sorted(path.glob("{:%Y%m%d}00/BC12_1h_grid_U_2D_{:%Y%m%d}_{:%Y%m%d}.nc".format(folderday[i], date_list[i], date_list[i]))) for i in range(len(date_list))]

    mydata = xr.open_mfdataset(files, drop_variables=drop_vars)
    ut_h = mydata['ubar']

    # load V
    drop_vars = (
        "sbv", "tauvo", "time_counter_bounds", "time_instant_bounds", "vos",
    )

    files = [sorted(path.glob("{:%Y%m%d}00/BC12_1h_grid_V_2D_{:%Y%m%d}_{:%Y%m%d}.nc".format(folderday[i], date_list[i], date_list[i]))) for i in range(len(date_list))]

    mydata = xr.open_mfdataset(files, drop_variables=drop_vars)
    vt_h = mydata['vbar']

    # 2D hourly outputs are alraedy in barotropic form so you dont need to do any conversion! 

    # Now get the required data from the daily files!

    drop_vars = (
        "depthu_bounds", "nav_lat", "nav_lon", 'time_counter_bounds', 'time_instant',
        'time_instant_bounds', 
    )

    files = [sorted(path.glob("{:%Y%m%d}00/BC12_1d_grid_U_{:%Y%m%d}_{:%Y%m%d}.nc".format(folderday[i], date_list[i], date_list[i]))) for i in range(len(date_list))]

    mydata = xr.open_mfdataset(files, drop_variables=drop_vars)
    u_d = mydata['uo']

    drop_vars = (
        "depthv_bounds", "nav_lat", "nav_lon", 'time_counter_bounds', 'time_instant',
        'time_instant_bounds', 
    )

    files = [sorted(path.glob("{:%Y%m%d}00/BC12_1d_grid_V_{:%Y%m%d}_{:%Y%m%d}.nc".format(folderday[i], date_list[i], date_list[i]))) for i in range(len(date_list))]

    mydata = xr.open_mfdataset(files, drop_variables=drop_vars)
    v_d = mydata['vo']

    # we WILL need to do some conversions here, so get e3t from the mesh_mask file
    mydata = xr.open_dataset("/ocean/mdunphy/CIOPSW-BC12/grid/mesh_mask_Bathymetry_NEP36_714x1020_SRTM30v11_NOAA3sec_WCTSS_JdeFSalSea.nc")
    e3t = mydata['e3t_0']

    # convert e3t to e3u and to e3v
    e3t_xshift = e3t.shift(x=-1,fill_value=0)
    e3u = e3t_xshift+e3t
    e3u = e3u*0.5
    e3u = e3u.rename({'z': 'depthu'})

    e3t_yshift = e3t.shift(y=-1,fill_value=0)
    e3v = e3t_yshift+e3t
    e3v = e3v*0.5
    e3v = e3v.rename({'z': 'depthv'})

    #calcuate barotropic component
    ut_d = (u_d*e3u[0,:,:,:]).sum(dim='depthu')/e3u[0,:,:,:].sum(dim='depthu')

    #subtract from u to get baroclinic component
    uc_d = u_d-ut_d #does this work even though their ut_d lacks the depth dimension?


    # interpolate + resample uc_d to get it in an hourly format
    offset = dt.timedelta(hours=-23) # daily avg timestamp as next day time 0, offset by -23 hours to line up with hourly data
    uc_h_interp = uc_d.resample(time_counter="1H", loffset=offset).interpolate("linear")
    #instead of taking 12hours off each side like when done with SSC, this method takes off 24hours on the end!

    # added together this should give your final u!!!!!
    u_new = ut_h  + uc_h_interp
    u_new = u_new.isel(time_counter = np.arange(0,24,1)) #remove extra hour 
    u_new = u_new.rename('vozocrtx') #name it what you want it named in final netcdf

    # And save!
    # np.save("u_new.npy",u_new)
    path = '/data/rbeutel/analysis/BC12_tidesback/'
    u_new.to_netcdf(str(path)+'U_new_{:%Y%m%d}.nc'.format(date_list[0]))
#     print('u_new_{:%d%b%y}_{:%d%b%y}.nc complete'.format(date_list[0],date_list[-1]))
    # Now for V

    #calcuate bartropic component
    vt_d = (v_d*e3v[0,:,:,:]).sum(dim='depthv')/e3v[0,:,:,:].sum(dim='depthv')


    #subtract from v to get baroclinic component
    vc_d = v_d-vt_d 

    # vc_d.load(scheduler="processes", num_workers=n)

    # interpolate + resample uc_d to get it in an hourly format
    offset = dt.timedelta(hours=-23) # daily avg timestamp as next day time 0, offset by -23 hours to line up with hourly data
    vc_h_interp = vc_d.resample(time_counter="1H", loffset=offset).interpolate("linear")
    #instead of taking 12hours off each side like when done with SSC, this method takes off 24hours on the end!

    v_new = vt_h  + vc_h_interp
    v_new = v_new.isel(time_counter = np.arange(0,24,1)) #remove extra hour 
    v_new = v_new.rename('vomecrty') #name it what you want it named in final netcdf

    # np.save("v_new.npy",v_new)
    v_new.to_netcdf(str(path)+'V_new_{:%Y%m%d}.nc'.format(date_list[0]))
#     print('v_new_{:%d%b%y}_{:%d%b%y}.nc complete'.format(date_list[0],date_list[-1]))