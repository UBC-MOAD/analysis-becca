# # Adding Back the Tides to CIOPS BC12 output
# An effort to make the daily files more accurate as they are currently lacking the tidal pumping that is so important to the flow of the Salish Sea

import xarray as xr
from pathlib import Path
import numpy as np
import datetime as dt

# Outside of loop: we WILL need to do some conversions here, so get e3t from the mesh_mask file
mydata = xr.open_dataset("/ocean/mdunphy/CIOPSW-BC12/grid/mesh_mask_Bathymetry_NEP36_714x1020_SRTM30v11_NOAA3sec_WCTSS_JdeFSalSea.nc")
e3t = mydata['e3t_0'][0,:,:,:]
H = np.sum(e3t,axis=0)

# also out of the loop we want to get our u-mask and v-mask for the fixins later
xmesh_u = mydata.rename({'z': 'depthu'})
xmesh_v = mydata.rename({'z': 'depthv'})
umask = xmesh_u.vmask[0,:,:,:]*xmesh_u.umask[0,:,:,:]
vmask = xmesh_v.vmask[0,:,:,:]*xmesh_v.umask[0,:,:,:]

# set running dates:
startday = [dt.datetime(2016,5,1)+dt.timedelta(days=i) for i in range(int(18*7))]
folders = [dt.datetime(2016,5,1)+dt.timedelta(days=7*(i+1)) for i in range(int(18))]
folders = np.repeat(folders,7)
print(startday[0])
print(startday[-1])

for i in range(len(startday)):
# all within a for loop so you dont have to restart the code every day for 4 years

    # dates for each run
    date_list_daily = [startday[i],startday[i+1]]
    folderday_daily = [folders[i], folders[i+1]]
    date_list_hourly = [startday[i+1],startday[i+2]]
#     print(date_list_hourly)
    folderday_hourly = [folders[i+1], folders[i+2]]
#     print(folderday_hourly)

    path = Path("/ocean/mdunphy/CIOPSW-BC12/")
    
    # load U
    drop_vars = (
        "sbu", "tauuo", "time_counter_bounds", "time_instant_bounds", "uos", "time_instant",
    )

    files = [sorted(path.glob("{:%Y%m%d}00/BC12_1h_grid_U_2D_{:%Y%m%d}_{:%Y%m%d}.nc".format(folderday_hourly[i], date_list_hourly[i], date_list_hourly[i]))) for i in range(len(date_list_hourly))]
#     print(files)

    mydata = xr.open_mfdataset(files, drop_variables=drop_vars)
    ut_h = mydata['ubar']

    # load V
    drop_vars = (
        "sbv", "tauvo", "time_counter_bounds", "time_instant_bounds", "vos", "time_instant",
    )

    files = [sorted(path.glob("{:%Y%m%d}00/BC12_1h_grid_V_2D_{:%Y%m%d}_{:%Y%m%d}.nc".format(folderday_hourly[i], date_list_hourly[i], date_list_hourly[i]))) for i in range(len(date_list_hourly))]

    mydata = xr.open_mfdataset(files, drop_variables=drop_vars)
    vt_h = mydata['vbar']

    # 2D hourly outputs are alraedy in barotropic form so you dont need to do any conversion! 

    # Now get the required data from the daily files!

    drop_vars = (
        "depthu_bounds", "nav_lat", "nav_lon", 'time_counter_bounds', 'time_instant',
        'time_instant_bounds', 
    )

    files = [sorted(path.glob("{:%Y%m%d}00/BC12_1d_grid_U_{:%Y%m%d}_{:%Y%m%d}.nc".format(folderday_daily[i], date_list_daily[i], date_list_daily[i]))) for i in range(len(date_list_daily))]

    mydata = xr.open_mfdataset(files, drop_variables=drop_vars)
    u_d = mydata['uo']

    drop_vars = (
        "depthv_bounds", "nav_lat", "nav_lon", 'time_counter_bounds', 'time_instant',
        'time_instant_bounds', 
    )

    files = [sorted(path.glob("{:%Y%m%d}00/BC12_1d_grid_V_{:%Y%m%d}_{:%Y%m%d}.nc".format(folderday_daily[i], date_list_daily[i], date_list_daily[i]))) for i in range(len(date_list_daily))]

    mydata = xr.open_mfdataset(files, drop_variables=drop_vars)
    v_d = mydata['vo']
    
    #got to resample the veloties right away so that we can add the impact of stretching and compressing the grid
    offset = dt.timedelta(hours=1)
    u_d = u_d.resample(time_counter="1H", loffset=offset).interpolate("linear")
    v_d = v_d.resample(time_counter="1H", loffset=offset).interpolate("linear")

#     u_d = u_d.where(u_d != 0)
#     v_d = v_d.where(v_d != 0)

    # bring in ssh for the same days
    drop_vars = (
    "mldkz5", 'mldr10_1', "nav_lat", "nav_lon", 'time_counter_bounds', 'time_instant',
    'time_instant_bounds', 'sbs', 'sbt', 'sos', 'ssh_ib', 'ssh_tide', 'tos'
    )
    
    files = [sorted(path.glob("{:%Y%m%d}00/BC12_1h_grid_T_2D_{:%Y%m%d}_{:%Y%m%d}.nc".format(folderday_daily[i], date_list_daily[i], date_list_daily[i]))) for i in range(len(date_list_daily))]
    
    mydata = xr.open_mfdataset(files, drop_variables=drop_vars)
    eta = mydata['zos']
    
    #calculate the e3t over this time period
    e3t_ssh = e3t+((eta*e3t)/H)
    
    # convert e3t to e3u and to e3v
    e3t_xshift = e3t_ssh.shift(x=-1,fill_value=0)
    e3u = e3t_xshift+e3t
    e3u = e3u*0.5 #check to make sure this worked
    e3u = e3u.rename({'z': 'depthu'})

    e3t_yshift = e3t_ssh.shift(y=-1,fill_value=0)
    e3v = e3t_yshift+e3t
    e3v = e3v*0.5
    e3v = e3v.rename({'z': 'depthv'})
    
    #calcuate barotropic component
    ut_d = (u_d*e3u*umask).sum(dim='depthu')/(e3u*umask).sum(dim='depthu')

    #subtract from u to get baroclinic component
    uc_d = u_d-ut_d

    # added together this should give your final u!!!!!
    u_new = ut_h  + uc_d 
    
    # now multiply by u-mask and v-max to get rid of the silly edge-effects
    u_new = u_new*umask
    
    #name it what you want it named in final netcdf
    u_new = u_new.rename('vozocrtx')
    
    # order the variables the way you need em
    u_new = u_new.transpose('time_counter','depthu','y','x')

    # And save!
    encoding={
          "vozocrtx": {"zlib": True, "complevel": 4, "_FillValue": 0}
    }
    
    path = '/ocean/rbeutel/data/'
    u_new.to_netcdf(str(path)+'{:%Y%m}/U_new_{:%Y%m%d}.nc'.format(date_list_hourly[0],date_list_hourly[0]), encoding=encoding)
    print('U_new_{:%d%b%y}'.format(date_list_hourly[0]))
    
    # Now for V!

    #calcuate bartropic component
    vt_d = (v_d*e3v*vmask).sum(dim='depthv')/(e3v*vmask).sum(dim='depthv')

    #subtract from v to get baroclinic component
    vc_d = v_d-vt_d 

    v_new = vt_h  + vc_d 
    
    # now multiply by u-mask and v-max to get rid of the silly edge-effects
    v_new = v_new*vmask
    
    #name it what you want it named in final netcdf
    v_new = v_new.rename('vomecrty') 
    
    v_new = v_new.transpose('time_counter','depthv','y','x')
    
    encoding={
          "vomecrty": {"zlib": True, "complevel": 4, "_FillValue": 0}
    }
    v_new.to_netcdf(str(path)+'{:%Y%m}/V_new_{:%Y%m%d}.nc'.format(date_list_hourly[0],date_list_hourly[0]), encoding=encoding)
    print('V_new_{:%d%b%y}'.format(date_list_hourly[0]))