# # Adding Back the Tides
# An effort to make the daily files more accurate as they are currently lacking the tidal pumping that is so important to the flow of the Salish Sea

import xarray as xr
from pathlib import Path
import numpy as np
import datetime as dt

#Outside the loop lets to the e3t shtuff
mydata= xr.open_dataset('/home/sallen/MEOPAR/grid/mesh_mask201702.nc') 
e3t = mydata['e3t_0'].squeeze()

# convert e3t to e3u and to e3v
e3t_xshift = e3t.shift(x=-1,fill_value=0)
e3u = e3t_xshift+e3t
e3u = e3u*0.5
e3u = e3u.rename({'z': 'depthu'})

e3t_yshift = e3t.shift(y=-1,fill_value=0)
e3v = e3t_yshift+e3t
e3v = e3v*0.5
e3v = e3v.rename({'z': 'depthv'})

# also out of the loop we want to get our u-mask and v-mask for the fixins later
xmesh_u = mydata.rename({'z': 'depthu','t': 'time_counter'})
xmesh_v = mydata.rename({'z': 'depthv','t': 'time_counter'})

# set running dates
startday = [dt.datetime(2019,2,28)+dt.timedelta(days=i) for i in range(int(62))]

#start the run
for start in startday:

# dates for each run
    numdays = 3 #16 for all except the last run
    date_list = [start + dt.timedelta(days=x) for x in range(numdays)]

    
    path = Path("/results/SalishSea/nowcast-green.201812/")


    # load u data
    drop_vars = (
        "nav_lon", "bounds_lon", "nav_lat", "bounds_lat", "area", "depthu_bounds", 
        "time_centered", "time_centered_bounds", "time_counter_bounds",
    )

    files = [sorted(path.glob("{:%d%b%y}".format(day).lower()+"/SalishSea_1h_*_grid_U.nc")) for day in date_list]
    # files = [sorted(path.glob("{:02}mar19/SalishSea_1h_*_grid_U.nc".format(day))) for day in days]

    mydata = xr.open_mfdataset(files, drop_variables=drop_vars)
    u = mydata['vozocrtx']


    # load v
    drop_vars = (
        "nav_lon", "bounds_lon", "nav_lat", "bounds_lat", "area", "depthv_bounds", 
        "time_centered", "time_centered_bounds", "time_counter_bounds",
    )

    files = [sorted(path.glob("{:%d%b%y}".format(day).lower()+"/SalishSea_1h_*_grid_V.nc")) for day in date_list]
    # files = [sorted(path.glob("{:02}mar19/SalishSea_1h_*_grid_V.nc".format(day))) for day in days]

    mydata = xr.open_mfdataset(files, drop_variables=drop_vars)
    v = mydata['vomecrty']

    #calcuate bartropic component of u
    ut_h = (u*e3u).sum(dim='depthu')/e3u.sum(dim='depthu')

    #calcuate bartropic component of v
    vt_h = (v*e3v).sum(dim='depthv')/e3v.sum(dim='depthv')

    # Now get the required data from the daily files
    # for these you must add an extra day at the front and the back so that when the interpolation happens we have the correct number of hours

    drop_vars = (
        "nav_lon", "bounds_lon", "nav_lat", "bounds_lat", "area", "depthu_bounds", 
        "time_centered", "time_centered_bounds", "time_counter_bounds",
    )

    files = [sorted(path.glob("{:%d%b%y}".format(day).lower()+"/SalishSea_1d_*_grid_U.nc")) for day in date_list]
    # files = [sorted(path.glob("{:02}mar19/SalishSea_1d_*_grid_U.nc".format(day))) for day in days]

    mydata = xr.open_mfdataset(files, drop_variables=drop_vars)
    u_d = mydata['vozocrtx']

    drop_vars = (
        "nav_lon", "bounds_lon", "nav_lat", "bounds_lat", "area", "depthv_bounds", 
        "time_centered", "time_centered_bounds", "time_counter_bounds",
    )

    files = [sorted(path.glob("{:%d%b%y}".format(day).lower()+"/SalishSea_1d_*_grid_V.nc")) for day in date_list]
    # files = [sorted(path.glob("{:02}mar19/SalishSea_1d_*_grid_V.nc".format(day))) for day in days]

    mydata = xr.open_mfdataset(files, drop_variables=drop_vars)
    v_d = mydata['vomecrty']

    #calcuate bartropic component
    ut_d = (u_d*e3u).sum(dim='depthu')/e3u.sum(dim='depthu')

    #subtract from u to get baroclinic component
    uc_d = u_d-ut_d 


    # interpolate + resample uc_d to get it in an hourly format
    uc_h_interp = uc_d.resample(time_counter="1H", loffset="30min").interpolate("linear")
    # the interp value should have extra hours BUT when you add together in the next step xarray will just be so nice to you and align the two for the corect hours!

    # added together this should give your final u!!!!!
    u_new = ut_h  + uc_h_interp
    
    # select for the hours you want in a single file
    u_new = u_new.isel(time_counter = np.arange(12,36,1))
    
    # now multiply by u-mask and v-max to get rid of the silly edge-effects
    u_new = u_new*xmesh_u.umask[0,:,:,:]*xmesh_u.vmask[0,:,:,:]
    
    # multiply by factor of 2
    u_new = u_new*2
    
    #name it what you want it named in final netcdf
    u_new = u_new.rename('vozocrtx')

    # And save!
    encoding={
          "vozocrtx": {"zlib": True, "complevel": 4, "_FillValue": 0}
    }
    path = '/data/rbeutel/analysis/ssc_tidesback/factor_test/'
    u_new.to_netcdf(str(path)+'U_new_{:%d%b%y}.nc'.format(date_list[1]))
    print('U_new_{:%d%b%y}.nc'.format(date_list[1]))

    # Now for V

    #calcuate bartropic component
    vt_d = (v_d*e3v).sum(dim='depthv')/e3v.sum(dim='depthv')

    #subtract from v to get baroclinic component
    vc_d = v_d-vt_d 

    vc_h_interp = vc_d.resample(time_counter="1H", loffset="30min").interpolate("linear")

    v_new = vt_h  + vc_h_interp
    
    v_new = v_new.isel(time_counter = np.arange(12,36,1))
    
    v_new = v_new*xmesh_v.umask[0,:,:,:]*xmesh_v.vmask[0,:,:,:]
    
    # multiply by factor of 2
    v_new = v_new*2
    
    v_new = v_new.rename('vomecrty') 
    
    encoding={
          "vomecrty": {"zlib": True, "complevel": 4, "_FillValue": 0}
    }

    v_new.to_netcdf(str(path)+'V_new_{:%d%b%y}.nc'.format(date_list[1]))
    print('V_new_{:%d%b%y}.nc'.format(date_list[1]))