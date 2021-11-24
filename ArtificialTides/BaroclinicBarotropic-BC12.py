# # Adding Back the Tides to CIOPS BC12 output
# An effort to make the daily files more accurate as they are currently lacking the tidal pumping that is so important to the flow of the Salish Sea

import xarray as xr
from pathlib import Path
import numpy as np
import datetime as dt

start = dt.datetime(2015,11,22) #start day of your run

# dates for each run
numdays = 15 #15 for all except the last run (14)
date_list = [start + dt.timedelta(days=x) for x in range(numdays)]
folderday = [start+ dt.timedelta(days=7),start+ dt.timedelta(days=7*2),start+ dt.timedelta(days=7*2)]
folderday = np.repeat(folderday,7)
folderday = folderday[:-6]

# In[3]:

path = Path("/ocean/mdunphy/CIOPSW-BC12/")

#load e3t
drop_vars = (
    "bounds_lon", "bounds_lat", "area", "deptht_bounds", "PAR",
    "time_centered", "time_centered_bounds", "time_counter_bounds", "dissolved_oxygen",
    "sigma_theta", "Fraser_tracer", "dissolved_inorganic_carbon", "total_alkalnity",
)

files = [sorted(path.glob("{:%d%b%y}".format(day).lower()+"/SalishSea_1h_*_carp_T.nc")) for day in date_list]
# files = [sorted(path.glob("{:02}mar19/SalishSea_1h_*_carp_T.nc".format(day))) for day in days]

mydata = xr.open_mfdataset(files, drop_variables=drop_vars)
e3t = mydata['e3t']

# convert e3t to e3u and to e3v
e3t_xshift = e3t.shift(x=-1,fill_value=0)
e3u = e3t_xshift+e3t
e3u = e3u*0.5
e3u = e3u.rename({'deptht': 'depthu'})

e3t_yshift = e3t.shift(y=-1,fill_value=0)
e3v = e3t_yshift+e3t
e3v = e3v*0.5
e3v = e3v.rename({'deptht': 'depthv'})


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


# u, v, and e3* get reused bellow, so to be safe force the loading of ut_h and vt_h
# ut_h.load()
# vt_h.load()

# Now get the required data from the daily files
# for these you must add an extra day at the front and the back so that when the interpolation happens we have the correct number of hours
drop_vars = (
    "nav_lon", "bounds_lon", "nav_lat", "bounds_lat", "area", "deptht_bounds", "PAR",
    "time_centered", "time_centered_bounds", "time_counter_bounds", "dissolved_oxygen",
    "sigma_theta", "Fraser_tracer", "dissolved_inorganic_carbon", "total_alkalnity",
)

files = [sorted(path.glob("{:%d%b%y}".format(day).lower()+"/SalishSea_1d_*_carp_T.nc")) for day in date_list]
# files = [sorted(path.glob("{:02}mar19/SalishSea_1d_*_carp_T.nc".format(day))) for day in days]

mydata = xr.open_mfdataset(files, drop_variables=drop_vars)
e3t_d = mydata['e3t']

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

# convert e3t to e3u and to e3v
e3t_xshift = e3t_d.shift(x=-1,fill_value=0)
e3u_d = e3t_xshift+e3t_d
e3u_d = e3u_d*0.5
e3u_d = e3u_d.rename({'deptht': 'depthu'})

e3t_yshift = e3t_d.shift(y=-1,fill_value=0)
e3v_d = e3t_yshift+e3t_d
e3v_d = e3v_d*0.5
e3v_d = e3v_d.rename({'deptht': 'depthv'})

#calcuate bartropic component
ut_d = (u_d*e3u_d).sum(dim='depthu')/e3u_d.sum(dim='depthu')

#subtract from u to get baroclinic component
uc_d = u_d-ut_d #does this work even though their ut_d lacks the depth dimension?


# interpolate + resample uc_d to get it in an hourly format
uc_h_interp = uc_d.resample(time_counter="1H", loffset="30min").interpolate("linear")
# the interp value should have extra hours BUT when you add together in the next step xarray will just be so nice to you and align the two for the corect hours!

# added together this should give your final u!!!!!
u_new = ut_h  + uc_h_interp

# And save!
# np.save("u_new.npy",u_new)
path = '/data/rbeutel/analysis/ssc_tidesback/'
u_new.to_netcdf(str(path)+'u_new_{:%d%b%y}_{:%d%b%y}.nc'.format(date_list[0],date_list[-1]))

# Now for V

#calcuate bartropic component
vt_d = (v_d*e3v_d).sum(dim='depthv')/e3v_d.sum(dim='depthv')


#subtract from v to get baroclinic component
vc_d = v_d-vt_d 

# vc_d.load(scheduler="processes", num_workers=n)

vc_h_interp = vc_d.resample(time_counter="1H", loffset="30min").interpolate("linear")

v_new = vt_h  + vc_h_interp

# np.save("v_new.npy",v_new)
v_new.to_netcdf(str(path)+'v_new_{:%d%b%y}_{:%d%b%y}.nc'.format(date_list[0],date_list[-1]))
