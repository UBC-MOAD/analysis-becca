#!/usr/bin/env python
# coding: utf-8

# # Adding Back the Tides
# An effort to make the daily files more accurate as they are currently lacking the tidal pumping that is so important to the flow of the Salish Sea

# In[1]:


import xarray as xr
from pathlib import Path
import numpy as np
import datetime as dt


# In[2]:


start_day, end_day = 1, 5
days = range(start_day, end_day+1)

start = dt.datetime(2019,3,1)
split = dt.datetime(2019,6,1)


# In[3]:


if start >= split:
    path = Path("/results2/SalishSea/nowcast-green.201812/")
else:
    path = Path("/results/SalishSea/nowcast-green.201812/")

drop_vars = (
    "bounds_lon", "bounds_lat", "area", "deptht_bounds", "PAR",
    "time_centered", "time_centered_bounds", "time_counter_bounds", "dissolved_oxygen",
    "sigma_theta", "Fraser_tracer", "dissolved_inorganic_carbon", "total_alkalnity",
)

files = [sorted(path.glob("{:02}mar19/SalishSea_1h_*_carp_T.nc".format(day))) for day in days]

mydata = xr.open_mfdataset(files, drop_variables=drop_vars)
e3t = mydata['e3t']


# In[4]:


# convert e3t to e3u and to e3v
e3t_xshift = e3t.shift(x=-1,fill_value=0)
e3u = e3t_xshift+e3t
e3u = e3u*0.5
e3u = e3u.rename({'deptht': 'depthu'})

e3t_yshift = e3t.shift(y=-1,fill_value=0)
e3v = e3t_yshift+e3t
e3v = e3v*0.5
e3v = e3v.rename({'deptht': 'depthv'})


# In[5]:


drop_vars = (
    "nav_lon", "bounds_lon", "nav_lat", "bounds_lat", "area", "depthu_bounds", 
    "time_centered", "time_centered_bounds", "time_counter_bounds",
)

files = [sorted(path.glob("{:02}mar19/SalishSea_1h_*_grid_U.nc".format(day))) for day in days]

mydata = xr.open_mfdataset(files, drop_variables=drop_vars)
u = mydata['vozocrtx']


# In[6]:


drop_vars = (
    "nav_lon", "bounds_lon", "nav_lat", "bounds_lat", "area", "depthv_bounds", 
    "time_centered", "time_centered_bounds", "time_counter_bounds",
)

files = [sorted(path.glob("{:02}mar19/SalishSea_1h_*_grid_V.nc".format(day))) for day in days]

mydata = xr.open_mfdataset(files, drop_variables=drop_vars)
v = mydata['vomecrty']


# In[7]:


#calcuate bartropic component of u
ut_h = (u*e3u).sum(dim='depthu')/e3u.sum(dim='depthu')


# In[8]:


#calcuate bartropic component of v
vt_h = (v*e3v).sum(dim='depthv')/e3v.sum(dim='depthv')


# u, v, and e3* get reused bellow, so to be safe force the loading of ut_h and vt_h

# In[9]:


n=4
ut_h.load(scheduler="processes", num_workers=n)


# In[ ]:


vt_h.load(scheduler="processes", num_workers=n)


# In[9]:


# Now get the required data from the daily files
drop_vars = (
    "nav_lon", "bounds_lon", "nav_lat", "bounds_lat", "area", "deptht_bounds", "PAR",
    "time_centered", "time_centered_bounds", "time_counter_bounds", "dissolved_oxygen",
    "sigma_theta", "Fraser_tracer", "dissolved_inorganic_carbon", "total_alkalnity",
)

files = [sorted(path.glob("{:02}mar19/SalishSea_1d_*_carp_T.nc".format(day))) for day in days]

mydata = xr.open_mfdataset(files, drop_variables=drop_vars)
e3t = mydata['e3t']

drop_vars = (
    "nav_lon", "bounds_lon", "nav_lat", "bounds_lat", "area", "depthu_bounds", 
    "time_centered", "time_centered_bounds", "time_counter_bounds",
)

files = [sorted(path.glob("{:02}mar19/SalishSea_1d_*_grid_U.nc".format(day))) for day in days]

mydata = xr.open_mfdataset(files, drop_variables=drop_vars)
u = mydata['vozocrtx']

drop_vars = (
    "nav_lon", "bounds_lon", "nav_lat", "bounds_lat", "area", "depthv_bounds", 
    "time_centered", "time_centered_bounds", "time_counter_bounds",
)

files = [sorted(path.glob("{:02}mar19/SalishSea_1d_*_grid_V.nc".format(day))) for day in days]

mydata = xr.open_mfdataset(files, drop_variables=drop_vars)
v = mydata['vomecrty']


# In[10]:


# convert e3t to e3u and to e3v
e3t_xshift = e3t.shift(x=-1,fill_value=0)
e3u = e3t_xshift+e3t
e3u = e3u*0.5
e3u = e3u.rename({'deptht': 'depthu'})

e3t_yshift = e3t.shift(y=-1,fill_value=0)
e3v = e3t_yshift+e3t
e3v = e3v*0.5
e3v = e3v.rename({'deptht': 'depthv'})


# In[11]:


#calcuate bartropic component
ut_d = (u*e3u).sum(dim='depthu')/e3u.sum(dim='depthu')


# In[12]:


#subtract from u to get baroclinic component
uc_d = u-ut_d #does this work even though their ut_d lacks the depth dimension?


# In[30]:


uc_d.load(scheduler="processes", num_workers=n)


# interpolate + resample uc_d to get it in an hourly format

# In[ ]:


uc_h_interp = uc_d.resample(time_counter="1H", loffset="30min").interpolate("linear")


# In[ ]:


u_new = ut_h  + uc_h_interp


# In[ ]:


np.save("u_new.npy",u_new)


# In[ ]:


#calcuate bartropic component
vt_d = (v*e3v).sum(dim='depthv')/e3v.sum(dim='depthv')


# In[ ]:


#subtract from v to get baroclinic component
vc_d = v-vt_d 


# In[ ]:


vc_d.load(scheduler="processes", num_workers=n)


# In[ ]:


vc_h_interp = vc_d.resample(time_counter="1H", loffset="30min").interpolate("linear")


# In[ ]:


v_new = vt_h  + vc_h_interp


# In[ ]:


np.save("v_new.npy",v_new)

