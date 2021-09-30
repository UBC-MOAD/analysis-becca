#!/usr/bin/env python
# coding: utf-8

# # Adding Back the Tides
# An effort to make the daily files more accurate as they are currently lacking the tidal pumping that is so important to the flow of the Salish Sea

# In[1]:


import xarray as xr
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np


# In[2]:


start_day, end_day = 1, 30
days = range(start_day, end_day+1)


# In[3]:


drop_vars = (
    "bounds_lon", "bounds_lat", "area", "deptht_bounds", "PAR",
    "time_centered_bounds", "time_counter_bounds", "dissolved_oxygen",
    "sigma_theta", "Fraser_tracer", "dissolved_inorganic_carbon", "total_alkalnity",
)

path = Path("/results2/SalishSea/nowcast-green.201905/")
files = [sorted(path.glob("{:02}mar19/SalishSea_1h_*_carp_T.nc".format(day))) for day in days]

mydata = xr.open_mfdataset(files, drop_variables=drop_vars)
e3t = mydata['e3t']


# In[4]:


# convert e3t to e3u and to e3v
e3u = np.zeros(np.shape(e3t))
e3v = np.zeros(np.shape(e3t))

for i in range(np.shape(e3t)[2]-1):
    e3u[:,:,i,:] = 1/2*(e3t[:,:,i,:]+e3t[:,:,i+1,:])
    
for i in range(np.shape(e3t)[3]-1):
    e3v[:,:,:,j] = 1/2*(e3t[:,:,:,j]+e3t[:,:,:,j+1])
    
e3u[:,:,-1,:] = e3u[:,:,-2,:]
e3v[:,:,:,-1] = e3u[:,:,:,-2]


# In[4]:


drop_vars = (
    "bounds_lon", "bounds_lat", "area", "depthu_bounds", 
    "time_centered_bounds", "time_counter_bounds",
)

path = Path("/results2/SalishSea/nowcast-green.201905/")
files = [sorted(path.glob("{:02}mar19/SalishSea_1h_*_grid_U.nc".format(day))) for day in days]

mydata = xr.open_mfdataset(files, drop_variables=drop_vars)
u = mydata['vozocrtx']


# In[12]:


drop_vars = (
    "bounds_lon", "bounds_lat", "area", "depthu_bounds", 
    "time_centered_bounds", "time_counter_bounds",
)

path = Path("/results2/SalishSea/nowcast-green.201905/")
files = [sorted(path.glob("{:02}mar19/SalishSea_1h_*_grid_V.nc".format(day))) for day in days]

mydata = xr.open_mfdataset(files, drop_variables=drop_vars)
v = mydata['vomecrty']


# In[ ]:


#calcuate bartropic component of u

ut_h = np.zeros((len(u.time_counter), len(u.y), len(u.x)))

for t in range(len(u.time_counter)):
    for y in range(len(u.y)):
        for x in range(len(u.x)):
            ut_h[t,y,x] = sum(u[t,:,y,x].values*e3u[t,:,y,x].values)/sum(e3u[t,:,y,x].values)


# In[ ]:


#calcuate bartropic component of v

vt_h = np.zeros((len(v.time_counter), len(v.y), len(v.x)))

for t in range(len(v.time_counter)):
    for y in range(len(v.y)):
        for x in range(len(v.x)):
            vt_h[t,y,x] = sum(v[t,:,y,x].values*e3v[t,:,y,x].values)/sum(e3v[t,:,y,x].values)


# In[ ]:


# Now get the required data from the daily files

drop_vars = (
    "bounds_lon", "bounds_lat", "area", "deptht_bounds", "PAR",
    "time_centered_bounds", "time_counter_bounds", "dissolved_oxygen",
    "sigma_theta", "Fraser_tracer", "dissolved_inorganic_carbon", "total_alkalnity",
)

path = Path("/results2/SalishSea/nowcast-green.201905/")
files = [sorted(path.glob("{:02}mar19/SalishSea_1d_*_carp_T.nc".format(day))) for day in days]

mydata = xr.open_mfdataset(files, drop_variables=drop_vars)
e3t = mydata['e3t']

drop_vars = (
    "bounds_lon", "bounds_lat", "area", "depthu_bounds", 
    "time_centered_bounds", "time_counter_bounds",
)

path = Path("/results2/SalishSea/nowcast-green.201905/")
files = [sorted(path.glob("{:02}mar19/SalishSea_1d_*_grid_U.nc".format(day))) for day in days]

mydata = xr.open_mfdataset(files, drop_variables=drop_vars)
u = mydata['vozocrtx']

drop_vars = (
    "bounds_lon", "bounds_lat", "area", "depthu_bounds", 
    "time_centered_bounds", "time_counter_bounds",
)

path = Path("/results2/SalishSea/nowcast-green.201905/")
files = [sorted(path.glob("{:02}mar19/SalishSea_1d_*_grid_V.nc".format(day))) for day in days]

mydata = xr.open_mfdataset(files, drop_variables=drop_vars)
v = mydata['vomecrty']


# In[ ]:


# convert e3t to e3u and to e3v
e3u = np.zeros(np.shape(e3t))
e3v = np.zeros(np.shape(e3t))

for i in range(np.shape(e3t)[2]-1):
    e3u[:,:,i,:] = 1/2*(e3t[:,:,i,:]+e3t[:,:,i+1,:])
    
for i in range(np.shape(e3t)[3]-1):
    e3v[:,:,:,j] = 1/2*(e3t[:,:,:,j]+e3t[:,:,:,j+1])
    
e3u[:,:,-1,:] = e3u[:,:,-2,:]
e3v[:,:,:,-1] = e3u[:,:,:,-2]


# In[ ]:


#calcuate bartropic component

ut_d = np.zeros((len(u.time_counter), len(u.y), len(u.x)))

for t in range(len(u.time_counter)):
    for y in range(len(u.y)):
        for x in range(len(u.x)):
            ut_d[t,y,x] = sum(u[t,:,y,x].values*e3u[t,:,y,x].values)/sum(e3u[t,:,y,x].values)


# In[ ]:


#subtract from u to get baroclinic component

uc_d = np.zeros(np.shape(u))

for t in range(len(u.time_counter)):
    for z in range(40):
        for y in range(len(u.y)):
            for x in range(len(u.x)):
                uc_d[t,z,y,x] = u[t,z,y,x]-ut_d[t,y,x]


# In[ ]:


# treat baroclinic daily as a constant hourly component and add barotropic hourly to it to add back tides

u_new = np.zeros(len(ut_h[:,0,0]),40,len(u.y), len(u.x))

for t in range(len(ut_h[:,0,0])):
    for z in range(40):
        for y in range(len(u.y)):
            for x in range(len(u.x)):
                u_new[t,z,y,x] = ut_h[t,y,x] + uc_d[t//24,z,y,x] 
                
np.save("u_new.npy",u_new)


# In[ ]:


#calcuate bartropic component

vt_d = np.zeros((len(v.time_counter), len(v.y), len(v.x)))

for t in range(len(v.time_counter)):
    for y in range(len(v.y)):
        for x in range(len(v.x)):
            vt_d[t,y,x] = sum(v[t,:,y,x].values*e3v[t,:,y,x].values)/sum(e3v[t,:,y,x].values)


# In[ ]:


#subtract from v to get baroclinic component

vc_d = np.zeros(np.shape(v))

for t in range(len(v.time_counter)):
    for z in range(40):
        for y in range(len(v.y)):
            for x in range(len(v.x)):
                vc_d[t,z,y,x] = v[t,z,y,x]-vt_d[t,y,x]


# In[ ]:


# treat baroclinic daily as a constant hourly component and add barotropic hourly to it to add back tides

v_new = np.zeros(len(vt_h[:,0,0]),40,len(v.y), len(v.x))

for t in range(len(vt_h[:,0,0])):
    for z in range(40):
        for y in range(len(v.y)):
            for x in range(len(v.x)):
                v_new[t,z,y,x] = vt_h[t,y,x] + vc_d[t//24,z,y,x] 
                
np.save("v_new.npy",v_new)

