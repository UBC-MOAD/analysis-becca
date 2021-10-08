#!/usr/bin/env python
# coding: utf-8

# # Adding Back the Tides
# An effort to make the daily files more accurate as they are currently lacking the tidal pumping that is so important to the flow of the Salish Sea

# In[1]:


import xarray as xr
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import datetime as dt


# In[2]:


start_day, end_day = 1, 30
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
    "time_centered_bounds", "time_counter_bounds", "dissolved_oxygen",
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

# e3u = np.zeros(np.shape(e3t))
# e3v = np.zeros(np.shape(e3t))

# for i in range(np.shape(e3t)[2]-1):
#     e3u[:,:,i,:] = 1/2*(e3t[:,:,i,:]+e3t[:,:,i+1,:])
    
# for i in range(np.shape(e3t)[3]-1):
#     e3v[:,:,:,j] = 1/2*(e3t[:,:,:,j]+e3t[:,:,:,j+1])
    
# e3u[:,:,-1,:] = e3u[:,:,-2,:]
# e3v[:,:,:,-1] = e3u[:,:,:,-2]


# In[5]:


drop_vars = (
    "bounds_lon", "bounds_lat", "area", "depthu_bounds", 
    "time_centered_bounds", "time_counter_bounds",
)

files = [sorted(path.glob("{:02}mar19/SalishSea_1h_*_grid_U.nc".format(day))) for day in days]

mydata = xr.open_mfdataset(files, drop_variables=drop_vars)
u = mydata['vozocrtx']


# In[6]:


drop_vars = (
    "bounds_lon", "bounds_lat", "area", "depthu_bounds", 
    "time_centered_bounds", "time_counter_bounds",
)

files = [sorted(path.glob("{:02}mar19/SalishSea_1h_*_grid_V.nc".format(day))) for day in days]

mydata = xr.open_mfdataset(files, drop_variables=drop_vars)
v = mydata['vomecrty']


# In[7]:


#calcuate bartropic component of u
ut_h = (u*e3u).sum(dim='depthu')/e3u.sum(dim='depthu')

# ut_h = np.zeros((len(u.time_counter), len(u.y), len(u.x)))

# for t in range(len(u.time_counter)):
#     for y in range(len(u.y)):
#         for x in range(len(u.x)):
#             ut_h[t,y,x] = sum(u[t,:,y,x].values*e3u[t,:,y,x].values)/sum(e3u[t,:,y,x].values)


# In[8]:


#calcuate bartropic component of v
vt_h = (v*e3v).sum(dim='depthv')/e3v.sum(dim='depthv')

# vt_h = np.zeros((len(v.time_counter), len(v.y), len(v.x)))

# for t in range(len(v.time_counter)):
#     for y in range(len(v.y)):
#         for x in range(len(v.x)):
#             vt_h[t,y,x] = sum(v[t,:,y,x].values*e3v[t,:,y,x].values)/sum(e3v[t,:,y,x].values)


# In[9]:


# Now get the required data from the daily files
drop_vars = (
    "bounds_lon", "bounds_lat", "area", "deptht_bounds", "PAR",
    "time_centered_bounds", "time_counter_bounds", "dissolved_oxygen",
    "sigma_theta", "Fraser_tracer", "dissolved_inorganic_carbon", "total_alkalnity",
)

files = [sorted(path.glob("{:02}mar19/SalishSea_1d_*_carp_T.nc".format(day))) for day in days]

mydata = xr.open_mfdataset(files, drop_variables=drop_vars)
e3t = mydata['e3t']

drop_vars = (
    "bounds_lon", "bounds_lat", "area", "depthu_bounds", 
    "time_centered_bounds", "time_counter_bounds",
)

files = [sorted(path.glob("{:02}mar19/SalishSea_1d_*_grid_U.nc".format(day))) for day in days]

mydata = xr.open_mfdataset(files, drop_variables=drop_vars)
u = mydata['vozocrtx']

drop_vars = (
    "bounds_lon", "bounds_lat", "area", "depthu_bounds", 
    "time_centered_bounds", "time_counter_bounds",
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

# e3u = np.zeros(np.shape(e3t))
# e3v = np.zeros(np.shape(e3t))

# for i in range(np.shape(e3t)[2]-1):
#     e3u[:,:,i,:] = 1/2*(e3t[:,:,i,:]+e3t[:,:,i+1,:])
    
# for i in range(np.shape(e3t)[3]-1):
#     e3v[:,:,:,j] = 1/2*(e3t[:,:,:,j]+e3t[:,:,:,j+1])
    
# e3u[:,:,-1,:] = e3u[:,:,-2,:]
# e3v[:,:,:,-1] = e3u[:,:,:,-2]


# In[11]:


#calcuate bartropic component
ut_d = (u*e3u).sum(dim='depthu')/e3u.sum(dim='depthu')

# ut_d = np.zeros((len(u.time_counter), len(u.y), len(u.x)))

# for t in range(len(u.time_counter)):
#     for y in range(len(u.y)):
#         for x in range(len(u.x)):
#             ut_d[t,y,x] = sum(u[t,:,y,x].values*e3u[t,:,y,x].values)/sum(e3u[t,:,y,x].values)


# In[ ]:


#subtract from u to get baroclinic component
uc_d = u-ut_d #does this work even though their ut_d lacks the depth dimension?

# uc_d = np.zeros(np.shape(u))

# for t in range(len(u.time_counter)):
#     for z in range(40):
#         for y in range(len(u.y)):
#             for x in range(len(u.x)):
#                 uc_d[t,z,y,x] = u[t,z,y,x]-ut_d[t,y,x]


# In[ ]:


# treat baroclinic daily as a constant hourly component and add barotropic hourly to it to add back tides

u_new = np.zeros([ut_h.shape[0],40,ut_h.shape[1], ut_h.shape[2]])

for t in range(len(ut_h[:,0,0])):
    for z in range(40):
        for y in range(len(u.y)):
            for x in range(len(u.x)):
                u_new[t,z,y,x] = ut_h[t,y,x] + uc_d[t//24,z,y,x] 
                
np.save("u_new.npy",u_new)


# In[ ]:


#calcuate bartropic component
vt_d = (v*e3v).sum(dim='depthv')/e3v.sum(dim='depthv')

# vt_d = np.zeros((len(v.time_counter), len(v.y), len(v.x)))

# for t in range(len(v.time_counter)):
#     for y in range(len(v.y)):
#         for x in range(len(v.x)):
#             vt_d[t,y,x] = sum(v[t,:,y,x].values*e3v[t,:,y,x].values)/sum(e3v[t,:,y,x].values)


# In[ ]:


#subtract from v to get baroclinic component
vc_d = v-vt_d #does this work even though vt_d lacks the depth dimension?

# vc_d = np.zeros(np.shape(v))

# for t in range(len(v.time_counter)):
#     for z in range(40):
#         for y in range(len(v.y)):
#             for x in range(len(v.x)):
#                 vc_d[t,z,y,x] = v[t,z,y,x]-vt_d[t,y,x]


# In[ ]:


# treat baroclinic daily as a constant hourly component and add barotropic hourly to it to add back tides

v_new = np.zeros([len(vt_h[:,0,0]),40,len(v.y), len(v.x)])

for t in range(len(vt_h[:,0,0])):
    for z in range(40):
        for y in range(len(v.y)):
            for x in range(len(v.x)):
                v_new[t,z,y,x] = vt_h[t,y,x] + vc_d[t//24,z,y,x] 
                
np.save("v_new.npy",v_new)

