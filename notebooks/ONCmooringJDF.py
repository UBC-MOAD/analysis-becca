#!/usr/bin/env python
# coding: utf-8

# # ONC Mouring at the Mouth of the JDF
# lat = 48.508133 <br>
# lon = -124.749117<br>
# depth = 226 m <br>
# station name = OceanNetworksCanada-Pacific-SalishSea-JuandeFucaStrait-JFCNMooring

# In[1]:


import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import datetime as dt
from dateutil.relativedelta import relativedelta
from salishsea_tools import evaltools as et, viz_tools
import cmocean as cmo
import netCDF4 as nc
from matplotlib.colors import LogNorm
from salishsea_tools import geo_tools
import gsw

# %matplotlib inline


# In[ ]:


def calc_stats(x, y):
    """
    """
    
    stats = {}
    MSE = np.mean((y - x)**2)
    stats['RMSE'] = np.sqrt(MSE)
    stats['bias'] = np.mean(y) - np.mean(x)
    stats['WSS'] = 1 - MSE / np.mean((abs(y - np.mean(x)) + abs(x - np.mean(x)))**2)
    
    return stats


def plot_panel(ax, x, y, lims, units):
    """
    """
    
    stats = calc_stats(x, y)

    statstext = f"RMSE = {stats['RMSE']:.3f} {units}\nbias = {stats['bias']:.3f} {units}\nWSS = {stats['WSS']:.3f}"
    
    props = dict(boxstyle='round', facecolor='w', alpha=0.9)
    c = ax.text(0.5, 0.04, statstext, bbox=props, transform=ax.transAxes)
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    
    return c


# ## ONC Mooring with More Data - JF2C
# lat = 48.360383 <br>
# lon = -124.213367 <br>
# depth = 175 m <br>
# station name =  OceanNetworksCanada-Pacific-SalishSea-JuandeFucaStrait-JF2CMooring

# In[2]:


JF2C = pd.read_csv('JuandeFucaStrait_JF2CMooring_ConductivityTemperatureDepth_20121101T000539Z_20200926T135042Z-NaN_clean.csv', header=113, skiprows=[114])


# In[3]:


JF2C['dtUTC'] = [ pd.to_datetime('{:%Y-%m-%d} {:%H:%M:%S}'.format(pd.to_datetime(ii), pd.to_datetime(ii))) for ii in JF2C['#"Time UTC (yyyy-mm-ddThh:mm:ss.fffZ)"']]


# In[4]:


JF2C[' "Practical Salinity (psu)"'] = pd.to_numeric(JF2C[' "Practical Salinity (psu)"'],errors='coerce')
JF2C[' "Temperature (C)"'] = pd.to_numeric(JF2C[' "Temperature (C)"'],errors='coerce')


# In[5]:


fig, ax = plt.subplots(figsize=(14, 7))

ax.plot(JF2C['dtUTC'],JF2C[' "Practical Salinity (psu)"'], c='r')
ax.set_ylabel("Salinity (psu)",fontsize=14, c='r')
ax.set_ylim(31.9,34.1)
ax.set_xlim(pd.to_datetime('2012-10-31'), pd.to_datetime('2020-09-30'))

ax2 = ax2=ax.twinx()
ax2.plot(JF2C['dtUTC'],JF2C[' "Temperature (C)"'], c='b')
ax2.set_ylabel("Temperature (C)",fontsize=14, c='b')

# #set summer upwelling as apr1-sep30
season = [pd.to_datetime(pd.to_datetime('2013-04-01') + relativedelta(months=6*i)) for i in range(16)]
for i in range(0, len(season),2):
    ax.fill_between(JF2C['dtUTC'], 30, 35, where=(JF2C['dtUTC']<season[i+1])&(JF2C['dtUTC']>=season[i]), color='grey', alpha=0.3)


# grey = summer upwelling (april 01 to september 30) <br>
# white = winter downwelling (october 01 to march 31)

# In[6]:


with nc.Dataset('/ocean/rbeutel/MEOPAR/grid/mesh_mask201702.nc') as mesh:
    tmask=np.copy(mesh.variables['tmask'][0,:,:,:])
    navlat=np.copy(mesh.variables['nav_lat'][:,:])
    navlon=np.copy(mesh.variables['nav_lon'][:,:])
    e3t_0=np.copy(mesh.variables['e3t_0'][0,:,:,:])
    gdepw=np.copy(mesh.variables['gdepw_0'][0,:,:,:])
    print(mesh.variables.keys())


# In[7]:


obslon=-124.213367
obslat=48.360383
j,i=geo_tools.find_closest_model_point(-124.213367,obslat,navlon,navlat)


# In[8]:


JF2C['i'] = i
JF2C['j'] = j


# In[115]:


# [(x,y) for (x,y) in enumerate(gdepw[:,j,i])]


# In[10]:


JF2C['k'] = 29


# In[11]:


JF2C.head()


# In[14]:


# path to model files:
PATH= '/results2/SalishSea/nowcast-green.201905/' #path to the results file

# start and end dates for analysis: #will cut off anything before and after
start_date = dt.datetime(2013,1,1)
end_date = dt.datetime(2020,9,26)

# number of days per model file:
flen=1 #the number of days per model output file (almost always 1 for susan, elise does 10)

# dictionary mapping desired model variables to the file types where they are found
filemap={'vosaline':'grid_T','votemper':'grid_T'} 

# dictionary mapping model file types to their time resolution in hours (1 is hourly files, 24 is daily)
fdict={'ptrc_T':1,'grid_T':1} #loading data from hourly files, 1 for 1 hour (daily replace with a 24)

# results format (naming convention and storage format of the files)
# -- nowcast: files like 01jan15/SalishSea_1h_20150101_20150101_ptrc_T.nc
# -- long: files like SalishSea_1h_20150206_20150804_ptrc_T_20150427-20150506.nc, all in one directory
namfmt='nowcast'


# In[15]:


# match model output to observations and return both in a dataframe
# the model variables will have their original names prefixed by mod_, eg mod_vosaline
# the observation file names are unchanged. 
data=et.matchData(data=JF2C,filemap=filemap, fdict=fdict, mod_start=start_date, mod_end=end_date, 
                  mod_nam_fmt=namfmt, mod_basedir=PATH, mod_flen=flen, preIndexed=True)
# 


# In[20]:


# convert observed from psu to g/kg
data['lat'] = 48.360383
data['lon'] =  -124.213367
data["Sea Pressure (db)"] = pd.to_numeric(JF2C[' "Pressure (decibar)"'],errors='coerce')-10.1325

data['Real Salinity (g/kg)'] = gsw.conversions.SA_from_SP(np.array(data[' "Practical Salinity (psu)"']), np.array(data["Sea Pressure (db)"]),np.array(data['lon']),np.array(data['lat']))


# In[21]:


# convert modelled from CT to in-situ T
data['modelled T (C)'] = gsw.conversions.t_from_CT(np.array(data['mod_vosaline']), np.array(data['mod_votemper']), np.array(data["Sea Pressure (db)"]))


# ### mod vs. obs plots

# In[22]:


fig, axs = plt.subplots(1, 2, figsize = (12, 6))

axs[0].plot((32,34.6),(32,34.6),'k-',alpha=.2)
axs[1].plot((5.4,10),(5.4,10),'k-',alpha=.2)

iiT=(~np.isnan(data[' "Temperature (C)"']))&(~np.isnan(data['modelled T (C)']))
iiS=(~np.isnan(data['Real Salinity (g/kg)']))&(~np.isnan(data['mod_vosaline']))
counts, xedges, yedges, m1=axs[1].hist2d(data.loc[iiT,[' "Temperature (C)"']].values.flatten(),
                                      data.loc[iiT,['modelled T (C)']].values.flatten(),bins=40,norm=LogNorm())
counts, xedges, yedges, m2=axs[0].hist2d(data.loc[iiS,['Real Salinity (g/kg)']].values.flatten(),
                                      data.loc[iiS,['mod_vosaline']].values.flatten(),bins=40,norm=LogNorm())

cb0=fig.colorbar(m2,ax=axs[0])
cb0.set_label('Count')
cb1=fig.colorbar(m1,ax=axs[1])
cb1.set_label('Count')

ntick=np.arange(32.2, 34.4, 0.4)
axs[0].set_xlim((32,34.6))
axs[0].set_ylim((32,34.6))
axs[0].set_xticks(ntick)
axs[0].set_yticks(ntick)
    
stick=np.arange(5.6,9.8, 0.6)
axs[1].set_xlim((5.4,10))
axs[1].set_ylim((5.4,10))
axs[1].set_xticks(stick)
axs[1].set_yticks(stick)
    
for ax in (axs[1],axs[0]):
    ax.set_aspect(1, adjustable='box')
    
axs[1].set_ylabel('Modeled',fontsize=12)
axs[0].set_ylabel('Modeled',fontsize=12)
axs[1].set_xlabel('Observed',fontsize=12)
axs[0].set_xlabel('Observed',fontsize=12)

axs[0].set_title('Salinity (g/kg)',fontsize=12)
axs[1].set_title('Temperature ($^{\circ}$C)',fontsize=12)

# plot the stats pannel
plot_panel(axs[0], data['Real Salinity (g/kg)'], data['mod_vosaline'], (32,34.6), 'g/kg')
plot_panel(axs[1], data[' "Temperature (C)"'], data['modelled T (C)'], (5.4,10), '$^{\circ}$C')

plt.savefig('SSCeval.png')


# ^ looks like there's a lsight tendency in the model to nunderestimate the salinity and temperature

# # do same comparison to CIOPS-BC12?

# In[23]:


# the question is - does the function above only work for SSC files? - they don't
# use the same process, try to make similar functions for CIOPS


# In[6]:


with nc.Dataset('/ocean/mdunphy/CIOPSW-BC12/grid/mesh_mask_Bathymetry_NEP36_714x1020_SRTM30v11_NOAA3sec_WCTSS_JdeFSalSea.nc') as mesh:
    tmask=np.copy(mesh.variables['tmask'][0,:,:,:])
    navlat=np.copy(mesh.variables['nav_lat'][:,:])
    navlon=np.copy(mesh.variables['nav_lon'][:,:])
    e3t_0=np.copy(mesh.variables['e3t_0'][0,:,:,:])
    gdepw=np.copy(mesh.variables['gdepw_0'][0,:,:,:])
    print(mesh.variables.keys())


# In[7]:


obslon=-124.213367
obslat=48.360383
j,i = np.unravel_index(np.argmin([np.abs(navlon - obslon)+np.abs(navlat - obslat)], axis=None), navlon.shape)


# In[8]:


JF2C['i'] = i
JF2C['j'] = j


# In[9]:


# [(x,y) for (x,y) in enumerate(gdepw[:,j,i])]


# In[40]:


k=29
JF2C['k'] = k


# In[61]:


# add column for "real" salinity instead of practical
# convert observed from psu to g/kg
JF2C['lat'] = 48.360383
JF2C['lon'] =  -124.213367
JF2C["Sea Pressure (db)"] = pd.to_numeric(JF2C[' "Pressure (decibar)"'],errors='coerce')-10.1325

JF2C['Real Salinity (g/kg)'] = gsw.conversions.SA_from_SP(np.array(JF2C[' "Practical Salinity (psu)"']), np.array(JF2C["Sea Pressure (db)"]),np.array(JF2C['lon']),np.array(JF2C['lat']))


# In[63]:


# do we take the daily average of the mooring since we only have the daily average at that location for CIOPS?
# - makes sense to me!
df = JF2C[JF2C.dtUTC>=dt.datetime(2015,11,22)][JF2C.dtUTC<dt.datetime(2019,12,29)].resample('D',on='dtUTC').mean()#.loc[JF2C.dtUTC].fillna(method='ffill')
df


# In[34]:


dates= pd.date_range(dt.datetime(2015,11,22), periods=1498)


# In[33]:


fold = [dt.datetime(2015,11,22)+dt.timedelta(days=7*(i+1)) for i in range(int(214))]
fold = np.repeat(fold,7)


# In[45]:


potT = [float(xr.open_dataset("/ocean/mdunphy/CIOPSW-BC12/{:%Y%m%d}00/BC12_1d_grid_T_{:%Y%m%d}_{:%Y%m%d}.nc".format(fold[d], dates[d], dates[d]))['thetao'][0,k,j,i].values) for d in range(len(dates))]


# In[46]:


S = [float(xr.open_dataset("/ocean/mdunphy/CIOPSW-BC12/{:%Y%m%d}00/BC12_1d_grid_T_{:%Y%m%d}_{:%Y%m%d}.nc".format(fold[d], dates[d], dates[d]))['so'][0,k,j,i].values) for d in range(len(dates))]


# In[58]:


P = np.zeros(len(S))
P.fill(gdepw[k][0][0])


# In[59]:


# convert modelled from potential T to in-situ T
CT = gsw.conversions.CT_from_pt(S,potT)
T = gsw.conversions.t_from_CT(S, CT, P)


# In[89]:


# make a pandas datafram with your new fun info
d = {'T': T, 'S': S, 'P': P}
ciops = pd.DataFrame(data=d, index=dates)


# In[114]:


fig, axs = plt.subplots(1, 2, figsize = (12, 6))

axs[0].plot((32,34.6),(32,34.6),'k-',alpha=.2)
axs[1].plot((5.4,10),(5.4,10),'k-',alpha=.2)

iiT=(~np.isnan(df[' "Temperature (C)"']))&(~np.isnan(ciops['T']))
iiS=(~np.isnan(df['Real Salinity (g/kg)']))&(~np.isnan(ciops['S']))
counts, xedges, yedges, m1=axs[1].hist2d(df.loc[iiT,[' "Temperature (C)"']].values.flatten(),
                                      ciops.loc[iiT,['T']].values.flatten(),bins=30,norm=LogNorm())
counts, xedges, yedges, m2=axs[0].hist2d(df.loc[iiS,['Real Salinity (g/kg)']].values.flatten(),
                                      ciops.loc[iiS,['S']].values.flatten(),bins=30,norm=LogNorm())

cb0=fig.colorbar(m2,ax=axs[0])
cb0.set_label('Count')
cb1=fig.colorbar(m1,ax=axs[1])
# cb1.set_label('Count')

ntick=np.arange(32.2, 34.4, 0.4)
axs[0].set_xlim((32,34.6))
axs[0].set_ylim((32,34.6))
axs[0].set_xticks(ntick)
axs[0].set_yticks(ntick)
    
stick=np.arange(5.6,9.8, 0.6)
axs[1].set_xlim((5.4,10))
axs[1].set_ylim((5.4,10))
axs[1].set_xticks(stick)
axs[1].set_yticks(stick)
    
for ax in (axs[1],axs[0]):
    ax.set_aspect(1, adjustable='box')
    
axs[1].set_ylabel('Modeled',fontsize=12)
axs[0].set_ylabel('Modeled',fontsize=12)
axs[1].set_xlabel('Observed',fontsize=12)
axs[0].set_xlabel('Observed',fontsize=12)

axs[0].set_title('Salinity (g/kg)',fontsize=12)
axs[1].set_title('Temperature ($^{\circ}$C)',fontsize=12)

# plot the stats pannel
plot_panel(axs[0], df['Real Salinity (g/kg)'], ciops['S'], (32,34.6), 'g/kg')
plot_panel(axs[1], df[' "Temperature (C)"'], ciops['T'], (5.4,10), '$^{\circ}$C')

plt.savefig('CIOPSeval.png')


# In[ ]:




