# make average salinity files for the whole dang season

import xarray as xr
import numpy as np
import pandas as pd

# maybe having a dataframe would simplify all of this excessive messy code??
starttime = 1
endtime = 720
yearjumps = [0,1,-2,1,0,1,0,1,1,0,1,0]

transport = np.array([])
salinity_start = np.array([])
temp_start = np.array([])
salinity_end = np.array([])
temp_end = np.array([])
dic_start = np.array([])
nit_start = np.array([])
dic_end = np.array([])
nit_end = np.array([])
sil_start = np.array([])
sil_end = np.array([])
section = np.array([])
age = np.array([])
depth = np.array([])

year = [16,17]

for num in year:
    file = '/ocean/rbeutel/MOAD/analysis-becca/Ariane/1yr_runs/201905_1hr/forward_01jan'+str(num)+'/ariane_positions_quantitative.nc'
    mydata = xr.open_dataset(file)
    
    file = '/ocean/rbeutel/MOAD/analysis-becca/Ariane/1yr_runs/201905_Car_Sal/forward_01jan'+str(num)+'/ariane_positions_quantitative.nc'
    mydata_dic = xr.open_dataset(file)
    
    file = '/ocean/rbeutel/MOAD/analysis-becca/Ariane/1yr_runs/201905_Nit_Sil/forward_01jan'+str(num)+'/ariane_positions_quantitative.nc'
    mydata_nut = xr.open_dataset(file)
    for i in range(12):
        
        Xtrans = (mydata.final_transp[(np.isnan(mydata.final_section)==False)])
        Xsali_start = (mydata.init_salt[(np.isnan(mydata.final_section)==False)])
        Xsali_end = (mydata.final_salt[(np.isnan(mydata.final_section)==False)])
        Xtemp_start = (mydata.init_temp[(np.isnan(mydata.final_section)==False)])
        Xtemp_end = (mydata.final_temp[(np.isnan(mydata.final_section)==False)])
        Xdic_start = (mydata_dic.init_temp[(np.isnan(mydata.final_section)==False)])
        Xdic_end = (mydata_dic.final_temp[(np.isnan(mydata.final_section)==False)])
        Xnit_start = (mydata_nut.init_temp[(np.isnan(mydata.final_section)==False)])
        Xnit_end = (mydata_nut.final_temp[(np.isnan(mydata.final_section)==False)])
        Xsil_start = (mydata_nut.init_salt[(np.isnan(mydata.final_section)==False)])
        Xsil_end = (mydata_nut.final_salt[(np.isnan(mydata.final_section)==False)])
        Xage = np.abs((mydata.init_t[(np.isnan(mydata.final_section)==False)])-(mydata.final_t[(np.isnan(mydata.final_section)==False)]))
        Xsection = (mydata_dic.final_section[(np.isnan(mydata.final_section)==False)])
        Xdepth = (mydata_dic.final_depth[(np.isnan(mydata.final_section)==False)])
        
        
        transport = np.append(transport,Xtrans)
        salinity_start = np.append(salinity_start,Xsali_start)
        temp_start = np.append(temp_start,Xtemp_start)
        salinity_end = np.append(salinity_end,Xsali_end)
        temp_end = np.append(temp_end,Xtemp_end)
        dic_start = np.append(dic_start,Xdic_end)
        nit_start = np.append(nit_start,Xnit_start)
        dic_end = np.append(dic_end,Xdic_end)
        nit_end = np.append(nit_end,Xnit_end)
        sil_start = np.append(sil_start,Xsil_start)
        sil_end = np.append(sil_end,Xsil_end)
        section = np.append(section,Xsection)
        age = np.append(age,Xage)
        depth = np.append(depth,Xdepth)
        
#2019 and 2018 (both done month by month)
year = [18,19]
str_mo = ['jan', 'feb', 'mar', 'apr', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec']
for num in year:
    for i in range(len(str_mo)):
        
        file = '/ocean/rbeutel/MOAD/analysis-becca/Ariane/201905_1hr/forward_01'+str_mo[i]+str(num)+'/ariane_positions_quantitative.nc'
        mydata = xr.open_dataset(file)
        
        file = '/ocean/rbeutel/MOAD/analysis-becca/Ariane/201905_Car_Sal/forward_01'+str_mo[i]+str(num)+'/ariane_positions_quantitative.nc'
        mydata_dic = xr.open_dataset(file)
        
        file = '/ocean/rbeutel/MOAD/analysis-becca/Ariane/201905_Nit_Sil/forward_01'+str_mo[i]+str(num)+'/ariane_positions_quantitative.nc'
        mydata_nut = xr.open_dataset(file)
        
        Xtrans = (mydata.final_transp[(np.isnan(mydata.final_section)==False)])
        Xsali_start = (mydata.init_salt[(np.isnan(mydata.final_section)==False)])
        Xsali_end = (mydata.final_salt[(np.isnan(mydata.final_section)==False)])
        Xtemp_start = (mydata.init_temp[(np.isnan(mydata.final_section)==False)])
        Xtemp_end = (mydata.final_temp[(np.isnan(mydata.final_section)==False)])
        Xdic_start = (mydata_dic.init_temp[(np.isnan(mydata.final_section)==False)])
        Xdic_end = (mydata_dic.final_temp[(np.isnan(mydata.final_section)==False)])
        Xnit_start = (mydata_nut.init_temp[(np.isnan(mydata.final_section)==False)])
        Xnit_end = (mydata_nut.final_temp[(np.isnan(mydata.final_section)==False)])
        Xsil_start = (mydata_nut.init_salt[(np.isnan(mydata.final_section)==False)])
        Xsil_end = (mydata_nut.final_salt[(np.isnan(mydata.final_section)==False)])
        Xage = np.abs((mydata.init_t[(np.isnan(mydata.final_section)==False)])-(mydata.final_t[(np.isnan(mydata.final_section)==False)]))
        Xsection = (mydata_dic.final_section[(np.isnan(mydata.final_section)==False)])
        Xdepth = (mydata_dic.final_depth[(np.isnan(mydata.final_section)==False)])
        
        
        transport = np.append(transport,Xtrans)
        salinity_start = np.append(salinity_start,Xsali_start)
        temp_start = np.append(temp_start,Xtemp_start)
        salinity_end = np.append(salinity_end,Xsali_end)
        temp_end = np.append(temp_end,Xtemp_end)
        dic_start = np.append(dic_start,Xdic_end)
        nit_start = np.append(nit_start,Xnit_start)
        dic_end = np.append(dic_end,Xdic_end)
        nit_end = np.append(nit_end,Xnit_end)
        sil_start = np.append(sil_start,Xsil_start)
        sil_end = np.append(sil_end,Xsil_end)
        section = np.append(section,Xsection)
        age = np.append(age,Xage)
        depth = np.append(depth,Xdepth)
        
year = [20]

for num in year:
    file = '/ocean/rbeutel/MOAD/analysis-becca/Ariane/1yr_runs/201905_1hr/forward_01jan'+str(num)+'/ariane_positions_quantitative.nc'
    mydata = xr.open_dataset(file)
    
    file = '/ocean/rbeutel/MOAD/analysis-becca/Ariane/1yr_runs/201905_Car_Sal/forward_01jan'+str(num)+'/ariane_positions_quantitative.nc'
    mydata_dic = xr.open_dataset(file)
    
    file = '/ocean/rbeutel/MOAD/analysis-becca/Ariane/1yr_runs/201905_Nit_Sil/forward_01jan'+str(num)+'/ariane_positions_quantitative.nc'
    mydata_nut = xr.open_dataset(file)
    for i in range(12):
        
        Xtrans = (mydata.final_transp[(np.isnan(mydata.final_section)==False)])
        Xsali_start = (mydata.init_salt[(np.isnan(mydata.final_section)==False)])
        Xsali_end = (mydata.final_salt[(np.isnan(mydata.final_section)==False)])
        Xtemp_start = (mydata.init_temp[(np.isnan(mydata.final_section)==False)])
        Xtemp_end = (mydata.final_temp[(np.isnan(mydata.final_section)==False)])
        Xdic_start = (mydata_dic.init_temp[(np.isnan(mydata.final_section)==False)])
        Xdic_end = (mydata_dic.final_temp[(np.isnan(mydata.final_section)==False)])
        Xnit_start = (mydata_nut.init_temp[(np.isnan(mydata.final_section)==False)])
        Xnit_end = (mydata_nut.final_temp[(np.isnan(mydata.final_section)==False)])
        Xsil_start = (mydata_nut.init_salt[(np.isnan(mydata.final_section)==False)])
        Xsil_end = (mydata_nut.final_salt[(np.isnan(mydata.final_section)==False)])
        Xage = np.abs((mydata.init_t[(np.isnan(mydata.final_section)==False)])-(mydata.final_t[(np.isnan(mydata.final_section)==False)]))
        Xsection = (mydata_dic.final_section[(np.isnan(mydata.final_section)==False)])
        Xdepth = (mydata_dic.final_depth[(np.isnan(mydata.final_section)==False)])
        
        
        transport = np.append(transport,Xtrans)
        salinity_start = np.append(salinity_start,Xsali_start)
        temp_start = np.append(temp_start,Xtemp_start)
        salinity_end = np.append(salinity_end,Xsali_end)
        temp_end = np.append(temp_end,Xtemp_end)
        dic_start = np.append(dic_start,Xdic_end)
        nit_start = np.append(nit_start,Xnit_start)
        dic_end = np.append(dic_end,Xdic_end)
        nit_end = np.append(nit_end,Xnit_end)
        sil_start = np.append(sil_start,Xsil_start)
        sil_end = np.append(sil_end,Xsil_end)
        section = np.append(section,Xsection)
        age = np.append(age,Xage)
        depth = np.append(depth,Xdepth)
        
data = {'transport':transport, 'salinity_start':salinity_start, 'salinity_end':salinity_end, 'temp_start':temp_start, 'temp_end':temp_end,
       'dic_start':dic_start, 'dic_end':dic_end, 'nit_start':nit_start, 'nit_end':nit_end, 'sil_start':sil_start, 'sil_end':sil_end,
       'section':section, 'age':age, 'depth':depth}
dataframe = pd.DataFrame(data=data)

dataframe.to_csv('FiveYearProperties.csv')