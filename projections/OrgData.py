# lets get the ariane output into a nice formate for easy analysis of source waters at different inflow depths
import datetime as dt
import xarray as xr
import pandas as pd
import numpy as np
import math

# upwelling runs
upendday = [dt.datetime(2014, 9, 3), 
            dt.datetime(2015, 9, 5), dt.datetime(2016, 9, 13), 
            dt.datetime(2017, 10, 12), dt.datetime(2018, 9, 6), 
            dt.datetime(2019, 11, 5),dt.datetime(2020, 10, 17),
            dt.datetime(2021, 9, 22), dt.datetime(2022, 10, 15), 
            dt.datetime(2023, 9, 22)]
uplen = [161, 144, 158, 156,128,193,244,189, 117,136]  # Number of days in the season for each year
upstart = [upendday[i]-dt.timedelta(days=uplen[i]) for i in range(len(upendday))]  # Start dates

# downwelling runs
dwendday = [dt.datetime(2014, 3, 6), dt.datetime(2015, 2, 12), 
            dt.datetime(2016, 3, 19),dt.datetime(2017, 4, 19), 
            dt.datetime(2018, 2, 1),dt.datetime(2019, 4, 6), 
            dt.datetime(2020, 1, 27),dt.datetime(2021, 2, 2), 
            dt.datetime(2022, 1, 25),dt.datetime(2023, 4, 19)]
dwlen = [113, 143, 149, 198, 91, 163, 53, 82, 105,176]  # Number of days in the season for each year
dwstart = [dwendday[i]-dt.timedelta(days=dwlen[i]) for i in range(len(dwendday))]  # Start dates

# spring runs
spendday = [dt.datetime(2014, 3, 25), 
            dt.datetime(2015, 4, 13), dt.datetime(2016, 4, 7), 
            dt.datetime(2017, 5, 8), dt.datetime(2018, 4, 30), 
            dt.datetime(2019, 4, 25),dt.datetime(2020, 2, 15),
            dt.datetime(2021, 3, 18), dt.datetime(2022, 6, 11), 
            dt.datetime(2023, 5, 8)]
slen = [18,60,19,19,88,19,19,42,137,19]
spstart = [spendday[i]-dt.timedelta(days=slen[i]) for i in range(len(spendday))]  # Start dates

# fall runs
flendday = [dt.datetime(2014, 9, 21), 
            dt.datetime(2015, 10, 21), dt.datetime(2016, 10, 2), 
            dt.datetime(2017, 11, 1), dt.datetime(2018, 10, 24), 
            dt.datetime(2019, 12, 4),dt.datetime(2020, 11, 11),
            dt.datetime(2021, 10, 6), dt.datetime(2022, 10,24), 
            dt.datetime(2023, 10, 13)]
flen = [18,45,19,20,48,28,25,21,20,20]
fstart = [flendday[i]-dt.timedelta(days=flen[i]) for i in range(len(flendday))]  # Start dates


##################################################################################
# Combine data for a whole year instead of separating into upwelling/downwelling #
##################################################################################

def get_tracers_all(dwend,spend,upend,flend):
    #function to get the temperature, and salinity + section info for each parcel
    #south div is to set up boolean for the three south water masses, 1=CUC, 2=south shelf/davidson, 3=Columbia, 0=NA 

    section = np.array([])
    Idepth = np.array([])
    Fdepth = np.array([])
    lon = np.array([])
    lat = np.array([])
    start = np.array([])
    end = np.array([])
    date = np.array([])
    trans = np.array([])
    salt = np.array([])


    dw_st_files = ['/ocean/rbeutel/MOAD/biogeo_paper/FRDR/model/ariane/results/down_cas7/S_T/{:%Y%m%d}/ariane_positions_quantitative.nc'.format(day) for day in dwend]
    sp_st_files = ['/ocean/rbeutel/MOAD/biogeo_paper/FRDR/model/ariane/results/spring_cas7/S_T/{:%Y%m%d}/ariane_positions_quantitative.nc'.format(day) for day in spend]
    up_st_files = ['/ocean/rbeutel/MOAD/biogeo_paper/FRDR/model/ariane/results/up_cas7/S_T/{:%Y%m%d}/ariane_positions_quantitative.nc'.format(day) for day in upend]
    fl_st_files = ['/ocean/rbeutel/MOAD/biogeo_paper/FRDR/model/ariane/results/fall_cas7/S_T/{:%Y%m%d}/ariane_positions_quantitative.nc'.format(day) for day in flend]

    for i in range(len(upend)):
        print(i)
        dws_t = xr.open_dataset(dw_st_files[i])
        sps_t = xr.open_dataset(sp_st_files[i])
        ups_t = xr.open_dataset(up_st_files[i])
        fls_t = xr.open_dataset(fl_st_files[i])

        dwtides = ((abs(dws_t.init_t-dws_t.final_t) > 24) & ~np.isnan(dws_t.final_section)) # boolean to ignore tidally pumped and lost parcels
        sptides = ((abs(sps_t.init_t-sps_t.final_t) > 24) & ~np.isnan(sps_t.final_section))
        uptides = ((abs(ups_t.init_t-ups_t.final_t) > 24) & ~np.isnan(ups_t.final_section))
        fltides = ((abs(fls_t.init_t-fls_t.final_t) > 24) & ~np.isnan(fls_t.final_section))

        section = np.append(section,np.append(dws_t.final_section[dwtides],
                                              np.append(sps_t.final_section[sptides],
                                                        np.append(ups_t.final_section[uptides],fls_t.final_section[fltides]))))
        trans = np.append(trans,np.append(dws_t.final_transp[dwtides],
                                              np.append(sps_t.final_transp[sptides],
                                                        np.append(ups_t.final_transp[uptides],fls_t.final_transp[fltides]))))
        Idepth = np.append(Idepth,np.append(dws_t.init_depth[dwtides],
                                              np.append(sps_t.init_depth[sptides],
                                                        np.append(ups_t.init_depth[uptides],fls_t.init_depth[fltides]))))
        Fdepth = np.append(Fdepth,np.append(dws_t.final_depth[dwtides],
                                              np.append(sps_t.final_depth[sptides],
                                                        np.append(ups_t.final_depth[uptides],fls_t.final_depth[fltides]))))
        lon = np.append(lon,np.append(dws_t.final_lon[dwtides],
                                              np.append(sps_t.final_lon[sptides],
                                                        np.append(ups_t.final_lon[uptides],fls_t.final_lon[fltides]))))
        lat = np.append(lat,np.append(dws_t.final_lat[dwtides],
                                              np.append(sps_t.final_lat[sptides],
                                                        np.append(ups_t.final_lat[uptides],fls_t.final_lat[fltides]))))
        salt = np.append(salt,np.append(dws_t.final_salt[dwtides],
                                              np.append(sps_t.final_salt[sptides],
                                                        np.append(ups_t.final_salt[uptides],fls_t.final_salt[fltides]))))
        
        # print("before dates")
        # write out startdate (init_t - so when they're actually at JdF) of each parcel
        dw_dates = dwstart[i] + np.array([
            dt.timedelta(days=int(np.floor(x / 24)) - 100)
            for x in dws_t.init_t[dwtides].values
        ])

        sp_dates = spstart[i] + np.array([
            dt.timedelta(days=int(np.floor(x / 24)) - 100)
            for x in sps_t.init_t[sptides].values
        ])

        up_dates = upstart[i] + np.array([
            dt.timedelta(days=int(np.floor(x / 24)) - 100)
            for x in ups_t.init_t[uptides].values
        ])

        fl_dates = fstart[i] + np.array([
            dt.timedelta(days=int(np.floor(x / 24)) - 100)
            for x in fls_t.init_t[fltides].values
        ])
        date = np.append(date, np.append(dw_dates, np.append(sp_dates, np.append(up_dates, fl_dates))))
        # print("after dates")

        start = np.append(start,np.append(dws_t.init_t[dwtides],
                                              np.append(sps_t.init_t[sptides],
                                                        np.append(ups_t.init_t[uptides],fls_t.init_t[fltides]))))
        end = np.append(end,np.append(dws_t.final_t[dwtides],
                                              np.append(sps_t.final_t[sptides],
                                                        np.append(ups_t.final_t[uptides],fls_t.final_t[fltides]))))


    d = {'date':date,  'start':start,'end':end,'section':section, 'Idepth':Idepth,'Fdepth':Fdepth, 
         'lon':lon, 'lat':lat, 'transport':trans, 'salt':salt}
    df = pd.DataFrame(d)
    return df

all = get_tracers_all(dwendday,spendday,upendday,flendday)
all.to_csv('/ocean/rbeutel/MOAD/analysis-becca/projections/combineddata.csv')
print('done')