# Alter and save boundary condition files based on upwelling change only
# a bit more complicated than the direct property change, this is based on a change of CUC versus shelf contribution
# and the average (based on Beutel2025) property differences between those water masses
# contribution differs month by month

# bring in packages
import datetime as dt
import pandas as pd
import xarray as xr
import numpy as np
import os
import gsw
from salishsea_tools.LiveOcean_BCs import convect, stabilize

#############
# FUNCTIONS #
#############
# get depth fraction profile of a certain water mass on a specific day 
def get_depth_profile(df, date_str):
    row = df.loc[df["day"] == date_str]
    prof = row.select_dtypes(include=np.number).to_numpy().squeeze()
    return prof  # shape (40,)


##########################
# ORGANIZE PRE-LOOP INFO # 
##########################

# want to do this for every day between 2018/01/01 and 2022/12/31
startday=dt.datetime(2018,1,1) 
endday=dt.datetime(2022,12,31)
dates = pd.date_range(start=startday,end=endday, freq="D")


# bring in water mass fraction data (already setup for every day in that span)
floop = pd.read_csv('./output/inflowfraction_D_loop.csv').drop(columns=['Unnamed: 0'])
ffresh = pd.read_csv('./output/inflowfraction_D_fresh.csv').drop(columns=['Unnamed: 0'])
fsouth = pd.read_csv('./output/inflowfraction_D_south.csv').drop(columns=['Unnamed: 0'])
fcuc = pd.read_csv('./output/inflowfraction_D_cuc.csv').drop(columns=['Unnamed: 0'])
foffD = pd.read_csv('./output/inflowfraction_D_offD.csv').drop(columns=['Unnamed: 0'])
foffS = pd.read_csv('./output/inflowfraction_D_offS.csv').drop(columns=['Unnamed: 0'])
fnorth = pd.read_csv('./output/inflowfraction_D_north.csv').drop(columns=['Unnamed: 0'])


# get the mean upwelling period property differences from your lagrangian results
S_diff_south = 33.90-32.97       # mean salinity difference between CUC and south shelf, g/kg
S_diff_north = 33.90-33.24       # mean salinity difference between CUC and north shelf, g/kg

T_diff_south = 6.25-8.92         # mean temperature difference between CUC and south shelf, degC
T_diff_north = 6.25-7.55         # mean temperature difference between CUC and north shelf, degC

O_diff_south = 62.26-200.76      # mean oxygen difference between CUC and south shelf, umol/kg
O_diff_north = 62.26-153.26      # mean oxygen difference between CUC and north shelf, umol/kg

TA_diff_south = 2262-2218        # mean alkalinity difference between CUC and south shelf, umol/kg
TA_diff_north = 2262-2227        # mean alkalinity difference between CUC and north shelf, umol/kg

DIC_diff_south = 2263-2124       # mean DIC difference between CUC and south shelf, umol/kg
DIC_diff_north = 2263-2167       # mean DIC difference between CUC and north shelf, umol/kg

N_diff_south = 37.27-19.09       # mean NO3 difference between CUC and south shelf, mmol/m3
N_diff_north = 37.27-24.26       # mean NO3 difference between CUC and north shelf, mmol/m3


# differences in units of umol/kg must be converted to uM to match the boundary condition files
# multiply by the density (kg/m3) and devide by 1000 (1000 L per m3)
# use the mean S and T from the two water masses, and take depth as 150 m because that's approximately where the mixing between the two will occur
rho = gsw.rho(np.mean([33.90,32.97,33.24]), np.mean([6.25,8.92,7.55]), 150) # kg/m3

O_diff_south = O_diff_south*rho/1000
O_diff_north = O_diff_north*rho/1000

TA_diff_south = TA_diff_south*rho/1000
TA_diff_north = TA_diff_north*rho/1000

DIC_diff_south = DIC_diff_south*rho/1000
DIC_diff_north = DIC_diff_north*rho/1000



########
# LOOP #
########
# loop though changing the boundary one day at a time
for day in dates:

    # SETUP

    # load boundary file
    data = xr.open_dataset('/results/forcing/LiveOcean/boundary_conditions/LiveOcean_v201905_y{:%Y}m{:%m}d{:%d}.nc'.format(day,day,day))

    # the fraction of change depends on the month we are in, based on Brady2017
    # skip process entirely if the change to upwelling is insignificant, save symbolic link for the sake of keeping things organized
    if day.month == 4: # april
        frac_change = 0.16
    elif day.month == 7: # july
        frac_change = -0.12
    elif day.month == 8: #august
        frac_change = -0.15
    elif day.month == 9: # september
        frac_change = -0.26
    elif day.month == 10: # october
        frac_change = -0.65
    else:
        os.symlink(os.path.abspath("/ocean/rbeutel/MOAD/analysis-becca/projections/output/bdy_UpwelllingChange/LiveOcean_v201905_y{:%Y}m{:%m}d{:%d}.nc".format(day,day,day)),
                   os.path.abspath("./output/bdy_UpwellingChange/LiveOcean_v201905_y{:%Y}m{:%m}d{:%d}.nc".format(day,day,day)))
        continue 

    # build an easily editable tracer dictionary 
    interps = {
        "salt":       data["vosaline"][0, :, :, :].values.copy(),   # (40, 1, 950)
        "temp":       data["votemper"][0, :, :, :].values.copy(),
        "NO3":        data["NO3"][0, :, :, :].values.copy(),
        "Si":         data["Si"][0, :, :, :].values.copy(),
        "oxygen":     data["OXY"][0, :, :, :].values.copy(),
        "TIC":        data["DIC"][0, :, :, :].values.copy(),
        "alkalinity": data["TA"][0, :, :, :].values.copy(),
    }

    # to use the profile function need to get date into the same format as the fraction dataframes
    date_str = "{:%Y}-{:%m}-{:%d}".format(day,day,day)

    # EDITS
    # for this case we are changing the CUC and south or north shelf properties
    south_prof = get_depth_profile(fsouth, date_str)
    north_prof = get_depth_profile(fnorth, date_str)
    cuc_prof   = get_depth_profile(fcuc, date_str)

    # change the CUC contribution by the fractional change
    requested_delta_f = frac_change * cuc_prof

    # balance out change with shelf water
    # If CUC increases, cap transfer by available shelf water.
    # If CUC decreases, cap transfer by available CUC water.
    if np.sum(south_prof)>np.sum(north_prof): # select transfer from south if larger than north
        delta_f = np.where(
            requested_delta_f >= 0,
            np.minimum(requested_delta_f, south_prof),
            np.maximum(requested_delta_f, -cuc_prof),
        )

        delta_S_profile = delta_f * S_diff_south
        delta_T_profile = delta_f * T_diff_south
        delta_O_profile = delta_f * O_diff_south
        delta_TA_profile = delta_f * TA_diff_south
        delta_DIC_profile = delta_f * DIC_diff_south
        delta_N_profile = delta_f * N_diff_south

    else:
        delta_f = np.where(
            requested_delta_f >= 0,
            np.minimum(requested_delta_f, north_prof),
            np.maximum(requested_delta_f, -cuc_prof),
        )

        delta_S_profile = delta_f * S_diff_north
        delta_T_profile = delta_f * T_diff_north
        delta_O_profile = delta_f * O_diff_north
        delta_TA_profile = delta_f * TA_diff_north
        delta_DIC_profile = delta_f * DIC_diff_north
        delta_N_profile = delta_f * N_diff_north

    interps["salt"] += delta_S_profile[:, None, None]
    interps["temp"] += delta_T_profile[:, None, None]
    interps["oxygen"] += delta_O_profile[:, None, None]
    interps["alkalinity"] += delta_TA_profile[:, None, None]
    interps["TIC"] += delta_DIC_profile[:, None, None]
    interps["NO3"] += delta_N_profile[:, None, None]


    # STABILIZE 

    # start by applying convect
    # this function simply convects vertically based on the density of cells on top of eachother
    sigma = gsw.sigma0(interps['salt'], interps['temp']) # need density first
    sigma, interps = convect(sigma, interps)

    # stabilize (small changes to salinity to stabilize marginally stable cells)
    interps = stabilize(sigma, interps)


    # SAVE EDITED FILE AS NETCDF
    out = data.copy(deep=True)

    out["vosaline"][0, :, :, :] = interps["salt"]
    out["votemper"][0, :, :, :] = interps["temp"]
    out["NO3"][0, :, :, :]      = interps["NO3"]
    out["Si"][0, :, :, :]       = interps["Si"]
    out["OXY"][0, :, :, :]      = interps["oxygen"]
    out["DIC"][0, :, :, :]      = interps["TIC"]
    out["TA"][0, :, :, :]       = interps["alkalinity"]

    # update metadata
    out.attrs["history"] = out.attrs.get("history", "") + " | edited for 2100, upwelling change only"

    # save!
    out.to_netcdf("/ocean/rbeutel/MOAD/analysis-becca/projections/output/bdy_UpwellingChange/LiveOcean_v201905_y{:%Y}m{:%m}d{:%d}.nc".format(day,day,day))
    print(date_str)