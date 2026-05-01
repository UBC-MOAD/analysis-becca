# Alter and save boundary condition files based on localchange only
# a bit more complicated than the direct property change, this is based on a change of fresh versus south shelf contribution based on Fickin2016
# the average (based on Beutel2025) property differences between those water masses
# and the projected property changes in south shelf, noth shelf, offshore surface, and south brackish based on Siedlecki2021

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

def oxygen_mlL_to_uM(oxygen_mlL):
    """
    Convert dissolved oxygen from mL/L to uM using molar volume.
    """

    # in seawater the 1 ml/L of oxygen is approximately equal to 44.66 umol/L
    # this is based on the molar volume of an ideal gas as STP of 22.39 L/mol
    C = 44.66
    
    # Convert mL/L to µmol/L using the in-situ molar volume
    uM = oxygen_mlL * C
    
    return uM



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


# get the mean property differences from your lagrangian results
# difference = fresh source property - south source property
S_diff_fresh_south   = 30.04-32.86 # g/kg
T_diff_fresh_south   = 10.07-9.40  # degC
rho = gsw.rho(np.mean([30.04,32.86]), np.mean([10.07,9.40]), 30) # kg/m3
O_diff_fresh_south   = (292.97-200.42)*rho/1000  # mmol/m3 from umol/kg
TA_diff_fresh_south  = (2111.7-2213.8)*rho/1000  # mmol/m3 from umol/kg
DIC_diff_fresh_south = (1943.4-2116.9)*rho/1000  # mmol/m3 from umol/kg
N_diff_fresh_south   = 7.48-17.84  # mmol/m3

# based on the difference between the CMIP5 enssemble mean and the historical columbia river hydrograph in Fickin2016 (figure 5) this is our monthly scaling of flow:
monthly_factor = np.array([
    2.5,  # Jan
    2.6,  # Feb
    2.3,  # Mar
    1.4,  # Apr
    0.8,  # May
    0.5,  # Jun
    0.5,  # Jul
    0.7,  # Aug
    0.8,  # Sep
    0.9,  # Oct
    1.2,  # Nov
    2.2   # Dec
])

########
# LOOP #
########
# loop though changing the boundary one day at a time
for day in dates:

    # SETUP

    # load boundary file
    data = xr.open_dataset('/results/forcing/LiveOcean/boundary_conditions/LiveOcean_v201905_y{:%Y}m{:%m}d{:%d}.nc'.format(day,day,day))

    # the fraction of change depends on the month we are in
    frac_change = monthly_factor[day.month - 1]

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
    # for this case we are changing shelf and surface water properties
    fresh_prof = get_depth_profile(ffresh, date_str)
    south_prof = get_depth_profile(fsouth, date_str)
    north_prof = get_depth_profile(fnorth, date_str)
    offS_prof  = get_depth_profile(foffS, date_str)


    # brackish versus south shelf fraction
    # -----------------------------------------------------

    # change fraction
    requested_delta_fresh = (frac_change - 1.0) * fresh_prof

    # If fresh increases, cap by available south water.
    # If fresh decreases, cap by available fresh water.
    delta_fresh = np.where(
        requested_delta_fresh >= 0,
        np.minimum(requested_delta_fresh, south_prof),
        np.maximum(requested_delta_fresh, -fresh_prof),
    )

    # Updated fractions after exchange
    new_fresh_prof = fresh_prof + delta_fresh
    new_south_prof = south_prof - delta_fresh

    # property change based on fraction change
    delta_S_profile   = delta_fresh * S_diff_fresh_south
    delta_T_profile   = delta_fresh * T_diff_fresh_south
    delta_O_profile   = delta_fresh * O_diff_fresh_south
    delta_TA_profile  = delta_fresh * TA_diff_fresh_south
    delta_DIC_profile = delta_fresh * DIC_diff_fresh_south
    delta_N_profile   = delta_fresh * N_diff_fresh_south

    interps["salt"]       += delta_S_profile[:, None, None]
    interps["temp"]       += delta_T_profile[:, None, None]
    interps["oxygen"]     += delta_O_profile[:, None, None]
    interps["alkalinity"] += delta_TA_profile[:, None, None]
    interps["TIC"]        += delta_DIC_profile[:, None, None]
    interps["NO3"]        += delta_N_profile[:, None, None]


    # property change based on projections, using new south shelf and brackish fractions
    # ----------------------------------------------------------------------------------------------
    delta_T_profile = (
        new_fresh_prof * 2.82
        + offS_prof * 2.89
        + north_prof * 0.87
        + new_south_prof * 0.87
    )

    delta_O_profile = (
        new_fresh_prof * oxygen_mlL_to_uM(-0.39)
        + offS_prof * oxygen_mlL_to_uM(-0.4)
        + north_prof * oxygen_mlL_to_uM(-0.22)
        + new_south_prof * oxygen_mlL_to_uM(-0.22)
    )

    delta_TA_profile = (
        new_fresh_prof * -10.53
        + offS_prof * -7.46
        + north_prof * -4.05
        + new_south_prof * -4.05
    )

    delta_DIC_profile = (
        new_fresh_prof * 73.57
        + offS_prof * 80.98
        + north_prof * 75.93
        + new_south_prof * 75.93
    )

    delta_N_profile = (
        new_fresh_prof * -0.37 
        + offS_prof * -0.32
        + north_prof * 0
        + new_south_prof * 0
    )

    interps["temp"]       += delta_T_profile[:, None, None]
    interps["oxygen"]     += delta_O_profile[:, None, None]
    interps["alkalinity"] += delta_TA_profile[:, None, None]
    interps["TIC"]        += delta_DIC_profile[:, None, None]
    interps["NO3"]        += delta_N_profile[:, None, None]


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
    out.attrs["history"] = out.attrs.get("history", "") + " | edited for 2100, south local changes only"

    # save!
    out.to_netcdf("/ocean/rbeutel/MOAD/analysis-becca/projections/output/bdy_LocalChange/LiveOcean_v201905_y{:%Y}m{:%m}d{:%d}.nc".format(day,day,day))
    print(date_str)