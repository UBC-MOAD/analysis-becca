# Combine all Pacific boundary
# this involves combining the fractioning from both upwelling and columbia river alterations
# the property changes are all wrapped into one number

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

# difference = CUC property - shelf source property
S_diff_south = 33.90-32.97       # g/kg
S_diff_north = 33.90-33.24       # g/kg
T_diff_south = 6.25-8.92         # degC
T_diff_north = 6.25-7.55         # degC
rho = gsw.rho(np.mean([33.90,32.97,33.24]), np.mean([6.25,8.92,7.55]), 150) # kg/m3
O_diff_south = (62.26-200.76)*rho/1000      # mmol/m3 from umol/kg
O_diff_north = (62.26-153.26)*rho/1000      # mmol/m3 from umol/kg
TA_diff_south = (2262-2218)*rho/1000        # mmol/m3 from umol/kg
TA_diff_north = (2262-2227)*rho/1000        # mmol/m3 from umol/kg
DIC_diff_south = (2263-2124)*rho/1000       # mmol/m3 from umol/kg
DIC_diff_north = (2263-2167)*rho/1000       # mmol/m3 from umol/kg
N_diff_south = 37.27-19.09       # mmol/m3
N_diff_north = 37.27-24.26       # mmol/m3


# columbia river flow factor based on the difference between the CMIP5 enssemble mean and the historical columbia river hydrograph in Fickin2016 (figure 5) this is our monthly scaling of flow:
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

    # bring in profiles for each source water (all are being edited)
    fresh_prof = get_depth_profile(ffresh, date_str)
    south_prof = get_depth_profile(fsouth, date_str)
    north_prof = get_depth_profile(fnorth, date_str)
    offS_prof  = get_depth_profile(foffS, date_str)
    offD_prof  = get_depth_profile(foffD, date_str)
    cuc_prof  = get_depth_profile(fcuc, date_str)


    # -----------------------
    # ------ Upwelling ------
    # -----------------------

    # the fraction of change depends on the month we are in
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
        frac_change = 0.0 # for all other months, no change

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

        # Updated fractions after exchange
        cuc_prof = cuc_prof + delta_f
        south_prof = south_prof - delta_f

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

        # Updated fractions after exchange
        cuc_prof = cuc_prof + delta_f
        north_prof = north_prof - delta_f


    interps["salt"] += delta_S_profile[:, None, None]
    interps["temp"] += delta_T_profile[:, None, None]
    interps["oxygen"] += delta_O_profile[:, None, None]
    interps["alkalinity"] += delta_TA_profile[:, None, None]
    interps["TIC"] += delta_DIC_profile[:, None, None]
    interps["NO3"] += delta_N_profile[:, None, None]


    # ----------------------------
    # ------ Columbia River ------
    # ----------------------------

    # the fraction of change depends on the month we are in
    frac_change = monthly_factor[day.month - 1]
    requested_delta_fresh = (frac_change - 1.0) * fresh_prof

    # If fresh increases, cap by available south water.
    # If fresh decreases, cap by available fresh water.
    delta_fresh = np.where(
        requested_delta_fresh >= 0,
        np.minimum(requested_delta_fresh, south_prof),
        np.maximum(requested_delta_fresh, -fresh_prof),
    )

    # Updated fractions after exchange
    fresh_prof = fresh_prof + delta_fresh
    south_prof = south_prof - delta_fresh

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


    # ------------------------------
    # ------ Other Properties ------
    # ------------------------------

    delta_T_profile = (
        fresh_prof * 2.82
        + offS_prof * 2.89
        + north_prof * 2.62
        + south_prof * 2.62
        + offD_prof * 2.05
        + cuc_prof * 2.05
    )

    delta_O_profile = (
        fresh_prof * oxygen_mlL_to_uM(-0.39)
        + offS_prof * oxygen_mlL_to_uM(-0.4)
        + north_prof * oxygen_mlL_to_uM(-0.5)
        + south_prof * oxygen_mlL_to_uM(-0.5) 
        + offD_prof * oxygen_mlL_to_uM(-0.6)
        + cuc_prof * oxygen_mlL_to_uM(-0.6)
    )

    delta_TA_profile = (
        fresh_prof * -10.53
        + offS_prof * -7.46
        + north_prof * -8.22
        + south_prof * -8.22
        + offD_prof * -1.66
        + cuc_prof * -1.66
    )

    delta_DIC_profile = (
        fresh_prof * 73.57
        + offS_prof * 80.98
        + north_prof * 80.93
        + south_prof * 80.93
        + offD_prof * 77.66
        + cuc_prof * 77.66
    )

    delta_N_profile = (
        fresh_prof * -0.37 
        + offS_prof * -0.32
        + north_prof * 0.84
        + south_prof * 0.84
        + offD_prof * 0.91
        + cuc_prof * 0.91
    )

    interps["temp"]       += delta_T_profile[:, None, None]
    interps["oxygen"]     += delta_O_profile[:, None, None]
    interps["alkalinity"] += delta_TA_profile[:, None, None]
    interps["TIC"]        += delta_DIC_profile[:, None, None]
    interps["NO3"]        += delta_N_profile[:, None, None]

    # ------------------------------------------------------------------------------------------------------------
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
    out.attrs["history"] = out.attrs.get("history", "") + " | edited for 2100, all Pacific boundary changes"

    # save!
    out.to_netcdf("/ocean/rbeutel/MOAD/analysis-becca/projections/output/bdy_AllChange/LiveOcean_v201905_y{:%Y}m{:%m}d{:%d}.nc".format(day,day,day))
    print(date_str)