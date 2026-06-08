# taking changes applied to the JdF boundary and applying them to the Johstone bdy
# based on salinity bins

import argparse
import xarray as xr
import pandas as pd 
import numpy as np
import datetime as dt
from pathlib import Path
from salishsea_tools.LiveOcean_BCs import convect, stabilize
import gsw

#allow folder (ie. change scenario applied) to be input in the command line
parser = argparse.ArgumentParser(description='Get folder/scenario name.')
parser.add_argument('folder', type=str,
                    help='Folder to get changed JdF boundary info from')
args = parser.parse_args()

folder = args.folder
# for example: "bdy_AllChange"

#######################
###### FUNCTIONS ######
#######################

# get depth fraction profile of a certain water mass on a specific day 
def get_depth_profile(df, date_str):
    row = df.loc[df["day"] == date_str]
    prof = row.select_dtypes(include=np.number).to_numpy().squeeze()
    return prof  # shape (40,)

def salinity_binned_delta_lookup(S_source, delta_source, salinity_edges):
    """
    Compute mean source anomaly as a function of source salinity.

    S_source: JdF salinity array
    delta_source: JdF anomaly array for any variable
    salinity_edges: 1-D salinity bin edges

    Returns
    -------
    bin_centres : 1-D array
    mean_delta_filled : 1-D array of anomaly vs salinity
    """

    bin_centres = 0.5 * (salinity_edges[:-1] + salinity_edges[1:])
    mean_delta = np.full(len(bin_centres), np.nan)

    for i in range(len(bin_centres)):
        Smin = salinity_edges[i]
        Smax = salinity_edges[i + 1]

        mask = (S_source >= Smin) & (S_source < Smax)

        if np.any(mask):
            mean_delta[i] = np.nanmean(delta_source[mask])

    valid = np.isfinite(mean_delta)

    if valid.sum() < 2:
        raise ValueError("Not enough populated salinity bins to interpolate anomalies.")

    mean_delta_filled = np.interp(
        bin_centres,
        bin_centres[valid],
        mean_delta[valid],
    )

    return bin_centres, mean_delta_filled


def apply_salinity_binned_lookup(S_target, bin_centres, mean_delta_filled):
    """
    Apply salinity-binned anomaly lookup to a target boundary.

    S_target: Johnstone salinity array
    bin_centres : 1-D array
    mean_delta_filled : 1-D array of anomaly vs salinity

    Returns
    -------
    delta_target: anomaly field with same shape as S_target
    """

    delta_target = np.interp(
        S_target.ravel(),
        bin_centres,
        mean_delta_filled,
        left=mean_delta_filled[0],
        right=mean_delta_filled[-1],
    ).reshape(S_target.shape)

    return delta_target



###################
###### SETUP ######
###################

# Johnstone strait data
data = xr.open_dataset('/ocean/rbeutel/MOAD/analysis-becca/projections/output/Original_Johnstone_bdy.nc')
original_df = {
    "salt":       data["vosaline"][:, :, :, :].values.copy(),   # [12,40,10,30]
    "temp":       data["votemper"][:, :, :, :].values.copy(),
    "oxygen":     data["OXY"][:, :, :, :].values.copy(),
    "Si":         data["Si"][:, :, :, :].values.copy(),
    "NO3":        data["NO3"][:, :, :, :].values.copy(),
    "TIC":        data["DIC"][:, :, :, :].values.copy(),
    "alkalinity": data["TA"][:, :, :, :].values.copy(),
}

df = {k: v.copy() for k, v in original_df.items()}


# JdF every day between 2018/01/01 and 2022/12/31
startday=dt.datetime(2018,1,1) 
endday=dt.datetime(2022,12,31)
dates = pd.date_range(start=startday,end=endday, freq="D")
# not loading here, do once in loop

# salinity range for johnstone strait
salinity_edges = np.arange(np.min(data.vosaline.values),np.max(data.vosaline.values)+0.1, 0.1)

# salinity based lookup table separated by month
monthly_source = {
    month: {
        "S": [],
        "salt": [],
        "temp": [],
        "oxygen": [],
        "NO3": [],
        "alkalinity": [],
        "TIC": [],
    }
    for month in range(1, 13)
}



########################
###### DAILY LOOP ######
########################

# build monthly JdF change
for day in dates:
    month = day.month

    # bring in and organise original and edited JdF boundary
    _original_jdf = xr.open_dataset('/results/forcing/LiveOcean/boundary_conditions/LiveOcean_v201905_y{:%Y}m{:%m}d{:%d}.nc'.format(day,day,day))
    original_jdf = {
        "salt":       _original_jdf["vosaline"][0, :, :, :].values.copy(),   # (40, 1, 950)
        "temp":       _original_jdf["votemper"][0, :, :, :].values.copy(),
        "NO3":        _original_jdf["NO3"][0, :, :, :].values.copy(),
        "Si":         _original_jdf["Si"][0, :, :, :].values.copy(),
        "oxygen":     _original_jdf["OXY"][0, :, :, :].values.copy(),
        "TIC":        _original_jdf["DIC"][0, :, :, :].values.copy(),
        "alkalinity": _original_jdf["TA"][0, :, :, :].values.copy(),
    }

    a = Path("/ocean/rbeutel/MOAD/analysis-becca/projections/output/{}/LiveOcean_v201905_y{:%Y}m{:%m}d{:%d}.nc".format(folder,day,day,day))

    # Check if the file exists
    if a.exists():
        _interps_jdf = xr.open_dataset(a)
    else:
        print(day)
        S_jdf = original_jdf["salt"]
        monthly_source[month]["S"].append(S_jdf.ravel())

        # calculate monthly deltas
        monthly_source[month]["salt"].append(
            (np.zeros(np.shape(S_jdf.ravel())))
        )
        monthly_source[month]["temp"].append(
            (np.zeros(np.shape(S_jdf.ravel())))
        )
        monthly_source[month]["oxygen"].append(
            (np.zeros(np.shape(S_jdf.ravel())))
        )
        monthly_source[month]["NO3"].append(
            (np.zeros(np.shape(S_jdf.ravel())))
        )
        monthly_source[month]["alkalinity"].append(
            (np.zeros(np.shape(S_jdf.ravel())))
        )
        monthly_source[month]["TIC"].append(
            (np.zeros(np.shape(S_jdf.ravel())))
        )
        continue
    
    interps_jdf = {
        "salt":       _interps_jdf["vosaline"][0, :, :, :].values.copy(),   # (40, 1, 950)
        "temp":       _interps_jdf["votemper"][0, :, :, :].values.copy(),
        "NO3":        _interps_jdf["NO3"][0, :, :, :].values.copy(),
        "Si":         _interps_jdf["Si"][0, :, :, :].values.copy(),
        "oxygen":     _interps_jdf["OXY"][0, :, :, :].values.copy(),
        "TIC":        _interps_jdf["DIC"][0, :, :, :].values.copy(),
        "alkalinity": _interps_jdf["TA"][0, :, :, :].values.copy(),
    }
    _original_jdf.close()
    _interps_jdf.close()

    # save original salt for use in comparing to Johnstone salt
    S_jdf = original_jdf["salt"]
    monthly_source[month]["S"].append(S_jdf.ravel())

    # calculate monthly deltas
    monthly_source[month]["salt"].append(
        (interps_jdf["salt"] - original_jdf["salt"]).ravel()
    )
    monthly_source[month]["temp"].append(
        (interps_jdf["temp"] - original_jdf["temp"]).ravel()
    )
    monthly_source[month]["oxygen"].append(
        (interps_jdf["oxygen"] - original_jdf["oxygen"]).ravel()
    )
    monthly_source[month]["NO3"].append(
        (interps_jdf["NO3"] - original_jdf["NO3"]).ravel()
    )
    monthly_source[month]["alkalinity"].append(
        (interps_jdf["alkalinity"] - original_jdf["alkalinity"]).ravel()
    )
    monthly_source[month]["TIC"].append(
        (interps_jdf["TIC"] - original_jdf["TIC"]).ravel()
    )


# bin by salinities relevant to Johnstone strait
monthly_lookup = {}

for month in range(1, 13):

    S_all = np.concatenate(monthly_source[month]["S"])

    monthly_lookup[month] = {}

    for var in ["salt", "temp", "oxygen", "NO3", "alkalinity", "TIC"]:

        delta_all = np.concatenate(monthly_source[month][var])

        bin_centres, mean_delta_filled = salinity_binned_delta_lookup(
            S_all,
            delta_all,
            salinity_edges,
        )

        monthly_lookup[month][var] = {
            "bin_centres": bin_centres,
            "mean_delta": mean_delta_filled,
        }

# apply information from the lookup to each month of the Johnstone boundary
for month in np.arange(1,13):
    S_north = df["salt"][month-1].copy()

    df["salt"][month-1] += apply_salinity_binned_lookup(
        S_north,
        monthly_lookup[month]["salt"]["bin_centres"],
        monthly_lookup[month]["salt"]["mean_delta"],
    )

    df["temp"][month-1] += apply_salinity_binned_lookup(
        S_north,
        monthly_lookup[month]["temp"]["bin_centres"],
        monthly_lookup[month]["temp"]["mean_delta"],
    )

    df["oxygen"][month-1] += apply_salinity_binned_lookup(
        S_north,
        monthly_lookup[month]["oxygen"]["bin_centres"],
        monthly_lookup[month]["oxygen"]["mean_delta"],
    )

    df["NO3"][month-1] += apply_salinity_binned_lookup(
        S_north,
        monthly_lookup[month]["NO3"]["bin_centres"],
        monthly_lookup[month]["NO3"]["mean_delta"],
    )

    df["alkalinity"][month-1] += apply_salinity_binned_lookup(
        S_north,
        monthly_lookup[month]["alkalinity"]["bin_centres"],
        monthly_lookup[month]["alkalinity"]["mean_delta"],
    )

    df["TIC"][month-1] += apply_salinity_binned_lookup(
        S_north,
        monthly_lookup[month]["TIC"]["bin_centres"],
        monthly_lookup[month]["TIC"]["mean_delta"],
    )

#######################
###### STABILIZE ######
#######################

# unfortunately the stability functions can't handle all the months at once so need to loop through
for m in range(12):

    # build 3-D monthly dictionary
    df_m = {
        "salt":       df["salt"][m].copy(),
        "temp":       df["temp"][m].copy(),
        "oxygen":     df["oxygen"][m].copy(),
        "NO3":        df["NO3"][m].copy(),
        "TIC":        df["TIC"][m].copy(),
        "alkalinity": df["alkalinity"][m].copy(),
    }

    # density for this month only
    sigma = gsw.sigma0(df_m["salt"], df_m["temp"])

    # convect this monthly boundary
    sigma, df_m = convect(sigma, df_m)

    # recompute sigma after convect
    sigma = gsw.sigma0(df_m["salt"], df_m["temp"])

    # stabilize this monthly boundary
    df_m = stabilize(sigma, df_m)

    # write stabilized 3-D fields back into 4-D climatology
    df["salt"][m]       = df_m["salt"]
    df["temp"][m]       = df_m["temp"]
    df["oxygen"][m]     = df_m["oxygen"]
    df["NO3"][m]        = df_m["NO3"]
    df["TIC"][m]        = df_m["TIC"]
    df["alkalinity"][m] = df_m["alkalinity"]


##################
###### SAVE ######
##################

# SAVE EDITED FILE AS NETCDF
out = data.copy(deep=True)

out["vosaline"].values[...] = df["salt"]
out["votemper"].values[...] = df["temp"]
out["NO3"].values[...]      = df["NO3"]
out["Si"].values[...]       = df["Si"]
out["OXY"].values[...]      = df["oxygen"]
out["DIC"].values[...]      = df["TIC"]
out["TA"].values[...]       = df["alkalinity"]

# update metadata
out.attrs["history"] = out.attrs.get("history", "") + " | edited for 2100"

# save!
out.to_netcdf("/ocean/rbeutel/MOAD/analysis-becca/projections/output/{}/Johnstone_bdy.nc".format(folder),unlimited_dims=('time_counter'))