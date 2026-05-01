# Alter and save boundary condition files based on southern source water change only

# bring in packages
import datetime as dt
import pandas as pd
import xarray as xr
import gsw
from salishsea_tools.LiveOcean_BCs import convect, stabilize

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


#############
# FUNCTIONS #
#############
# get depth fraction profile of a certain water mass on a specific day 
def get_depth_profile(df, date_str):
    row = df.loc[df["day"] == date_str]
    prof = row.select_dtypes(include=np.number).to_numpy().squeeze()
    return prof  # shape (40,)

# convert oxygen change from Siedlecki2021 to SalishSeaCast units
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


    # EDITS

    # for this case we are only changing CUC and south shelf properties
    south_prof = get_depth_profile(fsouth, date_str)
    cuc_prof   = get_depth_profile(fcuc, date_str)

    # edit temperature
    for k in range(interps["temp"].shape[0]):
        delta = (
            south_prof[k] * 1.75
            + cuc_prof[k] * 1.03
        )
        interps["temp"][k, :, :] += delta

    # edit oxygen, remember unit conversion
    for k in range(interps["oxygen"].shape[0]):
        delta = (
            south_prof[k] * oxygen_mlL_to_uM(-0.28)
            + cuc_prof[k] * oxygen_mlL_to_uM(-0.3)
        )
        interps["oxygen"][k, :, :] += delta

    # edit TA
    for k in range(interps["alkalinity"].shape[0]):
        delta = (
            south_prof[k] * -4.17
            + cuc_prof[k] * -0.83
        )
        interps["alkalinity"][k, :, :] += delta

    # edit DIC
    for k in range(interps["TIC"].shape[0]):
        delta = (
            south_prof[k] * 5
            + cuc_prof[k] * 38.83
        )
        interps["TIC"][k, :, :] += delta

    # edit NO3
    for k in range(interps["NO3"].shape[0]):
        delta = (
            south_prof[k] * 0.84
            + cuc_prof[k] * 0.46
        )
        interps["NO3"][k, :, :] += delta


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
    out.attrs["history"] = out.attrs.get("history", "") + " | edited for 2100, south source water change only"

    # save!
    out.to_netcdf("/ocean/rbeutel/MOAD/analysis-becca/projections/output/bdy_SouthChange/LiveOcean_v201905_y{:%Y}m{:%m}d{:%d}.nc".format(day,day,day))
    print(date_str)