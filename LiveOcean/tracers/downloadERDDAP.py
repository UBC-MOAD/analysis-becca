#######################################
# Downloading tracer data from ERDDAP #
#######################################

# NOTE - this is NOT for model evaluation, we are just getting as many observations as we can

# to run:
# python3 downloadERDDAP.py organization_acronym measurement

import argparse
from erddapy import ERDDAP
import xarray as xr
import pandas as pd
import numpy as np
import gsw

from datetime import timedelta

#allow station id to be input in the command line
parser = argparse.ArgumentParser(description='Get ERDAPP data.')
parser.add_argument('org', type=str,
                    help='Acronym for the organization in lower case (onc, ios, ooi)')
# parser.add_argument('type', type=str,
#                     help='Type of measurement taken (ctd, oxy, bottle)')
args = parser.parse_args()

org = args.org
# type = args.type

#######
# ONC #
#######
if org == 'onc':

    # removing ones that are deeper than 500: 'scalar_1190542','scalar_1194359','scalar_1192794','scalar_1200145', 'scalar_1206556', 'scalar_1209192', 'scalar_1210110', 'scalar_1200128', 'scalar_1207539',	
    # 'scalar_1211656','scalar_1214292',

    IDs = ['scalar_1196717', 'scalar_1200697', 'scalar_1189272', 'scalar_1192789', 'scalar_1200624', 'scalar_1206551', 'scalar_1211690',
    'scalar_1215572','scalar_117206','scalar_117641','scalar_118186', 'scalar_118900', 'scalar_119070', 'scalar_1189343', 'scalar_1190578', 'scalar_1194377', 'scalar_1200638',	
    'scalar_1210115','scalar_117354','scalar_117686','profile_118840', 'scalar_118809', 'scalar_1189201', 'profile_1190483', 'scalar_1190473', 'scalar_1192792', 'scalar_1194404',	
    'profile_1200310', 'scalar_1200612', 'scalar_1202976', 'scalar_1205219', 'profile_1209160', 'scalar_1209170', 'scalar_1214273', 'scalar_1203269', 'scalar_1203272', 'scalar_1203275',
    'scalar_1205756', 'scalar_1205780', 'scalar_1205786', 'scalar_1206033', 'scalar_1209404', 'scalar_1209815', 'scalar_1209821', 'scalar_1209825', 'scalar_1210591', 'scalar_1210595',
    'scalar_1210608', 'scalar_1210979', 'scalar_1210987', 'scalar_1210995', 'scalar_1211254', 'scalar_1211259', 'scalar_1211526', 'scalar_1211531', 'scalar_1211537', 'scalar_1212183',
    'scalar_1212187', 'scalar_1212194', 'scalar_1212237', 'scalar_1212254', 'scalar_1212274', 'scalar_1212816', 'scalar_1212823', 'scalar_1212833', 'scalar_1213400', 'scalar_1213483',
    'scalar_1214120', 'scalar_1214127', 'scalar_1215222', 'scalar_1215232', 'scalar_118135', 'scalar_1189647', 'scalar_1192652', 'scalar_1194530', 'scalar_1196674', 'scalar_1200300',
    'scalar_1207525', 'scalar_1213677', 'scalar_1215677', 'scalar_1211295', 'scalar_1211300', 'scalar_1211285', 'scalar_1211905', 'scalar_1214536', 'scalar_1211290', 'scalar_1193600',
    'scalar_1193604', 'scalar_1196090', 'scalar_1196120', 'scalar_1201632', 'scalar_1205254', 'scalar_1209589', 'scalar_1210442', 'scalar_1212040', 'scalar_1214399', 'scalar_1205971',
    'scalar_1206537', 'scalar_1210344', 'scalar_1214374', 'scalar_1200374', 'scalar_1194675', 'scalar_1195535', 'scalar_1195536', 'scalar_1198139', 'scalar_1201774', 'scalar_1207706',
    'scalar_1211013', 'scalar_1215125', 'scalar_1193580', 'scalar_1196049', 'scalar_1196050',	'scalar_1201558', 'scalar_1207155', 'scalar_1208296',	'scalar_1212013', 'scalar_1214324',
    'scalar_1214355', 'scalar_1211668', 'scalar_117788', 'scalar_117790', 'scalar_117792', 'scalar_117794', 'scalar_118279', 'scalar_118281', 'scalar_118283', 'scalar_118285',
    'scalar_119011', 'scalar_119013', 'scalar_119015', 'scalar_119017', 'scalar_1200219', 'scalar_1200224', 'scalar_1200226', 'scalar_1200228','scalar_1200246', 'scalar_1200248',
    'scalar_1200250', 'scalar_1200253', 'scalar_1210222', 'scalar_1210225', 'scalar_1210226', 'scalar_1210227', 'scalar_1210314', 'scalar_1210317','scalar_1210322', 'scalar_1210323',
    'scalar_1215625', 'scalar_119027', 'scalar_119029', 'scalar_119031', 'scalar_119033', 'scalar_119036', 'scalar_119038', 'scalar_119040', 'scalar_119042', 'scalar_1202958',
    'scalar_1202960', 'scalar_1202962', 'scalar_1202964', 'scalar_1205520', 'scalar_1205522', 'scalar_1205524', 'scalar_1205526', 'scalar_1213746','scalar_1213748', 'scalar_117491',
    'scalar_118213', 'scalar_118947', 'scalar_1189320', 'scalar_1190513', 'scalar_1192588',	'scalar_1196277', 'scalar_1200329', 'scalar_1206517','scalar_1213661', 'scalar_1214337',
    'scalar_119211', 'scalar_1189528', 'scalar_1190903', 'scalar_1192239', 'scalar_1193781', 'scalar_1203388', 'scalar_1209305','scalar_1211396', 'scalar_1213346', 'scalar_119330',
    'scalar_1189541', 'scalar_1189892', 'scalar_1190727', 'scalar_1191182', 'scalar_117406', 'scalar_118131', 'scalar_1189474','scalar_1190647', 'scalar_1196755', 'scalar_1208038',
    'scalar_1213768', 'scalar_1195410', 'scalar_1197448', 'scalar_1201585', 'scalar_1209616', 'scalar_1212055', 'scalar_1215309', 'scalar_117767', 'scalar_118870', 'scalar_1210082',
    'scalar_117498']

    d = { "time (UTC)": np.array([]),
          "salinity (1e-3)": np.array([]),
          "temperature (K)": np.array([]),
          "pressure (dbar)": np.array([]),
          "longitude (degrees_east)": np.array([]),
          "latitude (degrees_north)": np.array([]),
          "depth (m)": np.array([])}
    df = pd.DataFrame(d)
    df.set_index('time (UTC)',inplace=True)
    df_TP = df
    df_tp = df

    # upper case T and P
    for id in IDs:
        print(f"Processing ID: {id}")
        e = ERDDAP(
            server="http://dap.onc.uvic.ca/erddap",
            protocol="tabledap",
            response="nc"
        )
        e.dataset_id = id
        e.constraints = {'depth<=': 500}

        try:
            # Retrieve metadata for variables that have a standard_name attribute
            metadata = e.get_var_by_attr(dataset_id=id, standard_name=lambda v: v is not None)

            # Initialize variable names
            temperature_name = None
            pressure_name = None
            
            # Search through metadata to check for correct variable names for temperature and salinity
            for var in metadata:
                if 'Temperature' in var:
                    temperature_name = var
                if 'Pressure' in var:
                    pressure_name = var

            if not temperature_name or not pressure_name:
                print(f"Required variables missing in dataset {id}, skipping...")
                continue

            e.variables = [
                "time",
                "salinity",
                "Temperature",
                "Pressure",
                "longitude",
                "latitude",
                "depth"
            ]

            # Download the data
            d = e.to_pandas()
            d.rename(columns={temperature_name: 'temperature (K)', pressure_name: 'pressure (dbar)'}, inplace=True)
            d['temperature (K)'] = pd.to_numeric(d['temperature (K)'], errors='coerce')
            d['pressure (dbar)'] = pd.to_numeric(d['pressure (dbar)'], errors='coerce')

            # resample such that its hourly (time resolution much higher than necessary for our purposes in ONC data)
            d['time (UTC)'] = pd.to_datetime(d['time (UTC)'])
            d.set_index('time (UTC)',inplace=True)
            d = d.resample('H').mean()

            if not d.empty:
                df_TP = pd.concat([df_TP, d])
                print(f"Data added for ID: {id}")
            else:
                print(f"No data found within constraints for ID: {id}")

        except Exception as ex:
            print(f"Failed to download data for ID: {id}. Error: {ex}")

    # lower case T and P
    for id in IDs:
        print(f"Processing ID: {id}")
        e = ERDDAP(
            server="http://dap.onc.uvic.ca/erddap",
            protocol="tabledap",
            response="nc"
        )
        e.dataset_id = id
        e.constraints = {'depth<=': 500}

        try:
            # Retrieve metadata for variables that have a standard_name attribute
            metadata = e.get_var_by_attr(dataset_id=id, standard_name=lambda v: v is not None)

            # Initialize variable names
            temperature_name = None
            pressure_name = None
            
            # Search through metadata to check for correct variable names for temperature and salinity
            for var in metadata:
                if 'temperature' in var:
                    temperature_name = var
                if 'pressure' in var:
                    pressure_name = var

            if not temperature_name or not pressure_name:
                print(f"Required variables missing in dataset {id}, skipping...")
                continue

            e.variables = [
                "time",
                "salinity",
                "temperature",
                "pressure",
                "longitude",
                "latitude",
                "depth"
            ]

            # Download the data
            d = e.to_pandas(index_col='time', parse_dates=['time'])
            print(d.columns)
            d['temperature (K)'] = pd.to_numeric(d['temperature (K)'], errors='coerce')
            d['pressure (dbar)'] = pd.to_numeric(d['pressure (dbar)'], errors='coerce')

            # resample such that its hourly (time resolution much higher than necessary for our purposes in ONC data)
            d['time (UTC)'] = pd.to_datetime(d['time (UTC)'])
            d.set_index('time (UTC)',inplace=True)
            d = d.resample('H').mean()

            if not d.empty:
                df_tp = pd.concat([df_tp, d])
                print(f"Data added for ID: {id}")
            else:
                print(f"No data found within constraints for ID: {id}")

        except Exception as ex:
            print(f"Failed to download data for ID: {id}. Error: {ex}")

    # Final DataFrame processing
    print("erddap worked")
    df = pd.concat[(df_tp, df_TP)]
    df['datetime'] = np.array(df.index)
    index = pd.Index(range(len(df)))
    df.set_index(index,inplace=True)

    name = 'ONC.p'
    df.to_pickle(name)
    print('ONC data download and processing complete.')


#######
# IOS #
#######
if org == 'ios':

    # first the mooring data
    df = pd.DataFrame()
    # download obs from ERDDAP
    e = ERDDAP(
    server="https://data.cioospacific.ca/erddap",
    protocol="tabledap",
    )

    e.response = "nc"
    e.dataset_id = "IOS_CTD_Moorings"
    e.constraints = {'depth<=':500} 

    e.variables = [  
        "time",
        "longitude",
        "latitude",
        "sea_water_pressure",
        "depth",
        "sea_water_temperature",
        "sea_water_practical_salinity",
        "DOXYZZ01"
    ]

    df = e.to_pandas()
    print('erddap worked')

    # lets be careful with the resampling since many data points (at different locations) may have been taken in the same hour
    df['time (UTC)'] = pd.to_datetime(df['time (UTC)'])
    df.set_index('time (UTC)',inplace=True)
    df = df.groupby(['latitude (degrees_north)', 'longitude (degrees_east)', 'depth (m)']).resample('H').mean()
    df.reset_index(inplace=True)

    name = 'IOS_ctd_moor.p'
    df.to_pickle(name)

    
    
    # second, the profile data
    df = pd.DataFrame()
    # download obs from ERDDAP
    e = ERDDAP(
    server="https://data.cioospacific.ca/erddap",
    protocol="tabledap",
    )

    e.response = "nc"
    e.dataset_id = "IOS_CTD_Profiles"
    e.constraints = {'depth<=':500} 

    e.variables = [  
        "time",
        "longitude",
        "latitude",
        "sea_water_pressure",
        "depth",
        "sea_water_temperature",
        "sea_water_practical_salinity",
        "DOXYZZ01"
    ]

    df = e.to_pandas()
    df['time (UTC)'] = pd.to_datetime(df['time (UTC)'])
    df.set_index('time (UTC)',inplace=True)
    df = df.groupby(['latitude (degrees_north)', 'longitude (degrees_east)', 'depth (m)']).resample('H').mean()
    df.reset_index(inplace=True)

    print('erddap worked')

    name = 'IOS_ctd_prof.p'
    df.to_pickle(name)

#######
# OOI #
#######
if org == 'ooi':

    #** Nitrate **#
    IDs =['ooi-ce01issm-rid16-07-nutnrb000','ooi-ce04ossm-rid26-07-nutnrb000','ooi-ce02shsm-rid26-07-nutnrb000','ooi-ce06issm-rid16-07-nutnrb000','ooi-ce09ossm-rid26-07-nutnrb000','ooi-ce07shsm-rid26-07-nutnrb000']

    d = {"time (UTC)":np.array([]),
    "z (m)":np.array([]),
    "latitude (degrees_north)":np.array([]),
    "longitude (degrees_east)":np.array([]),
    "mole_concentration_of_nitrate_in_sea_water (micromol.L-1)":np.array([]),
    "mole_concentration_of_nitrate_in_sea_water_suna_qc_agg":np.array([])}
    df = pd.DataFrame(d)
    df.set_index('time (UTC)',inplace=True)

    for id in IDs:
        e = ERDDAP(
        server="http://erddap.dataexplorer.oceanobservatories.org/erddap",
        protocol="tabledap",
        )
        e.response = "nc"
        e.dataset_id = id
        e.constraints = None
        e.variables = [  
            "time",
            "z",
            "latitude",
            "longitude",
            "mole_concentration_of_nitrate_in_sea_water",
            "mole_concentration_of_nitrate_in_sea_water_suna_qc_agg"
        ]
        d = e.to_pandas()

        # convert to hourly here so that data from different locations isn't accidentally grouped
        d['time (UTC)'] = pd.to_datetime(d['time (UTC)'])
        d.set_index('time (UTC)',inplace=True)
        d = d.resample('h',axis=0).mean()

        df = pd.concat([df,d])

    print('nitrate erddap worked')
    
    # reset index
    df['datetime'] = np.array(df.index)
    index = pd.Index(range(len(df)))
    df.set_index(index,inplace=True)

    name = 'OOI_nitrate.p'
    df.to_pickle(name)

    #** Chlorophyll **#
    IDs =['ooi-ce01issm-rid16-02-flortd000','ooi-ce01issm-sbd17-06-flortd000','ooi-ce04ospd-dp01b-04-flntua103','ooi-ce04ossm-rid27-02-flortd000','ooi-ce02shsm-rid27-02-flortd000',
          'ooi-ce06issm-rid16-02-flortd000','ooi-ce06issm-sbd17-06-flortd000','ooi-ce09ossm-rid27-02-flortd000','ooi-ce07shsm-rid27-02-flortd000','ooi-gp03flma-ris01-05-flortd000',
          'ooi-gp03flmb-ris01-05-flortd000','ooi-rs01sbpd-dp01a-04-flntua102']#,'ooi-rs03axps-pc03a-4c-flordd303','ooi-rs01sbps-pc01a-4c-flordd103']

    d = { "time (UTC)":np.array([]),
    "z (m)":np.array([]),
    "latitude (degrees_north)":np.array([]),
    "longitude (degrees_east)":np.array([]),
    "mass_concentration_of_chlorophyll_a_in_sea_water (microg.L-1)":np.array([]),
    "sea_water_temperature (degree_Celsius)":np.array([]),
    "sea_water_practical_salinity (1e-3)":np.array([]),
    "mass_concentration_of_chlorophyll_a_in_sea_water_qc_agg":np.array([])}
    df = pd.DataFrame(d)
    df.set_index('time (UTC)',inplace=True)


    for id in IDs:
        # print(id)
        e = ERDDAP(
        server="http://erddap.dataexplorer.oceanobservatories.org/erddap",
        protocol="tabledap",
        )
        e.response = "nc"
        e.dataset_id = id
        e.constraints = None
        e.variables = [  
            "time",
            "z",
            "latitude",
            "longitude",
            "mass_concentration_of_chlorophyll_a_in_sea_water",
            "mass_concentration_of_chlorophyll_a_in_sea_water_qc_agg",
            "sea_water_temperature",
            "sea_water_practical_salinity"
        ]
        d = e.to_pandas()

        # convert to hourly here so that data from different locations isn't accidentally grouped
        d['time (UTC)'] = pd.to_datetime(d['time (UTC)'])
        d.set_index('time (UTC)',inplace=True)
        d = d.resample('h',axis=0).mean()

        df = pd.concat([df,d])

    print('chlorophyll erddap worked')
    
    # reset index
    df['datetime'] = np.array(df.index)
    index = pd.Index(range(len(df)))
    df.set_index(index,inplace=True)

    name = 'OOI_chlorophyll.p'
    df.to_pickle(name)



    #** Oxygen **#
    # some datasets have different naming conventions.... 
    # trying to get them all at once  is really not working so download separately 
    
    IDs=['ooi-ce01issm-rid16-03-dostad000','ooi-ce01issm-mfd37-03-dostad000','ooi-ce04osbp-lj01c-06-ctdbpo108','ooi-ce04osbp-lj01c-06-dostad108','ooi-ce04osps-pc01b-4a-ctdpfa109','ooi-ce04osps-pc01b-4a-dostad109','ooi-ce04ossm-rid27-04-dostad000','ooi-ce02shbp-lj01d-06-ctdbpn106','ooi-ce02shbp-lj01d-06-dostad106','ooi-ce02shsm-rid27-04-dostad000','ooi-ce06issm-rid16-03-dostad000','ooi-ce06issm-mfd37-03-dostad000','ooi-ce09ossm-rid27-04-dostad000','ooi-ce09ossm-mfd37-03-dostad000','ooi-ce07shsm-rid27-04-dostad000','ooi-ce07shsm-mfd37-03-dostad000','ooi-gp03flma-ris01-03-dostad000','ooi-gp03flmb-ris01-03-dostad000','ooi-rs03ashs-mj03b-10-ctdpfb304','ooi-rs03axpd-dp03a-06-dostad304','ooi-rs03axbs-lj03a-12-ctdpfb301','ooi-rs03axbs-lj03a-12-dostad301','ooi-rs03axps-pc03a-4a-ctdpfa303','ooi-rs03axps-pc03a-4a-dostad303','ooi-rs03ccal-mj03f-12-ctdpfb305','ooi-rs03ecal-mj03e-12-ctdpfb306','ooi-rs01sbpd-dp01a-06-dostad104','ooi-rs01slbs-lj01a-12-ctdpfb101','ooi-rs01slbs-lj01a-12-dostad101','ooi-rs01sbps-pc01a-4a-ctdpfa103','ooi-rs01sbps-pc01a-4a-dostad103']

    # Define a base DataFrame to store all the data
    d = {
        "time (UTC)": np.array([]),
        "z (m)": np.array([]),
        "latitude (degrees_north)": np.array([]),
        "longitude (degrees_east)": np.array([]),
        "mole_concentration_of_dissolved_molecular_oxygen_in_sea_water (micromol.L-1)": np.array([]),
        "mole_concentration_of_dissolved_molecular_oxygen_in_sea_water_qc_agg":np.array([]),
        "sea_water_temperature (degree_Celsius)": np.array([]),
        "sea_water_practical_salinity (1e-3)": np.array([])
    }
    df = pd.DataFrame(d)
    df.set_index('time (UTC)',inplace=True)
    df_moor = df
    df_prof = df

    ## Mooring data first
    for id in IDs:
        print(f"Processing ID: {id}")
        e = ERDDAP(
            server="http://erddap.dataexplorer.oceanobservatories.org/erddap",
            protocol="tabledap"
        )
        e.response = "nc"
        e.dataset_id = id
        e.constraints = None

        try:
            # Retrieve metadata for variables that have a standard_name attribute
            metadata = e.get_var_by_attr(dataset_id=id, standard_name=lambda v: v is not None)

            # Initialize variable names
            temperature_name = None
            sainity_name = None
            
            # Search through metadata to check for correct variable names for temperature and salinity
            for var in metadata:
                if 'sea_water_temperature' in var:
                    temperature_name = var
                if 'sea_water_practical_salinity' in var:
                    salinity_name = var

            if not temperature_name or not salinity_name:
                print(f"Required variables missing in dataset {id}, skipping...")
                continue

            # Configure variables for download
            e.variables = [
                "time",
                "z",
                "latitude",
                "longitude",
                'mole_concentration_of_dissolved_molecular_oxygen_in_sea_water',
                'mole_concentration_of_dissolved_molecular_oxygen_in_sea_water_qc_agg',
                'sea_water_temperature',
                'sea_water_practical_salinity'
            ]
                
            # Download the data
            d = e.to_pandas()

            # convert to hourly here so that data from different locations isn't accidentally grouped
            d['time (UTC)'] = pd.to_datetime(d['time (UTC)'])
            d.set_index('time (UTC)',inplace=True)
            d = d.resample('h',axis=0).mean()

            # Check and convert data types
            d['sea_water_temperature (degree_Celsius)'] = pd.to_numeric(d['sea_water_temperature (degree_Celsius)'], errors='coerce')
            d['sea_water_practical_salinity (1e-3)'] = pd.to_numeric(d['sea_water_practical_salinity (1e-3)'], errors='coerce')

            df_moor = pd.concat([df_moor, d])

        except Exception as ex:
            print(f"Failed to download data for ID: {id}. Error: {ex}")


    # and the profile data
    for id in IDs:
        print(f"Processing ID: {id}")
        e = ERDDAP(
            server="http://erddap.dataexplorer.oceanobservatories.org/erddap",
            protocol="tabledap"
        )

        e.response = "nc"
        e.dataset_id = id
        e.constraints = None

        try:
            # Retrieve metadata for variables that have a standard_name attribute
            metadata = e.get_var_by_attr(dataset_id=id, standard_name=lambda v: v is not None)

            # Initialize variable names
            temperature_name = None
            sainity_name = None
            
            # Search through metadata to check for correct variable names for temperature and salinity
            for var in metadata:
                if 'sea_water_temperature_profiler_depth_enabled' in var:
                    temperature_name = var
                if 'sea_water_practical_salinity_profiler_depth_enabled' in var:
                    salinity_name = var

            if not temperature_name or not salinity_name:
                print(f"Required variables missing in dataset {id}, skipping...")
                continue

            # Attempt to fetch variable metadata
            # Configure variables for download
            e.variables = [
                "time",
                "z",
                "latitude",
                "longitude",
                'mole_concentration_of_dissolved_molecular_oxygen_in_sea_water',
                'mole_concentration_of_dissolved_molecular_oxygen_in_sea_water_qc_agg',
                'sea_water_temperature_profiler_depth_enabled',
                'sea_water_practical_salinity_profiler_depth_enabled'
            ]
            
            # Download the data
            d = e.to_pandas()

            # convert to hourly here so that data from different locations isn't accidentally grouped
            d['time (UTC)'] = pd.to_datetime(d['time (UTC)'])
            d.set_index('time (UTC)',inplace=True)
            d = d.resample('h',axis=0).mean()

            # Check and convert data types
            d['sea_water_temperature_profiler_depth_enabled (degree_Celsius)'] = pd.to_numeric(d['sea_water_temperature_profiler_depth_enabled (degree_Celsius)'], errors='coerce')
            d['sea_water_practical_salinity_profiler_depth_enabled (1e-3)'] = pd.to_numeric(d['sea_water_practical_salinity_profiler_depth_enabled (1e-3)'], errors='coerce')
            d.rename(columns={'sea_water_temperature_profiler_depth_enabled (degree_Celsius)': 'sea_water_temperature (degree_Celsius)', 'sea_water_practical_salinity_profiler_depth_enabled (1e-3)': 'sea_water_practical_salinity (1e-3)'}, inplace=True)

            df_prof = pd.concat([df_prof, d])
        
        except Exception as ex:
            print(f"Failed to download data for ID: {id}. Error: {ex}")

    print('oxygen erddap worked')
    df = pd.concat([df_moor, df_prof])
    
    # reset index
    df['datetime'] = np.array(df.index)
    index = pd.Index(range(len(df)))
    df.set_index(index,inplace=True)

    name = 'OOI_oxygen.p'
    df.to_pickle(name)


    #** Carbon Dioxide **# 
    IDs=['ooi-ce01issm-rid16-05-pco2wb000', 'ooi-ce01issm-mfd35-05-pco2wb000', 'ooi-ce04osbp-lj01c-09-pco2wb104', 'ooi-ce04osps-pc01b-4d-pco2wa105', 'ooi-ce02shbp-lj01d-09-pco2wb103', 'ooi-ce06issm-rid16-05-pco2wb000', 'ooi-ce06issm-mfd35-05-pco2wb000', 'ooi-ce09ossm-mfd35-05-pco2wb000', 'ooi-ce07shsm-mfd35-05-pco2wb000']

    d = { "time (UTC)":np.array([]),
    "z (m)":np.array([]),
    "latitude (degrees_north)":np.array([]),
    "longitude (degrees_east)":np.array([]),
    "partial_pressure_of_carbon_dioxide_in_sea_water (microatm)":np.array([]),
    "partial_pressure_of_carbon_dioxide_in_sea_water_qc_agg":np.array([])}
    df = pd.DataFrame(d)
    df.set_index('time (UTC)',inplace=True)

    for id in IDs:
        e = ERDDAP(
        server="http://erddap.dataexplorer.oceanobservatories.org/erddap",
        protocol="tabledap",
        )
        e.response = "nc"
        e.dataset_id = id
        e.constraints = None
        e.variables = [  
            "time",
            "z",
            "latitude",
            "longitude",
            "partial_pressure_of_carbon_dioxide_in_sea_water",
            "partial_pressure_of_carbon_dioxide_in_sea_water_qc_agg"
        ]
        d = e.to_pandas()

        # convert to hourly here so that data from different locations isn't accidentally grouped
        d['time (UTC)'] = pd.to_datetime(d['time (UTC)'])
        d.set_index('time (UTC)',inplace=True)
        d = d.resample('h',axis=0).mean()

        df = pd.concat([df,d])

    print('carbon dioxide erddap worked')
    
    # reset index
    df['datetime'] = np.array(df.index)
    index = pd.Index(range(len(df)))
    df.set_index(index,inplace=True)

    name = 'OOI_co2.p'
    df.to_pickle(name)