import pandas as pd
import numpy as np
import os
from collections import OrderedDict
import gzip

import geocoder
from ish_parser import ish_report, ish_reportException

DIR_THIS_SCRIPT = os.path.dirname(os.path.realpath(__file__))


# Function to download raw
def ftp_to_raw_entry_list(file_name_noaa, file_name_local, ftp):
    from ftplib import FTP
    '''This function download raw weather data entries from noaa
    '''
    # Retrieve binary gzip file from NOAA
    try:
        file = open(file_name_local, 'wb')
        ftp.retrbinary('RETR '+ file_name_noaa, file.write)
        file.close()
    except:
        os.remove(file_name_local)
        print(f"Failed to download the file: {file_name_noaa}")
        return
    
    # Read binary gzip file into a list
    v_noaa_raw_entries = []
    gzip_file = gzip.open(file_name_local, 'rb')
    
    for i, line in enumerate(gzip_file):
        temp_bytes = line[0:-1]
        temp_str = str(temp_bytes)[2:-1]
        v_noaa_raw_entries.append(temp_str)

    return v_noaa_raw_entries

# Functions to get data from parsed ish report.
get_datetime_from_rpt = lambda x: ish_report().loads(x).datetime
get_T_C_from_rpt = lambda x: ish_report().loads(x).air_temperature.get_numeric()
get_T_F_from_rpt = lambda x: ish_report().loads(x).air_temperature.get_fahrenheit().get_numeric()


def download_noaa_weather_element(station_ID, year, ftp_instance):
    ftp_path = '/pub/data/noaa/' + str(year)
    file_name_noaa = station_ID + '-' + str(year) + '.gz'
    raw_file_out_dir = ''
    file_name_local = os.path.join(raw_file_out_dir, station_ID + '-' + str(year) + '.gz')
    ftp_instance.cwd(ftp_path)
    v_noaa_raw_elements = ftp_to_raw_entry_list(file_name_noaa, file_name_local, ftp_instance)
    v_noaa_datetime = list(map(get_datetime_from_rpt, v_noaa_raw_elements))
    v_noaa_T_F = list(map(get_T_C_from_rpt, v_noaa_raw_elements))
    df_out = pd.DataFrame(OrderedDict({
        'Datetime': v_noaa_datetime,
        'Temperature': pd.to_numeric(v_noaa_T_F, errors='coerce'),
    }))
    # Clean the raw zip file
    os.remove(file_name_local)
    return(df_out)


# Function to wrap all
def download_noaa_weather(station_list_csv_path, years=[2019], work_dir='./'):
    """
    { item_description }
    """


    from ftplib import FTP
    try:
        df_stations = pd.read_csv(station_list_csv_path)
    except:
        print(f"Fail to read station list from {station_list_csv_path}")
        return
    
    # Log in to NOAA FTP
    ftp=FTP('ftp.ncdc.noaa.gov')
    ftp.login()
    print('NOAA FTP login succeeded.')

    n_total_stations = len(df_stations['StationID'].index)
    
    for year in years:
        # Create the year sub-folder if not exists.
        year_sub_dir = os.path.join(work_dir, str(year))
        v_missing = []
        if not os.path.exists(year_sub_dir):
            os.mkdir(year_sub_dir)
        for i, station_ID in enumerate(df_stations['StationID']):
            try:
                print(f"Downloading {year} weather data for station: {station_ID} -- ({i}/{n_total_stations}) -- {round(i/n_total_stations*100, 2)}%")
                df_out = download_noaa_weather_element(station_ID, year, ftp)
                df_out.to_csv(os.path.join(year_sub_dir, str(year) + "_" + station_ID + ".csv"), index = False)
            except:
                print(f"Failed to download weather data for station: {station_ID}")
                v_missing.append(station_ID)
        df_missing = pd.DataFrame({'Missing Station ID': v_missing})
        df_missing.to_csv(str(year) + "_" + 'missing.csv', index=False)
        
    ftp.quit()
    print('Weather download finished. Logout FTP.')


def add_weather_to_ts(df_ts, noaa_station_id, dir_all_weather=None, new_cols=None, ftp_instance=None):
    '''
    This function add weather data to a dataframe of time-series data. The weather data is from the specified NOAA weather station. 
    If the directory of pre-downloaded weather files is not provided, it will download the weather data from NOAA on the fly.   
    '''
    if dir_all_weather==None:
        print('Will download from NOAA')
        if ftp_instance == None:
            from ftplib import FTP
            ftp_instance=FTP('ftp.ncdc.noaa.gov')
            ftp_instance.login()
            print('NOAA FTP login succeeded.')
        for i, year in enumerate(set(pd.to_datetime(df_ts['Datetime']).dt.year)):
            df_weather_temp = download_noaa_weather_element(noaa_station_id, year, ftp_instance)
            df_weather_temp['Datetime'] = pd.to_datetime(df_weather_temp['Datetime'], format='%Y-%m-%d %H:%M:%S').dt.tz_localize(None)
            if i == 0:
                df_weather = df_weather_temp
            else:
                df_weather = df_weather.append(df_weather_temp)
        ftp_instance.quit()
        print('---> NOAA FTP Logout succeeded.')
    else:
        print('Use pre-downloaded weather data')
        # Get the raw weather data
        for i, year in enumerate(set(pd.to_datetime(df_ts['Datetime']).dt.year)):
            dir_weather_csv = os.path.join(dir_all_weather, str(year), f"{year}_{noaa_station_id}.csv")
            df_weather_temp = pd.read_csv(dir_weather_csv, parse_dates=True)
            df_weather_temp['Datetime'] = pd.to_datetime(df_weather_temp['Datetime'], format='%Y-%m-%d %H:%M:%S').dt.tz_localize(None)
            if i == 0:
                df_weather = df_weather_temp
            else:
                df_weather = df_weather.append(df_weather_temp)

    # Merge the weather data to the time-series data
    df_ts_weather = pd.merge_asof(df_ts, df_weather.dropna().sort_values('Datetime'), on='Datetime', tolerance=pd.Timedelta('3600s'))
    return df_ts_weather if new_cols==None else df_ts_weather[new_cols]


################################################################################
# Some utility functions
################################################################################
def geocode_address(input_address):
    """
    Geocode an address to (lat, lon)
    """
    try:
        latlng = geocoder.arcgis(input_address).latlng
        return latlng
    except:
        pass

    try:
        latlng = geocoder.osm(input_address).latlng
        return latlng
    except:
        pass

    try:
        latlng = geocoder.ottawa(input_address).latlng
        return latlng
    except:
        pass

    try:
        latlng = [float(s) for s in geocoder.yandex(input_address).latlng]
        return latlng
    except:
        pass

def haversine_distance(lat1, lon1, lat2, lon2):
    # Get radians from decimals
    r_lat1, r_lon1, r_lat2, r_lon2 = map(np.radians, [lat1, lon1, lat2, lon2])

    # Calculate the distance between the two locations
    temp = np.sin((r_lat2 - r_lat1) / 2) ** 2 + np.cos(r_lat1) * np.cos(r_lat2) * np.sin((r_lon2 - r_lon1) / 2) ** 2
    distance = 2 * 6371 * np.arcsin(np.sqrt(temp))
    return (distance)

def find_closest_weather_station(tuple_lat_lon, 
    df_weather_station_list=pd.read_csv(os.path.join(DIR_THIS_SCRIPT, 'station_list.csv'))
    ):
    """
    Finds a closest weather station.
    
    :param      tuple_lat_lon:            The tuple (lat, lon)
    :type       tuple_lat_lon:            { type_description }
    :param      df_weather_station_list:  The df weather station list
    :type       df_weather_station_list:  { pandas dataframe }
    """
    latitude, longitude = float(tuple_lat_lon[0]), float(tuple_lat_lon[1])

    v_coord = np.asarray(df_weather_station_list[['LAT', 'LON']].values)
    # Find the closest and second closest weather station (backup if the closest doesn't work)
    v_distance = [haversine_distance(latitude, longitude, coord[0], coord[1]) for coord in v_coord]
    closest_index = np.argmin(v_distance)
    second_closest_index = np.argpartition(v_distance, 2)[2]
    third_closest_index = np.argpartition(v_distance, 3)[3]

    closest_weather_station_ID = df_weather_station_list.loc[closest_index, 'StationID']
    closest_weather_station_name = df_weather_station_list.loc[closest_index, 'STATION NAME']
    second_closest_weather_station_ID = df_weather_station_list.loc[second_closest_index, 'StationID']
    second_closest_weather_station_name = df_weather_station_list.loc[second_closest_index, 'STATION NAME']
    third_closest_weather_station_ID = df_weather_station_list.loc[third_closest_index, 'StationID']
    third_closest_weather_station_name = df_weather_station_list.loc[third_closest_index, 'STATION NAME']
    
    return closest_weather_station_ID, closest_weather_station_name

if __name__ == '__main__':
    work_dir = os.path.dirname(os.path.abspath( __file__ ))

    # e.g.1, Utility function: find the geographically closest station for a (lat, lon) coordinate
    df_stations = pd.read_csv('station_list.csv')
    str_address = 'Lawrence Berkeley National Lab'
    station_id, station_name = find_closest_weather_station(geocode_address(str_address), df_stations)
    print(f'The closest NOAA weather station for "{str_address}" is ID: {station_id}, Name: {station_name}')


    # e.g.2, download weather data for a single station and year 
    # from ftplib import FTP
    # ftp=FTP('ftp.ncdc.noaa.gov')
    # ftp.login()
    # print('---> NOAA FTP login succeeded.')
    # print(download_noaa_weather_element('722780-23183', 2018, ftp).describe())
    # ftp.quit()
    # print('---> NOAA FTP Logout succeeded.')
    
    # e.g.3, Download weather file for all weather stations in the station_list.csv for the specified years 
    # download_noaa_weather('station_list.csv', [2019, 2020], work_dir)


    
