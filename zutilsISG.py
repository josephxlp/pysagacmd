import os 
import rasterio
from rasterio.transform import from_origin
import re
import numpy as np
from pprint import pprint

def extract_metadata_and_data(filename):
    with open(filename, 'r', encoding='utf-8') as file:
        content = file.read()
    
    # Find the header section
    match = re.search(r'begin_of_head[\s=]*([\s\S]*?)end_of_head', content)
    if not match:
        raise ValueError("Header section not found in the file.")
    
    header_content = match.group(1).strip()
    
    # Extract key-value pairs
    metadata = {}
    for line in header_content.splitlines():
        line = line.strip()
        if not line:
            continue
        
        # Handle key-value pairs with colon or equals sign
        match = re.match(r'([\w\s]+)[:=]\s*(.*)', line)
        if match:
            key = match.group(1).strip()
            value = match.group(2).strip()
            metadata[key] = value
    
    # Extract numerical data after the header
    data_start = match.end()
    data_content = content[data_start:].strip()
    
    # Process data while ignoring non-numeric lines
    data_rows = []
    for row in data_content.splitlines():
        row = row.strip()
        if not row:
            continue
        try:
            data_rows.append(list(map(float, row.split())))
        except ValueError:
            continue  # Ignore lines that can't be converted to float
    
    data_array = np.array(data_rows)
    
    return metadata, data_array

def dms_to_dd(dms):
    """Convert degrees, minutes, seconds to decimal degrees."""
    parts = dms.split('°')
    degrees = float(parts[0])
    minutes = float(parts[1].split('\'')[0]) / 60.0
    seconds = float(parts[1].split('\'')[1].strip('"')) / 3600.0
    return degrees + minutes + seconds if degrees >= 0 else degrees - minutes - seconds



# try:isg_fpath and txt_fpath
isg_fpath = r"c:\Users\Joseph\Downloads\Topo_DEM_Mekong_delta_excl_rivers_and_bedrock\igs\Brazil_2015_MAPGEO2015_gravG_20201222.isg"
tif_fpath = isg_fpath.replace('.isg','_A.tif')
# metadata, data_array = extract_metadata_and_data(isg_fpath)
# print(metadata)
# print(data_array.shape)

def isg2tif_f1(isg_fpath,output_file):
    metadata, data_array = extract_metadata_and_data(isg_fpath)
    delta_lat = float(metadata['delta lat'])  # Resolution in latitude (decimal degrees)
    delta_lon = float(metadata['delta lon'])  # Resolution in longitude (decimal degrees)
    lat_min = float(metadata['lat min'])      # Minimum latitude (decimal degrees)
    lat_max = float(metadata['lat max'])      # Maximum latitude (decimal degrees)
    lon_min = float(metadata['lon min'])      # Minimum longitude (decimal degrees)
    lon_max = float(metadata['lon max'])      # Maximum longitude (decimal degrees)
    nrows = int(metadata['nrows'])            # Number of rows
    ncols = int(metadata['ncols'])            # Number of columns
    nodata = float(metadata['nodata'])        # No-data value

    #if metadata['reference'] == 'WGS84':
    crs = "EPSG:4326"                         # CRS for WGS84

    transform = rasterio.transform.from_origin(
                    west=lon_min,
                    north=lat_max,
                    xsize=delta_lon,
                    ysize=delta_lat
                )
    with rasterio.open(
                        output_file,
                        'w',
                        driver='GTiff',
                        height=nrows,
                        width=ncols,
                        count=1,  # Single band
                        dtype=data_array.dtype,
                        crs=crs,
                        transform=transform,
                        nodata=nodata
                    ) as dst:
                        dst.write(data_array, 1)

    print(f"isg2tif_f1 write to {output_file}")


def isg2tif_f2(input_file,output_file):

    metadata, data_array = extract_metadata_and_data(input_file)
    #pprint(metadata)
    print(data_array.shape)

    lat_min = dms_to_dd(metadata['lat min'])
    lat_max = dms_to_dd(metadata['lat max'])
    lon_min = dms_to_dd(metadata['lon min'])
    lon_max = dms_to_dd(metadata['lon max'])

    # Step 3: Define resolution in decimal degrees
    delta_lat = 5 / 60.0  # 5 minutes in degrees
    delta_lon = 5 / 60.0  # 5 minutes in degrees

    transform = rasterio.transform.from_origin(
        west=lon_min,
        north=lat_max,
        xsize=delta_lon,
        ysize=delta_lat
    )

    with rasterio.open(
        output_file,
        'w',
        driver='GTiff',
        height=int(metadata['nrows']),
        width=int(metadata['ncols']),
        count=1,  # Single band
        dtype=data_array.dtype,
        crs=f"EPSG:{metadata['EPSG code']}",
        transform=transform,
        nodata=float(metadata['nodata'])
    ) as dst:
        dst.write(data_array, 1)

    print(f"isg2tif_f2 write to {output_file}")