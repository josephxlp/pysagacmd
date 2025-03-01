import os

import geopandas as gpd
import rasterio
import requests
from rasterio.warp import transform_bounds

"""
- SRTMGL3 (SRTM GL3 90m)
- SRTMGL1 (SRTM GL1 30m)
- SRTMGL1_E (SRTM GL1 Ellipsoidal 30m)
- AW3D30 (ALOS World 3D 30m)
- AW3D30_E (ALOS World 3D Ellipsoidal, 30m)
- SRTM15Plus (Global Bathymetry SRTM15+ V2.1 500m)
- NASADEM (NASADEM Global DEM)
- COP30 (Copernicus Global DSM 30m)
- COP90 (Copernicus Global DSM 90m)
- EU_DTM (DTM 30m)
- GEDI_L3 (DTM 1000m)
- GEBCOIceTopo (Global Bathymetry 500m)
- GEBCOSubIceTopo (Global Bathymetry 500m)
"""

dem_types = {
    # "SRTM GL3 90m": "SRTMGL3",
    "SRTM GL1 30m": "SRTMGL1",
    "SRTM GL1 Ellipsoidal 30m": "SRTMGL1_E",
    "ALOS World 3D 30m": "AW3D30",
    "ALOS World 3D Ellipsoidal, 30m": "AW3D30_E",
    "Global Bathymetry SRTM15+ V2.1 500m": "SRTM15Plus",
    "NASADEM Global DEM": "NASADEM",
    "Copernicus Global DSM 30m": "COP30",
    # "Copernicus Global DSM 90m": "COP90",
    "DTM 30m": "EU_DTM",
    "DTM 1000m": "GEDI_L3",
    "Global Bathymetry 500m": "GEBCOIceTopo",
    "Global Bathymetry Sub Ice 500m": "GEBCOSubIceTopo",
}
print("available datasets on globaldem API")
print(dem_types)


def get_WGS84_bounds_vector(vpath):
    """
    Get the bounding box of a vector file in WGS84 (EPSG:4326).

    Parameters:
        vpath (str): Path to the vector file.

    Returns:
        tuple: (Xmin, Ymin, Xmax, Ymax) in WGS84.
    """
    # Read vector file using geopandas
    gdf = gpd.read_file(vpath)

    # Get bounding box and CRS
    bounds = gdf.total_bounds  # (Xmin, Ymin, Xmax, Ymax)
    src_crs = gdf.crs

    # Reproject bounds if the CRS is not WGS84
    if src_crs and src_crs.to_string() != "EPSG:4326":
        Xmin, Ymin, Xmax, Ymax = transform_bounds(
            src_crs, "EPSG:4326", bounds[0], bounds[1], bounds[2], bounds[3]
        )
    else:
        Xmin, Ymin, Xmax, Ymax = bounds

    return Xmin, Ymin, Xmax, Ymax


def get_WGS84_bounds_raster(rpath):
    """
    Get the bounding box of a raster file in WGS84 (EPSG:4326).

    Parameters:
        rpath (str): Path to the raster file.

    Returns:
        tuple: (Xmin, Ymin, Xmax, Ymax) in WGS84.
    """
    with rasterio.open(rpath) as src:
        bounds = src.bounds
        src_crs = src.crs

        if src_crs and src_crs.to_string() != "EPSG:4326":
            bounds = transform_bounds(
                src_crs,
                "EPSG:4326",
                bounds.left,
                bounds.bottom,
                bounds.right,
                bounds.top,
            )

    return bounds.left, bounds.bottom, bounds.right, bounds.top


def download_globaldem(outdir, bbox, varname="SRTMGL3", api_key="api"):
    """
    Downloads DEM data from the OpenTopography API for the given bounding box.

    Parameters:
        outdir (str): Directory to save the downloaded DEM file.
        bbox (tuple): Bounding box (Xmin, Ymin, Xmax, Ymax) in WGS84.
        varname (str): DEM type (default is "SRTMGL3").
        api_key (str): OpenTopography API key.

    Returns:
        None
    """
    Xmin, Ymin, Xmax, Ymax = bbox
    print(
        f"Extracted Coordinates ::: Xmin = {Xmin}\t Ymin = {Ymin}\t Xmax = {Xmax}\t Ymax = {Ymax}"
    )

    url = (
        f"https://portal.opentopography.org/API/globaldem?"
        f"demtype={varname}&south={Ymin}&north={Ymax}&west={Xmin}&east={Xmax}"
        f"&outputFormat=GTiff&API_Key={api_key}"
    )

    output_filename = f"{varname}_{Xmin:.4f}_{Ymin:.4f}_{Xmax:.4f}_{Ymax:.4f}.tif"
    output_path = os.path.join(outdir, output_filename)


    if not os.path.isfile(output_path):
        try:
            response = requests.get(url)

            if response.status_code == 200:
                with open(output_path, "wb") as f:
                    f.write(response.content)
                print(f"Data downloaded and saved to {output_path}")
            else:
                print(
                    f"Error: HTTP {response.status_code}\nResponse: {response.text}\nURL: {url}"
                )
        except requests.RequestException as e:
            print(f"Request failed: {e}")
        except Exception as e:
            print(f"An error occurred: {e}")
    else:
        print(f"File already exists: {output_path}")

def is_vector_file(filepath):
    """Check if the file is a vector file based on its extension."""
    vector_extensions = [".shp", ".geojson", ".gpkg"]
    return any(filepath.lower().endswith(ext) for ext in vector_extensions)


def is_raster_file(filepath):
    """Check if the file is a raster file based on its extension."""
    raster_extensions = [".tif", ".tiff", ".img"]
    return any(filepath.lower().endswith(ext) for ext in raster_extensions)


def extract_tile_name_from_filename(filepath):
    """
    Extract the tile name from the file name.
    Assumes the tile name is in the format like 'N58E025'.
    """
    filename = os.path.basename(filepath)  # Get the file name (e.g., 'N58E025.gpkg')
    tile_name = os.path.splitext(filename)[0]  # Remove the extension (e.g., 'N58E025')
    return tile_name


if __name__ == "__main__":
    from uvars import ot_dpath, gpkg_pattern
    from dotenv import load_dotenv
    import numpy as np
    from glob import glob

    # Load environment variables
    load_dotenv()
    api_key = os.getenv("OT_API_KEY")
    #print(f"API Key Loaded: {api_key}")
    OPEN_TOPOGRAPHY_DPATH = ot_dpath

    varnames = np.array(list(dem_types.values()))
    varnames = np.array(['COP90','GEDI_L3'])
    file_paths = glob(gpkg_pattern)
    print(f"Found {len(file_paths)} vector files")

    for rpath in file_paths:
        # Extract tile name from the file name
        tile = extract_tile_name_from_filename(rpath)
        tile_dpath = os.path.join(OPEN_TOPOGRAPHY_DPATH, tile)
        os.makedirs(tile_dpath, exist_ok=True)

        if is_vector_file(rpath):
            print(f"Processing vector file: {rpath}")
            bbox = get_WGS84_bounds_vector(rpath)
        elif is_raster_file(rpath):
            print(f"Processing raster file: {rpath}")
            bbox = get_WGS84_bounds_raster(rpath)
        else:
            print(f"Unsupported file type: {rpath}")
            continue

        # Print extracted bbox
        Xmin, Ymin, Xmax, Ymax = bbox
        for varname in varnames:
            tile_var_dpath = os.path.join(tile_dpath, varname)
            os.makedirs(tile_var_dpath, exist_ok=True)
            download_globaldem(
                outdir=tile_var_dpath, bbox=bbox, varname=varname, api_key=api_key
            )
            print(f"DEM downloaded for {varname} at {tile_var_dpath}")

    