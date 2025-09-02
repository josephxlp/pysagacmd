import os
import rasterio
import geopandas as gpd
from rasterio.mask import mask
from shapely.geometry import box
from tqdm import tqdm
# skip if already exisit ?

def ensure_crs(gdf, target_crs):
    """Reproject GeoDataFrame if needed."""
    if gdf.crs != target_crs:
        return gdf.to_crs(target_crs)
    return gdf

def clip_raster_to_bbox(raster_path, bbox_geom, output_path):
    """Clip raster using bounding box and save to output path."""
    with rasterio.open(raster_path) as src:
        bbox = [bbox_geom]
        clipped, transform = mask(src, bbox, crop=True)

        out_meta = src.meta.copy()
        out_meta.update({
            "driver": "GTiff",
            "height": clipped.shape[1],
            "width": clipped.shape[2],
            "transform": transform
        })

        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        with rasterio.open(output_path, "w", **out_meta) as dest:
            dest.write(clipped)

def process_single_polygon(id, sg_bounds, sg_crs, tif_files, tilename, vboxes_dir, show_progress=True):
    """Process and clip all rasters for a single polygon."""
    iterator = tqdm(tif_files, desc=f'ID{id}', leave=False, disable=not show_progress)
    for tif_path in iterator:
        with rasterio.open(tif_path) as src:
            if sg_crs != src.crs:
                sg_bbox = gpd.GeoSeries([box(*sg_bounds)], crs=sg_crs).to_crs(src.crs).geometry.values[0]
            else:
                sg_bbox = box(*sg_bounds)

            out_dir = os.path.join(vboxes_dir, tilename, f'ID{id}')
            out_filename = os.path.join(out_dir, os.path.basename(tif_path))
            clip_raster_to_bbox(tif_path, sg_bbox, out_filename)

def clip_tifs_by_vbox_rois(vpath, tif_files, tilename, vboxes_dir):
    """Clip all TIFs using bounding boxes from vector file."""
    os.makedirs(vboxes_dir, exist_ok=True)
    gdf = gpd.read_file(vpath)
    gdf['id'] = gdf.index

    for id, row in tqdm(gdf.iterrows(), total=len(gdf), desc=f'Processing {tilename}', unit='ROI'):
        sg_bounds = row.geometry.bounds
        sg_crs = gdf.crs
        process_single_polygon(id, sg_bounds, sg_crs, tif_files, tilename, vboxes_dir)

def generate_tile_file_paths(base_path, tilename, suffixes):
    """Construct full paths to each TIF based on suffix list."""
    return [os.path.join(base_path, tilename, f"{tilename}_{suffix}.tif") for suffix in suffixes]

if __name__ == '__main__':
    #from uvars import TILES12_DPATH, VBOX_DPATH
    from uvars import tiles12_dir,vec_atroi_dir,tif_atroi_dir
    os.makedirs(vec_atroi_dir, exist_ok=True)
    os.makedirs(tif_atroi_dir, exist_ok=True)

    tilenames = ["N10E105", "S01W063", "N13E103"]

    # List of suffixes to avoid repetition
    suffixes = [
        "s1", "s2","etchm","esawc_x","esawc","tdem_hem","egm08",
        "tdem_dem_egm_v_3_gap","tdem_dem_egm","edem_egm"
       
    ]

    for tilename in tilenames:
        tif_files = generate_tile_file_paths(tiles12_dir, tilename, suffixes)
        vpath = f"{vec_atroi_dir}/{tilename}_ROIs.gpkg"
        
        clip_tifs_by_vbox_rois(vpath=vpath, tif_files=tif_files,
                               tilename=tilename, vboxes_dir=tif_atroi_dir)
        




