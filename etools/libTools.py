
import os
import rasterio

def split_raster_stack_to_bands(stack_path, outdir, basename, overwrite=False):
    #make this parallel
    """
    Splits a raster stack into individual bands and saves each as a separate file.

    Args:
        stack_path (str): Path to the input raster stack.
        outdir (str): Directory to save the output bands.
        basename (str): Base name to use in output filenames.
        overwrite (bool, optional): If True, overwrite existing files. Default is False.
    """
    os.makedirs(outdir, exist_ok=True)

    with rasterio.open(stack_path) as src:
        for i in range(1, src.count + 1):
            band = src.read(i)
            out_path = os.path.join(outdir, f"{basename}_B{i}.tif")

            if not overwrite and os.path.exists(out_path):
                print(f"Skipping existing file: {out_path}")
                continue

            profile = src.profile.copy()
            profile.update(count=1)

            with rasterio.open(out_path, 'w', **profile) as dst:
                dst.write(band, 1)
            print(f"Saved: {out_path}")


import os
import rasterio
from rasterio.warp import reproject, Resampling
import numpy as np

def _reproject_align_raster(src_path, ref_profile):
    """
    Reproject and align a raster to match the reference profile.
    """
    with rasterio.open(src_path) as src:
        dst_array = np.empty((ref_profile['height'], ref_profile['width']), dtype=ref_profile['dtype'])

        reproject(
            source=rasterio.band(src, 1),
            destination=dst_array,
            src_transform=src.transform,
            src_crs=src.crs,
            dst_transform=ref_profile['transform'],
            dst_crs=ref_profile['crs'],
            resampling=Resampling.bilinear
        )

    return dst_array

def create_raster_stack_from_list(raster_paths, outdir, basename, align_path=None, overwrite=False):
    """
    Stacks a list of single-band rasters into a multi-band raster after aligning them.

    Args:
        raster_paths (list of str): Paths to input rasters.
        outdir (str): Directory to save the output stack.
        basename (str): Base name for the output stack.
        align_path (str or None): Path to the raster used as alignment reference.
        overwrite (bool): If True, overwrite existing file. Default is False.
    """
    os.makedirs(outdir, exist_ok=True)
    out_path = os.path.join(outdir, f"{basename}_stack.tif")

    if os.path.exists(out_path) and not overwrite:
        print(f"Stack file already exists: {out_path}. Skipping.")
        return

    if align_path is None:
        align_path = raster_paths[0]

    with rasterio.open(align_path) as ref:
        ref_profile = ref.profile.copy()
        ref_profile.update(count=len(raster_paths))

        aligned_bands = []
        for path in raster_paths:
            aligned = _reproject_align_raster(path, ref_profile)
            aligned_bands.append(aligned)

    stack_array = np.stack(aligned_bands, axis=0)

    with rasterio.open(out_path, 'w', **ref_profile) as dst:
        dst.write(stack_array)
    
    print(f"Saved stacked raster: {out_path}")



""" 
Library Tools
https://saga-gis.sourceforge.io/saga_tool_doc/8.3.0/grid_tools.html

-------------------------------

Aggregate
Change Cell Values (interactive)
Change Data Storage
Change Grid Values
Change Grid Values - Flood Fill [+]
Change Grid Values - Flood Fill (interactive) [+]
Change a Grid's No-Data Value
Clip Grids
Clip Grids (interactive)
Close Gaps [x]
Close Gaps with Spline [x]
Close Gaps with Stepwise Resampling [x]
Close One Cell Gaps
Combine Classes
Combine Grids
Constant Grid
Copy Grid
Create Grid System
Crop to Data
Grid Buffer
Grid Cell Index
Grid Masking
Grid Proximity Buffer
Grid Value Request (interactive)
Grids from Classified Grid and Table
Interpolate Values Along Line (interactive)
Invert Data/No-Data
Invert Grid
Mirror Grid
Mosaicking
Mosaicking (Grid Collections)
Patching
Proximity Grid
Reclassify Grid Values
Resampling
Select Grid from List
Shrink and Expand
Threshold Buffer
Tiling [x]
Transpose Grids



#benchmark agains rasterio+parallel
"""

import time
from typing import List, Optional
import subprocess

""" 
give me an example using the functino above
add one other argument call outdir : make sure a folder is crated on this path 
all the ouputs should be saved in the output directory
all the output name should be modified inside the function to take all the other arguments so its easier to track each experiment
check if the output exist, if not run, else skip. add overwrite default to False, but if True runs regardless 
# keep the comments inside the fuction in cmd but non outside 
returns sdat path
"""

def close_gaps_stepwise_resampling(
    input_grid: str,
    output_grid: str = None,
    mask_grid: str = None,
    resampling: int = 3,
    grow_factor: float = 2.0
):
    """
    Fills NoData values in a raster by iteratively resampling and growing known values.

    This stepwise interpolation approach expands known values into unknown areas
    using progressive upscaling and resampling methods (e.g. B-spline or bilinear).

    Use case:
    - Useful for DEMs or grids with large or complex voids.
    - Suitable when gaps are not easily filled by local interpolation alone.

    Parameters
    ----------
    input_grid : str
        Path to the input SAGA grid file (.sgrd) with gaps.

    output_grid : str, optional
        Path to output grid with gaps filled. If not set, original is modified.

    mask_grid : str, optional
        Optional mask grid to define fill regions.

    resampling : int, optional
        Resampling method used during upscaling:
        0 = Nearest Neighbour,
        1 = Bilinear,
        2 = Bicubic Spline,
        3 = B-Spline (default)

    grow_factor : float, optional
        Factor by which known-value areas grow at each step (≥ 1.0). Default = 2.0
    """
    cmd = [
        "saga_cmd", "grid_tools", "29",
        "-INPUT", input_grid,
        "-RESAMPLING", str(resampling),
        "-GROW", str(grow_factor)
    ]

    if mask_grid:
        cmd += ["-MASK", mask_grid]
    if output_grid:
        cmd += ["-RESULT", output_grid]

    start = time.time()
    subprocess.run(cmd, check=True)
    end = time.time()
    print(f"Stepwise gap filling completed in {end - start:.2f} seconds.")



def close_gaps_with_spline(
    input_grid: str,
    output_grid: str = None,
    mask_grid: str = None,
    max_gap_cells: int = 0,
    max_points: int = 1000,
    local_points: int = 20,
    extended_neighborhood: bool = False,
    neighbours: int = 0,
    radius: int = 0,
    relaxation: float = 0.0
):
    """
    Fills NoData (gap) cells in a grid using spline-based local interpolation in SAGA-GIS.

    Use case:
    More refined and locally adapted interpolation than 'Close Gaps'.
    Ideal for elevation models or continuous raster surfaces with irregular missing areas.

    Parameters:
    ----------
    input_grid : str
        Path to the input SAGA grid file (.sgrd) with gaps to close.

    output_grid : str, optional
        Path to save the output gap-filled grid. If not provided, modifies in place.

    mask_grid : str, optional
        Optional mask grid to limit filling only to selected areas.

    max_gap_cells : int, optional
        Only process gaps smaller than this number of cells. 0 = all gaps.

    max_points : int, optional
        Maximum number of nearby known points used in interpolation.

    local_points : int, optional
        Number of points used in local interpolation.

    extended_neighborhood : bool, optional
        Use extended neighborhood search. Default is False.

    neighbours : int, optional
        Neighborhood type: 0 = Neumann, 1 = Moore.

    radius : int, optional
        Neighborhood radius in cells.

    relaxation : float, optional
        Relaxation factor for interpolation smoothness. 0 = none.
    """
    cmd = [
        "saga_cmd", "grid_tools", "25",
        "-GRID", input_grid,
        "-MAXGAPCELLS", str(max_gap_cells),
        "-MAXPOINTS", str(max_points),
        "-LOCALPOINTS", str(local_points),
        "-EXTENDED", str(int(extended_neighborhood)),
        "-NEIGHBOURS", str(neighbours),
        "-RADIUS", str(radius),
        "-RELAXATION", str(relaxation)
    ]

    if mask_grid:
        cmd.extend(["-MASK", mask_grid])
    if output_grid:
        cmd.extend(["-CLOSED", output_grid])

    start = time.time()
    subprocess.run(cmd, check=True)
    end = time.time()
    print(f"Spline gap filling completed in {end - start:.2f} seconds.")


def close_gaps_saga(
    input_grid: str,
    output_grid: str = None,
    mask_grid: str = None,
    threshold: float = 0.1
):
    """
    Fills NoData (gap) cells in a grid using SAGA-GIS's 'Close Gaps' tool.

    Use case:
    Use this to clean or preprocess elevation models or any gridded data 
    containing missing values (NoData) before analysis, modeling, or visualization.

    Parameters:
    ----------
    input_grid : str
        Path to the input grid file (.sgrd or .tif) with gaps to fill.

    output_grid : str, optional
        Path to save the output grid with filled gaps.
        If not specified, input grid will be modified in place.

    mask_grid : str, optional
        Optional mask grid to restrict where gaps should be filled.

    threshold : float, optional
        Tension threshold for surface fitting interpolation.
        Lower values mean tighter (smoother) interpolation.
        Default: 0.1
    """
    cmd = [
        "saga_cmd", "grid_tools", "7",
        "-INPUT", input_grid,
        "-THRESHOLD", str(threshold)
    ]

    if mask_grid:
        cmd.extend(["-MASK", mask_grid])
    if output_grid:
        cmd.extend(["-RESULT", output_grid])

    start = time.time()
    subprocess.run(cmd, check=True)
    end = time.time()
    print(f"Gap filling completed in {end - start:.2f} seconds.")



def saga_mosaic_grids(
    input_grids: Optional[List[str]] = None,
    file_list: Optional[str] = None,
    target_grid: Optional[str] = None,
    target_template: Optional[str] = None,
    name: str = "Mosaic",
    output_path: Optional[str] = None,
    cellsize: float = 1.0,
    data_type: int = 9,
    resampling: int = 3,
    overlap: int = 1,
    blend_distance: float = 10.0,
    blend_boundary: int = 0,
    match: int = 0,
    target_def: int = 0,
    extent: Optional[dict] = None,
    fit: int = 0,
):
    """
    Merge multiple SAGA grids into one mosaic grid.

    Parameters
    ----------
    input_grids : list of str, optional
        List of .sgrd files to mosaic. Ignored if `file_list` is provided.

    file_list : str, optional
        Path to a text file listing .sgrd paths line-by-line. Use when merging many files.

    target_grid : str, optional
        Output mosaic grid file (.sgrd). If not given, will print process only.

    target_template : str, optional
        Path to template grid for setting output grid system (takes precedence over extent/cellsize).

    name : str, optional
        Output grid base name. Default = 'Mosaic'.

    output_path : str, optional
        Output directory to write the final mosaic grid.

    cellsize : float
        Target cell size. Used if `target_def = 0`. Default = 1.0.

    data_type : int
        Output data type (see SAGA docs for codes, e.g., 9 = same as first grid).

    resampling : int
        Resampling method: 0 = NN, 1 = Bilinear, 2 = Bicubic, 3 = B-Spline.

    overlap : int
        Overlap resolution: 0 = first, 1 = last, ..., 6 = feathering.

    blend_distance : float
        Distance (in map units) used for blending/feathering.

    blend_boundary : int
        Blending boundary definition.

    match : int
        Histogram matching: 0 = none, 1 = to first grid, etc.

    target_def : int
        Target grid system: 0 = user-defined, 1 = use template/grid system.

    extent : dict, optional
        Bounding box extent: {'xmin': float, 'xmax': float, 'ymin': float, 'ymax': float}

    fit : int
        Grid fit: 0 = nodes, 1 = cells.
    """
    cmd = ["saga_cmd", "grid_tools", "3"]

    if input_grids:
        cmd += ["-GRIDS", ";".join(input_grids)]
    elif file_list:
        cmd += ["-FILE_LIST", file_list]

    cmd += [
        "-NAME", name,
        "-TYPE", str(data_type),
        "-RESAMPLING", str(resampling),
        "-OVERLAP", str(overlap),
        "-BLEND_DIST", str(blend_distance),
        "-BLEND_BND", str(blend_boundary),
        "-MATCH", str(match),
        "-TARGET_DEFINITION", str(target_def)
    ]

    if target_template:
        cmd += ["-TARGET_TEMPLATE", target_template]

    if target_def == 0 and extent:
        cmd += [
            "-TARGET_USER_SIZE", str(cellsize),
            "-TARGET_USER_XMIN", str(extent["xmin"]),
            "-TARGET_USER_XMAX", str(extent["xmax"]),
            "-TARGET_USER_YMIN", str(extent["ymin"]),
            "-TARGET_USER_YMAX", str(extent["ymax"]),
            "-TARGET_USER_FITS", str(fit)
        ]

    if target_grid:
        cmd += ["-TARGET_OUT_GRID", target_grid]

    print("Running SAGA mosaic tool...")
    start = time.time()
    subprocess.run(cmd, check=True)
    print(f"Mosaic complete in {time.time() - start:.2f} seconds.")


def saga_patch_grid(
    original_grid: str,
    patch_grid: str,
    output_grid: str = None,
    resampling: int = 3
):
    """
    Patch NoData gaps in a grid using another grid.

    Parameters
    ----------
    original_grid : str
        Path to the input grid (.sgrd) with gaps to be filled.

    patch_grid : str
        Path to the patch grid (.sgrd) to use as source for filling gaps.

    output_grid : str, optional
        Output grid file (.sgrd) with filled values. If not set, modifies original.

    resampling : int, optional
        Resampling method if grids are not aligned:
        0 = Nearest Neighbour, 1 = Bilinear, 2 = Bicubic Spline, 3 = B-Spline (default).
    """
    cmd = [
        "saga_cmd", "grid_tools", "5",
        "-ORIGINAL", original_grid,
        "-ADDITIONAL", patch_grid,
        "-RESAMPLING", str(resampling)
    ]

    if output_grid:
        cmd += ["-COMPLETED", output_grid]

    print("Running SAGA grid patching...")
    start = time.time()
    subprocess.run(cmd, check=True)
    print(f"Patch completed in {time.time() - start:.2f} seconds.")


def saga_resample_grids(
    input_grids: list,
    output_grids: str = None,
    target_template: str = None,
    keep_type: bool = False,
    scale_up: int = 5,
    scale_down: int = 3,
    user_defined: bool = False,
    cellsize: float = 1.0,
    xmin: float = 0.0,
    xmax: float = 100.0,
    ymin: float = 0.0,
    ymax: float = 100.0,
    cols: int = 101,
    rows: int = 101,
    fit: int = 0
):
    """
    Resample a list of SAGA grids using specified target system or template.

    Parameters
    ----------
    input_grids : list
        List of input grid paths (.sgrd) to be resampled.

    output_grids : str, optional
        Output path for resampled grids (grid list output).

    target_template : str, optional
        Grid file (.sgrd) to define the target grid system.

    keep_type : bool, default=False
        Whether to preserve original data type.

    scale_up : int, default=5
        Method for upscaling:
        [0]=NN, [1]=Bilinear, [2]=Bicubic, [3]=B-Spline, [4]=Mean,
        [5]=Mean Area Weighted (default), [6]=Min, [7]=Max, [8]=Majority.

    scale_down : int, default=3
        Method for downscaling:
        [0]=NN, [1]=Bilinear, [2]=Bicubic, [3]=B-Spline.

    user_defined : bool, default=False
        Whether to define the grid system manually (True) or via template (False).

    cellsize, xmin, xmax, ymin, ymax : float
        Define target extent when using `user_defined=True`.

    cols, rows : int
        Number of cells (only used if `user_defined=True`).

    fit : int, default=0
        Grid alignment: [0]=nodes, [1]=cells
    """

    cmd = ["saga_cmd", "grid_tools", "0"]

    # Input
    cmd += ["-INPUT", ",".join(input_grids)]

    # Optional Output
    if output_grids:
        cmd += ["-OUTPUT", output_grids]

    # Preserve type
    cmd += ["-KEEP_TYPE", str(int(keep_type))]

    # Resampling methods
    cmd += ["-SCALE_UP", str(scale_up)]
    cmd += ["-SCALE_DOWN", str(scale_down)]

    # Target system
    if user_defined:
        cmd += [
            "-TARGET_DEFINITION", "0",
            "-TARGET_USER_SIZE", str(cellsize),
            "-TARGET_USER_XMIN", str(xmin),
            "-TARGET_USER_XMAX", str(xmax),
            "-TARGET_USER_YMIN", str(ymin),
            "-TARGET_USER_YMAX", str(ymax),
            "-TARGET_USER_COLS", str(cols),
            "-TARGET_USER_ROWS", str(rows),
            "-TARGET_USER_FITS", str(fit)
        ]
    elif target_template:
        cmd += [
            "-TARGET_DEFINITION", "1",
            "-TARGET_TEMPLATE", target_template
        ]

    print("Running SAGA resampling...")
    start = time.time()
    subprocess.run(cmd, check=True)
    print(f"Resampling completed in {time.time() - start:.2f} seconds.")


def saga_tile_grid(
    grid_input: str,
    tiles_output: str = None,
    save_tiles: bool = False,
    output_dir: str = None,
    base_name: str = "tile",
    overlap: int = 0,
    overlap_sym: int = 0,
    method: int = 0,
    nx: int = 100,
    ny: int = 100,
    xrange_min: float = None,
    xrange_max: float = None,
    yrange_min: float = None,
    yrange_max: float = None,
    dcell: float = 1.0,
    dx: float = 100.0,
    dy: float = 100.0
):
    """
    Tile a large grid into smaller subgrids using SAGA GIS.

    Parameters
    ----------
    grid_input : str
        Path to the input grid (.sgrd).

    tiles_output : str, optional
        Output grid list file to hold tile objects (not used if `save_tiles=True`).

    save_tiles : bool, default=False
        Whether to save each tile individually to disk.

    output_dir : str, optional
        Output folder for saving tiles if `save_tiles=True`.

    base_name : str, default='tile'
        Base filename used when saving individual tiles.

    overlap : int, default=0
        Number of overlapping cells between adjacent tiles.

    overlap_sym : int, default=0
        Where to apply overlapping cells:
        [0]=symmetric, [1]=bottom/left, [2]=top/right.

    method : int, default=0
        Tiling method:
        [0]=number of grid cells per tile (NX/NY),
        [1]=coordinate-based (XRANGE, DY, etc.)

    nx, ny : int
        Number of columns and rows per tile (used if method=0).

    xrange_min, xrange_max, yrange_min, yrange_max : float
        Used only if method=1. Defines offset/range coordinates.

    dcell : float
        Cell size (only if method=1).

    dx, dy : float
        Tile size in map units (used if method=1).
    """

    cmd = ["saga_cmd", "grid_tools", "27"]

    # Required input
    cmd += ["-GRID", grid_input]

    # Optional output grid list (if not saving to disk)
    if tiles_output:
        cmd += ["-TILES", tiles_output]

    # Save tiles to disk
    cmd += ["-TILES_SAVE", str(int(save_tiles))]
    if save_tiles and output_dir:
        cmd += ["-TILES_PATH", output_dir]
        cmd += ["-TILES_NAME", base_name]

    # Overlap
    cmd += ["-OVERLAP", str(overlap)]
    cmd += ["-OVERLAP_SYM", str(overlap_sym)]

    # Tiling method
    cmd += ["-METHOD", str(method)]

    if method == 0:
        cmd += ["-NX", str(nx)]
        cmd += ["-NY", str(ny)]
    else:
        if None in [xrange_min, xrange_max, yrange_min, yrange_max]:
            raise ValueError("XRANGE and YRANGE values are required for method=1")
        cmd += [
            "-XRANGE_MIN", str(xrange_min),
            "-XRANGE_MAX", str(xrange_max),
            "-YRANGE_MIN", str(yrange_min),
            "-YRANGE_MAX", str(yrange_max),
            "-DCELL", str(dcell),
            "-DX", str(dx),
            "-DY", str(dy),
        ]

    print("Running SAGA Tiling...")
    start = time.time()
    subprocess.run(cmd, check=True)
    print(f"Tiling completed in {time.time() - start:.2f} seconds.")


# ############################################################################# 
# #############################################################################
# #############################################################################

#latter to find use for it
# ############################################################################# 
# #############################################################################
# #############################################################################

def tile_grid_saga(
    grid_input: str,
    grid_output: str,
    method: int = 0,
    nx: int = 100,
    ny: int = 100,
    x_range: tuple = None,
    y_range: tuple = None,
    dcell: float = 1.0,
    dx: float = 100.0,
    dy: float = 100.0,
    overlap: int = 0,
    overlap_sym: int = 0,
    save_tiles: bool = False,
    tiles_path: str = None,
    tiles_name: str = "tile"
):
    """
    Run the SAGA-GIS Tiling tool on a grid to split it into tiles.

    Use case:
    Use this function to split a large raster into smaller grid tiles for 
    parallel processing or memory-efficient analysis.

    Parameters:
    ----------
    grid_input : str
        Path to input grid file (.sgrd or .tif)
    grid_output : str
        Path to output grid list or directory for generated tiles
    method : int, optional
        0: tile by number of cells (nx, ny),
        1: tile by coordinate ranges and sizes
    nx : int, optional
        Number of cells in x-direction (used if method=0)
    ny : int, optional
        Number of cells in y-direction (used if method=0)
    x_range : tuple, optional
        (min_x, max_x), used if method=1
    y_range : tuple, optional
        (min_y, max_y), used if method=1
    dcell : float, optional
        Cell size (used if method=1)
    dx : float, optional
        Tile size in x (map units, used if method=1)
    dy : float, optional
        Tile size in y (map units, used if method=1)
    overlap : int, optional
        Number of overlapping cells between tiles
    overlap_sym : int, optional
        0: symmetric, 1: bottom/left, 2: top/right
    save_tiles : bool, optional
        If True, saves each tile to disk
    tiles_path : str, optional
        Output directory where tiles are saved (required if save_tiles is True)
    tiles_name : str, optional
        Base filename for output tiles
    """
    cmd = [
        "saga_cmd", "grid_tools", "27",
        "-GRID", grid_input,
        "-TILES", grid_output,
        "-METHOD", str(method),
        "-OVERLAP", str(overlap),
        "-OVERLAP_SYM", str(overlap_sym),
        "-TILES_SAVE", str(int(save_tiles)),
        "-TILES_NAME", tiles_name,
    ]

    if save_tiles and tiles_path:
        cmd.extend(["-TILES_PATH", tiles_path])

    if method == 0:
        # Method 0: tile by number of cells
        cmd.extend(["-NX", str(nx), "-NY", str(ny)])
    elif method == 1:
        # Method 1: tile by coordinates and map units
        if not (x_range and y_range):
            raise ValueError("x_range and y_range must be provided for method 1")
        cmd.extend([
            "-XRANGE_MIN", str(x_range[0]),
            "-XRANGE_MAX", str(x_range[1]),
            "-YRANGE_MIN", str(y_range[0]),
            "-YRANGE_MAX", str(y_range[1]),
            "-DCELL", str(dcell),
            "-DX", str(dx),
            "-DY", str(dy),
        ])

    start = time.time()
    subprocess.run(cmd, check=True)
    end = time.time()
    print(f"Tiling completed in {end - start:.2f} seconds.")
