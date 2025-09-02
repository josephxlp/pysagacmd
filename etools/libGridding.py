""" 
libInterp or libGridding
"""
import subprocess
import os

def run_saga_adw_interpolation(
    points_path, field, output_grid, outdir,
    cellsize=10.0, search_radius=1000, 
    search_range=1, search_points_all=1,
    search_points_min=1, search_points_max=20,
    dw_weighting=1, dw_idw_power=2.0, dw_bandwidth=1.0,
    cv_method=0, cv_samples=10,
    target_definition=0, target_extent=None,
    target_template=None,
    overwrite=False
):
    """
    Python wrapper to run SAGA-GIS Angular Distance Weighted (ADW) interpolation.

    Required arguments:
        points_path (str): Path to input points (SHP or compatible).
        field (str): Attribute to interpolate.
        output_grid (str): Output .sgrd or .tif file path.
        outdir (str): Output folder (created if missing).
    
    Optional parameters (with common DEM settings):
        cellsize (float): Grid cell size (e.g. 10 for 10m DEM).
        search_radius (float): Max search distance (default: 1000).
        search_range (int): 0 = local, 1 = global.
        search_points_all (int): 0 = max N points, 1 = all in distance.
        search_points_min (int): Minimum number of neighbors (default 1).
        search_points_max (int): Maximum number of neighbors (default 20).
        dw_weighting (int): 0=None, 1=IDW, 2=Exp, 3=Gauss. Recommend 1.
        dw_idw_power (float): IDW power, e.g., 2.0.
        dw_bandwidth (float): For exp/gauss weighting.
        cv_method (int): 0=None, 1=LOO, 2=2-fold, 3=k-fold.
        cv_samples (int): Used only if cv_method=3.
        target_definition (int): 0=User defined extent, 1=From template.
        target_extent (dict): Required if target_definition==0:
            {'xmin': float, 'xmax': float, 'ymin': float, 'ymax': float}
        target_template (str): Path to template grid (if used).
        overwrite (bool): If True, overwrite output.

    Returns:
        None
    """
    os.makedirs(outdir, exist_ok=True)
    if os.path.exists(output_grid) and not overwrite:
        print(f"Output {output_grid} already exists. Skipping.")
        return

    cmd = [
        "saga_cmd", "grid_gridding", "7",
        "-POINTS", points_path,
        "-FIELD", field,
        "-TARGET_OUT_GRID", output_grid,
        "-TARGET_DEFINITION", str(target_definition),
        "-TARGET_USER_SIZE", str(cellsize),
        "-SEARCH_RANGE", str(search_range),
        "-SEARCH_RADIUS", str(search_radius),
        "-SEARCH_POINTS_ALL", str(search_points_all),
        "-SEARCH_POINTS_MIN", str(search_points_min),
        "-SEARCH_POINTS_MAX", str(search_points_max),
        "-DW_WEIGHTING", str(dw_weighting),
        "-DW_IDW_POWER", str(dw_idw_power),
        "-DW_BANDWIDTH", str(dw_bandwidth),
        "-CV_METHOD", str(cv_method),
        "-CV_SAMPLES", str(cv_samples)
    ]

    # Add grid extent if user-defined
    if target_definition == 0:
        cmd += [
            "-TARGET_USER_XMIN", str(target_extent['xmin']),
            "-TARGET_USER_XMAX", str(target_extent['xmax']),
            "-TARGET_USER_YMIN", str(target_extent['ymin']),
            "-TARGET_USER_YMAX", str(target_extent['ymax']),
        ]
        cols = int((target_extent['xmax'] - target_extent['xmin']) / cellsize)
        rows = int((target_extent['ymax'] - target_extent['ymin']) / cellsize)
        cmd += ["-TARGET_USER_COLS", str(cols), "-TARGET_USER_ROWS", str(rows)]
    elif target_template:
        cmd += ["-TARGET_TEMPLATE", target_template]

    print("Running command:")
    print(" ".join(cmd))

    subprocess.run(cmd, check=True)

