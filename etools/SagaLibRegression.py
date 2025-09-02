"""
https://saga-gis.sourceforge.io/saga_tool_doc/8.3.0/statistics_regression.html

Library Regression
-------------------------------------
Cellwise Trend for Grids [x]
GWR for Grid Downscaling [x]
GWR for Multiple Predictor Grids [x]
GWR for Multiple Predictors (Gridded Model Output) [x]
GWR for Multiple Predictors (Shapes) [-]
GWR for Single Predictor (Gridded Model Output) [-]
GWR for Single Predictor Grid [x]
Multiple Linear Regression Analysis (Shapes) [-]
Multiple Linear Regression Analysis (Table) [-]
Multiple Regression Analysis (Grid and Predictor Grids) [x]
Multiple Regression Analysis (Points and Predictor Grids) [-]
Polynomial Regression [-]
Regression Analysis (Points and Predictor Grid) [-]
Trend Analysis [-]
Trend Analysis (Shapes) [-]
Zonal Multiple Regression Analysis (Points and Predictor Grids) [-]
"""
# test this with a sample ROI with gaps and no gaps
#- TLS and MKD # gaps of different sizes 
# add the preprocessing and postprocessing steps [saga # geotif]
# add prprocessing for shapes, like tif to gdf, to vector 
import os
import time
import subprocess

""" 
give me an example using the functino above
add one other argument call outdir : make sure a folder is crated on this path 
all the ouputs should be saved in the output directory
all the output name should be modified inside the function to take all the other arguments so its easier to track each experiment
check if the output exist, if not run, else skip. add overwrite default to False, but if True runs regardless 
# keep the comments inside the fuction in cmd
returns sdat path
"""


# ############################################################################# 
# #############################################################################
# #############################################################################

#latter to find use for it
# ############################################################################# 
# #############################################################################
# #############################################################################

def run_gwr_grid_downscaling(
        # for this one add also a function to make prediction using the model 
    predictors, dependent, regression, quality, residuals,
    log_reg=0,
    model_out=0,
    search_range=0,
    search_radius=10,
    dw_weighting=0,
    dw_idw_power=2.0,
    dw_bandwidth=1.0,
    reg_rescorr=None,
    model=None,
):
    """
    Runs SAGA-GIS Geographically Weighted Regression (GWR) for grid downscaling.

    Use case:
    Perform GWR using grid predictors and dependent variable grids.
    Outputs regression, residuals, coefficient of determination, and optionally model parameters.
    Supports logistic regression and various spatial weighting schemes.
    """
    cmd = [
        "saga_cmd", "statistics_regression", "14",

        # Required inputs/outputs
        "-PREDICTORS", predictors,          # Grid list of predictor variables; required
        "-DEPENDENT", dependent,            # Grid for dependent variable; required
        "-REGRESSION", regression,          # Output regression grid; required
        "-QUALITY", quality,                # Output coefficient of determination grid; required
        "-RESIDUALS", residuals,            # Output residuals grid; required

        # Optional flags and outputs
        "-LOGISTIC", str(log_reg),          # Logistic regression: 0=no(default), 1=yes
        "-MODEL_OUT", str(model_out),        # Output model parameters grids: 0=no(default), 1=yes
        "-SEARCH_RANGE", str(search_range),  # Search range: 0=local(default), 1=global
        "-SEARCH_RADIUS", str(search_radius), # Search radius in cells; min=1; default=10
        "-DW_WEIGHTING", str(dw_weighting),  # Distance weighting: 0=none(default), 1=IDW, 2=exponential, 3=gaussian
        "-DW_IDW_POWER", str(dw_idw_power),  # Power for IDW weighting; min=0; default=2.0
        "-DW_BANDWIDTH", str(dw_bandwidth),  # Bandwidth for exponential/gaussian; min=0; default=1.0
    ]
    if reg_rescorr:
        cmd.extend(["-REG_RESCORR", reg_rescorr])  # Optional regression residual correction output
    if model:
        cmd.extend(["-MODEL", model])              # Optional regression parameters output grids

    start = time.time()
    subprocess.run(cmd, check=True)
    end = time.time()
    print(f"GWR Grid Downscaling completed in {end - start:.2f} seconds")


def run_cellwise_trend(
    y_grids,                # [required] Path to input grid list (dependent variable)
    coeff,                  # [required] Path to output grid list (polynomial coefficients)

    r2=None,                # [optional] Path to output grid (R² values)
    r2adj=None,             # [optional] Path to output grid (adjusted R² values)
    stderr=None,            # [optional] Path to output grid (standard error)
    p=None,                 # [optional] Path to output grid (significance level)

    linear=1,               # [optional] Use linear trend? (1=True, 0=False)
    order=2,                # [optional] Polynomial order (default=2)
    xsource=0,              # [optional] Independent variable source (0=list index, 1=table, 2=grids)

    x_table=None,           # [optional] Path to table with independent variable values (if xsource=1)
    x_grids=None,           # [optional] Grid list of independent variable (if xsource=2)

    outdir=".",             # [optional] Output directory to save result files
):
    """
    Run SAGA-GIS Cellwise Trend fitting polynomial trends per cell.

    Parameters:
        y_grids   : Input dependent variable grid list (.sgrd/.sdat)
        coeff     : Output grid list for polynomial coefficients
        r2        : Output R² grid
        r2adj     : Output adjusted R² grid
        stderr    : Output standard error grid
        p         : Output significance level grid
        linear    : Whether to use linear trend (1) or not (0)
        order     : Polynomial order (e.g., 1 for linear, 2 for quadratic)
        xsource   : Independent variable source: 0=list index, 1=table, 2=grid list
        x_table   : Table for X if xsource = 1
        x_grids   : Grid list for X if xsource = 2
        outdir    : Directory to save all output files
    """

    # Ensure output directory exists
    os.makedirs(outdir, exist_ok=True)

    cmd = [
        "saga_cmd", "statistics_regression", "9",
        "-Y_GRIDS", y_grids,
        "-COEFF", os.path.join(outdir, coeff),
        "-LINEAR", str(linear),
        "-ORDER", str(order),
        "-XSOURCE", str(xsource),
    ]

    if r2:
        cmd.extend(["-R2", os.path.join(outdir, r2)])
    if r2adj:
        cmd.extend(["-R2ADJ", os.path.join(outdir, r2adj)])
    if stderr:
        cmd.extend(["-STDERR", os.path.join(outdir, stderr)])
    if p:
        cmd.extend(["-P", os.path.join(outdir, p)])
    if x_table:
        cmd.extend(["-X_TABLE", x_table])
    if x_grids:
        cmd.extend(["-X_GRIDS", x_grids])

    start = time.time()
    subprocess.run(cmd, check=True)
    print(f"✅ Cellwise Trend completed in {time.time() - start:.2f} seconds")


def run_gwr_multiple_predictor_grids(
    points, dependent, predictors, regression, quality,
    residuals=None,
    model=None,
    log_reg=0,
    model_out=0,
    resolution=1,
    resolution_val=1.0,
    dw_weighting=0,
    dw_idw_power=2.0,
    dw_bandwidth=1.0,
    search_range=1,
    search_radius=1000.0,
    search_points_all=1,
    search_points_min=16,
    search_points_max=20,
    search_direction=0,
):
    """
    Runs SAGA-GIS GWR for multiple predictor grids using points input.
    
    Use case:
    Spatial regression modeling with multiple predictor grids and points data.
    Outputs regression grid, quality grid, optionally residuals and model parameters.
    Supports logistic regression, spatial weighting, and user-defined resolution.
    """
    cmd = [
        "saga_cmd", "statistics_regression", "6",

        # Required inputs/outputs
        "-POINTS", points,                  # Input points shapefile; required
        "-DEPENDENT", dependent,            # Dependent variable field; required
        "-PREDICTORS", predictors,          # Predictor grids list; required
        "-REGRESSION", regression,          # Output regression grid; required
        "-QUALITY", quality,                # Output coefficient of determination grid; required

        # Optional outputs
    ]
    if residuals:
        cmd.extend(["-RESIDUALS", residuals])  # Residuals output shapes; optional
    if model:
        cmd.extend(["-MODEL", model])          # Regression parameters grids list; optional

    # Options and flags
    cmd.extend([
        "-LOGISTIC", str(log_reg),              # Logistic regression: 0=no(default), 1=yes
        "-MODEL_OUT", str(model_out),           # Output regression parameters: 0=no(default), 1=yes
        "-RESOLUTION", str(resolution),         # Model resolution: 0=same as predictors, 1=user defined(default)
        "-RESOLUTION_VAL", str(resolution_val), # Resolution value in map units; min=0, default=1.0
        "-DW_WEIGHTING", str(dw_weighting),     # Weighting function: 0=none(default), 1=IDW, 2=exponential, 3=gaussian
        "-DW_IDW_POWER", str(dw_idw_power),     # Power for IDW; min=0, default=2.0
        "-DW_BANDWIDTH", str(dw_bandwidth),     # Bandwidth for exp/gauss; min=0, default=1.0
        "-SEARCH_RANGE", str(search_range),     # Search range: 0=local, 1=global(default)
        "-SEARCH_RADIUS", str(search_radius),   # Max search distance in map units; min=0, default=1000.0
        "-SEARCH_POINTS_ALL", str(search_points_all), # Number of points: 0=max nearest, 1=all within distance(default)
        "-SEARCH_POINTS_MIN", str(search_points_min), # Minimum points; min=1, default=16
        "-SEARCH_POINTS_MAX", str(search_points_max), # Maximum points; min=1, default=20
        "-SEARCH_DIRECTION", str(search_direction)    # Direction: 0=all(default), 1=quadrants
    ])

    start = time.time()
    subprocess.run(cmd, check=True)
    end = time.time()
    print(f"GWR Multiple Predictor Grids completed in {end - start:.2f} seconds")


def run_gwr_single_predictor_grid(
    points, dependent, predictor, regression,
    residuals=None, quality=None, intercept=None, slope=None,
    log_reg=0,
    dw_weighting=0,
    dw_idw_power=2.0,
    dw_bandwidth=1.0,
    search_range=1,
    search_radius=1000.0,
    search_points_all=1,
    search_points_min=16,
    search_points_max=20,
    search_direction=0,
):
    """
    Run SAGA-GIS Geographically Weighted Regression (GWR) for a single predictor grid.
    
    Use case:
    Applies GWR using point data and one predictor grid, producing regression output grid,
    optionally residuals, quality (R²), intercept, and slope grids.
    Supports logistic regression and spatial weighting functions.
    """
    cmd = [
        "saga_cmd", "statistics_regression", "4",

        # Required inputs/outputs
        "-POINTS", points,              # Input points shapefile; required
        "-DEPENDENT", dependent,        # Dependent variable field; required
        "-PREDICTOR", predictor,        # Predictor grid; required
        "-REGRESSION", regression,      # Output regression grid; required
    ]
    if residuals:
        cmd.extend(["-RESIDUALS", residuals])   # Residuals output shapes; optional
    if quality:
        cmd.extend(["-QUALITY", quality])       # Coefficient of determination grid; optional
    if intercept:
        cmd.extend(["-INTERCEPT", intercept])   # Intercept grid; optional
    if slope:
        cmd.extend(["-SLOPE", slope])           # Slope grid; optional

    # Options and flags
    cmd.extend([
        "-LOGISTIC", str(log_reg),          # Logistic regression: 0=no(default), 1=yes
        "-DW_WEIGHTING", str(dw_weighting),# Weighting function: 0=none(default), 1=IDW, 2=exp, 3=gauss
        "-DW_IDW_POWER", str(dw_idw_power),# Power for IDW; min=0, default=2.0
        "-DW_BANDWIDTH", str(dw_bandwidth),# Bandwidth for exponential/gaussian; min=0, default=1.0
        "-SEARCH_RANGE", str(search_range),# Search range: 0=local, 1=global(default)
        "-SEARCH_RADIUS", str(search_radius),# Max search distance (map units); min=0, default=1000.0
        "-SEARCH_POINTS_ALL", str(search_points_all), # Number of points: 0=max nearest, 1=all within distance(default)
        "-SEARCH_POINTS_MIN", str(search_points_min), # Minimum points; min=1, default=16
        "-SEARCH_POINTS_MAX", str(search_points_max), # Maximum points; min=1, default=20
        "-SEARCH_DIRECTION", str(search_direction),   # Direction: 0=all(default), 1=quadrants
    ])

    start = time.time()
    subprocess.run(cmd, check=True)
    end = time.time()
    print(f"GWR Single Predictor Grid completed in {end - start:.2f} seconds")


def run_multiple_regression_analysis(
    dependent, predictors, regression,
    residuals=None,
    info_coeff=None,
    info_model=None,
    info_steps=None,
    resampling=3,
    coord_x=0,
    coord_y=0,
    method=3,
    p_value=5.0,
    crossval=0,
    crossval_k=10,
):
    """
    Run SAGA-GIS Multiple Regression Analysis (Grid and Predictor Grids).
    
    Use case:
    Performs linear regression of one dependent grid on multiple predictor grids.
    Supports forward selection, stepwise methods, and cross-validation.
    Outputs regression grid, residuals, and detailed coefficient/model/step tables.
    """
    cmd = [
        "saga_cmd", "statistics_regression", "8",

        # Required inputs/outputs
        "-DEPENDENT", dependent,      # Dependent variable grid; required
        "-PREDICTORS", predictors,    # Predictor grids list; required
        "-REGRESSION", regression,    # Output regression grid; required
    ]
    if residuals:
        cmd.extend(["-RESIDUALS", residuals])       # Residuals grid; optional
    if info_coeff:
        cmd.extend(["-INFO_COEFF", info_coeff])     # Coefficients table; optional
    if info_model:
        cmd.extend(["-INFO_MODEL", info_model])     # Model info table; optional
    if info_steps:
        cmd.extend(["-INFO_STEPS", info_steps])     # Steps info table; optional

    # Options
    cmd.extend([
        "-RESAMPLING", str(resampling),   # Resampling: 0=Nearest,1=Bilinear,2=Bicubic,3=B-Spline (default=3)
        "-COORD_X", str(coord_x),         # Include X coordinate (0=no, default=0)
        "-COORD_Y", str(coord_y),         # Include Y coordinate (0=no, default=0)
        "-METHOD", str(method),            # Method: 0=all,1=forward,2=backward,3=stepwise(default=3)
        "-P_VALUE", str(p_value),          # Significance level (p-value %) for predictor selection (default=5.0)
        "-CROSSVAL", str(crossval),        # Cross validation: 0=none,1=leave-one-out,2=2-fold,3=k-fold(default=0)
        "-CROSSVAL_K", str(crossval_k),    # Subsamples for k-fold CV; min=2 (default=10)
    ])

    start = time.time()
    subprocess.run(cmd, check=True)
    end = time.time()
    print(f"Multiple Regression Analysis completed in {end - start:.2f} seconds")

