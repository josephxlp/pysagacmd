import subprocess
import os

def run_saga_dtm_filter_custom(input_grid, radius=2, terrain_slope=30.0, use_confidence_interval=False, 
                              output_dir="output", overwrite=False):
    """
    Run SAGA-GIS DTM Filter (slope-based) with custom parameters and dynamic output filenames.
    
    Parameters:
    - input_grid (str): Path to the input grid file (e.g., DSM in SAGA grid format).
    - radius (int): Search radius of kernel in cells (minimum: 1). Default: 2.
    - terrain_slope (float): Approximate terrain slope in percent (minimum: 0.0). Default: 30.0.
    - use_confidence_interval (bool): Use 5 percent confidence interval. Default: False.
    - output_dir (str): Directory for output files. Default: 'output'.
    - overwrite (bool): If True, overwrite existing output files. Default: False.
    
    Returns:
    - dict: Result containing status, message, and output file paths.
    """
    try:
        # Validate parameters
        if radius < 1:
            return {"status": "error", "message": "Radius must be at least 1."}
        if terrain_slope < 0.0:
            return {"status": "error", "message": "Terrain slope must be non-negative."}
        if not os.path.exists(input_grid):
            return {"status": "error", "message": f"Input grid file {input_grid} does not exist."}

        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)

        # Generate dynamic output filenames
        input_filename = os.path.splitext(os.path.basename(input_grid))[0]
        confidence_suffix = "_conf" if use_confidence_interval else ""
        ground_output = os.path.join(
            output_dir, 
            f"{input_filename}_ground_r{radius}_s{terrain_slope:.1f}{confidence_suffix}.sgrd"
        )
        nonground_output = os.path.join(
            output_dir, 
            f"{input_filename}_nonground_r{radius}_s{terrain_slope:.1f}{confidence_suffix}.sgrd"
        )

        # Check if output files exist
        if not overwrite and (os.path.exists(ground_output) or os.path.exists(nonground_output)):
            return {
                "status": "skipped",
                "message": f"Output files already exist: {ground_output}, {nonground_output}. Set overwrite=True to regenerate.",
                "ground_output": ground_output,
                "nonground_output": nonground_output
            }

        # Construct the saga_cmd command
        cmd = [
            "saga_cmd", "grid_filter", "7",
            "-INPUT", input_grid,
            "-RADIUS", str(radius),
            "-TERRAINSLOPE", str(terrain_slope),
            "-STDDEV", "1" if use_confidence_interval else "0",
            "-GROUND", ground_output,
            "-NONGROUND", nonground_output
        ]

        # Execute the command
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
        
        return {
            "status": "success",
            "message": f"DTM Filter completed successfully. Outputs saved to {ground_output} and {nonground_output}.",
            "ground_output": ground_output,
            "nonground_output": nonground_output,
            "stdout": result.stdout
        }
    
    except subprocess.CalledProcessError as e:
        return {
            "status": "error",
            "message": f"Error running SAGA-GIS: {e.stderr}",
            "ground_output": ground_output,
            "nonground_output": nonground_output
        }
    except FileNotFoundError:
        return {
            "status": "error",
            "message": "saga_cmd not found. Ensure SAGA-GIS is installed and added to the system PATH.",
            "ground_output": ground_output,
            "nonground_output": nonground_output
        }
    except Exception as e:
        return {
            "status": "error",
            "message": f"Unexpected error: {str(e)}",
            "ground_output": ground_output,
            "nonground_output": nonground_output
        }
    
""" 
rad_list = [12,21,32,64]#[1,2,4,8]
slp_list = [1,8, 15,30]
cit_list = [False, True]
c = 0
for radius in rad_list:
    for slope in slp_list:
        for ci in cit_list:
            print(f'ruining id{c} radius{radius} slope {slope} ci {ci}')
            run_saga_dtm_filter_custom(input_grid=edem_fn, 
                           radius=radius, terrain_slope=slope, 
                           use_confidence_interval=ci, 
                            output_dir=f1dir, overwrite=False)
            c +=1

print(c)
"""
