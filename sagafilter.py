import os 
import subprocess 

def saga_filter_dtm_by_slope_threshold(saga_cmd_path,input_dem, output_ground, output_nonground=None, 
                 radius=5, terrain_slope=30.0, filter_mod=0, stddev=0.1):
    # create mine where the slope map is provided to make decision  on the filter 
    # play with slope and other parameters to get the best result here 
    # https://saga-gis.sourceforge.io/saga_tool_doc/9.7.1/grid_filter_7.html 
    """
    Run SAGA GIS grid_filter module 7 using saga_cmd via Python subprocess.

    Parameters:
        input_dem (str): Path to the input DEM file.
        output_ground (str): Path to save the bare earth output.
        output_nonground (str, optional): Path to save the removed objects output. Defaults to None.
        radius (int, optional): Kernel radius. Defaults to 5.
        terrain_slope (float, optional): Terrain slope percentage. Defaults to 30.0.
        filter_mod (int, optional): Filter modification mode (0 = none, 1 = relax, 2 = amplify). Defaults to 0.
        stddev (float, optional): Standard deviation value. Defaults to 0.1.
    """

    # Construct the command as a list of arguments
    command = [
        saga_cmd_path,
        "grid_filter", "7",
        "-INPUT", input_dem,
        "-GROUND", output_ground,
    ]

    # Add optional parameters if provided
    if output_nonground:
        command.extend(["-NONGROUND", output_nonground])
    command.extend([
        "-RADIUS", str(radius),
        "-TERRAINSLOPE", str(terrain_slope),
        "-FILTERMOD", str(filter_mod),
        "-STDDEV", str(stddev)
    ])

    # Execute the command using subprocess
    try:
        print("Running SAGA GIS command...")
        result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        # Print the output and error (if any)
        print("SAGA GIS Command Output:")
        print(result.stdout)
        if result.stderr:
            print("SAGA GIS Command Errors:")
            print(result.stderr)

    except subprocess.CalledProcessError as e:
        print("An error occurred while running the SAGA GIS command.")
        print(f"Return Code: {e.returncode}")
        print(f"Error Output: {e.stderr}")


def saga_calculate_slope_statistics(saga_cmd_path, input_dem, method="avg"):
    """
    Calculate the average or median slope of the input DEM using SAGA GIS.
    
    Parameters:
        saga_cmd_path (str): Path to the saga_cmd executable.
        input_dem (str): Path to the input DEM file.
        method (str): Method for calculating slope statistic ('avg' or 'median').
    
    Returns:
        float: The calculated slope statistic (average or median).
    """
    # Temporary output file for slope grid
    slope_output = "temp_slope.sgrd"
    
    # Command to calculate slope using SAGA GIS
    command = [
        saga_cmd_path,
        "ta_morphometry", "0",  # Slope calculation tool
        "-ELEVATION", input_dem,
        "-SLOPE", slope_output
    ]
    
    try:
        # Run the slope calculation command
        print("Calculating slope...")
        result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if result.stderr:
            print(f"Slope Calculation Errors: {result.stderr}")
        
        # Read the slope grid data
        print("Reading slope grid...")
        slope_data_command = [
            saga_cmd_path,
            "io_gdal", "2",  # Export to GeoTIFF
            "-GRIDS", slope_output,
            "-FILE", "temp_slope.tif"
        ]
        subprocess.run(slope_data_command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        
        # Load the slope data into a NumPy array
        import rasterio
        with rasterio.open("temp_slope.tif") as dataset:
            slope_array = dataset.read(1)
        
        # Remove no-data values
        #slope_array = slope_array[slope_array > 0]
        
        # Calculate the desired statistic
        if method == "avg":
            return np.mean(slope_array)
        elif method == "median":
            return np.median(slope_array)
        else:
            raise ValueError("Invalid method. Use 'avg' or 'median'.")
    
    except subprocess.CalledProcessError as e:
        print("An error occurred while calculating slope.")
        print(f"Return Code: {e.returncode}")
        print(f"Error Output: {e.stderr}")
        return None
    finally:
        # Clean up temporary files
        import os
        for temp_file in ["temp_slope.sgrd", "temp_slope.tif"]:
            if os.path.exists(temp_file):
                os.remove(temp_file)

def saga_filter_dtm_by_slope_map_agg(saga_cmd_path, input_dem, output_ground, output_nonground=None, 
                 radius=5, terrain_slope_method="avg", filter_mod=0, stddev=0.1):
    """
    Run SAGA GIS grid_filter module 7 using saga_cmd via Python subprocess.
    
    Parameters:
        saga_cmd_path (str): Path to the saga_cmd executable.
        input_dem (str): Path to the input DEM file.
        output_ground (str): Path to save the bare earth output.
        output_nonground (str, optional): Path to save the removed objects output. Defaults to None.
        radius (int, optional): Kernel radius. Defaults to 5.
        terrain_slope_method (str, optional): Method to calculate terrain slope ('avg' or 'median'). Defaults to 'avg'.
        filter_mod (int, optional): Filter modification mode (0 = none, 1 = relax, 2 = amplify). Defaults to 0.
        stddev (float, optional): Standard deviation value. Defaults to 0.1.
    """
    # Calculate terrain slope dynamically
    terrain_slope = saga_calculate_slope_statistics(saga_cmd_path, input_dem, method=terrain_slope_method)
    if terrain_slope is None:
        print("Failed to calculate terrain slope. Using default value of 30.0.")
        terrain_slope = 30.0
    
    # Construct the command as a list of arguments
    command = [
        saga_cmd_path,
        "grid_filter", "7",
        "-INPUT", input_dem,
        "-GROUND", output_ground,
    ]
    # Add optional parameters if provided
    if output_nonground:
        command.extend(["-NONGROUND", output_nonground])
    command.extend([
        "-RADIUS", str(radius),
        "-TERRAINSLOPE", str(terrain_slope),
        "-FILTERMOD", str(filter_mod),
        "-STDDEV", str(stddev)
    ])
    
    # Execute the command using subprocess
    try:
        print("Running SAGA GIS command...")
        result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        # Print the output and error (if any)
        print("SAGA GIS Command Output:")
        print(result.stdout)
        if result.stderr:
            print("SAGA GIS Command Errors:")
            print(result.stderr)
    except subprocess.CalledProcessError as e:
        print("An error occurred while running the SAGA GIS command.")
        print(f"Return Code: {e.returncode}")
        print(f"Error Output: {e.stderr}")