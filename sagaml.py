import os 
import subprocess

def run_gwr_downscaling(saga_cmd_path, predictors, dependent_variable, regression_output,
                        quality_output, residuals_output, model_output=None, logistic_regression=False,
                        model_out=False, search_range=0, search_radius=10, weighting_function=3,
                        idw_power=2.0, bandwidth=7.0):
    """
    Run SAGA GIS GWR for Grid Downscaling tool using saga_cmd via Python subprocess.

    Parameters:
        saga_cmd_path (str): Path to the saga_cmd executable.
        predictors (list of str): List of paths to predictor grids.
        dependent_variable (str): Path to the dependent variable grid.
        regression_output (str): Path to save the regression output grid.
        quality_output (str): Path to save the coefficient of determination (R²) grid.
        residuals_output (str): Path to save the residuals grid.
        model_output (list of str, optional): List of paths to save the regression parameters grids.
        logistic_regression (bool, optional): Enable logistic regression. Defaults to False.
        model_out (bool, optional): Output model parameters. Defaults to False.
        search_range (int, optional): Search range type (0: local, 1: global). Defaults to 0.
        search_radius (int, optional): Search distance in cells. Defaults to 10.
        weighting_function (int, optional): Weighting function type. Defaults to 3 (gaussian).
        idw_power (float, optional): Power for inverse distance weighting. Defaults to 2.0.
        bandwidth (float, optional): Bandwidth for exponential/Gaussian weighting. Defaults to 7.0.
    """
    # Construct the command as a list of arguments
    command = [
        saga_cmd_path,
        "statistics_regression", "14",  # Tool ID for GWR for Grid Downscaling
        "-PREDICTORS", ";".join(predictors),  # Predictors are passed as a semicolon-separated list
        "-REGRESSION", regression_output,
        "-DEPENDENT", dependent_variable,
        "-QUALITY", quality_output,
        "-RESIDUALS", residuals_output,
        "-LOGISTIC", "1" if logistic_regression else "0",
        "-MODEL_OUT", "1" if model_out else "0",
        "-SEARCH_RANGE", str(search_range),
        "-SEARCH_RADIUS", str(search_radius),
        "-DW_WEIGHTING", str(weighting_function),
        "-DW_IDW_POWER", str(idw_power),
        "-DW_BANDWIDTH", str(bandwidth)
    ]

    # Add optional model output if provided
    if model_output:
        command.extend(["-MODEL", ";".join(model_output)])

    # Execute the command using subprocess
    try:
        print("Running SAGA GIS GWR for Grid Downscaling command...")
        result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        # Print the output and error (if any)
        print("SAGA GIS Command Output:")
        print(result.stdout)
        if result.stderr:
            print("SAGA GIS Command Errors:")
            print(result.stderr)

    except subprocess.CalledProcessError as e:
        print("An error occurred while running the SAGA GIS GWR for Grid Downscaling command.")
        print(f"Return Code: {e.returncode}")
        print(f"Error Output: {e.stderr}")

def run_terrain_clustering(saga_cmd_path, input_dem, output_clusters, num_classes=5, max_iterations=25):
    """
    Run SAGA GIS Terrain Clustering tool using saga_cmd via Python subprocess.

    Parameters:
        input_dem (str): Path to the input DEM file.
        output_clusters (str): Path to save the clusters output.
        num_classes (int, optional): Number of terrain classes. Defaults to 5.
        max_iterations (int, optional): Maximum iterations for clustering. Defaults to 25.
    """
    # Construct the command as a list of arguments
    command = [
        saga_cmd_path,
        "ta_morphometry", "clustering",
        "-ELEVATION", input_dem,
        "-CLUSTER", output_clusters,
        "-NCLUSTER", str(num_classes),
        "-MAXITER", str(max_iterations)
    ]

    # Execute the command using subprocess
    try:
        print("Running SAGA GIS Terrain Clustering command...")
        result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        # Print the output and error (if any)
        print("SAGA GIS Command Output:")
        print(result.stdout)
        if result.stderr:
            print("SAGA GIS Command Errors:")
            print(result.stderr)

    except subprocess.CalledProcessError as e:
        print("An error occurred while running the SAGA GIS Terrain Clustering command.")
        print(f"Return Code: {e.returncode}")
        print(f"Error Output: {e.stderr}")