import os
import subprocess

def saga_close_gaps(
        saga_cmd_path, 
        input_grid, 
        mask_grid=None, 
        output_grid=None, 
        tension_threshold=0.1, 
        overwrite=False):
    """
    Run SAGA GIS Close Gaps tool using saga_cmd via Python subprocess.

    Parameters:
        saga_cmd_path (str): Path to the saga_cmd executable.
        input_grid (str): Path to the input grid file with gaps (no data values).
        mask_grid (str, optional): Path to the mask grid file. Defaults to None.
        output_grid (str, optional): Path to save the output grid with gaps closed. If not provided,
                                     changes will be applied to the original grid.
        tension_threshold (float, optional): Tension threshold for interpolation. Defaults to 0.1.
        overwrite (bool, optional): Whether to overwrite the output file if it already exists. Defaults to False.
    """
    # Check if the output file exists and handle overwrite logic
    if output_grid and os.path.isfile(output_grid):
        if overwrite:
            print(f"Output file '{output_grid}' exists but overwrite is enabled. Proceeding...")
        else:
            print(f"Output file '{output_grid}' already exists and overwrite is disabled. Skipping execution.")
            return

    # Construct the command as a list of arguments
    command = [
        saga_cmd_path,
        "grid_tools", "7",  # Tool ID for Close Gaps
        "-INPUT", input_grid,
        "-THRESHOLD", str(tension_threshold)
    ]

    # Add optional parameters if provided
    if mask_grid:
        command.extend(["-MASK", mask_grid])
    if output_grid:
        command.extend(["-RESULT", output_grid])

    # Execute the command using subprocess
    try:
        print("Running SAGA GIS Close Gaps command...")
        result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        # Print the output and error (if any)
        print("SAGA GIS Command Output:")
        print(result.stdout)
        if result.stderr:
            print("SAGA GIS Command Errors:")
            print(result.stderr)

    except subprocess.CalledProcessError as e:
        print("An error occurred while running the SAGA GIS Close Gaps command.")
        print(f"Return Code: {e.returncode}")
        print(f"Error Output: {e.stderr}")


def saga_close_gaps_with_spline(
    saga_cmd_path,
    input_grid,
    mask_grid=None,
    output_grid=None,
    max_gap_cells=0,
    max_points=1000,
    local_points=20,
    extended_neighbourhood=False,
    neighbours="Neumann",
    radius=0,
    relaxation=0.0,
    overwrite=False
):
    """
    Run SAGA GIS Close Gaps with Spline tool using saga_cmd via Python subprocess.

    Parameters:
        saga_cmd_path (str): Path to the saga_cmd executable.
        input_grid (str): Path to the input grid file with gaps (no data values).
        mask_grid (str, optional): Path to the mask grid file. Defaults to None.
        output_grid (str, optional): Path to save the output grid with gaps closed. If not provided,
                                     changes will be applied to the original grid.
        max_gap_cells (int, optional): Only process gaps with fewer cells than this value. Ignored if set to 0. Defaults to 0.
        max_points (int, optional): Maximum points used for interpolation. Defaults to 1000.
        local_points (int, optional): Number of points for local interpolation. Defaults to 20.
        extended_neighbourhood (bool, optional): Whether to use an extended neighbourhood. Defaults to False.
        neighbours (str, optional): Neighbourhood type. Options: "Neumann" or "Moore". Defaults to "Neumann".
        radius (int, optional): Radius (in cells) for neighbourhood consideration. Defaults to 0.
        relaxation (float, optional): Relaxation parameter for spline interpolation. Defaults to 0.0.
        overwrite (bool, optional): Whether to overwrite the output file if it already exists. Defaults to False.
    """
    # Check if the output file exists and handle overwrite logic
    if output_grid and os.path.isfile(output_grid):
        if overwrite:
            print(f"Output file '{output_grid}' exists but overwrite is enabled. Proceeding...")
        else:
            print(f"Output file '{output_grid}' already exists and overwrite is disabled. Skipping execution.")
            return

    # Map neighbour options to their corresponding integer values
    neighbour_options = {"Neumann": 0, "Moore": 1}
    if neighbours not in neighbour_options:
        raise ValueError(f"Invalid neighbours option: {neighbours}. Choose from 'Neumann' or 'Moore'.")
    neighbours_value = neighbour_options[neighbours]

    # Construct the command as a list of arguments
    command = [
        saga_cmd_path,
        "grid_tools", "25",  # Tool ID for Close Gaps with Spline
        "-GRID", input_grid,
        "-MAXGAPCELLS", str(max_gap_cells),
        "-MAXPOINTS", str(max_points),
        "-LOCALPOINTS", str(local_points),
        "-EXTENDED", "1" if extended_neighbourhood else "0",
        "-NEIGHBOURS", str(neighbours_value),
        "-RADIUS", str(radius),
        "-RELAXATION", str(relaxation)
    ]

    # Add optional parameters if provided
    if mask_grid:
        command.extend(["-MASK", mask_grid])
    if output_grid:
        command.extend(["-CLOSED", output_grid])

    # Execute the command using subprocess
    try:
        print("Running SAGA GIS Close Gaps with Spline command...")
        result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        # Print the output and error (if any)
        print("SAGA GIS Command Output:")
        print(result.stdout)
        if result.stderr:
            print("SAGA GIS Command Errors:")
            print(result.stderr)

    except subprocess.CalledProcessError as e:
        print("An error occurred while running the SAGA GIS Close Gaps with Spline command.")
        print(f"Return Code: {e.returncode}")
        print(f"Error Output: {e.stderr}")


def saga_close_gaps_with_resampling(
    saga_cmd_path,
    input_grid,
    mask_grid=None,
    output_grid=None,
    resampling="B-Spline Interpolation",
    grow_factor=2.0,
    overwrite=False
):
    """
    Run SAGA GIS Close Gaps with Stepwise Resampling tool using saga_cmd via Python subprocess.

    Parameters:
        saga_cmd_path (str): Path to the saga_cmd executable.
        input_grid (str): Path to the input grid file with gaps (no data values).
        mask_grid (str, optional): Path to the mask grid file. Defaults to None.
        output_grid (str, optional): Path to save the output grid with gaps closed. If not provided,
                                     changes will be applied to the original grid.
        resampling (str, optional): Resampling method. Options: 
                                    "Nearest Neighbour", "Bilinear Interpolation", 
                                    "Bicubic Spline Interpolation", "B-Spline Interpolation".
                                    Defaults to "B-Spline Interpolation".
        grow_factor (float, optional): Grow factor for stepwise resampling. Minimum: 1.0. Defaults to 2.0.
        overwrite (bool, optional): Whether to overwrite the output file if it already exists. Defaults to False.
    """
    # Check if the output file exists and handle overwrite logic
    if output_grid and os.path.isfile(output_grid):
        if overwrite:
            print(f"Output file '{output_grid}' exists but overwrite is enabled. Proceeding...")
        else:
            print(f"Output file '{output_grid}' already exists and overwrite is disabled. Skipping execution.")
            return

    # Map resampling options to their corresponding integer values
    resampling_options = {
        "Nearest Neighbour": 0,
        "Bilinear Interpolation": 1,
        "Bicubic Spline Interpolation": 2,
        "B-Spline Interpolation": 3
    }
    if resampling not in resampling_options:
        raise ValueError(f"Invalid resampling option: {resampling}. Choose from {list(resampling_options.keys())}.")
    resampling_value = resampling_options[resampling]

    # Validate grow factor
    if grow_factor < 1.0:
        raise ValueError("Grow factor must be at least 1.0.")

    # Construct the command as a list of arguments
    command = [
        saga_cmd_path,
        "grid_tools", "29",  # Tool ID for Close Gaps with Stepwise Resampling
        "-INPUT", input_grid,
        "-RESAMPLING", str(resampling_value),
        "-GROW", str(grow_factor)
    ]

    # Add optional parameters if provided
    if mask_grid:
        command.extend(["-MASK", mask_grid])
    if output_grid:
        command.extend(["-RESULT", output_grid])

    # Execute the command using subprocess
    try:
        print("Running SAGA GIS Close Gaps with Stepwise Resampling command...")
        result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        # Print the output and error (if any)
        print("SAGA GIS Command Output:")
        print(result.stdout)
        if result.stderr:
            print("SAGA GIS Command Errors:")
            print(result.stderr)

    except subprocess.CalledProcessError as e:
        print("An error occurred while running the SAGA GIS Close Gaps with Stepwise Resampling command.")
        print(f"Return Code: {e.returncode}")
        print(f"Error Output: {e.stderr}")


def saga_close_one_cell_gaps(
    saga_cmd_path,
    input_grid,
    output_grid=None,
    neighbourhood="Moore",
    method="arithmetic mean",
    overwrite=False
):
    """
    Run SAGA GIS Close One Cell Gaps tool using saga_cmd via Python subprocess.

    Parameters:
        saga_cmd_path (str): Path to the saga_cmd executable.
        input_grid (str): Path to the input grid file with gaps (no data values).
        output_grid (str, optional): Path to save the output grid with gaps closed. If not provided,
                                     changes will be applied to the original grid.
        neighbourhood (str, optional): Neighbourhood type. Options: "Neumann" or "Moore". Defaults to "Moore".
        method (str, optional): Method for filling gaps. Options: 
                                "arithmetic mean", "median", "majority", "minority". Defaults to "arithmetic mean".
        overwrite (bool, optional): Whether to overwrite the output file if it already exists. Defaults to False.
    """
    # Check if the output file exists and handle overwrite logic
    if output_grid and os.path.isfile(output_grid):
        if overwrite:
            print(f"Output file '{output_grid}' exists but overwrite is enabled. Proceeding...")
        else:
            print(f"Output file '{output_grid}' already exists and overwrite is disabled. Skipping execution.")
            return

    # Map neighbourhood options to their corresponding integer values
    neighbourhood_options = {"Neumann": 0, "Moore": 1}
    if neighbourhood not in neighbourhood_options:
        raise ValueError(f"Invalid neighbourhood option: {neighbourhood}. Choose from 'Neumann' or 'Moore'.")
    neighbourhood_value = neighbourhood_options[neighbourhood]

    # Map method options to their corresponding integer values
    method_options = {
        "arithmetic mean": 0,
        "median": 1,
        "majority": 2,
        "minority": 3
    }
    if method not in method_options:
        raise ValueError(f"Invalid method option: {method}. Choose from {list(method_options.keys())}.")
    method_value = method_options[method]

    # Construct the command as a list of arguments
    command = [
        saga_cmd_path,
        "grid_tools", "6",  # Tool ID for Close One Cell Gaps
        "-INPUT", input_grid,
        "-MODE", str(neighbourhood_value),
        "-METHOD", str(method_value)
    ]

    # Add optional parameters if provided
    if output_grid:
        command.extend(["-RESULT", output_grid])

    # Execute the command using subprocess
    try:
        print("Running SAGA GIS Close One Cell Gaps command...")
        result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        # Print the output and error (if any)
        print("SAGA GIS Command Output:")
        print(result.stdout)
        if result.stderr:
            print("SAGA GIS Command Errors:")
            print(result.stderr)

    except subprocess.CalledProcessError as e:
        print("An error occurred while running the SAGA GIS Close One Cell Gaps command.")
        print(f"Return Code: {e.returncode}")
        print(f"Error Output: {e.stderr}")


def saga_patching(
    saga_cmd_path,
    original_grid,
    patch_grid,
    output_grid=None,
    resampling="B-Spline Interpolation",
    overwrite=False
):
    """
    Run SAGA GIS Patching tool using saga_cmd via Python subprocess.

    Parameters:
        saga_cmd_path (str): Path to the saga_cmd executable.
        original_grid (str): Path to the input grid file with gaps (no data values).
        patch_grid (str): Path to the patch grid file used to fill gaps in the original grid.
        output_grid (str, optional): Path to save the patched grid. If not provided,
                                     changes will be applied to the original grid.
        resampling (str, optional): Resampling method. Options: 
                                    "Nearest Neighbour", "Bilinear Interpolation", 
                                    "Bicubic Spline Interpolation", "B-Spline Interpolation".
                                    Defaults to "B-Spline Interpolation".
        overwrite (bool, optional): Whether to overwrite the output file if it already exists. Defaults to False.
    """
    # Check if the output file exists and handle overwrite logic
    if output_grid and os.path.isfile(output_grid):
        if overwrite:
            print(f"Output file '{output_grid}' exists but overwrite is enabled. Proceeding...")
        else:
            print(f"Output file '{output_grid}' already exists and overwrite is disabled. Skipping execution.")
            return

    # Map resampling options to their corresponding integer values
    resampling_options = {
        "Nearest Neighbour": 0,
        "Bilinear Interpolation": 1,
        "Bicubic Spline Interpolation": 2,
        "B-Spline Interpolation": 3
    }
    if resampling not in resampling_options:
        raise ValueError(f"Invalid resampling option: {resampling}. Choose from {list(resampling_options.keys())}.")
    resampling_value = resampling_options[resampling]

    # Construct the command as a list of arguments
    command = [
        saga_cmd_path,
        "grid_tools", "5",  # Tool ID for Patching
        "-ORIGINAL", original_grid,
        "-ADDITIONAL", patch_grid,
        "-RESAMPLING", str(resampling_value)
    ]

    # Add optional parameters if provided
    if output_grid:
        command.extend(["-COMPLETED", output_grid])

    # Execute the command using subprocess
    try:
        print("Running SAGA GIS Patching command...")
        result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        # Print the output and error (if any)
        print("SAGA GIS Command Output:")
        print(result.stdout)
        if result.stderr:
            print("SAGA GIS Command Errors:")
            print(result.stderr)

    except subprocess.CalledProcessError as e:
        print("An error occurred while running the SAGA GIS Patching command.")
        print(f"Return Code: {e.returncode}")
        print(f"Error Output: {e.stderr}")