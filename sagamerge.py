import os 
import subprocess

def saga_mosaicking(
    saga_cmd_path,
    input_grids,
    output_mosaic,
    target_template=None,
    data_type="same as first grid in list",
    resampling="B-Spline Interpolation",
    overlap="last",
    blend_distance=10.0,
    blend_boundary="valid data cells",
    match="none",
    target_definition="grid or grid system",
    cellsize=1.0,
    west=0.0,
    east=100.0,
    south=0.0,
    north=100.0,
    columns=101,
    rows=101,
    rounding=True,
    fit="nodes",
    overwrite=False
):
    """
    Run SAGA GIS Mosaicking (Grid Collections) tool using saga_cmd via Python subprocess.

    Parameters:
        saga_cmd_path (str): Path to the saga_cmd executable.
        input_grids (list): List of paths to input grid collections to be merged.
        output_mosaic (str): Path to save the output mosaic grid collection.
        target_template (str, optional): Path to a grid file whose system will be used for output grids. Defaults to None.
        data_type (str, optional): Data storage type. Options: 
                                   "bit", "unsigned 1 byte integer", "signed 1 byte integer", etc.
                                   Defaults to "same as first grid in list".
        resampling (str, optional): Resampling method. Options: 
                                    "Nearest Neighbour", "Bilinear Interpolation", 
                                    "Bicubic Spline Interpolation", "B-Spline Interpolation".
                                    Defaults to "B-Spline Interpolation".
        overlap (str, optional): Overlapping areas handling method. Options: 
                                 "first", "last", "minimum", "maximum", "mean", "blend boundary", "feathering".
                                 Defaults to "last".
        blend_distance (float, optional): Blending distance in map units. Minimum: 0.0. Defaults to 10.0.
        blend_boundary (str, optional): Blending boundary for distance calculation. Options: 
                                        "valid data cells", "grid boundaries", "vertical grid boundaries", "horizontal grid boundaries".
                                        Defaults to "valid data cells".
        match (str, optional): Matching method. Options: 
                               "none", "match histogram of first grid in list", "match histogram of overlapping area", "regression".
                               Defaults to "none".
        target_definition (str, optional): Target grid system definition. Options: 
                                           "user defined", "grid or grid system".
                                           Defaults to "grid or grid system".
        cellsize (float, optional): Cell size for user-defined grid system. Minimum: 0.0. Defaults to 1.0.
        west (float, optional): West boundary for user-defined grid system. Defaults to 0.0.
        east (float, optional): East boundary for user-defined grid system. Defaults to 100.0.
        south (float, optional): South boundary for user-defined grid system. Defaults to 0.0.
        north (float, optional): North boundary for user-defined grid system. Defaults to 100.0.
        columns (int, optional): Number of columns for user-defined grid system. Minimum: 1. Defaults to 101.
        rows (int, optional): Number of rows for user-defined grid system. Minimum: 1. Defaults to 101.
        rounding (bool, optional): Whether to round bounding coordinates to multiples of cell size. Defaults to True.
        fit (str, optional): Fit method for user-defined grid system. Options: "nodes", "cells". Defaults to "nodes".
        overwrite (bool, optional): Whether to overwrite the output file if it already exists. Defaults to False.
    """
    # Check if the output file exists and handle overwrite logic
    if os.path.isfile(output_mosaic):
        if overwrite:
            print(f"Output file '{output_mosaic}' exists but overwrite is enabled. Proceeding...")
        else:
            print(f"Output file '{output_mosaic}' already exists and overwrite is disabled. Skipping execution.")
            return

    # Map data type options to their corresponding integer values
    data_type_options = {
        "bit": 0,
        "unsigned 1 byte integer": 1,
        "signed 1 byte integer": 2,
        "unsigned 2 byte integer": 3,
        "signed 2 byte integer": 4,
        "unsigned 4 byte integer": 5,
        "signed 4 byte integer": 6,
        "unsigned 8 byte integer": 7,
        "signed 8 byte integer": 8,
        "4 byte floating point number": 9,
        "8 byte floating point number": 10,
        "same as first grid in list": 11
    }
    if data_type not in data_type_options:
        raise ValueError(f"Invalid data type option: {data_type}. Choose from {list(data_type_options.keys())}.")
    data_type_value = data_type_options[data_type]

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

    # Map overlap options to their corresponding integer values
    overlap_options = {
        "first": 0,
        "last": 1,
        "minimum": 2,
        "maximum": 3,
        "mean": 4,
        "blend boundary": 5,
        "feathering": 6
    }
    if overlap not in overlap_options:
        raise ValueError(f"Invalid overlap option: {overlap}. Choose from {list(overlap_options.keys())}.")
    overlap_value = overlap_options[overlap]

    # Map blend boundary options to their corresponding integer values
    blend_boundary_options = {
        "valid data cells": 0,
        "grid boundaries": 1,
        "vertical grid boundaries": 2,
        "horizontal grid boundaries": 3
    }
    if blend_boundary not in blend_boundary_options:
        raise ValueError(f"Invalid blend boundary option: {blend_boundary}. Choose from {list(blend_boundary_options.keys())}.")
    blend_boundary_value = blend_boundary_options[blend_boundary]

    # Map match options to their corresponding integer values
    match_options = {
        "none": 0,
        "match histogram of first grid in list": 1,
        "match histogram of overlapping area": 2,
        "regression": 3
    }
    if match not in match_options:
        raise ValueError(f"Invalid match option: {match}. Choose from {list(match_options.keys())}.")
    match_value = match_options[match]

    # Map target definition options to their corresponding integer values
    target_definition_options = {
        "user defined": 0,
        "grid or grid system": 1
    }
    if target_definition not in target_definition_options:
        raise ValueError(f"Invalid target definition option: {target_definition}. Choose from {list(target_definition_options.keys())}.")
    target_definition_value = target_definition_options[target_definition]

    # Map fit options to their corresponding integer values
    fit_options = {
        "nodes": 0,
        "cells": 1
    }
    if fit not in fit_options:
        raise ValueError(f"Invalid fit option: {fit}. Choose from {list(fit_options.keys())}.")
    fit_value = fit_options[fit]

    # Construct the command as a list of arguments
    command = [
        saga_cmd_path,
        "grid_tools", "38",  # Tool ID for Mosaicking (Grid Collections)
        "-GRIDS", ";".join(input_grids),
        "-TYPE", str(data_type_value),
        "-RESAMPLING", str(resampling_value),
        "-OVERLAP", str(overlap_value),
        "-BLEND_DIST", str(blend_distance),
        "-BLEND_BND", str(blend_boundary_value),
        "-MATCH", str(match_value),
        "-TARGET_DEFINITION", str(target_definition_value),
        "-TARGET_USER_SIZE", str(cellsize),
        "-TARGET_USER_XMIN", str(west),
        "-TARGET_USER_XMAX", str(east),
        "-TARGET_USER_YMIN", str(south),
        "-TARGET_USER_YMAX", str(north),
        "-TARGET_USER_COLS", str(columns),
        "-TARGET_USER_ROWS", str(rows),
        "-TARGET_USER_FLAT", "1" if rounding else "0",
        "-TARGET_USER_FITS", str(fit_value),
        "-MOSAIC", output_mosaic
    ]

    # Add optional parameters if provided
    if target_template:
        command.extend(["-TARGET_TEMPLATE", target_template])

    # Execute the command using subprocess
    try:
        print("Running SAGA GIS Mosaicking (Grid Collections) command...")
        result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        # Print the output and error (if any)
        print("SAGA GIS Command Output:")
        print(result.stdout)
        if result.stderr:
            print("SAGA GIS Command Errors:")
            print(result.stderr)

    except subprocess.CalledProcessError as e:
        print("An error occurred while running the SAGA GIS Mosaicking (Grid Collections) command.")
        print(f"Return Code: {e.returncode}")
        print(f"Error Output: {e.stderr}")

