import os 
import subprocess 


def saga_terrain_flooding(
    saga_cmd_path,
    dem,
    seed_points,
    water_level,
    output_water_body,
    output_flooded_dem=None,
    water_level_default=0.5,
    level_reference=0,
    constant_level=True
):
    """
    Run the SAGA GIS Terrain Flooding tool using saga_cmd via Python subprocess.

    Parameters:
        saga_cmd_path (str): Path to the saga_cmd executable.
        dem (str): Path to the input digital elevation model grid.
        seed_points (str): Path to the input seed points shapefile.
        water_level (str): Attribute field with the local water level.
        output_water_body (str): Path to save the water body grid output.
        output_flooded_dem (str, optional): Path to save the flooded DEM output.
        water_level_default (float, optional): Default water level if no attribute is selected. Defaults to 0.5.
        level_reference (int, optional): Reference for water level [0: relative to DEM, 1: absolute height]. Defaults to 0.
        constant_level (bool, optional): Whether to model the water level as constant. Defaults to True.
    """
    # Construct the command
    command = [
        saga_cmd_path,
        "ta_hydrology", "32",  # Tool ID for Terrain Flooding
        "-DEM", dem,
        "-SEEDS", seed_points,
        "-WATER_LEVEL", water_level,
        "-WATER_LEVEL_DEFAULT", str(water_level_default),
        "-LEVEL_REFERENCE", str(level_reference),
        "-CONSTANT_LEVEL", str(int(constant_level)),
        "-WATER_BODY", output_water_body
    ]
    if output_flooded_dem:
        command.extend(["-DEM_FLOODED", output_flooded_dem])
    
    # Execute command
    try:
        print("Running SAGA GIS Terrain Flooding command...")
        result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        
        # Output results
        print("SAGA GIS Command Output:")
        print(result.stdout)
        if result.stderr:
            print("SAGA GIS Command Errors:")
            print(result.stderr)
    except subprocess.CalledProcessError as e:
        print("An error occurred while running the SAGA GIS Terrain Flooding command.")
        print(f"Return Code: {e.returncode}")
        print(f"Error Output: {e.stderr}")


# New block for Quasi-Dynamic Flow Accumulation
def saga_quasi_dynamic_flow(
    saga_cmd_path,
    dem,
    output_flow_accumulation,
    manning_k=None,
    output_flow_through=None,
    output_travel_time=None,
    output_concentration=None,
    output_velocity=None,
    flow_k_default=20.0,
    time=5.0,
    routing=1,
    flow_depth=1,
    flow_const=10.0,
    reset=True,
    rain=10.0
):
    """
    Run the SAGA GIS Quasi-Dynamic Flow Accumulation tool using saga_cmd via Python subprocess.

    Parameters:
        saga_cmd_path (str): Path to the saga_cmd executable.
        dem (str): Path to the input digital elevation model grid.
        output_flow_accumulation (str): Path to save the flow accumulation grid output.
        manning_k (str, optional): Path to the Manning-Strickler coefficient grid.
        output_flow_through (str, optional): Path to save the flow-through output.
        output_travel_time (str, optional): Path to save the flow travel time output.
        output_concentration (str, optional): Path to save the flow concentration output.
        output_velocity (str, optional): Path to save the flow velocity output.
        flow_k_default (float, optional): Default Manning-Strickler coefficient if no grid is selected. Defaults to 20.0.
        time (float, optional): Simulation time in minutes. Defaults to 5.0.
        routing (int, optional): Flow routing method [0: D8, 1: MFD]. Defaults to 1.
        flow_depth (int, optional): Flow depth model [0: constant, 1: dynamic]. Defaults to 1.
        flow_const (float, optional): Constant flow depth in mm. Defaults to 10.0.
        reset (bool, optional): Whether to reset the flow accumulation raster. Defaults to True.
        rain (float, optional): Flow portion added to each raster cell before simulation starts (mm). Defaults to 10.0.
    """
    # Construct the command
    command = [
        saga_cmd_path,
        "sim_hydrology", "8",  # Tool ID for Quasi-Dynamic Flow Accumulation
        "-DEM", dem,
        "-FLOW_ACC", output_flow_accumulation,
        "-FLOW_K_DEFAULT", str(flow_k_default),
        "-TIME", str(time),
        "-ROUTING", str(routing),
        "-FLOW_DEPTH", str(flow_depth),
        "-FLOW_CONST", str(flow_const),
        "-RESET", str(int(reset)),
        "-RAIN", str(rain)
    ]
    if manning_k:
        command.extend(["-FLOW_K", manning_k])
    if output_flow_through:
        command.extend(["-FLOW", output_flow_through])
    if output_travel_time:
        command.extend(["-TIME_MEAN", output_travel_time])
    if output_concentration:
        command.extend(["-TIME_CONC", output_concentration])
    if output_velocity:
        command.extend(["-VELOCITY", output_velocity])
    
    # Execute command
    try:
        print("Running SAGA GIS Quasi-Dynamic Flow Accumulation command...")
        result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        
        # Output results
        print("SAGA GIS Command Output:")
        print(result.stdout)
        if result.stderr:
            print("SAGA GIS Command Errors:")
            print(result.stderr)
    except subprocess.CalledProcessError as e:
        print("An error occurred while running the SAGA GIS Quasi-Dynamic Flow Accumulation command.")
        print(f"Return Code: {e.returncode}")
        print(f"Error Output: {e.stderr}")


def run_river_basin(
    saga_cmd_path,
    dtm,
    hg_grid,
    static_use_grid=None,
    output_grad=None,
    output_direc=None,
    output_hggrad=None,
    output_rivspeed=None,
    output_coordinates=None,
    output_basinshare=None,
    output_statwuse=None,
    output_numinflowcells=None,
    w_cons=False,
    w_cons2=0,
    p_cr=0.0035,
    n_cr=1,
    enf_vmax=True,
    v_tresh=4.0
):
    """
    Run SAGA GIS RiverBasin tool via Python subprocess.

    Parameters:
        saga_cmd_path (str): Path to the saga_cmd executable.
        dtm (str): Path to the input digital terrain model (DTM) grid.
        hg_grid (str): Path to the main river channel grid.
        static_use_grid (str, optional): Path to the static water use raster.
        output_grad (str, optional): Output file path for flow gradients.
        output_direc (str, optional): Output file path for flow directions.
        output_hggrad (str, optional): Output file path for main channel gradients.
        output_rivspeed (str, optional): Output file path for river speeds.
        output_coordinates (str, optional): Output file path for cell coordinates.
        output_basinshare (str, optional): Output file path for basin share raster.
        output_statwuse (str, optional): Output file path for static water use per cell.
        output_numinflowcells (str, optional): Output file path for number of inflow cells.
        w_cons (bool, optional): Enable proportional area water extraction. Defaults to False.
        w_cons2 (int, optional): Type of water use extraction [0: river cells, 1: sub-basin cells]. Defaults to 0.
        p_cr (float, optional): Parameter for main channel lag time calculation. Defaults to 0.0035.
        n_cr (int, optional): Number of storage units in river cascade. Defaults to 1.
        enf_vmax (bool, optional): Consider maximum river velocity. Defaults to True.
        v_tresh (float, optional): Maximum river velocity in km/h. Defaults to 4.0.
    """
    command = [
        saga_cmd_path,
        "sim_rivflow", "0",
        "-INPUT", dtm,
        "-INPUT2", hg_grid,
        "-WCons", str(int(w_cons)),
        "-WCons2", str(w_cons2),
        "-pCr", str(p_cr),
        "-nCr", str(n_cr),
        "-EnfVmax", str(int(enf_vmax)),
        "-VTresh", str(v_tresh)
    ]
    
    if static_use_grid:
        command.extend(["-INPUT3", static_use_grid])
    if output_grad:
        command.extend(["-OUTPUT2", output_grad])
    if output_direc:
        command.extend(["-OUTPUT3", output_direc])
    if output_hggrad:
        command.extend(["-OUTPUT4", output_hggrad])
    if output_rivspeed:
        command.extend(["-OUTPUT5", output_rivspeed])
    if output_coordinates:
        command.extend(["-OUTPUT6", output_coordinates])
    if output_basinshare:
        command.extend(["-OUTPUT7", output_basinshare])
    if output_statwuse:
        command.extend(["-OUTPUT8", output_statwuse])
    if output_numinflowcells:
        command.extend(["-OUTPUT9", output_numinflowcells])
    
    try:
        print("Running SAGA GIS RiverBasin command...")
        result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        print("SAGA GIS Command Output:")
        print(result.stdout)
        if result.stderr:
            print("SAGA GIS Command Errors:")
            print(result.stderr)
    except subprocess.CalledProcessError as e:
        print("An error occurred while running the SAGA GIS RiverBasin command.")
        print(f"Return Code: {e.returncode}")
        print(f"Error Output: {e.stderr}")


def generate_river_grid(dtm_path, output_path, source_x, source_y, mouth_x, mouth_y, overwrite=False):
    """
    Runs the SAGA GIS RiverGridGeneration tool.
    
    Parameters:
    - dtm_path (str): Path to the input Digital Elevation Model (DEM)
    - output_path (str): Path for the output river grid raster
    - source_x (int): X-coordinate of the source cell
    - source_y (int): Y-coordinate of the source cell
    - mouth_x (int): X-coordinate of the mouth cell
    - mouth_y (int): Y-coordinate of the mouth cell
    - overwrite (bool): Whether to overwrite existing river grid cells (default: False)
    
    Returns:
    - None (Runs the SAGA GIS command-line tool)
    """
    overwrite_flag = 1 if overwrite else 0
    
    command = [
        "saga_cmd", "sim_rivflow", "3",
        "-INPUT", dtm_path,
        "-OUTPUT", output_path,
        "-SX", str(source_x),
        "-SY", str(source_y),
        "-MX", str(mouth_x),
        "-MY", str(mouth_y),
        "-Owrite", str(overwrite_flag)
    ]
    
    try:
        subprocess.run(command, check=True)
        print("River grid generation completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error running SAGA GIS tool: {e}")


def run_compound_basic_terrain_analysis(
    elevation, shade, slope, aspect, hcurv, vcurv, convergence, sinks, flow,
    wetness, lsfactor, channels, basins, chnl_base, chnl_dist, vall_depth, rsp,
    threshold=5
):
    """
    Runs the Compound Basic Terrain Analysis tool from SAGA GIS.
    
    Parameters:
        elevation (str): Path to input DEM grid.
        shade (str): Path to output analytical hillshading grid.
        slope (str): Path to output slope grid.
        aspect (str): Path to output aspect grid.
        hcurv (str): Path to output plan curvature grid.
        vcurv (str): Path to output profile curvature grid.
        convergence (str): Path to output convergence index grid.
        sinks (str): Path to output closed depressions grid.
        flow (str): Path to output total catchment area grid.
        wetness (str): Path to output topographic wetness index grid.
        lsfactor (str): Path to output LS-factor grid.
        channels (str): Path to output channel network shapefile.
        basins (str): Path to output drainage basins shapefile.
        chnl_base (str): Path to output channel network base level grid.
        chnl_dist (str): Path to output channel network distance grid.
        vall_depth (str): Path to output valley depth grid.
        rsp (str): Path to output relative slope position grid.
        threshold (int, optional): Strahler order to begin a channel. Default is 5.
    """
    command = [
        "saga_cmd", "ta_compound", "0",
        "-ELEVATION", elevation,
        "-SHADE", shade,
        "-SLOPE", slope,
        "-ASPECT", aspect,
        "-HCURV", hcurv,
        "-VCURV", vcurv,
        "-CONVERGENCE", convergence,
        "-SINKS", sinks,
        "-FLOW", flow,
        "-WETNESS", wetness,
        "-LSFACTOR", lsfactor,
        "-CHANNELS", channels,
        "-BASINS", basins,
        "-CHNL_BASE", chnl_base,
        "-CHNL_DIST", chnl_dist,
        "-VALL_DEPTH", vall_depth,
        "-RSP", rsp,
        "-THRESHOLD", str(threshold)
    ]
    
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode == 0:
        print("Compound Basic Terrain Analysis completed successfully.")
    else:
        print("Error:", result.stderr)