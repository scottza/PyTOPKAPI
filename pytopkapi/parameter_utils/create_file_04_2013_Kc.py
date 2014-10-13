"""
*** OBJECTIVE

1. creating the TOPKAPI parameter file from GIS binary files (binary
grid format).  -main routine 'creat_param_file'

*** COMMENT

In the present program, to write the parameter file,
generally the grid are tranformed in an array by rearranging the cell
order from the West to East and North to South, as in the following
example:

The grid format is like:

GIS bingrid file
-9999  -9999  -9999  -9999  -9999  -9999  -9999
-9999      0      1  -9999  -9999  -9999  -9999
-9999      2      3     4   -9999  -9999  -9999
    5      6      7     8       9  -9999  -9999
-9999     10     11    12      13     14  -9999
-9999  -9999     15    16      17  -9999  -9999
-9999  -9999  -9999  -9999     18  -9999  -9999
-9999  -9999  -9999 -9999   -9999  -9999  -9999

The corresponding array extracted and ordered from West to East, North
to South is: 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18

"""
import sys
from warnings import warn
from ConfigParser import SafeConfigParser, NoOptionError

import numpy as np
from numpy import ma
import networkx as nx
#from osgeo import gdal # that was on python win 32 bit RUBBISH
from osgeo import gdal

from pytopkapi import utils as ut
from pytopkapi import pretreatment as pm

def compute_cell_coordinates(mask_fname):
    dset = gdal.Open(mask_fname)
    mask = dset.ReadAsArray()

    x0, dx, fy, y0, fx, dy = dset.GetGeoTransform()

    # Adjust x0 and y0 to give pixel centres. GDAL specifies (x0, y0)
    # as the top left corner of the top left pixel. PyTOPKAPI expects
    # pixel centres. At the centre of the first pixel the following
    # holds: (Xpixel, Yline) == (0.5, 0.5)
    x0 = x0 + dx/2.0 + fy/2.0
    y0 = y0 + fx/2.0 + dy/2.0

    Yline, Xpixel = np.nonzero(mask == 1)

    Xgeo = x0 + Xpixel*dx + Yline*fy
    Ygeo = y0 + Xpixel*fx + Yline*dy

    return Xgeo, Ygeo

def read_raster(rast_fname, file_format='GTiff'):
    """Read the data in a raster file

    Parameters
    ----------
    rast_fname : string
        The path to the raster file.
    file_format : string
        The file format of the raster file, currently this can only be
        GeoTIFF ('GTiff' - default).

    Returns
    -------
    data : Numpy ndarray
        The 2D raster data in a Numpy array. The dtype of the returned
        array is the same as the data type stored in the raster file.

    """
    if file_format != 'GTiff':
        err_str = 'Reading %s files not implemented.' % file_format
        raise NotImplementedError(err_str)
    else:
        dset = gdal.Open(rast_fname)
        data = dset.ReadAsArray()

    return data

def generate_param_file(ini_fname, isolated_cells=False):
    """Create a PyTOPKAPI parameter file

    Generates a parameter file from the catchment data provided in a
    set of georeferenced raster files. The input files are currently
    limited to GeoTIFF. Support for 32-bit raster files with ArcGIS
    style headers is planned.

    Parameters
    ----------
    ini_fname : string

        The path to an ini-style config file specifying the input data
        file locations and format. The ini file should contain the
        following sections and parameters:

        [raster_files]
        dem_fname = <path to DEM file>
        mask_fname = <path to catchment mask file>
        soil_depth_fname = <path to soil depth file>
        conductivity_fname = <path to saturated conductivity file>
        hillslope_fname = <path to hill slope file>
        sat_moisture_content_fname = <path to saturated moisture content file>
        resid_moisture_content_fname = <path to residual moisture content file>
        bubbling_pressure_fname = <path to bubbling pressure file>
        pore_size_dist_fname = <path to pore size index file>
        overland_manning_fname = <path to overland Manning roughness file>
        channel_network_fname = <path to channel network file>
        flowdir_fname = <path to flow direction file>
        flowdir_source = <source of flowdir file. Can be `GRASS` or `ARCGIS`>

        [output]
        param_fname = <path to output parameter file>

        [numerical_values]
        pVs_t0 = <initial percent saturation of soil stores>
        Vo_t0 = <initial volume of overland stores>
        Qc_t0 = <initial flow rate in channels>
        Kc = <crop factor>

    isolated_cells : bool
        A flag indicating whether the generated parameter file should
        consist of a set of isolated and unconnected cells. The
        default value of `False` generates a parameter file for
        modelling the catchment as a network of inter-connected cells.

    Returns
    -------
    Nothing

    Notes
    -----

    Note that the initial soil store and overland volumes, as well as
    the initial channel flow are constant for all cells. Spatially
    varying values can be assigned by directly manipulating the
    paramter file, or by using the routines in
    `pytopkapi.parameter_utils.modify_file`.

    """
    config             = SafeConfigParser()
    config.read(ini_fname)

    dem_fname          = config.get('raster_files', 'dem_fname')
    mask_fname         = config.get('raster_files', 'mask_fname')
    soil_depth_fname   = config.get('raster_files', 'soil_depth_fname')
    conductivity_fname = config.get('raster_files', 'conductivity_fname')
    hillslope_fname    = config.get('raster_files', 'hillslope_fname')
    theta_sat_fname    = config.get('raster_files', 'sat_moisture_content_fname')
    theta_r_fname      = config.get('raster_files', 'resid_moisture_content_fname')
    psi_b_fname        = config.get('raster_files', 'bubbling_pressure_fname')
    lamda_fname        = config.get('raster_files', 'pore_size_dist_fname')
    n_o_fname          = config.get('raster_files', 'overland_manning_fname')
    network_fname      = config.get('raster_files', 'channel_network_fname')
    flowdir_fname      = config.get('raster_files', 'flowdir_fname')
    fdir_source        = config.get('raster_files', 'flowdir_source')
    KC                 = 1
#    FEDDES             = 1
    WWB                = 1
    try:
        kc_ras_fname   = config.get('raster_files', 'kc_ras_fname')
    except NoOptionError:
        KC             = 0    
#    try:
#        brks_fname     = config.get('raster_files', 'brooks_lamda_fname')
#    except NoOptionError:
#        FEDDES         = 0
    try:
        wwb_ras_fname  = config.get('raster_files', 'wwb_fname')
    except NoOptionError:
        WWB            = 0

    pVs_t0 = config.get('numerical_values', 'pVs_t0')
    Vo_t0  = config.get('numerical_values', 'Vo_t0')
    Qc_t0  = config.get('numerical_values', 'Qc_t0')
    Kc     = config.getfloat('numerical_values', 'Kc')


    param_fname = config.get('output', 'param_fname')
    
    # Check if any of the raster files is missing:
    rasters_fn = [dem_fname,mask_fname,soil_depth_fname,conductivity_fname,\
                  hillslope_fname,theta_sat_fname,theta_r_fname,psi_b_fname,\
                  lamda_fname,n_o_fname,network_fname,flowdir_fname]
    if KC:
        rasters_fn.append(kc_ras_fname)
    if WWB:
        rasters_fn.append(wwb_ras_fname)
#    if FEDDES:
#        rasters_fn.append(brks_fname)
    
    
    # Check if all required raster files that were listed in the ini file exist
    for fn in rasters_fn:
        try:
           with open(fn): pass
        except IOError:
           raise ValueError('Raster file \'%s\' does not exist.'%(fn.split('/')[-1]))

    # Read the input rasters
    dem             = read_raster(dem_fname)
    mask            = read_raster(mask_fname)
    hillslope       = read_raster(hillslope_fname)
    depth           = read_raster(soil_depth_fname)
    theta_sat       = read_raster(theta_sat_fname)
    theta_r         = read_raster(theta_r_fname)
    conductivity    = read_raster(conductivity_fname)
    psi_b           = read_raster(psi_b_fname)
    lamda           = read_raster(lamda_fname)
    n_o             = read_raster(n_o_fname)
    channel_network = read_raster(network_fname)
    flowdir         = read_raster(flowdir_fname)
    if KC:
        Kc_ras      = read_raster(kc_ras_fname)
        print 'KC is on..'
    if WWB:
        print 'WWB are on..'
        wwb_ras     = read_raster(wwb_ras_fname)
#    if FEDDES:
#        brooks_l    = read_raster(brks_fname)

    # Calculate parameters
    ncells      = mask[mask == 1].size
    nparams     = 22
    cell_labels = np.arange(ncells)
    tan_beta    = np.tan((np.pi/180.0)*hillslope)
    X, Y        = compute_cell_coordinates(mask_fname)

    channel_network[channel_network < 255]  = 1
    channel_network[channel_network == 255] = 0

    if isolated_cells == True:
        cell_down = -999
        channel_length = 0
        n_c = 0
        tan_beta_channel = 0
    else:
        # Calculate the network connections and channel lengths.
        cell_down = cell_connectivity(flowdir, mask, fdir_source)

        channel_length, tan_beta_channel = channel_properties(cell_labels,
                                                     channel_network[mask == 1],
                                                     X, Y, cell_down,
                                                     dem[mask == 1])

        n_c = strahler_to_channel_manning(cell_labels,
                                          channel_network[mask == 1],
                                          cell_down)

    # Write parameter file
    param_table       = np.zeros((ncells, nparams))
    param_table[:,0]  = cell_labels
    param_table[:,1]  = X
    param_table[:,2]  = Y
    param_table[:,3]  = channel_network[mask == 1]
    param_table[:,4]  = channel_length
    if WWB:
        param_table[:,5]  = wwb_ras[mask == 1]
    param_table[:,6]  = tan_beta[mask == 1]
    param_table[:,7]  = tan_beta_channel
    param_table[:,8]  = depth[mask == 1]
    param_table[:,9]  = conductivity[mask == 1]
    param_table[:,10] = theta_r[mask == 1]
    param_table[:,11] = theta_sat[mask == 1]
    param_table[:,12] = n_o[mask == 1]
    param_table[:,13] = n_c
    param_table[:,14] = cell_down
    param_table[:,15] = pVs_t0
    param_table[:,16] = Vo_t0
    param_table[:,17] = Qc_t0
    if KC:
        param_table[:,18] = Kc_ras[mask == 1] * Kc
    else:
        param_table[:,18] = Kc
    param_table[:,19] = psi_b[mask == 1]
    param_table[:,20] = lamda[mask == 1]
    
    # Check Kc-raster file for any invalid values
    if KC:
        kc_tmp = Kc_ras[mask == 1] * Kc
        blanks = kc_tmp[kc_tmp < 0.].shape[0]
        if blanks:
            err_msg = 'The crop factor raster file has %i blank cells of data value(s) within the mask' \
            %(blanks)
            raise ValueError(err_msg)
#    if FEDDES:
#        param_table[:,21] = brooks_l[mask == 1]

    format = '%d %f %f %d %f %d %f %f %f %0.10f %f %f %f %f %d %f %f %f %f %f %f %f'
    np.savetxt(param_fname, param_table, fmt=format)

def generate_param_file_DL(ini_fname, isolated_cells=False):
    """Create a PyTOPKAPI parameter file

    Generates a parameter file from the catchment data provided in a
    set of georeferenced raster files. The input files are currently
    limited to GeoTIFF. Support for 32-bit raster files with ArcGIS
    style headers is planned.

    Parameters
    ----------
    ini_fname : string

        The path to an ini-style config file specifying the input data
        file locations and format. The ini file should contain the
        following sections and parameters:

        [raster_files]
        dem_fname = <path to DEM file>
        mask_fname = <path to catchment mask file>
        soil_depth_fname = <path to soil depth file>
        conductivity_fname = <path to saturated conductivity file>
        hillslope_fname = <path to hill slope file>
        sat_moisture_content_fname = <path to saturated moisture content file>
        resid_moisture_content_fname = <path to residual moisture content file>
        bubbling_pressure_fname = <path to bubbling pressure file>
        pore_size_dist_fname = <path to pore size index file>
        overland_manning_fname = <path to overland Manning roughness file>
        channel_network_fname = <path to channel network file>
        flowdir_fname = <path to flow direction file>
        flowdir_source = <source of flowdir file. Can be `GRASS` or `ARCGIS`>

        [output]
        param_fname = <path to output parameter file>

        [numerical_values]
        pVs_t0 = <initial percent saturation of soil stores>
        Vo_t0 = <initial volume of overland stores>
        Qc_t0 = <initial flow rate in channels>
        Kc = <crop factor>

    isolated_cells : bool
        A flag indicating whether the generated parameter file should
        consist of a set of isolated and unconnected cells. The
        default value of `False` generates a parameter file for
        modelling the catchment as a network of inter-connected cells.

    Returns
    -------
    Nothing

    Notes
    -----

    Note that the initial soil store and overland volumes, as well as
    the initial channel flow are constant for all cells. Spatially
    varying values can be assigned by directly manipulating the
    paramter file, or by using the routines in
    `pytopkapi.parameter_utils.modify_file`.

    """
    config             = SafeConfigParser()
    config.read(ini_fname)

    dem_fname           = config.get('raster_files', 'dem_fname')
    mask_fname          = config.get('raster_files', 'mask_fname')
    soil_depth_fname    = config.get('raster_files', 'soil_depth_fname')
    conductivity_fname  = config.get('raster_files', 'conductivity_fname')
    hillslope_fname     = config.get('raster_files', 'hillslope_fname')
    theta_sat_fname     = config.get('raster_files', 'sat_moisture_content_fname')
    theta_r_fname       = config.get('raster_files', 'resid_moisture_content_fname')
    psi_b_fname         = config.get('raster_files', 'bubbling_pressure_fname')
    lamda_fname         = config.get('raster_files', 'pore_size_dist_fname')
    n_o_fname           = config.get('raster_files', 'overland_manning_fname')
    network_fname       = config.get('raster_files', 'channel_network_fname')
    flowdir_fname       = config.get('raster_files', 'flowdir_fname')
    fdir_source         = config.get('raster_files', 'flowdir_source')
    # raster files required for the second soil layer:
    soil_depth_fname_   = config.get('raster_files', 'soil_depth_fname_II')
    conductivity_fname_ = config.get('raster_files', 'conductivity_fname_II')
    theta_sat_fname_    = config.get('raster_files', 'sat_moisture_content_fname_II')
    theta_r_fname_      = config.get('raster_files', 'resid_moisture_content_fname_II')
    psi_b_fname_        = config.get('raster_files', 'bubbling_pressure_fname_II')
    lamda_fname_        = config.get('raster_files', 'pore_size_dist_fname_II')

    KC                 = 1
#    FEDDES             = 1
    WWB                = 1
    try:
        kc_ras_fname   = config.get('raster_files', 'kc_ras_fname')
    except NoOptionError:
        KC             = 0    
#    try:
#        brks_fname     = config.get('raster_files', 'brooks_lamda_fname')
#    except NoOptionError:
#        FEDDES         = 0
    try:
        wwb_ras_fname  = config.get('raster_files', 'wwb_fname')
    except NoOptionError:
        WWB            = 0

    pVs_t0  = config.get('numerical_values', 'pVs_t0')
    pVs_t0_ = config.get('numerical_values', 'pVs_t0_II')
    Vo_t0   = config.get('numerical_values', 'Vo_t0')
    Qc_t0   = config.get('numerical_values', 'Qc_t0')
    Kc      = config.getfloat('numerical_values', 'Kc')


    param_fname = config.get('output', 'param_fname')
    
    # Check if any of the raster files is missing:
    rasters_fn = [dem_fname,mask_fname,soil_depth_fname,conductivity_fname,\
                  hillslope_fname,theta_sat_fname,theta_r_fname,psi_b_fname,\
                  lamda_fname,n_o_fname,network_fname,flowdir_fname]
    if KC:
        rasters_fn.append(kc_ras_fname)
    if WWB:
        rasters_fn.append(wwb_ras_fname)
#    if FEDDES:
#        rasters_fn.append(brks_fname)

    # Add the raster names for the second soil layer
    for fn in [soil_depth_fname_,conductivity_fname_,theta_sat_fname_,
               theta_r_fname_,psi_b_fname_,lamda_fname_]:
        rasters_fn.append(fn)

    # Check if all required raster files that were listed in the ini file exist
    for fn in rasters_fn:
        try:
           with open(fn): pass
        except IOError:
           raise ValueError('Raster file \'%s\' does not exist.'%(fn.split('/')[-1]))

    # Read the input rasters
    dem             = read_raster(dem_fname)
    mask            = read_raster(mask_fname)
    hillslope       = read_raster(hillslope_fname)
    depth           = read_raster(soil_depth_fname)
    depth_          = read_raster(soil_depth_fname_)
    theta_sat       = read_raster(theta_sat_fname)
    theta_sat_      = read_raster(theta_sat_fname)
    theta_r         = read_raster(theta_r_fname)
    theta_r_        = read_raster(theta_r_fname)
    conductivity    = read_raster(conductivity_fname)
    conductivity_   = read_raster(conductivity_fname)
    psi_b           = read_raster(psi_b_fname)
    psi_b_          = read_raster(psi_b_fname)
    lamda           = read_raster(lamda_fname)
    lamda_          = read_raster(lamda_fname)
    n_o             = read_raster(n_o_fname)
    channel_network = read_raster(network_fname)
    flowdir         = read_raster(flowdir_fname)
    if KC:
        Kc_ras      = read_raster(kc_ras_fname)
        print 'KC is on..'
    if WWB:
        print 'WWB are on..'
        wwb_ras     = read_raster(wwb_ras_fname)

    # Calculate parameters
    ncells      = mask[mask == 1].size
    nparams     = 28
    cell_labels = np.arange(ncells)
    tan_beta    = np.tan((np.pi/180.0)*hillslope)
    X, Y        = compute_cell_coordinates(mask_fname)

    channel_network[channel_network < 255]  = 1
    channel_network[channel_network == 255] = 0

    if isolated_cells == True:
        cell_down = -999
        channel_length = 0
        n_c = 0
        tan_beta_channel = 0
    else:
        # Calculate the network connections and channel lengths.
        cell_down = cell_connectivity(flowdir, mask, fdir_source)

        channel_length, tan_beta_channel = channel_properties(cell_labels,
                                                     channel_network[mask == 1],
                                                     X, Y, cell_down,
                                                     dem[mask == 1])

        n_c = strahler_to_channel_manning(cell_labels,
                                          channel_network[mask == 1],
                                          cell_down)

    # Write parameter file
    param_table       = np.zeros((ncells, nparams))
    param_table[:,0]  = cell_labels
    param_table[:,1]  = X
    param_table[:,2]  = Y
    param_table[:,3]  = channel_network[mask == 1]
    param_table[:,4]  = channel_length
    if WWB:
        param_table[:,5]  = wwb_ras[mask == 1]
    param_table[:,6]  = tan_beta[mask == 1]
    param_table[:,7]  = tan_beta_channel
    param_table[:,8]  = depth[mask == 1]
    param_table[:,9]  = conductivity[mask == 1]
    param_table[:,10] = theta_r[mask == 1]
    param_table[:,11] = theta_sat[mask == 1]
    param_table[:,12] = n_o[mask == 1]
    param_table[:,13] = n_c
    param_table[:,14] = cell_down
    param_table[:,15] = pVs_t0
    param_table[:,16] = Vo_t0
    param_table[:,17] = Qc_t0
    if KC:
        param_table[:,18] = Kc_ras[mask == 1] * Kc
    else:
        param_table[:,18] = Kc
    param_table[:,19] = psi_b[mask == 1]
    param_table[:,20] = lamda[mask == 1]

    # Now add parameters from second soil layer
    param_table[:,21] = depth_[mask == 1]
    param_table[:,22] = conductivity_[mask == 1]
    param_table[:,23] = theta_r_[mask == 1]
    param_table[:,24] = theta_sat_[mask == 1]
    param_table[:,25] = pVs_t0_
    param_table[:,26] = psi_b_[mask == 1]
    param_table[:,27] = lamda_[mask == 1]
    
    # Check Kc-raster file for any invalid values
    if KC:
        kc_tmp = Kc_ras[mask == 1] * Kc
        blanks = kc_tmp[kc_tmp < 0.].shape[0]
        if blanks:
            err_msg = 'The crop factor raster file has %i blank cells of data value(s) within the mask' \
            %(blanks)
            raise ValueError(err_msg)

    format = '%d %f %f %d %f %d %f %f %f %0.10f %f %f %f %f %d %f %f %f %f %f %f %f %f %f %f %f %f %f'
    np.savetxt(param_fname, param_table, fmt=format)

##############################################
###  SUBROUTINE USED IN "creat_param_file" ###
##############################################

####### READING THE BINARY FILES ###########

def read_bin_data(fname):
    """Read data from a float32 binary file into an array.

    """

    f = open(fname, "rb")
    raw = f.read()
    f.close()
    data = np.fromstring(raw, 'f')
    if sys.byteorder == 'big':
        data = data.byteswap()

    return data

def read_arc_bin(bingrid_name):
    """Read data from an ArcGIS float32 binary file into an array.

    """

    li_headers=read_headers_arc_bin(bingrid_name)

    rows = li_headers[1] # fixed for UM grid
    cols = li_headers[0] # fixed for UM grid

    bin_name = bingrid_name + '.flt'

    a = read_bin_data(bin_name)

    a = a.reshape(rows, cols)

    return a

def read_headers_arc_bin(bingrid_name):
    """Read the ascii headers of the binary grid file:
    ncols         62
    nrows         121
    xllcorner     -288595.47161281
    yllcorner     -3158065.5722693
    cellsize      1000
    NODATA_value  -9999
    byteorder     LSBFIRST
    """

    hdr_name = bingrid_name + '.hdr'
    f=open(hdr_name,'r')
    tab_read=f.readlines()
    f.close()

    li_headers=[]
    i=-1
    for line in tab_read:
        i=i+1
        donnees=line.split()
        if i<6:
            li_headers.append(float(donnees[1]))
        else:
            li_headers.append(donnees[1])

    return li_headers

def _make_strahler_dicts(G):
    """Prepare dictionaries for the Strahler algorithm"""
    nodes_per_arc = {}
    arcs_per_node = {}

    for edge_id, edge in enumerate(G.edges()):
        nodes_per_arc[edge_id] = edge

        for node in G.nodes():
            if node in edge:
                if node in arcs_per_node:
                    arcs_per_node[node].append(edge_id)
                else:
                    arcs_per_node[node] = [edge_id]

    return nodes_per_arc, arcs_per_node

def strahler_stream_order(start_arc_id, start_up_node,
                          nodes_per_arc, arcs_per_node, stream_orders):
    """Calculate the Strahler stream order

    This function recursively computes the Strahler stream order using
    the algorithm described by Gleyzer et al. (2004). The sequence of
    stream orders for the starting arc and each upstream arc is
    returned in the dictionary `stream_orders`. To compute the
    Strahler order for the entire network, `start_arc_id` should be
    the arc ID for the stream arc closest to the catchment outlet and
    `start_up_node` should be the node ID at the upstream end of
    `start_arc_id`.

    Parameters
    ----------
    start_arc_id : int
        The integer ID of the current stream arc as defined in
        `nodes_per_arc`.
    start_up_node : int
        The integer ID of the upstream node for the current stream
        arc.
    nodes_per_arc : dict
        A dictionary containing an ordered tuple representing the
        upstream and downstream node IDs for each stream arc
        ID. e.g. {0 : (upstream_node, downstream_node)}
    arcs_per_node : dict
        A dictionary containing a list of the stream arc IDs for
        stream arcs adjacent to each node in the network.
    stream_orders : dict
        A dictionary with the (key, value) pairs representing the
        stream arc ID and associated Strahler order.

    Returns
    -------
    order : int
        The stream order of the stream arc described by
        `start_arc_id`.

    References
    ----------
    Alexander Gleyzer, Michael Denisyuk, Alon Rimmer and Yigal
    Salingar, 2004. A Fast Recursive GIS Algorithm for Computing
    Strahler Stream Order in Braided and Nonbraided Networks. Journal
    of the American Water Resources Association (JAWRA) 40(4):937-946.

    """
    if len(arcs_per_node[start_up_node]) == 1:
        stream_orders[start_arc_id] = 1
    else:
        upstream_orders = {}

        for arc_id in arcs_per_node[start_up_node]:
            if arc_id != start_arc_id:
                up_node, down_node = nodes_per_arc[arc_id]
                if up_node != start_up_node:
                    upstream_orders[arc_id] = strahler_stream_order(arc_id,
                                                                    up_node,
                                                                  nodes_per_arc,
                                                                  arcs_per_node,
                                                                  stream_orders)
                else:
                    upstream_orders[arc_id] = strahler_stream_order(arc_id,
                                                                    down_node,
                                                                  nodes_per_arc,
                                                                  arcs_per_node,
                                                                  stream_orders)

        max_order = 0
        max_order_count = 0
        up_orders = upstream_orders.values()
        up_orders.sort(reverse=True)

        for order in up_orders:
            if order > max_order:
                max_order = order
                max_order_count += 1
            elif order == max_order:
                max_order_count += 1

        if max_order_count > 1:
            stream_orders[start_arc_id] = max_order + 1
        else:
            stream_orders[start_arc_id] = max_order

    return stream_orders[start_arc_id]

def strahler_to_channel_manning(cell_labels, channel_network, cell_down):
    """Calculate the Manning roughness for channel cells

    Computes the Strahler order for the channel in each channel cell
    and assigns a Manning roughness to each using a table of the
    correspondance between the Strahler order and the values of
    Manning roughness, as proposed in Liu and Todini (2002).

    Parameters
    ----------
    cell_labels : 1D Numpy ndarray
        An array of the labels associated with each cell in the
        catchment
    channel_network : 1D Numpy ndarray
        An ordered array with each channel cell indicated by a value
        of one, zero otherwise.
    cell_down : 1D Numpy ndarray
        An ordered array giving the label of the downstream cell in
        the catchment network. The outlet of the catchment is
        indicated by a negative number.

    Returns
    -------
    n_c : 1D array
        An array containing the values of the manning coefficient for
        each channel cell and zero for non-channel cells.

    """
    strahler_manning = {1 : 0.050,
                        2 : 0.040,
                        3 : 0.035,
                        4 : 0.030,
                        5 : 0.030,
                        6 : 0.025,
                        7 : 0.025}

    # ensure input arrays are integers
    cell_labels = np.asarray(cell_labels, dtype=np.int)
    channel_network = np.asarray(channel_network, dtype=np.int)
    cell_down = np.asarray(cell_down, dtype=np.int)

    # compute strahler order
    nodes = cell_labels[channel_network == 1]

    edges = []
    for k in cell_labels:
        if (channel_network[k] == 1) and (cell_down[k] >= 0):
            edges.append((k, cell_down[k]))

    G = nx.DiGraph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)

    # determine the outlet stream arc ID and it's upstream node using
    # obfuscated list comprehension. Strictly for speed of course ;-)
    outlet_node = nx.topological_sort(G)[-1]

    outlet_info = [(edge_id, edge) for edge_id, edge in enumerate(G.edges())
           if outlet_node in edge]

    outlet_edge_id = outlet_info[0][0]
    outlet_up_node = outlet_info[0][1][0]

    stream_orders = {}

    nodes_per_arc, arcs_per_node = _make_strahler_dicts(G)
    strahler_stream_order(outlet_edge_id, outlet_up_node,
                          nodes_per_arc, arcs_per_node, stream_orders)

    # assign strahler value of stream arcs to cells
    strahler_per_node = {}

    for key in nodes_per_arc.keys():
        edge = nodes_per_arc[key]

        strahler_per_node[edge[0]] = stream_orders[key]

        if outlet_node in edge:
            strahler_per_node[outlet_node] = stream_orders[key]

    # assign manning based on strahler order
    n_c = np.zeros(cell_labels.size)
    for node in nx.topological_sort(G):
        key = strahler_per_node[node]
        n_c[node] = strahler_manning[key]

    return n_c

def cell_connectivity(flowdir, mask, source='GRASS'):
    """Compute the connectivity between cells in the catchment

    Associate each cell in the catchment with the label of it's
    downstream neighbour. This defines the directed network of
    connections between the cells which make up a catchment in the
    model.

    Parameters
    ----------
    flowdir : Numpy ndarray
        8D flow direction codes as defined by your favourite GIS
        toolbox.
    mask : Numpy ndarray
        A 2D array with the same shape as `flowdir`. Cells which
        comprise the catchment should be marked by a value of 1.
    source : string
        A string describing the source of the flow direction
        codes. Current options are 'ArcGIS' or 'GRASS' (default)

    Returns
    -------
    cell_down : Numpy ndarray
        An ordered array containing the label of the immediate
        downstream cell for each cell in the catchment network. A 1D
        array with length equal to the number of cells. The catchment
        outlet is indiacted by a value of -999.

    """
    sentinel_value = -999
    # Specify flow direction code from GRASS GIS r.watershed or ArcGIS
    # Hydrology toolbox flow-direction tool
    if source == 'GRASS':
        ddict = {1 : (-1,  1),
                 2 : (-1,  0),
                 3 : (-1, -1),
                 4 : ( 0, -1),
                 5 : ( 1, -1),
                 6 : ( 1,  0),
                 7 : ( 1,  1),
                 8 : ( 0,  1)}
    elif source == 'ArcGIS':
        ddict = {128 : (-1,  1),
                 64  : (-1,  0),
                 32  : (-1, -1),
                 16  : ( 0, -1),
                 8   : ( 1, -1),
                 4   : ( 1,  0),
                 2   : ( 1,  1),
                 1   : ( 0,  1)}
    else:
        raise ValueError('Unknown flow direction scheme: %s' % source)

    ncells = mask[mask == 1].size
#    print 'ncells', ncells
    int_min = np.iinfo(np.int).min
#    print 'int_min', int_min
    cell_id = np.ones(mask.shape, dtype=np.int)*int_min
    cell_id[mask == 1] = np.arange(ncells)

    cell_down = np.ones(ncells, dtype=np.int)*int_min

    outlet_found = False
    nrows, ncols = mask.shape
#    print 'mask.shape', mask.shape
    for i in range(nrows):
        for j in range(ncols):
            fdir = flowdir[i, j]

            if fdir in ddict.keys():
                r, c = ddict[fdir]
                m, n = i+r, j+c

                # If the raster is closely cropped to the catchment
                # boundary, (m, n) computed for the catchment outlet
                # will fall outside the array bounds.

                if (m < 0) or (m > nrows-1) or \
                   (n < 0) or (n > ncols-1):
                    if outlet_found:
                        error_txt = """More than one catchment outlet detected.
Check that the flow direction raster and catchment mask are compatible.

First outlet detected at:
cell_id = %s
row = %s
col = %s

Second outlet detected at:
cell_id = %s
row = %s
col = %s""" % (outlet_loc[0], outlet_loc[1], outlet_loc[2], cell_id[i, j], i, j)
                        raise ValueError(error_txt)
                    else:
                        outlet_found = True
                        outlet_loc = (cell_id[i, j], i, j)

                        cell_down[cell_id[i, j]] = sentinel_value
                else:
#                    print i, j
#                    print m, n
#                    print cell_id[i,j]
#                    print cell_down[cell_id[i,j]]
                    cell_down[cell_id[i, j]] = cell_id[m, n]

    if cell_down[cell_down == int_min].size > 1:
        warn_txt = """There are %d catchment cells without a downstream link.
Check the validity of the flow-direction raster."""  \
        % cell_down[cell_down == int_min].size

        warn(warn_txt)

    # Handle case where outlet isn't at the raster boundary.
    cell_down[cell_down == int_min] = sentinel_value
    if cell_down[cell_down == sentinel_value].size > 1:
        error_txt = """More than one catchment outlet detected.
Check that the flow direction raster and catchment mask are compatible.
"""
        raise ValueError(error_txt)

    return cell_down

def channel_properties(cell_labels, channel_network, X, Y, cell_down, dem):
    """Compute the length and slope of the channels

    Cells draining diagonally have a different channel length from
    cells draining North, South, East or West. This function computes
    the channel length as a function of the drainage direction (based
    on the catchment's cell connectivity). The slope is calculated as
    the height difference over the length, in a downstream direction
    (i.e. negative slopes indicate an inconsistency in the input DEM).

    Parameters
    ----------
    cell_labels : 1D Numpy ndarray
        An array of the labels associated with each cell in the
        catchment
    channel_network : 1D Numpy ndarray
        An ordered array with each channel cell indicated by a value
        of one, zero otherwise.
    X : 1D Numpy ndarray
        An ordered array of the X coordinate of the centre of each
        cell.
    Y : 1D Numpy ndarray
        An ordered array of the Y coordinate of the centre of each
        cell.
    cell_down : 1D Numpy ndarray
        An ordered array giving the label of the downstream cell in
        the catchment network. The outlet of the catchment is
        indicated by a negative number.
    dem : 1D Numpy ndarray
        An ordered array of cell elevations.

    Returns
    -------
    Xc : 1D Numpy ndarray
        An array containing the length of the channel in each channel
        cell, zero otherwise.
    tan_beta_channel : 1D Numpy ndarray
        An array containing the slope of the channel in each channel
        cell, zero otherwise.

    """
    # Ensure input arrays are numpy arrays of the correct dtype.
    cell_labels = np.asarray(cell_labels, dtype=np.int)
    channel_network = np.asarray(channel_network, dtype=np.int)
    X = np.asarray(X, dtype=np.float)
    Y = np.asarray(Y, dtype=np.float)
    cell_down = np.asarray(cell_down, dtype=np.int)
    dem = np.asarray(dem, dtype=np.float)

    Xc = np.zeros(cell_labels.shape, dtype=np.float)
    tan_beta_channel = np.zeros(cell_labels.shape, dtype=np.float)

    for i in cell_labels[channel_network == 1]:
        indx = cell_down[i]
        if indx >= 0:
            # Channel cells upstream of the catchment outlet.
            Xcell = X[i]
            Ycell = Y[i]

            Xcell_down = X[cell_labels == indx]
            Ycell_down = Y[cell_labels == indx]

            Xc[i] = ut.distance(Xcell, Ycell, Xcell_down, Ycell_down)
            tan_beta_channel[i] = (dem[i]
                                   - dem[cell_labels == indx][0])/Xc[i]

    # Assign sensible values to the catchment outlet cell, since it
    # has no downstream neighbour.
    outlet_indx = np.nonzero(cell_down < 0)
    cond = (channel_network == 1) & (cell_down == outlet_indx[0][0])
    upstream_indx = cell_labels[cond][0]

    Xc[outlet_indx] = Xc[upstream_indx]
    tan_beta_channel[outlet_indx] = tan_beta_channel[upstream_indx]

    return Xc, tan_beta_channel


def texture_lookup(soil_codes, soil_lookup_table):
    """
    * Objective:
      Extraction of the parameters L (soil depth) and theta_s (porosity or 
      humidity at saturation) for each catchment cell from the SIRI map
    * Input
      - file_bin_WRC90 is the binary grid file containing the WRC90 soil 
        property codes (Here only three 3 for Loamy Sand, 2 for Sandy Loam, 
        1 for Clay)
      - file_table_WRC90_soil is an ASCII file containing a table of 
        correspondance between the WRC90 codes and the values of Ks 
        (permeability) and theta_r (residual soil moisture)
    * Ouput
      This routine returns two 1D array (ar_theta_r, ar_theta_s) containing 
      respectively the values of Ks and theta_r for each cell. Cells are 
      ordered from West to East, North to South.
    """
    #Read the binary grid file of GLCC land use type
    tab=read_arc_bin(soil_codes)
    nrows=np.shape(tab)[0]
    ncols=np.shape(tab)[1]
    tab=np.reshape(tab,ncols*nrows)
    ind=np.where(tab>-99)
    ar_WRC90=tab[ind]

    #Read the Table file within a header line
    tab=pm.read_column_input(file_table_WRC90_soil,9)
    ar_code=tab[:,0]
    ar_theta_r_moy=tab[:,3]
    ar_theta_r_ect=tab[:,4]
    ar_conduct=tab[:,8]

    ar_theta_r=np.array(ar_WRC90)
    ar_Ks=np.array(ar_WRC90)
    for i in ar_code:
        ind=np.where(ar_WRC90==i)
        #!!!! TO BE CHANGED FOR PARAMETER ADJUSTMENT !!!#
        ar_theta_r[ind]=ar_theta_r_moy[np.where(ar_code==i)][0]
        ar_Ks[ind]=ar_conduct[np.where(ar_code==i)][0]
         
    return ar_theta_r, ar_Ks