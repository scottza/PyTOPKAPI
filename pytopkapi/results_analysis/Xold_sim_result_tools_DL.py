__all__ = ['extract_ssi', 'extract_ssi_to_file',
           'extract_Q_channel', 'extract_Q_down',
           'extract_overland_volume', 'extract_overland_volume_to_file',
           'run_version_log', 'Extract_flow_at_outlet']

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':'10'})
rc('text', usetex=True)

from datetime import timedelta
from ConfigParser import SafeConfigParser, NoOptionError, NoSectionError
import ConfigParser
import xlsxwriter as xxw
import matplotlib.pyplot as plt
import os

import h5py
import numpy as np
import numpy.ma as ma

import pytopkapi.utils as ut
import pytopkapi.pretreatment as pm

import datetime as dtm
from math import log10
import gc
import psutil
import random
import shutil


class MemoryWillRunOutError(Exception):
    pass

def memcheck(MC, MT):
    if MC > MT:
        print 'we will run out of memory.. Currently np arrays are \
holding %d MB. Threshold is %d MB.' %(MC, MT)
        raise MemoryWillRunOutError()

def memcheck_process(mp):
    '''
    mp: memory of python process [MB]
    (mc: memory of current np.array [MB])
    mt: memory treshold of python process [MB]
    mb: memory buffer - a memory threshold that will trigger an exeption
    '''
    mt = 4096
    mb = 500
    if mp > mt-mb:
        raise MemoryWillRunOutError()
    else:
        print 'Current total memory is %dMB out of available %dMB' %(mp,mt)

def memory_usage_psutil():
    # return the memory usage in MB
    process = psutil.Process(os.getpid())
    mem = process.get_memory_info()[0]/float(2**20)
    return mem

# gzip compression flag
comp = 6

def extract_Q_down(control_fname):
    """Extract combined soil and overland out flow rates.

    Read a PyTOPKAPI simulation file and return the combined overland
    andsoil store outflows in a Numpy array.

    Parameters
    ----------
    control_fname : string
        The file name of a PyTOPKAPI simulation control file. The name
        should contain the full path relative to the current
        directory.

    Returns
    -------
    Qdown : Numpy array
        A Numpy array containing the simulated outflow flow rates from
        the overland and soil store of each cell.

    """
    config = SafeConfigParser()
    config.read(control_fname)

    sim_fname = config.get('output_files', 'file_out')

    tkpi_file = h5py.File(sim_fname)
    Qdown = tkpi_file['/Q_down'][...]
    tkpi_file.close()

    return Qdown

def extract_Q_channel(control_fname):
    """Extract channel flow rates from a PyTOPKAPI simulation file.

    Read a PyTOPKAPI simulation file and return the simulated channel
    flows in a Numpy masked array.

    Parameters
    ----------
    control_fname : string
        The file name of a PyTOPKAPI simulation control file. The name
        should contain the full path relative to the current
        directory.

    Returns
    -------
    Qc : Numpy masked array
        A Numpy masked array containing the simulated flow rates for
        channel cells.

    """
    config = SafeConfigParser()
    config.read(control_fname)

    param_fname = config.get('input_files', 'file_cell_param')
    sim_fname = config.get('output_files', 'file_out')

    params = np.loadtxt(param_fname)

    tkpi_file = h5py.File(sim_fname)
    Qc = tkpi_file['/Channel/Qc_out'][...]
    tkpi_file.close()

    channel_mask = params[:, 3]
    cond = params[:, 3]*np.ones(Qc.shape, dtype=np.int) != 1

    Qc = np.ma.masked_where(cond, Qc)

    return Qc

def extract_overland_volume(control_fname):
    """Extract the volumes in the overland stores.

    Read a PyTOPKAPI simulation file and return the combined overland
    and store volumes in a Numpy array.

    Parameters
    ----------
    control_fname : string
        The file name of a PyTOPKAPI simulation control file. The name
        should contain the full path relative to the current directory
        (or the root of the file system).

    Returns
    -------
    Vo : Numpy array
        A Numpy array containing the simulated storage volume in the
        overland store of each cell.

    """
    config = SafeConfigParser()
    config.read(control_fname)

    sim_fname = config.get('output_files', 'file_out')

    tkpi_file = h5py.File(sim_fname)
    Vo = tkpi_file['/Overland/V_o'][...]
    tkpi_file.close()

    return Vo

def extract_overland_volume_to_file(sim_fname, param_fname,
                                    result_fname, start_dt, timestep):
    """Extract the volumes in the overland stores to a file.

    Read a TOPKAPI simulation file and it's associated parameter file
    and extract the overland store volumes for each timestep. Store
    the results in a new HDF5 file, grouped by date and containing
    datasets of latitude, longitude and storage volume.

    Parameters
    ----------
    sim_fname : string
        The name of a PyTOPKAPI simulation file. This should include
        the full or relative path.
    param_fname : string
        The name of a parameter file describing the catchment. This
        should include the full or relative path.
    result_fname : string
        The name of an HDF5 file to store the output. This should
        include the full or relative path.
    start_dt : datetime.datetime
        The starting date and time of the simulated results in
        `sim_fname`.
    timestep : int
        The length of each model time-step in seconds.

    Returns
    -------
    Nothing

    """
    params = np.loadtxt(param_fname)
    x = params[:, 1]
    y = params[:, 2]
    soil_depth = params[:, 8]

    soil_depth = ma.masked_values(soil_depth, 0.0)
    x = ma.array(x, mask=soil_depth.mask).compressed()
    y = ma.array(y, mask=soil_depth.mask).compressed()

    tkpi_file = h5py.File(sim_fname)
    result_file = h5py.File(result_fname, 'w')

    overland_vol = tkpi_file['/Overland/V_o'][...]
    tkpi_file.close()
    rows, cols = overland_vol.shape

    # y
    dset = result_file.require_dataset('y', shape=y.shape,
                                       dtype=np.float32, compression=comp)
    dset[...] = y

    dset.attrs['name'] = 'y coordinate'
    dset.attrs['units'] = 'Projection dependent (Metres or Decimal degrees)'

    # x
    dset = result_file.require_dataset('x', shape=x.shape,
                                       dtype=np.float32, compression=comp)
    dset[...] = x

    dset.attrs['name'] = 'x coordinate'
    dset.attrs['units'] = 'Projection dependent (Metres or Decimal degrees)'

    curr_dt = start_dt
    for k in range(rows):
        print curr_dt

        ov = ma.array(overland_vol[k], mask=soil_depth.mask).compressed()

        dset = result_file.require_dataset(curr_dt.strftime('%Y%m%d%H00'),
                                           shape=ov.shape,
                                           dtype=np.float32, compression=comp)
        dset[...] = ov

        dset.attrs['name'] = 'TOPKAPI overland store volume'
        dset.attrs['units'] = 'm^3'

        curr_dt += timedelta(seconds=timestep)

    tkpi_file.close()
    result_file.close()

def extract_ssi(control_fname):
    """Extract SSI from a PyTOPKAPI simulation file.

    Read a PyTOPKAPI simulation file and it's associated parameter
    file and compute the Soil Saturation Index (SSI) for each model
    cell and timestep. The results are returned as a Numpy array.

    Parameters
    ----------
    control_fname : string
        The file name of a PyTOPKAPI simulation control file. The name
        should contain the full path relative to the current
        directory.

    Returns
    -------
    ssi : Numpy ndarray
        A Numpy array containing the calculated SSI values.

    """
    config = SafeConfigParser()
    config.read(control_fname)

    global_param_fname = config.get('input_files', 'file_global_param')
    param_fname        = config.get('input_files', 'file_cell_param')
    sim_fname          = config.get('output_files', 'file_out')
    fac_L              = config.getfloat('calib_params', 'fac_L')
    fac_th_s           = config.getfloat('calib_params', 'fac_th_s')

    params = np.loadtxt(param_fname)
    glob_params = np.genfromtxt(global_param_fname, names=True)

    soil_depth = fac_L*params[:, 8]
    factor = params[:, 11] - params[:, 10]*fac_th_s         # theta r, theta s
    cell_area = glob_params['X']**2 # m^2

    soil_depth = ma.masked_values(soil_depth, 0.0)
    factor = ma.array(factor, mask=soil_depth.mask)
    div = factor*soil_depth*cell_area

    tkpi_file = h5py.File(sim_fname)
    soil_vol = tkpi_file['/Soil/V_s'][...]
    tkpi_file.close()

    # ssi = (Vs/cell_vol)*100
    # cell_vol = (theta_s - theta_r)*soil_depth*cell_area
    sv = ma.array(soil_vol, mask=soil_depth.mask)
    ssi = (sv/(div))*100.0

    return ssi

def extract_ssi_to_file(sim_fname, param_fname,
                        result_fname, start_dt, timestep):
    """Extract percentage saturation to a file

    Read a TOPKAPI simulation file and it's associated parameter file
    and compute the SSI for each timestep. Store the results in a new
    HDF5 file, grouped by date and containing datasets of latitude,
    longitude and SSI value.

    Parameters
    ----------
    sim_fname : string
        The name of a PyTOPKAPI simulation file. This should include
        the full or relative path.
    param_fname : string
        The name of a parameter file describing the catchment. This
        should include the full or relative path.
    result_fname : string
        The name of an HDF5 file to store the output. This should
        include the full or relative path.
    start_dt : datetime.datetime
        The starting date and time of the simulated results in
        `sim_fname`.
    timestep : int
        The length of each model time-step in seconds.

    Returns
    -------
    Nothing

    """
    params = np.loadtxt(param_fname)
    x = params[:, 1]
    y = params[:, 2]

    soil_depth = params[:, 8] # need multiplier for fac_l
    factor = params[:, 11] - params[:, 10] # Need multiplier for fac_th_s
    #cell_area = 1000.0**2 # m^2 need to read from global param

    soil_depth = ma.masked_values(soil_depth, 0.0)
    factor = ma.array(factor, mask=soil_depth.mask)
    x = ma.array(x, mask=soil_depth.mask).compressed()
    y = ma.array(y, mask=soil_depth.mask).compressed()

    div = factor*soil_depth*cell_area

    tkpi_file = h5py.File(sim_fname)
    result_file = h5py.File(result_fname, 'w')

    soil_vol = tkpi_file['/Soil/V_s'][...]
    tkpi_file.close()
    rows, cols = soil_vol.shape

    # y
    dset = result_file.require_dataset('y', shape=y.shape,
                                       dtype=np.float32, compression=comp)
    dset[...] = y

    dset.attrs['name'] = 'y coordinate'
    dset.attrs['units'] = 'Projection dependent (Metres or Decimal degrees)'

    # x
    dset = result_file.require_dataset('x', shape=x.shape,
                                       dtype=np.float32, compression=comp)
    dset[...] = x

    dset.attrs['name'] = 'x coordinate'
    dset.attrs['units'] = 'Projection dependent (Metres or Decimal degrees)'

    curr_dt = start_dt
    for k in range(rows):
        #print curr_dt
        # ssi = (Vs/cell_vol)*100
        # cell_vol = (theta_s - theta_r)*soil_depth*cell_area
        sv = ma.array(soil_vol[k], mask=soil_depth.mask)
        ssi = (sv/(div))*100.0

        ssi = ssi.compressed()

        # ssi
        dset = result_file.require_dataset(curr_dt.strftime('%Y%m%d%H00'),
                                           shape=ssi.shape,
                                           dtype=np.float32, compression=comp)
        dset[...] = ssi

        dset.attrs['name'] = 'TOPKAPI soil saturation index'
        dset.attrs['units'] = '% saturation'

        curr_dt += timedelta(seconds=timestep)

    tkpi_file.close()
    result_file.close()
    
def run_version_log(control_fname, run_version_fname):
    """Write some parameters from TOPKAIP.ini file to a version
    file including the version number. This helps to keep track
    of changes made during calibration.

    Parameters
    ----------
    control_fname : string
        The file name of a PyTOPKAPI simulation control file.
    run_version_fname : string
        The file name of the version log file.
    """
    
    config = SafeConfigParser()
    config.read(control_fname)
    
    # Read the directory and file name for the change log-file
    fn_change_log_out = config.get('output_files', 'file_change_log_out')
    
    # Write the current calibration parameters to the log
    f_out = open(fn_change_log_out, 'a')
    
    f_out.write('______________________________\n')
    for name in ConfigParser.ConfigParser.items(config, 'calib_params'):
        f_out.write(name[0] + ' = ' + name[1] + '\n')
        
    # Read the version number from the config file
    config.read(run_version_fname)
    ver_n = config.getint('version_number', 'version')
    
    # Now update the version number for the next run and overwrite the config file
    ver_n += 1
    
    config = ConfigParser.RawConfigParser()
    
    config.add_section('version_number')
    config.set('version_number', 'version', ver_n)
    
    with open(run_version_fname, 'w') as configfile:
        config.write(configfile)
    
    f_out.write('    Run Version = ' + str(ver_n) + '\n')
    f_out.close()

    
def Extract_flow_at_outlet(plot_ini='plot-flow-precip.ini',
                           start_stamp='', split_month=10,
                           topkapi_ini = 'TOPKAPI.ini',
                           ver_ini = 'run_version.ini',
                           WWBU_transfer_n2i = False):
    """ The following simulated variables are written to an excel file:
    - Rainfall (mm)
    - Runoff at outlet (L/T)
    Temporal resolution of outputs are daily, monthly, yearly and a mean
    annually value.
    - Percentiles (0.1 to .99) of the monthly flows
    - Month by month averages
    - Water Balance values (P = ET + Qo + Qs + dS) for mean catchment
      conditions in daily resolution
    - Water Balance values (P = ET + Qo + Qs + dS) for all Wetland Water
      Balance (WWB) units, all in m3/T and in mm respectively
    
    A flow duration curve based on monthly values is created.

    Parameters
    ----------
    plot_ini : string
        The file name of a PyTOPKAPI simulation control file.
    start_stamp : string
        The date when the simulation starts
    split_month : integer
        The month which defines the aveaging period for one year, e.g. hydro-
        logical year commences at Oktober 1.
    topkapi_ini : string
        TOPKAPI_ini file
    version_ini : string
        Versin number for the current run
    
    Returns: nothing
    Produces: Flow duration graph and one XLS-file.
    """
    mp = memory_usage_psutil()
#    print 'Process memory at %dMB (after loading all modules)' %mp
        
    config = SafeConfigParser()
    
    # Read the version number of the run to be processed
    config.read(ver_ini)
    
    ver_num = config.get('version_number','version')
    
    # Read the plot_ini file 
    config.read(plot_ini)
    
    file_Qsim  = config.get('files','file_Qsim')
    file_rain  = config.get('files','file_rain')

    group_name = config.get('groups','group_name')
    
    outlet_ID  = config.getint('parameters', 'outlet_ID')
    
    run_name   = os.getcwd().split('\\')[-2].split('_')[1]

    fn_Qout    = 'results\\Qsim_outlet_'+run_name+'_'+str(ver_num)+'.xlsx'
    fn_Qout2   = 'results\\Qsim_outlet_Static.xlsx'
    
#==============================================================================
#     Read some parameters from the TOPKAPI.ini file
#==============================================================================
    config = SafeConfigParser()
    config.read(topkapi_ini)
    
#    print 'Read the file ', topkapi_ini
    #print 'config', config.read(topkapi_ini)

    file_global_param = config.get('input_files', 'file_global_param')
    param_fname       = config.get('input_files', 'file_cell_param')
    rain_dist         = config.getboolean('forcing_options', 'rain_dist')
    fac_L             = config.getfloat('calib_params', 'fac_l')
    fac_th_s          = config.getfloat('calib_params', 'fac_th_s')

    try:
        Dams_on       = config.getboolean('forcing_options', 'Dams_on')
        WL_dam_avg    = config.getfloat('forcing_options', 'Average_WL_dam')
    except NoOptionError:
        Dams_on = 0
    Piezos_on = 1
    try:
        piezos_tmp = config.get('piezometer_locations', 'piezos')
        if len(piezos_tmp.split(',')) > 1:
            piezos = []
            for v in piezos_tmp.split(','):
                piezos.append(int(v))
            nb_piezos = len(piezos)
        elif piezos_tmp == '':
            Piezos_on = 0
            print 'No piezometer locations specified...'
        else:
            piezos = [int(piezos_tmp)]
            nb_piezos = 1
    except NoSectionError:
        Piezos_on = 0
        print 'No piezometer locations specified...'
    
    if Piezos_on:
        try:
            config.read(plot_ini)
            fn_in_pzo = config.get('files','piezom_in')
            pzo_plot_on = 1
            # Read the piezo data file
            pzo_data = []
            #len_pzo_file = ut.file_len(fn_in_pzo)
            cnt = 0
            with open(fn_in_pzo) as f_pzo:
                for l in f_pzo:
                    if cnt:
                        l_tmp1 = l.replace('/', '\t').split('\t')
                        l_tmp2 = [float(x) for x in l_tmp1]
                        pzo_data.append(l_tmp2)
                        
                    else:
                        print 'first line', l
                        pzo_hdr  = [int(x) for x in l.split('\t')[1:]]
                    cnt += 1
            print 'pzo_hdr', pzo_hdr
            pzo_data = np.array(pzo_data)
            print 'pzo_data.shape', pzo_data.shape
        except NoOptionError:
            pzo_plot_on = 0
            print 'No piezometer data file specified in %s...' %plot_ini
        if pzo_plot_on:
            for p in piezos:
                if not p in pzo_hdr:
                    pzo_hdr_s = ''
                    piezos_s  = ''
                    for i in pzo_hdr:
                        pzo_hdr_s += str(i)+', '
                    for i in piezos:
                        piezos_s += str(i)+', '
                    
                    raise ValueError('At least one of the piezometer names %s \
specified in the piezometer file %s mismatch with the the piezometer cell \
lables are spcified %s specified in %s ' \
%(pzo_hdr_s[:-2], fn_in_pzo, piezos_s[:-2], topkapi_ini))

    X,Dt,alpha_s,alpha_o,alpha_c,A_thres,W_min,W_max\
      = pm.read_global_parameters(file_global_param)
    
    ar_cell_label, ar_coorx, ar_coory, ar_lambda, ar_Xc, ar_wwb, \
    ar_tan_beta, ar_tan_beta_channel, ar_L0, ar_Ks0, ar_theta_r, \
    ar_theta_s, ar_n_o0, ar_n_c0, ar_cell_down, ar_pVs_t0, ar_Vo_t0, \
    ar_Qc_t0, ar_kc, psi_b, lamda = pm.read_cell_parameters(param_fname)
    
    mp = memory_usage_psutil()
#    print 'Process memory at %dMB (after reading cell parameter file)' %mp
#==============================================================================
#     READING DATA FROM THE RESULT OR FORCING FILE
#==============================================================================
    ftype   = 'float32'   # float type
    # SSI
    ar_SSI_s = np.average(extract_ssi(topkapi_ini),axis=1)
    
    # Rainfall
    with h5py.File(file_rain, 'r') as f_r:
        ndar_rainfl = f_r['/'+group_name+'/']['rainfall'].value.astype(dtype=ftype)
    ar_rain = np.average(ndar_rainfl, axis=1)    # in mm/TMSP
    
    mp = memory_usage_psutil()
#    print 'Process memory at %dMB (after reading rainfall data)' %mp
#    h5file_in = h5.openFile(file_rain, mode='r')
#    group='/'+group_name+'/'
#    node = h5file_in.getNode(group+'rainfall')
#    ndar_rainfl = node.read()
#    h5file_in.close()
    #Compute the mean catchment rainfall
#    ar_rain = np.average(ndar_rainfl, axis=1)
#    print 'ar_rain ', ar_rain.shape
    #print 'ndar_rainfl.shape', ndar_rainfl.shape
    
    #Read the simulated Q channel data
#    ndar_Qc_out = ut.read_one_array_hdf(file_Qsim,'/Channel/','Qc_out')
#    ar_Qsim = ndar_Qc_out[1:,outlet_ID]
#   Change it from individual array to the same array to overwrite space
#   otherwise there is too much garbage collection.

#    ndar_tmp = ut.read_one_array_hdf(file_Qsim,'/Channel/','Qc_out')
#    ar_Qsim = ndar_tmp[1:,outlet_ID]
#    print 'Read Qsim data..'

    mt      = 250      # Memory treshold in MB
    b2mb    = float(2**20)  # Bytes to MB devision factor
    #shapexy = (7218, 8877)
    
    #mp = memory_usage_psutil()
    #print 'Process memory at %dMB' %mp

    try:
        with h5py.File(file_Qsim, 'r') as f:
            mp = memory_usage_psutil()
#            print 'Process memory at %dMB - opened h5py result file' %mp
            #raise MemoryWillRunOutError
            ar_Qsim = f['Channel']['Qc_out'].value[1:,outlet_ID].astype(dtype=ftype)
            mc = ar_Qsim.nbytes/b2mb
            ndar_ET_out = f['/']['ET_out'].value.astype(dtype=ftype)
            shapexy = ndar_ET_out.shape
            m_ndar_ET_out = ndar_ET_out.nbytes/b2mb
            mc += m_ndar_ET_out
#            print 'Process memory at %dMB - ndar_ET_out' %mp
            ar_ET_a = np.average(ndar_ET_out[1:,], axis=1)   # mm/TMSP
            # m3/TSTP for all channel cells
            ar_Echn = np.sum(f['Channel']['Ec_out'].value[1:,], axis=1).astype(dtype=ftype)
            mc += ar_Echn.nbytes/b2mb   # Units in m3/TMSP
#            memcheck_process(mp)
#            if ndar_ET_out.nbytes/b2mb > 1024:                
#                raise MemoryWillRunOutError
#            memcheck_process(mp)
            mp = memory_usage_psutil()
            mc += ar_ET_a.nbytes/b2mb
#            memcheck_process(mp)

            ndar_Qs_out = f['Soil']['Qs_out'].value.astype(dtype=ftype)
            mc += ndar_Qs_out.nbytes/b2mb
#            memcheck_process(mp)
            mp = memory_usage_psutil()
#            print 'Process memory at %dMB - ndar_Qs_out' %mp
            
            ndar_Qs_out_ = f['Soil']['Qs_outII'].value.astype(dtype=ftype)
            mc += ndar_Qs_out_.nbytes/b2mb
#            memcheck_process(mp)
            mp = memory_usage_psutil()
            
            ndar_Vs_cel = f['Soil']['V_s'].value.astype(dtype=ftype)
            ar_Vs = np.sum(ndar_Vs_cel[1:,], axis=1)   # in m3/catchment
            mc += ndar_Vs_cel.nbytes/b2mb
#            memcheck_process(mp)
            mp = memory_usage_psutil()
            
            ndar_Vs_cel_ = f['Soil']['V_sII'].value.astype(dtype=ftype)
            ar_Vs = np.sum(ndar_Vs_cel_[1:,], axis=1)   # in m3/catchment
            mc += ndar_Vs_cel_.nbytes/b2mb
#            memcheck_process(mp)
            mp = memory_usage_psutil()
#            print 'Process memory at %dMB - ndar_Vs_cel' %mp
            
            ndar_Qo_out = f['Overland']['Qo_out'].value.astype(dtype=ftype)
            mc += ndar_Qo_out.nbytes/b2mb
#            memcheck_process(mp)
            mp = memory_usage_psutil()
#            print 'Process memory at %dMB - ndar_Qo_out' %mp
            
            if Dams_on:
                ndar_Vd_cel = f['Dam']['V_d'].value.astype(dtype=ftype)
                mc += ndar_Vd_cel.nbytes/b2mb
#                memcheck_process(mp)
                mp = memory_usage_psutil()
#                print 'Process memory at %dMB - ndar_Vd_cel' %mp
                
                ndar_Qd_inp = f['Dam']['Qd_in'].value.astype(dtype=ftype)
                mc += ndar_Qd_inp.nbytes/b2mb
#                memcheck_process(mp)
                mp = memory_usage_psutil()
#                print 'Process memory at %dMB - ndar_Qd_inp' %mp
                
                ndar_Qd_out = f['Dam']['Qd_out'].value.astype(dtype=ftype)
                mc += ndar_Qd_out.nbytes/b2mb
#                memcheck_process(mp)
                mp = memory_usage_psutil()
#                print 'Process memory at %dMB - ndar_Qd_out' %mp
    except MemoryWillRunOutError or MemoryError:
        varnames = ['ndar_ET_out','ndar_Qs_out','ndar_Qs_out_','ndar_Vs_cel',
                    'ndar_Vs_cel_','ndar_Qo_out','ndar_Vd_cel','ndar_Qd_inp',
                    'ndar_Qd_out']
        print 'Too much memory was used up loading all data into memory. Large\
 arrays will be mapped on harddrive - long processing times can be expected..'
        for i in xrange(len(varnames)):
            if varnames[i] in locals():
                MBsize = eval(varnames[i] + '.nbytes/(1024.**2)')
                print varnames[i], '=', '%0.2fMB.' %MBsize,'Deleting', varnames[i], '..'
                #print varnames[i],':', 'shape:', locals()[varnames[i]].shape
                exec 'locals()[varnames[i]] = 0' in globals(), locals()
                #print 'lets see now\n:', locals()[varnames[i]]
                #exec 'del locals()[v]' does not work for some reason
                #exec '%s = np.resize(%s,(1,1))' %(v,v) not working!
                #print 'now lets see the shape\n', locals()[v].shape
                #del locals()[v]
        #print ndar_ET_out
        mp = memory_usage_psutil()
        print 'Process memory at %dMB - after deleting ndarrays' %mp
        print 'starting to run garbage collector..'
        for r in xrange(3):
            collected = gc.collect()
            print r+1, ". run garbage collector: collected %d objects." % (collected)
        mp = memory_usage_psutil()
        print 'Process memory at %dMB - after GC' %mp

        varnms = ['ndar_ET_out','ndar_Qs_out','ndar_Qs_out_','ndar_Vs_cel',
                  'ndar_Vs_cel','ndar_Qo_out','ndar_Vd_cel','ndar_Qd_inp',
                  'ndar_Qd_out']
        groups = ['/',     '/Soil/','/Soil/',  'Soil','Soil', 'Overland','Dam','Dam',  'Dam']
        entris = ['ET_out','Qs_out','Qs_outII','V_s', 'V_sII','Qo_out',  'V_d','Qd_in','Qd_out']
        
        i = 0
        with h5py.File(file_Qsim, 'r') as f:
            for n in varnms:
                exec 'locals()[n] = np.memmap(\'Results/TMP/\'+n, dtype=ftype, mode=\'w+\', shape=shapexy)'
    #            with h5py.File(file_Qsim, 'r') as f:
                exec 'locals()[n][:] = f[groups[i]][entris[i]].value.astype(dtype=ftype)'
                exec 'locals()[n][:].flush()'
                exec 'locals()[n] = []'
                i += 1
                mp = memory_usage_psutil()
                print 'Process memory at %dMB - mmap ndar_ET_out' %mp
                gc.collect()

        for r in xrange(3):
            collected = gc.collect()
            print r+1, ". run garbage collector: collected %d objects." % (collected)
        mp = memory_usage_psutil()
        print 'Process memory at %dMB - after GC' %mp
        
        for n in varnms:
            exec 'locals()[n] = np.memmap(\'Results/TMP/\'+n, dtype=ftype, mode=\'r\', shape=shapexy)'
        
        '''
        ndar_ET_out    = np.memmap('Results\\TMP\\ndar_ET_out', dtype=ftype, mode='w+', shape=shapexy)
        ndar_ET_out.flush()
        mp = memory_usage_psutil()
        print 'Process memory at %dMB - mmap ndar_ET_out' %mp
        ndar_Qs_out    = np.memmap('Results\\TMP\\ndar_Qs_out', dtype=ftype, mode='w+', shape=shapexy)
        ndar_Qs_out.flush()
        mp = memory_usage_psutil()
        print 'Process memory at %dMB - mmap ndar_Qs_out' %mp
        ndar_Vs_cel    = np.memmap('Results\\TMP\\ndar_Vs_cel', dtype=ftype, mode='w+', shape=shapexy)
        ndar_Vs_cel.flush()
        mp = memory_usage_psutil()
        print 'Process memory at %dMB - mmap ndar_Vs_cel' %mp
        ndar_Qo_out    = np.memmap('Results\\TMP\\ndar_Qo_out', dtype=ftype, mode='w+', shape=shapexy)
        ndar_Qo_out.flush()
        mp = memory_usage_psutil()
        print 'Process memory at %dMB - mmap ndar_Qo_out' %mp
        if Dams_on:
            ndar_Vd_cel    = np.memmap('Results\\TMP\\ndar_Vd_cel', dtype=ftype, mode='w+', shape=shapexy)
            ndar_Vd_cel.flush()
            mp = memory_usage_psutil()
            print 'Process memory at %dMB - mmap ndar_Vd_cel' %mp
            ndar_Qd_inp    = np.memmap('Results\\TMP\\ndar_Qd_inp', dtype=ftype, mode='w+', shape=shapexy)
            ndar_Qd_inp.flush()
            mp = memory_usage_psutil()
            print 'Process memory at %dMB - mmap ndar_Qd_inp' %mp
            ndar_Qd_out    = np.memmap('Results\\TMP\\ndar_Qd_out', dtype=ftype, mode='w+', shape=shapexy)
            ndar_Qd_out.flush()
            mp = memory_usage_psutil()
            print 'Process memory at %dMB - mmap ndar_Qd_out' %mp
        for r in xrange(2):
            collected = gc.collect()
            mp = memory_usage_psutil()
            print r+1, ". run garbage collector: collected %d objects." % (collected)
            print 'Process memory at %dMB - after GC' %mp

        with h5py.File(file_Qsim, 'r') as f:
            mp = memory_usage_psutil()
            print 'Process memory at %dMB - after opening h5py result file' %mp
            # find out what the dimensions of the ndarray are
            #ar_Qsim        = f['Channel']['Qc_out'].value[1:,outlet_ID].astype(dtype=ftype)
            ndar_ET_out[:] = f['/']['ET_out'].value.astype(dtype=ftype)
            mp = memory_usage_psutil()
            print 'Process memory at %dMB - after writing first ndarray into mmap array' %mp
            #ar_ET_a        = np.average(ndar_ET_out[1:,], axis=1).astype(dtype=ftype)
            #ar_Echn        = np.average(f['Channel']['Ec_out'].value[1:,], axis=1).astype(dtype=ftype)
            #ar_Vs          = np.average(f['/Soil/']['V_s'].value[1:,], axis=1).astype(dtype=ftype)
            ndar_Qs_out[:] = f['Soil']['Qs_out'].value.astype(dtype=ftype)
            ndar_Vs_cel[:] = f['Soil']['V_s'].value.astype(dtype=ftype)
            ndar_Qo_out[:] = f['Overland']['Qo_out'].value.astype(dtype=ftype)
            if Dams_on:
                ndar_Vd_cel[:] = f['Dam']['V_d'].value.astype(dtype=ftype)
                ndar_Qd_inp[:] = f['Dam']['Qd_in'].value.astype(dtype=ftype)
                ndar_Qd_out[:] = f['Dam']['Qd_out'].value.astype(dtype=ftype)


    #print 'ar_Qsim ', ar_Qsim.shape

    # Read Evaporation data and average it for the catchment
    #h5file_in = h5.openFile(file_Qsim, mode='r')

#    #group='/'+group_name+'/'
#    node = h5file_in.getNode('/','ET_out')
#    ndar_ET_out = node.read()
#    #Compute the mean catchment ETa for each time step
#    ar_ET_a = np.average(ndar_ET_out[1:,], axis=1)
    #print 'ar_ET_a.shape = ', ar_ET_a.shape

    # New approach using read one array method
#    ndar_tmp = ut.read_one_array_hdf(file_Qsim,'/','ET_out')
#    ar_Qsim = ndar_tmp[1:,]
#    ar_ET_a = np.average(ndar_tmp[1:,], axis=1)
#    print 'Reading ET_a data complete..'
    
    with h5py.File(file_Qsim, 'r') as f:
        ndar_ET_out = f['/']['ET_out'].value
        ar_ET_a = np.average(ndar_ET_out[1:,], axis=1)
        
#    # Read evaporation from channel cells (m3/d)
#    node = h5file_in.getNode('/Channel/','Ec_out')
#    ndar_Ec_out = node.read()
#    ar_Echn = np.average(ndar_Ec_out[1:,], axis=1)
#    print 'Read Ec data..'
    
    # New approach using read one array method
#    ndar_tmp = ut.read_one_array_hdf(file_Qsim,'/Channel/','Ec_out')
#    ar_Echn = np.average(ndar_tmp[1:,], axis=1)
#    print 'Reading ET_a data complete..'
    
    with h5py.File(file_Qsim, 'r') as f:
        ar_Echn = np.average(f['Channel']['Ec_out'].value[1:,], axis=1)

    # Read Soil Qs_out
#    node = h5file_in.getNode('/Soil/','Qs_out')
#    ndar_Qs_out = node.read()
#    print 'Read Qs data..'

    with h5py.File(file_Qsim, 'r') as f:
        ndar_Qs_out = f['Soil']['Qs_out'].value

    # Read Soil Vs_cel
#    node = h5file_in.getNode('/Soil/','V_s')
#    ndar_Vs_cel = node.read()
#    print 'Read Vs data..'

    with h5py.File(file_Qsim, 'r') as f:
        ndar_Vs_cel = f['Soil']['V_s'].value

    # Read Overland Qo_out
#    node = h5file_in.getNode('/Overland/','Qo_out')
#    ndar_Qo_out = node.read()
#    print 'Read Qo data..'

    with h5py.File(file_Qsim, 'r') as f:
        ndar_Qo_out = f['Overland']['Qo_out'].value

    if Dams_on:
        # Read Dam volume data
#        node = h5file_in.getNode('/Dam/','V_d')
#        ndar_Vd_cel = node.read()
        with h5py.File(file_Qsim, 'r') as f:
            ndar_Vd_cel = f['Dam']['V_d'].value
    
        # Read the dam inflow
#        node = h5file_in.getNode('/Dam/','Qd_in')
#        ndar_Qd_inp = node.read()
        with h5py.File(file_Qsim, 'r') as f:
            ndar_Qd_inp = f['Dam']['Qd_in'].value
    
#        node = h5file_in.getNode('/Dam/','Qd_out')
#        ndar_Qd_out = node.read()
        with h5py.File(file_Qsim, 'r') as f:
            ndar_Qd_out = f['Dam']['Qd_out'].value

#   h5file_in.close()'''

#==============================================================================
#     Compute SSI for each cell
#==============================================================================
    # create ndarray for SSI
    nb_cells    = ndar_ET_out.shape[1]
    n_rec       = ar_Qsim.shape[0]
    _1cell_area = X**2

    ndar_SSI_cl  = np.zeros((n_rec, nb_cells)).astype(dtype=ftype)
    ndar_SSI_cl_ = np.zeros((n_rec, nb_cells)).astype(dtype=ftype)
    
    div = (ar_theta_s*fac_th_s - ar_theta_r)*ar_L0*fac_L*_1cell_area

    for r in xrange(n_rec):
        ndar_SSI_cl[r,:]  = (ndar_Vs_cel[r,:] /div)*100.
        ndar_SSI_cl_[r,:] = (ndar_Vs_cel_[r,:]/div)*100.

#    print 'shape of ndar_SSI_cl:', ndar_SSI_cl.shape

    # Check if outflow cell reports runoff
    nb_Qflow = ar_Qsim[ar_Qsim > 0.]
    
    if len(nb_Qflow) < ar_Qsim.shape[0]:
        print 'probably wrong cell number for outflow node:'
        print 'len(nb_Qflow) < nb_cells', len(nb_Qflow), ar_Qsim.shape[0]

#==============================================================================
# DAM VOLUMES
#==============================================================================

    if Dams_on:
        ar_dam_u = np.unique(ar_wwb[ar_wwb>199])
        
        dam_idx = []   # a list unique cell labels for all dams
        
        for d in ar_dam_u:
            d_tmp = np.where(ar_wwb == d)[0]
            dam_idx.append(d_tmp)
            #dam_idx.append(d_tmp[0]) for some reason does not work.
        
        Vd_cel_d = []
        Qd_input = []    
        Qd_outpt = []
    
        nb_dam = len(ar_dam_u)
        print 'Number of Dams', nb_dam
        print 'Cell labels dams', dam_idx
        #print 'Cell labels dams', dam_idx[0]
        #print 'size of first dam array:', ndar_Vd_cel[1:,dam_idx[0]].size
        #print ndar_Vd_cel[1:,dam_idx[0]]
    
        for d in xrange(nb_dam):
            #print ndar_Vd_cel[1:,dam_idx[d]]
            #print dam_idx[d]
            Vd_cel_d.append(np.sum(ndar_Vd_cel[1:,dam_idx[d]], axis=1))
            Qd_input.append(np.sum(ndar_Qd_inp[1:,dam_idx[d]], axis=1))
            Qd_outpt.append(np.sum(ndar_Qd_out[1:,dam_idx[d]], axis=1))
        Vd_cel_d = np.array(Vd_cel_d)
        Qd_input = np.array(Qd_input)
        Qd_outpt = np.array(Qd_outpt)


#==============================================================================
# WETLAND UNIT WATER BALANE DATA QUERIES BASED ON WETLAND UNIT LABELS
#==============================================================================

    # Filter rivers out of of WWB cells
    ar_wwb[ar_lambda==1] = -9999.

    ar_wwb_u = np.unique(ar_wwb[ar_wwb>0])
    ar_wwb_u = ar_wwb_u[ar_wwb_u<200]
    
    tmp = [int(x) for x in ar_wwb_u]
    
    #print 'Extracting information from WWB-units', tmp
    
    wwb_idx = []

    for u in ar_wwb_u:
        u_tmp = np.where(ar_wwb == u)[0]
        wwb_idx.append(u_tmp)

    nb_wwb = len(ar_wwb_u)
    
    ET_out_w  = []  # units were mm/d -> m/d. m/d*m2 = m3/d
    Vs_cel_w  = []  # units already in m3
    Vs_cel_w_ = []  # units already in m3
    Rainfl_w  = []  # units were mm/d -> m/d. m/d*m2 = m3/d
    SSIndx_w  = []  # Dimensionless
    SSIndx_w_ = []  # Dimensionless
    wwbu_A_w  = []  # Area of each wetland unit
    
    for w in xrange(nb_wwb):
        ET_out_w.append(np.sum(ndar_ET_out[1:,wwb_idx[w]], axis=1)*1e-3*X**2)
        Vs_cel_w.append(np.sum(ndar_Vs_cel[1:,wwb_idx[w]], axis=1))
        Vs_cel_w_.append(np.sum(ndar_Vs_cel_[1:,wwb_idx[w]], axis=1))
        SSIndx_w.append( np.average(ndar_SSI_cl[ :,wwb_idx[w]], axis=1))
        SSIndx_w_.append(np.average(ndar_SSI_cl_[:,wwb_idx[w]], axis=1))
        wwbu_A_w.append(len(wwb_idx[w])*X**2)
        if rain_dist:
            Rainfl_w.append(np.sum(ndar_rainfl[:,wwb_idx[w]], axis=1)*1e-3*X**2)
        else:
            Rainfl_w.append(ar_rain*len(wwb_idx[w])*1e-3*X**2)
    
    ET_out_w  = np.array(ET_out_w)    # m3/d
    Vs_cel_w  = np.array(Vs_cel_w)    # m3
    Vs_cel_w_ = np.array(Vs_cel_w)    # m3
    Rainfl_w  = np.array(Rainfl_w)    # m3/d
    SSIndx_w  = np.array(SSIndx_w)    # Dimensionless
    SSIndx_w_ = np.array(SSIndx_w)    # Dimensionless
    wwbu_A_w  = np.array(wwbu_A_w)    # m2
    
#    print 'shape of SSIndx_w:', SSIndx_w.shape
#    print 'shape of Vs_cel_w:', Vs_cel_w.shape

    for r in Rainfl_w:
        for rr in r:
            if type(rr) != np.float32:
                print 'Rain AAAAAAAAALAAAAAAAAAAAAARM\n\nALAAAAAAAAAAAAAAAAAAAAAARM\n'
    for e in ET_out_w:
        for ee in e:
            if type(ee) != np.float32:
                print 'ET AAAAAAAAALAAAAAAAAAAAAARM\n\nALAAAAAAAAAAAAAAAAAAAAAARM\n'
    
    #print 'wwb_idx', wwb_idx
    # print 'len ET_out_w', len(ET_out_w)
    # print 'shape of each unit:', [x.shape[0] for x in ET_out_w]
    #print ET_out_w
#    ET_out_w = np.array(ET_out_w)
#    print ET_out_w
    # print 'ET_out_w.shape', ET_out_w.shape

#==============================================================================
# CATCHMENT - WETLAND BOUNDARY CAPTURE Qs_out FLOW
#==============================================================================
    # all cell downs within a catchment, not in the wetland
    ar_cell_down_C = ar_cell_down[ar_wwb < 1]
    # all cell downs within a WWB unit, not in the catchment
    ar_cell_down_W = ar_cell_down[ar_wwb > 0]
    # cell labels for cell without WWB units
    ar_cell_labl_C = ar_cell_label[ar_wwb < 1]
    # cell labels for cell with WWB units
    ar_cell_labl_W = ar_cell_label[ar_wwb > 0]
    #ar_cell_catch  = ar_cell_label[ar_lambda < 1] # only cells without rivers
    ar_cell_labl_R = ar_cell_label[ar_lambda==1]
    
    # Make list of indices wehre ar_wwb is a balacne area
    # BCW = Boundary Catchment -> Wetland
    # BWC = Boundary Wetland   -> Catchment
    
    ar_bcw = []
    ar_bwc = []
    ar_bwr = []
    
    for u in xrange(nb_wwb):
        bcw_tmp = []
        bwc_tmp = []
        #ar_cell_labl_C_u = ar_cell_label[wwb_idx < 1]
        for i in xrange(len(ar_cell_label)):
            if ar_cell_down[i] in wwb_idx[u] and ar_cell_label[i] not in wwb_idx[u]:
                bcw_tmp.append(ar_cell_label[i])
            if ar_cell_down[i] not in wwb_idx[u] and ar_cell_label[i] in wwb_idx[u]:
                bwc_tmp.append(ar_cell_label[i])
        ar_bcw.append(bcw_tmp)
        ar_bwc.append(bwc_tmp)
        
        bwr_tmp = []    # Boundary Wetland to River
        #ar_cell_labl_C_u = ar_cell_label[wwb_idx < 1]
        for i in xrange(len(ar_cell_label)):
#            if ar_cell_down[i] in wwb_idx[u] and ar_cell_down[i] in ar_cell_labl_R:
            if ar_cell_label[i] in wwb_idx[u] and ar_cell_down[i] in ar_cell_labl_R:
                bwr_tmp.append(ar_cell_label[i])
        #print 'WWB:',u,' WBR:',bwr_tmp
        ar_bwr.append(bwr_tmp)

#==============================================================================
# WETLAND SOIL WATER BALANCE IN - OUT FLOWS
#==============================================================================

    print 'Number of time steps =', n_rec

    # sum up soil flow across the boundaries
    Qs_in_w  = np.zeros((nb_wwb, n_rec))   # from m3/s to m3/d
    Qs_ot_w  = np.zeros((nb_wwb, n_rec))   # from m3/s to m3/d
    Qs_in_w_ = np.zeros((nb_wwb, n_rec))   # from m3/s to m3/d
    Qs_ot_w_ = np.zeros((nb_wwb, n_rec))   # from m3/s to m3/d
    Qo_in_w  = np.zeros((nb_wwb, n_rec))   # from m3/s to m3/d
    Qo_ot_w  = np.zeros((nb_wwb, n_rec))   # from m3/s to m3/d
    Qs_in_r  = np.zeros((nb_wwb, n_rec))   # from m3/s to m3/d
    Qs_in_r_ = np.zeros((nb_wwb, n_rec))   # from m3/s to m3/d
    Qo_in_r  = np.zeros((nb_wwb, n_rec))   # from m3/s to m3/d
    
    for u in xrange(nb_wwb):
        Qs_in_w[u]  = np.sum(ndar_Qs_out[1:,ar_bcw[u]],  axis=1)*Dt
        Qs_ot_w[u]  = np.sum(ndar_Qs_out[1:,ar_bwc[u]],  axis=1)*Dt
        Qs_in_w_[u] = np.sum(ndar_Qs_out_[1:,ar_bcw[u]], axis=1)*Dt
        Qs_ot_w_[u] = np.sum(ndar_Qs_out_[1:,ar_bwc[u]], axis=1)*Dt
        Qo_in_w[u]  = np.sum(ndar_Qo_out[1:,ar_bcw[u]],  axis=1)*Dt
        Qo_ot_w[u]  = np.sum(ndar_Qo_out[1:,ar_bwc[u]],  axis=1)*Dt
        Qs_in_r[u]  = np.sum(ndar_Qs_out[1:,ar_bwr[u]],  axis=1)*Dt
        Qs_in_r_[u] = np.sum(ndar_Qs_out_[1:,ar_bwr[u]], axis=1)*Dt
        Qo_in_r[u]  = np.sum(ndar_Qo_out[1:,ar_bwr[u]],  axis=1)*Dt
        # Units in m3/d
    
    # Make lists for coordinates
    ar_bcw_x,ar_bcw_y,ar_bwc_x,ar_bwc_y,bcw_wwb_l,bwc_wwb_l = [],[],[],[],[],[]
    ar_wwb_x, ar_wwb_y, ar_wwb_l       = [], [], []
    ar_wwb_C_x, ar_wwb_C_y, ar_wwb_C_l = [], [], []
    ar_wwb_W_x, ar_wwb_W_y, ar_wwb_W_l = [], [], []
    ar_cdC_x, ar_cdC_y, ar_cdC_l       = [], [], []
    ar_cd_x, ar_cd_y, ar_cd_l          = [], [], []
    ar_cl_x, ar_cl_y, ar_cl_l          = [], [], []
    
    
    for u in xrange(nb_wwb):
        for c in ar_bcw[u]:
            x_t, y_t = ut.show_cell_cords(ar_cell_label, c, ar_coorx, ar_coory)
            ar_bcw_x.append(x_t[0])
            ar_bcw_y.append(y_t[0])
            bcw_wwb_l.append(ar_wwb_u[u])
#            pt6.x, pt6.y = x_t, y_t
#            ptGeoms.append(arcpy.PointGeometry(pt))
        for cc in ar_bwc[u]:
            x_t, y_t = ut.show_cell_cords(ar_cell_label, cc, ar_coorx, ar_coory)
            ar_bwc_x.append(x_t[0])
            ar_bwc_y.append(y_t[0])
            bwc_wwb_l.append(ar_wwb_u[u])
        for ccc in wwb_idx[u]:
            x_t, y_t = ut.show_cell_cords(ar_cell_label, ccc, ar_coorx, ar_coory)
            ar_wwb_x.append(x_t[0])
            ar_wwb_y.append(y_t[0])
            ar_wwb_l.append(ar_wwb_u[u])
    
    for cccc in ar_cell_labl_C:
        x_t, y_t = ut.show_cell_cords(ar_cell_label, cccc, ar_coorx, ar_coory)
        ar_wwb_C_x.append(x_t[0])
        ar_wwb_C_y.append(y_t[0])
        ar_wwb_C_l.append(cccc)
    
    for ccccc in ar_cell_labl_W:
        x_t, y_t = ut.show_cell_cords(ar_cell_label, ccccc, ar_coorx, ar_coory)
        ar_wwb_W_x.append(x_t[0])
        ar_wwb_W_y.append(y_t[0])
        ar_wwb_W_l.append(ccccc)
    
    for cccccc in ar_cell_down_C:
        label_tmp = ar_cell_label[np.where(ar_cell_down==cccccc)[0][0]]
        x_t, y_t = ut.show_cell_cords(ar_cell_label, label_tmp, ar_coorx, ar_coory)
        ar_cdC_x.append(x_t[0])
        ar_cdC_y.append(y_t[0])
        ar_cdC_l.append(cccccc)
        
    for ccccccc in ar_cell_down:
        label_tmp = ar_cell_label[np.where(ar_cell_down==ccccccc)[0][0]]
        x_t, y_t = ut.show_cell_cords(ar_cell_label, label_tmp, ar_coorx, ar_coory)
        ar_cd_x.append(x_t[0])
        ar_cd_y.append(y_t[0])
        ar_cd_l.append(ccccccc)
    
    for cccccccc in ar_cell_label:
        x_t, y_t = ut.show_cell_cords(ar_cell_label, cccccccc, ar_coorx, ar_coory)
        ar_cl_x.append(x_t[0])
        ar_cl_y.append(y_t[0])
        ar_cl_l.append(cccccccc)
        
    #print ar_bcw_x, ar_bcw_y, ar_bwc_x, ar_bwc_y
    
    #ar_BCW = ar_cell_down_C[ar_cell_down_C in ]
    
#    for cd in ar_cell_down_C:
#        if cd in ar_wwb

#==============================================================================
# CREATE POINT SHAPE FILES FROM THE FILES
#==============================================================================

#    pt6 = arcpy.Point()
#    ptGeoms = []
#    for p in ptList:
#        pt.x = p[0]
#        pt.Y = p[1]
#        ptGeoms.append(arcpy.PointGeometry(pt))
#    arcpy.CopyFeatures_management(ptGeoms, r"C:\Temp\test.shp")

#==============================================================================
#     WATER BALANCE SECTION
#==============================================================================

    ar_Qs_m3  = np.average(ndar_Qs_out[1:,])*Dt       # m3/cell/second -> m3/cell/TSTP
    ar_Qs_m3_ = np.average(ndar_Qs_out_[1:,])*Dt       # m3/cell/second -> m3/cell/TSTP
    ar_Qo_m3  = np.average(ndar_Qo_out[1:,])*Dt       # m3
    ar_Qs_mm  = (ar_Qs_m3/X**2)*1e3              # m3 -> mm
    ar_Qo_mm  = (ar_Qo_m3/X**2)*1e3              # m3 -> mm
    print 'ar_Qs_mm, ar_Qo_mm', ar_Qs_mm, ar_Qo_mm

    A_catch    = X**2 * nb_cells                # m2
    #print 'A_catch: ', A_catch
    ar_Qsim_m3 = ar_Qsim *Dt                    # m3
    #print 'Dt: ', Dt
    #print 'sum ar_Qsim_m3: ', np.sum(ar_Qsim_m3)
    ar_Qsim_mm = (ar_Qsim_m3/A_catch)*1e3       # mm/s
    Qsim_mm    = np.sum(ar_Qsim_mm)
    
    #print 'Qsim_mm: ', Qsim_mm
    P_mm       = np.sum(ar_rain)
    ETa_mm     = np.sum(ar_ET_a)
#    print 'ar_Echn', ar_Echn
    Ech_m3     = np.sum(ar_Echn)          # mm/TMSP
    print 'Ech_m3', Ech_m3
    Ech_mm     = (Ech_m3/A_catch)*1e3     # m3 -> mm
    print 'Ech_mm', Ech_mm
    
    Vs_0       = ar_Vs[0]
    Vs_1       = ar_Vs[-1]
#    print 'len(ar_Vs)', len(ar_Vs)
#    print 'type(ar_Vs)', type(ar_Vs)
#    print 'type(Vs_0)', type(Vs_0)
#    print 'ar_Vs[-1]', ar_Vs[-1]
#    print 'ar_Vs[n_rec-1]', ar_Vs[n_rec-1]
#    print 'ar_Vs[len(ar_Vs)-1]', ar_Vs[len(ar_Vs)-1]
    print 'Area of the catchment is', A_catch, 'm2'
    dVs_m3     = Vs_1 - Vs_0        # for the entrire catchment
    dVs_mm  = (dVs_m3/A_catch)*1e3        # m3 -> mm
#    print 'nb_cells, Vs_0, Vs_1, dVs_m3, dVs_mm'
#    print nb_cells, Vs_0, Vs_1, dVs_m3, dVs_mm
#    for i in xrange(n_rec):
#        print ar_Vs[i]        
    
#==============================================================================
#   WRITING RESULTS TO XLS FILE
#==============================================================================

    wb          = xxw.Workbook(fn_Qout)

    # Add an Excel date format.
    date_format = wb.add_format({'num_format': 'yyyy/mm/dd'})
    st_sc_frmat = wb.add_format({'num_format': '0.00E+00'})
    text_wrappr = wb.add_format()
    text_wrappr.set_text_wrap()
    bold_format = wb.add_format()
    bold_format.set_bold()

    ws_0  = wb.add_worksheet('META')
    ws_1  = wb.add_worksheet('P_Qsim')
    ws_2  = wb.add_worksheet('Qdur_m')
    ws_3  = wb.add_worksheet('WWB_d')
    ws_4  = wb.add_worksheet('WWB_m')
    ws_5  = wb.add_worksheet('WWB_y')
    ws_51 = wb.add_worksheet('WWB')
    ws_52 = wb.add_worksheet('WWB_mm')
    ws_55 = wb.add_worksheet('WWB_d_mm')
    ws_56 = wb.add_worksheet('WWB_m_mm')
    ws_57 = wb.add_worksheet('WWB_y_mm')
    ws_6  = wb.add_worksheet('cBWin')
    ws_7  = wb.add_worksheet('cBWout')
    ws_8  = wb.add_worksheet('cWWB')
    ws_9  = wb.add_worksheet('cTRD')
    #ws_10 = wb.add_worksheet('cWWB_Ws') # Dublicate.. not needed.
    #ws_11 = wb.add_worksheet('cd_C')  # Cell down catchment
    #ws_12 = wb.add_worksheet('cd')    # Cell down all
    ws_13 = wb.add_worksheet('cCL')    # Cell labels
    if Piezos_on:
        ws_15 = wb.add_worksheet('PZM')
    if Dams_on:
        ws_14 = wb.add_worksheet('DAMs')
        ws_14.set_column(1,1,13)
        # Writing dam data in sheet DAMs
        ws_h = ['FSL', 'Vol_d', 'Qd_in', 'Qd_in', 'Qd_out', 'Qd_out']
        ws_u = ['(m3)', '(m3)', '(m3/s)', '(m3/TSTP)', '(m3/s)', '(m3/TSTP)']
        nb_var = len(ws_h)
        ws_14.write(3,1,'Date')

        Qd_in_sum = 0    # in m3 integrated over model period
        Qd_ot_sum = 0
        for h in xrange(nb_dam):
            Qd_in_sum += np.sum(Qd_input[h]) * Dt
            Qd_ot_sum += np.sum(Qd_outpt[h]) * Dt
            FSL = WL_dam_avg * ar_dam_u[h]
            ws_14.write(2, nb_var*h+1+0, FSL)
            for i in xrange(nb_var):
                ws_14.write(0, nb_var*h+i+1, ws_h[i]+'_'+str(ar_dam_u[h]))
                ws_14.write(1, nb_var*h+i+1, ws_u[i])
            # write daily data into sheet 3
            for l in xrange(n_rec):
                if h == 0:
                    ws_14.write(l+2, 0, l+1)
                try:
                    ws_14.write(l+2, nb_var*h+1+1, float(Vd_cel_d[h,l]))
                except TypeError:
                    print h, l
                    print Vd_cel_d[h,l]
                    print type(Vd_cel_d[h,l])
                    raise ValueError('here you go.. that is wrong')
                ws_14.write(l+2, nb_var*h+1+2, float(Qd_input[h,l]))
                ws_14.write(l+2, nb_var*h+1+3, float(Qd_input[h,l])*Dt)
                ws_14.write(l+2, nb_var*h+1+4, float(Qd_outpt[h,l]))
                ws_14.write(l+2, nb_var*h+1+5, float(Qd_outpt[h,l])*Dt)
        
        # DAM losses:
        D_dam_loss_mm = ((Qd_in_sum-Qd_ot_sum)/A_catch)*1e3

#==============================================================================
# Extract volumes in the cells of the piezometers and obs pzm data from file
#==============================================================================
#        plt.rc('text', usetex=True, fontsize=10)
#        plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'], 'size':'10'})



    if Piezos_on:
        print 'piezo cell label(s):',piezos
        piezos_Vs = ndar_Vs_cel[1:,piezos]
        piezosSSI = ndar_SSI_cl[:,piezos]
        ws_15.write(0, 0, 'TSTP')
        ws_15.write(0, 1, 'Date')
        
        
        print 'piezos_Vs.shape',piezos_Vs.shape
        print 'piezosSSI.shape',piezosSSI.shape
        
        #include SSI
        #ndar_SSI_cl
        # x = column
        # y = row
        for x in xrange(nb_piezos):
            ws_15.write(0, nb_piezos*x+2, 'Piezo_'+str(piezos[x]))
            ws_15.write(0, nb_piezos*x+3, 'Piezo_'+str(piezos[x]))
            ws_15.write(1, nb_piezos*x+2, 'VS (m3)')
            ws_15.write(1, nb_piezos*x+3, 'SSI (-)')
            for y in xrange(piezos_Vs.shape[0]):
                ws_15.write(y+2, nb_piezos*x+2, piezos_Vs[y,x])
                ws_15.write(y+2, nb_piezos*x+3, piezosSSI[y,x])
                if x < 1:
                    ws_15.write(y+2, 0, y+1)

        if pzo_plot_on:
#            pzo_hdr
#            pzo_data
            
            date_s   = [int(s) for s in start_stamp.split('/')]
            date_s   = dtm.datetime(date_s[2],date_s[1],date_s[0])
            dat      = date_s
            date_sim = []
            date_mon = []

            for i in xrange(n_rec):
                dat = dat + dtm.timedelta(seconds=Dt)
                date_sim.append(dat)
            
            ar_date_mon = pzo_data[:,0:3]
            
            for t in xrange(pzo_data.shape[0]):
#                print ar_date_mon[t,0], ar_date_mon[t,1], ar_date_mon[t,2]
                date_mon.append(dtm.datetime(int(ar_date_mon[t,0]),\
                                             int(ar_date_mon[t,1]),\
                                             int(ar_date_mon[t,2])))

            plt.clf()
            fig     = plt.figure()
            ax1     = fig.add_subplot(111)
#            ax2     = ax1.twinx()
            lines   = []
            t_leg   = []
            col_sim = ['orange','red','m','yellow']
            col_mon = ['black','blue','gray','brown']

            for i in xrange(nb_piezos):
                i_sim    = piezos.index(piezos[i])
                i_mon    = pzo_hdr.index(piezos[i])+3
                sim_pzo_SSI = piezosSSI[:,i_sim]
                mon_pzo_SSI =  pzo_data[:,i_mon]

                lines += ax1.plot(date_sim,sim_pzo_SSI,col_sim[i],linewidth=2, linestyle='dotted')
                lines += ax1.plot(date_mon,mon_pzo_SSI,col_mon[i],linewidth=1)
                
                t_leg.append(r'$Simulated~SSI~%s$' %str(piezos[i]))
                t_leg.append(r'$Monitored~SSI~%s$' %str(piezos[i]))

#            for v in [0.75]:
#                ax2.axhline(v, color='k', linestyle='dotted')
            ax1.set_ylim((0,120))
            ax1.legend(lines, t_leg, loc='best', fancybox=True)
            leg = ax1.get_legend()
            leg.get_frame().set_alpha(0.75)
#            plt.show()
            
            fig_fn = 'Results\\PZM_vs_SSI_V_'+str(ver_num)+'.png'
            fig.savefig(fig_fn, dpi=200)


#==============================================================================
# WRITING THE META SHEET
#==============================================================================

    # first make a list of all variables which need explanation
    
    # Tabs
    tabs = ['META','P_Qsim','Qdur_m','WWB_d','WWB_m','WWB_y','WWB','WWB_mm',\
            'WWB_d_mm','WWB_m_mm','WWB_y_mm','cBWin','cBWout','cWWB','cTRD',\
            'cCL','DAMs']
    txts = ['Information about this spreadsheet.','Precipitation (P) and Runo\
ff (Q) simulation data and mass  balance for catchment: Time series data of P\
, Q and volumetric soil moisture and catchment water balance components, mass\
 balance error for the modelling period','Monthly flow duration (Qdur_m) data\
 table for flow at catchment outlet','Wetland Water Balance (WWB) components \
per WWB-unit (daily, volumetric)','Wetland Water Balance (WWB) components per\
 WWB-unit (monthly, volumetric)','Wetland Water Balance (WWB) components per \
WWB-unit (yearly, volumetric)','Wetland Water Balance (WWB) components per WW\
B-unit (yearly, volumetric) for the entire modelling period','Wetland Water B\
alance (WWB) components per WWB-unit (yearly, mm) for the entire modelling pe\
riod','Wetland Water Balance (WWB) components per WWB-unit (daily, mm)','Wetl\
and Water Balance (WWB) components per WWB-unit (monthly, mm)','Wetland Water\
 Balance (WWB) components per WWB-unit (yearly, mm)','Coordinates for the wet\
land boundary of the inflowing edge and WWB-unit number','Coordinates for the\
 wetland boundary of the outflowing edge and WWB-unit number','Coordinates fo\
r each cell of WWB-units including unit number','Coordinates for each cell of\
 Terestrial, River or Dam (TRD) including cell label number','Coordinates of \
Cell Label (cCL) numbers including cell labels','Time series data for dam inp\
ut, output, volume and full supply level for all dams']



    ws_0.set_column(0,0,13)
    ws_0.set_column(1,1,100)
    ws_0.write(0, 0, 'Worksheet', bold_format)
    ws_0.write(0, 1, 'Content Description', bold_format)
    for r in xrange(len(tabs)):
        ws_0.write(r+1, 0, tabs[r])
        ws_0.write(r+1, 1, txts[r], text_wrappr)

#    for r2 in xrange(len()):
#        ws_0.write(r2+r+1, 0, tabs[0])



# Writing daily basic catchment data in sheet P_Qsim
    XLS_h = ['TimeStep','Date','P(mm/D)','Q(m3/s)','V_soil(m3)','SSI',\
             'ETa(mm/D)','','Date','P(mm/M)','Q(m/M)', 'Date','P(mm/Y)',\
             'Q(mm/Y)','P(mm)','Q(mm)','ETa(mm)', 'Ech(mm)','dVsoil(mm)',\
             'DAM_loss(mm)','Error(mm)','Qs(mm)','Qo(mm)']

    # Write a header line into xls worksheet
    for h in xrange(len(XLS_h)):
        ws_1.write(0, h, XLS_h[h])

    for r in xrange(n_rec):
        ws_1.write(r+1, 0, r+1)
        ws_1.write(r+1, 2, float(ar_rain[r]))    # mm/TMSP
        ws_1.write(r+1, 3, float(ar_Qsim[r]))    # m3/s
        ws_1.write(r+1, 4, float(ar_Vs[r]))      # m3 for entire catchment
        ws_1.write(r+1, 5, float(ar_SSI_s[r]))
        ws_1.write(r+1, 6, float(ar_ET_a[r]))
    
    #np.savetxt('test_Vs', ar_Vs, fmt='%.4f')
    
    c_o = 7    # Column offsetter +1

    #===========================================================================
    # Now make monthly averages from the flow
    #===========================================================================

    # [dd mm yyyy]

    #start_stamp = '22/8/2013'
    
    split_mon = int(split_month)
    
    stmp = []
    for s in start_stamp.split('/'):
        stmp.append(int(s))
        #print 'stmp', stmp
        
    start_date = dtm.datetime(stmp[2],stmp[1],stmp[0])
    
    print 'start date:', start_date.strftime("%d %b %Y")
    
    # Create datetime object for each simulated day
    tme = start_date
    D_day = []
    D_mon = []
    P_mon = []
    Q_mon = []
    P_mc  = []
    Q_mc  = []
    P_yc  = []
    Q_yc  = []

    Pt_mc_w     = np.zeros(nb_wwb)
    Pt_yc_w     = np.zeros(nb_wwb)
    Ea_mc_w     = np.zeros(nb_wwb)
    Ea_yc_w     = np.zeros(nb_wwb)
#    Vs_mc_w    = np.zeros(nb_wwb)
#    Vs_yc_w    = np.zeros(nb_wwb)
    Qs_in_mc_w  = np.zeros(nb_wwb)
    Qs_in_mc_w_ = np.zeros(nb_wwb)
    Qs_ot_mc_w  = np.zeros(nb_wwb)
    Qs_ot_mc_w_ = np.zeros(nb_wwb)
    Qs_in_yc_w  = np.zeros(nb_wwb)
    Qs_in_yc_w_ = np.zeros(nb_wwb)
    Qs_ot_yc_w  = np.zeros(nb_wwb)
    Qs_ot_yc_w_ = np.zeros(nb_wwb)
    Qo_in_mc_w  = np.zeros(nb_wwb)
    Qo_ot_mc_w  = np.zeros(nb_wwb)
    Qo_in_yc_w  = np.zeros(nb_wwb)
    Qo_ot_yc_w  = np.zeros(nb_wwb)
    Qs_in_mc_r  = np.zeros(nb_wwb)
    Qs_in_mc_r_ = np.zeros(nb_wwb)
    Qs_in_yc_r  = np.zeros(nb_wwb)
    Qs_in_yc_r_ = np.zeros(nb_wwb)
    Qo_in_mc_r  = np.zeros(nb_wwb)
    Qo_in_yc_r  = np.zeros(nb_wwb)
    SSInx_mc_w  = np.zeros(nb_wwb)
    SSInx_mc_w_ = np.zeros(nb_wwb)
    SSInx_yc_w  = np.zeros(nb_wwb)
    SSInx_yc_w_ = np.zeros(nb_wwb)

    D_yer  = []
    P_yer  = []
    Q_yer  = []

    M_mon = []    # Array of the months only i
    Y_mon = []    # Array of all years in monthly resolution
    #check how many months there are going to be..
    n_mon = 0
    n_yer = 0    
    cnt_y = 0
    cnt_d = 0
    cnt_m = 0
    mon_t = 0
    years = []
    
    # Take time increment (time step) in seconds from global parameter file (Dt)
    
    for i in xrange(n_rec):
        cnt_d += 1
        tme = tme + dtm.timedelta(seconds=Dt)
        mon = tme.month
        if not i:
            mon_p = mon
        if mon_p != mon and cnt_d >= 28:
            n_mon += 1
            cnt_m += 1
#            print 'cnt_m', cnt_m
            mon_p = mon
            cnt_d = 0
            if not mon_t:
                y_first = tme.year
            mon_t += 1
            if mon == split_mon and cnt_m >= 12:
                years.append(tme.year)
                cnt_y += 1
                n_yer += 1
                cnt_m  = 0
    #print 'n_yer', n_yer
    #print 'y_first', y_first
    Q_mon2D = np.zeros((n_yer, 12))*-9999.
    yrs = []
    
    #print np.min(years), np.max(years)+1
    
    if len(years) > 0:
        for y in xrange(np.min(years), np.max(years)+1):
            yrs.append(y)
        
    
    yrs = np.array(yrs)
    
    # 3D Conteiner for month-averages (12 x WWB x vars)
    # for temporary and permanent
    
    nb_vars = 9
    
    M_WWB_V = np.ones((12,nb_wwb,nb_vars))
    M_count = np.zeros((12))
    
    
    Q_mon_p = np.zeros((12, 10))*-9999.

    cnt_d_m   = 0
    cnt_d_y   = 0
    cnt_m     = 0
    cnt_m_all = 0
    cnt_y_all = 0
    tme       = start_date
    t_e       = tme + dtm.timedelta(seconds=Dt*n_rec)    # End time
    
    print "...last time stamp:", t_e.strftime("%d %b %Y")
    
    t_e_mon   = t_e.month
    t_e_day   = t_e.day
    t_e_yer   = t_e.year
    t_e_LDOF  = ut.last_day_of_month(t_e)    #t_e_LDOF = end time last day of month
    
    if t_e_day < t_e_LDOF:
        exclude_last_mon = 1
    else:
        exclude_last_mon = 0
    
    avg_m_start = 0
    
    Pt_m_w    = np.zeros((nb_wwb, mon_t), dtype=ftype)
    Pt_y_w    = np.zeros((nb_wwb, cnt_y), dtype=ftype)
    Ea_m_w    = np.zeros((nb_wwb, mon_t), dtype=ftype)
    Ea_y_w    = np.zeros((nb_wwb, cnt_y), dtype=ftype)
    Vs_m_w    = np.zeros((nb_wwb, mon_t), dtype=ftype)
    Vs_m_w_   = np.zeros((nb_wwb, mon_t), dtype=ftype)
    Vs_y_w    = np.zeros((nb_wwb, cnt_y), dtype=ftype)
    Vs_y_w_   = np.zeros((nb_wwb, cnt_y), dtype=ftype)
    
    Qs_in_m_w  = np.zeros((nb_wwb, mon_t), dtype=ftype)
    Qs_in_m_w_ = np.zeros((nb_wwb, mon_t), dtype=ftype)
    Qs_ot_m_w  = np.zeros((nb_wwb, mon_t), dtype=ftype)
    Qs_ot_m_w_ = np.zeros((nb_wwb, mon_t), dtype=ftype)
    Qs_in_y_w  = np.zeros((nb_wwb, cnt_y), dtype=ftype)
    Qs_in_y_w_ = np.zeros((nb_wwb, cnt_y), dtype=ftype)
    Qs_ot_y_w  = np.zeros((nb_wwb, cnt_y), dtype=ftype)
    Qs_ot_y_w_ = np.zeros((nb_wwb, cnt_y), dtype=ftype)
    Qs_in_m_r  = np.zeros((nb_wwb, mon_t), dtype=ftype)
    Qs_in_m_r_ = np.zeros((nb_wwb, mon_t), dtype=ftype)
    Qs_in_y_r  = np.zeros((nb_wwb, cnt_y), dtype=ftype)
    Qs_in_y_r_ = np.zeros((nb_wwb, cnt_y), dtype=ftype)
    Qo_in_m_w  = np.zeros((nb_wwb, mon_t), dtype=ftype)
    Qo_ot_m_w  = np.zeros((nb_wwb, mon_t), dtype=ftype)
    Qo_in_y_w  = np.zeros((nb_wwb, cnt_y), dtype=ftype)
    Qo_ot_y_w  = np.zeros((nb_wwb, cnt_y), dtype=ftype)
    Qo_in_m_r  = np.zeros((nb_wwb, mon_t), dtype=ftype)
    Qo_in_y_r  = np.zeros((nb_wwb, cnt_y), dtype=ftype)
    SSInx_m_w  = np.zeros((nb_wwb, mon_t), dtype=ftype)
    SSInx_m_w_ = np.zeros((nb_wwb, mon_t), dtype=ftype)
    SSInx_y_w  = np.zeros((nb_wwb, cnt_y), dtype=ftype)
    SSInx_y_w_ = np.zeros((nb_wwb, cnt_y), dtype=ftype)
    
    #print 'nb_wwb:', nb_wwb 
    #print 'mon_t', mon_t
    
    for i in xrange(n_rec):
        cnt_d_m += 1
        cnt_d_y += 1
        tme = tme + dtm.timedelta(seconds=Dt)
        D_day.append(tme)
        day = tme.day
        mon = tme.month
        yer = tme.year
        if not i:
            mon_p = mon
            yer_p = yer
        # only start putting values into month-averages if day 1 was passed
            
#            make another condition here to not include the last values unless the month will be complete.. use end time t_e for that.
        if day == 1:# or nb_days_left > days_last_month...must be simpler:
            avg_m_start = 1
        if exclude_last_mon and mon == t_e_mon and yer == t_e_yer:
#            print 'these days get excluded', tme.strftime('%d %b %Y')
            avg_m_start = 0

        P_mc.append(float(ar_rain[i]))
        Q_mc.append(float(ar_Qsim[i]))
        for w in xrange(nb_wwb):
            if avg_m_start:
#                if i < 40 and w == 0:
#                    print tme.strftime('%d %b %Y')
                M_WWB_V[mon-1,w, 0] += Rainfl_w[w,i]
                M_WWB_V[mon-1,w, 1] += ET_out_w[w,i]
                M_WWB_V[mon-1,w, 2] += Qs_in_w[w,i]
                M_WWB_V[mon-1,w, 3] += Qo_in_w[w,i]
                M_WWB_V[mon-1,w, 4] += Qs_ot_w[w,i]
                M_WWB_V[mon-1,w, 5] += Qo_ot_w[w,i]
                M_WWB_V[mon-1,w, 6] += Qs_in_r[w,i]
                M_WWB_V[mon-1,w, 7] += Qo_in_r[w,i]
                M_WWB_V[mon-1,w, 8] += Qs_in_w_[w,i]
                M_WWB_V[mon-1,w, 9] += Qs_ot_w_[w,i]
                M_WWB_V[mon-1,w,10] += Qs_in_r_[w,i]
#                M_WWB_V[mon-1,w,8] += SSIndx_w[w,i] Need to find a way to average!
            Pt_mc_w[w]    += Rainfl_w[w,i]
            Pt_yc_w[w]    += Rainfl_w[w,i]
            Ea_mc_w[w]    += ET_out_w[w,i]      # units in m3 method for monthly avg
            Ea_yc_w[w]    += ET_out_w[w,i]      # units in m3 method for yearly  avg
#            Vs_mc_w[w]     = Vs_cel_w[w,i]     # units in m3 method for monthly avg
#            Vs_yc_w[w]     = Vs_cel_w[w,i]     # units in m3
            Qs_in_mc_w[w]  += Qs_in_w[w,i]      # units in m3
            Qs_in_mc_w_[w] += Qs_in_w_[w,i]     # units in m3
            Qs_in_yc_w[w]  += Qs_in_w[w,i]      # units in m3
            Qs_in_yc_w_[w] += Qs_in_w_[w,i]     # units in m3
            Qo_in_mc_w[w]  += Qo_in_w[w,i]      # units in m3
            Qo_in_yc_w[w]  += Qo_in_w[w,i]      # units in m3
            Qs_ot_mc_w[w]  += Qs_ot_w[w,i]      # units in m3
            Qs_ot_mc_w_[w] += Qs_ot_w_[w,i]     # units in m3
            Qs_ot_yc_w[w]  += Qs_ot_w[w,i]      # units in m3
            Qs_ot_yc_w_[w] += Qs_ot_w_[w,i]     # units in m3
            Qo_ot_mc_w[w]  += Qo_ot_w[w,i]      # units in m3
            Qo_ot_yc_w[w]  += Qo_ot_w[w,i]      # units in m3
            Qs_in_mc_r[w]  += Qs_in_r[w,i]      # units in m3
            Qs_in_mc_r_[w] += Qs_in_r_[w,i]     # units in m3
            Qs_in_yc_r[w]  += Qs_in_r[w,i]      # units in m3
            Qs_in_yc_r_[w] += Qs_in_r_[w,i]     # units in m3
            Qo_in_mc_r[w]  += Qo_in_r[w,i]      # units in m3
            Qo_in_yc_r[w]  += Qo_in_r[w,i]      # units in m3
            SSInx_mc_w[w]  += SSIndx_w[w,i]     # Dimensionless
            SSInx_mc_w_[w] += SSIndx_w_[w,i]    # Dimensionless
            SSInx_yc_w[w]  += SSIndx_w[w,i]     # Dimensionless
            SSInx_yc_w_[w] += SSIndx_w_[w,i]    # Dimensionless
        P_yc.append(float(ar_rain[i]))
        Q_yc.append(float(ar_Qsim[i]))
        if mon_p != mon and cnt_d_m >= 28:
            if avg_m_start:
                M_count[mon-1] += 1
            cnt_m += 1
            D_mon.append(tme)
            P_mon.append(sum(P_mc))
            Q_mon.append(np.average(Q_mc))
            for w in xrange(nb_wwb):
#                Pt_m_w[w,cnt_m_all]    = Pt_mc_w[w]/cnt_d_m     # monthly avg
#                Ea_m_w[w,cnt_m_all]    = Ea_mc_w[w]/cnt_d_m     # monthly avg
#                Vs_m_w[w,cnt_m_all]    = Vs_cel_w[w,i]          # current vol
#                Qs_in_m_w[w,cnt_m_all] = Qs_in_mc_w[w]/cnt_d_m  # monthly avg
#                Qs_ot_m_w[w,cnt_m_all] = Qs_ot_mc_w[w]/cnt_d_m  # monthly avg
#                Qo_in_m_w[w,cnt_m_all] = Qo_in_mc_w[w]/cnt_d_m  # monthly avg
#                Qo_ot_m_w[w,cnt_m_all] = Qo_ot_mc_w[w]/cnt_d_m  # monthly avg
#                Qs_in_m_r[w,cnt_m_all] = Qs_in_mc_r[w]/cnt_d_m  # monthly avg
#                Qo_in_m_r[w,cnt_m_all] = Qo_in_mc_r[w]/cnt_d_m  # monthly avg
#                SSInx_m_w[w,cnt_m_all] = SSInx_mc_w[w]/cnt_d_m  # monthly avg
                Pt_m_w[w,cnt_m_all]     = Pt_mc_w[w]       # monthly sum
                Ea_m_w[w,cnt_m_all]     = Ea_mc_w[w]       # monthly sum
                Vs_m_w[w,cnt_m_all]     = Vs_cel_w[w,i]    # current vol
                Vs_m_w_[w,cnt_m_all]    = Vs_cel_w_[w,i]    # current vol
                Qs_in_m_w[w,cnt_m_all]  = Qs_in_mc_w[w]    # monthly sum
                Qs_in_m_w_[w,cnt_m_all] = Qs_in_mc_w_[w]    # monthly sum
                Qs_ot_m_w[w,cnt_m_all]  = Qs_ot_mc_w[w]    # monthly sum
                Qs_ot_m_w_[w,cnt_m_all] = Qs_ot_mc_w_[w]    # monthly sum
                Qo_in_m_w[w,cnt_m_all]  = Qo_in_mc_w[w]    # monthly sum
                Qo_ot_m_w[w,cnt_m_all]  = Qo_ot_mc_w[w]    # monthly sum
                Qs_in_m_r[w,cnt_m_all]  = Qs_in_mc_r[w]    # monthly sum
                Qs_in_m_r_[w,cnt_m_all] = Qs_in_mc_r_[w]    # monthly sum
                Qo_in_m_r[w,cnt_m_all]  = Qo_in_mc_r[w]    # monthly sum
                SSInx_m_w[w,cnt_m_all]  = SSInx_mc_w[w]/cnt_d_m    # m-AVG dimensionless
                SSInx_m_w_[w,cnt_m_all] = SSInx_mc_w_[w]/cnt_d_m    # m-AVG dimensionless
            P_mc        = []
            Q_mc        = []
            Pt_mc_w     = np.zeros(nb_wwb)
            Ea_mc_w     = np.zeros(nb_wwb)
            Vs_mc_w     = np.zeros(nb_wwb)
            Qs_in_mc_w  = np.zeros(nb_wwb)
            Qs_in_mc_w_ = np.zeros(nb_wwb)
            Qs_ot_mc_w  = np.zeros(nb_wwb)
            Qs_ot_mc_w_ = np.zeros(nb_wwb)
            Qo_in_mc_w  = np.zeros(nb_wwb)
            Qo_ot_mc_w  = np.zeros(nb_wwb)
            Qs_in_mc_r  = np.zeros(nb_wwb)
            Qs_in_mc_r_ = np.zeros(nb_wwb)
            Qo_in_mc_r  = np.zeros(nb_wwb)
            SSInx_mc_w  = np.zeros(nb_wwb)
            SSInx_mc_w_ = np.zeros(nb_wwb)
            M_mon.append(mon)
            Y_mon.append(yer)
            #print tme.strftime("%d %b %Y")
            cnt_d_m = 0
            if mon == split_mon and cnt_m >= 12:
                cnt_m = 0
                D_yer.append(tme)
                P_yer.append(sum(P_yc))
                Q_yer.append(np.average(Q_yc))
                for w in xrange(nb_wwb):
#                    Pt_y_w[w,cnt_y_all]    = Pt_yc_w[w]/cnt_d_y     # yearly avg
#                    Ea_y_w[w,cnt_y_all]    = Ea_yc_w[w]/cnt_d_y     # yearly avg
#                    Vs_y_w[w,cnt_y_all]    = Vs_cel_w[w,i]          # current vol
#                    Qs_in_y_w[w,cnt_y_all] = Qs_in_yc_w[w]/cnt_d_y  # yearly avg
#                    Qs_ot_y_w[w,cnt_y_all] = Qs_ot_yc_w[w]/cnt_d_y  # yearly avg
#                    Qo_in_y_w[w,cnt_y_all] = Qo_in_yc_w[w]/cnt_d_y  # yearly avg
#                    Qo_ot_y_w[w,cnt_y_all] = Qo_ot_yc_w[w]/cnt_d_y  # yearly avg
#                    Qs_in_y_r[w,cnt_y_all] = Qs_in_yc_r[w]/cnt_d_y  # yearly avg
#                    Qo_in_y_r[w,cnt_y_all] = Qo_in_yc_r[w]/cnt_d_y  # yearly avg
#                    SSInx_y_w[w,cnt_y_all] = SSInx_yc_w[w]/cnt_d_y  # Dimensionless
                    Pt_y_w[w,cnt_y_all]     = Pt_yc_w[w]       # yearly sum
                    Ea_y_w[w,cnt_y_all]     = Ea_yc_w[w]       # yearly sum
                    Vs_y_w[w,cnt_y_all]     = Vs_cel_w[w,i]    # current vol
                    Vs_y_w_[w,cnt_y_all]    = Vs_cel_w_[w,i]    # current vol
                    Qs_in_y_w[w,cnt_y_all]  = Qs_in_yc_w[w]    # yearly sum
                    Qs_in_y_w_[w,cnt_y_all] = Qs_in_yc_w_[w]    # yearly sum
                    Qs_ot_y_w[w,cnt_y_all]  = Qs_ot_yc_w[w]    # yearly sum
                    Qs_ot_y_w_[w,cnt_y_all] = Qs_ot_yc_w_[w]    # yearly sum
                    Qo_in_y_w[w,cnt_y_all]  = Qo_in_yc_w[w]    # yearly sum
                    Qo_ot_y_w[w,cnt_y_all]  = Qo_ot_yc_w[w]    # yearly sum
                    Qs_in_y_r[w,cnt_y_all]  = Qs_in_yc_r[w]    # yearly sum
                    Qs_in_y_r_[w,cnt_y_all] = Qs_in_yc_r_[w]    # yearly sum
                    Qo_in_y_r[w,cnt_y_all]  = Qo_in_yc_r[w]    # yearly sum
                    SSInx_y_w[w,cnt_y_all]  = SSInx_yc_w[w]/cnt_d_y  # y-AVG dimensionless
                    SSInx_y_w_[w,cnt_y_all] = SSInx_yc_w_[w]/cnt_d_y  # y-AVG dimensionless
                P_yc  = []
                Q_yc  = []
                Pt_yc_w    = np.zeros(nb_wwb)
                Ea_yc_w    = np.zeros(nb_wwb)
                Vs_yc_w    = np.zeros(nb_wwb)
                Qs_in_yc_w = np.zeros(nb_wwb)
                Qs_ot_yc_w = np.zeros(nb_wwb)
                Qo_in_yc_w = np.zeros(nb_wwb)
                Qo_ot_yc_w = np.zeros(nb_wwb)
                Qs_in_yc_r = np.zeros(nb_wwb)
                Qo_in_yc_r = np.zeros(nb_wwb)
                SSInx_yc_w = np.zeros(nb_wwb)
                cnt_y_all += 1
                cnt_d_y    = 0
            cnt_m_all += 1
            mon_p      = mon
            yer_p      = yer

    M_mon = np.array(M_mon)
    Y_mon = np.array(Y_mon)


    for m in xrange(12):
        for w in xrange(nb_wwb):
            for v in xrange(nb_vars):
                M_WWB_V[m,w,v] = M_WWB_V[m,w,v]/M_count[m]

#    print M_count

    #print 'max year Y_mon:', np.max(Y_mon)
    #print 'max year yers:', np.max(yrs)
    
    for m in xrange(n_mon):
        if Y_mon[m] in yrs:
            try:
                Q_mon2D[np.where(yrs==Y_mon[m])[0][0],M_mon[m]-1] = Q_mon[m]
            except IndexError:
                print Q_mon2D.shape
                print np.where(yrs==Y_mon[m])
                print M_mon[m]
                print Y_mon[m]
                raise ValueError()
    if P_yer:
        P_MA = np.average(P_yer)    #MAP = Mean annual precipitation
        if np.isnan(P_MA):
            P_MA = 'nan'
    else:
        P_MA = 'nan'
    if Q_yer:
        Q_MA = np.average(Q_yer)    #MAQ = Mean annual runoff
        if np.isnan(Q_MA):
            Q_MA = 'nan'
    else:
        Q_MA = 'nan'

    #print 'D_day:', len(D_day)
    #print 'D_mon:', len(D_mon)
    #print 'P_mon:', len(P_mon)
    #print 'Q_mon:', len(Q_mon)
    
    ws_1.set_column(1,1,10)
    for rr in xrange(n_rec):
        ws_1.write(rr+1, 1, D_day[rr], date_format)

    ws_1.set_column(1,c_o+1,10)
    for rrr in xrange(len(D_mon)):
        ws_1.write(rrr+1, c_o+1, D_mon[rrr], date_format)
        ws_1.write(rrr+1, c_o+2, float(P_mon[rrr]))
        ws_1.write(rrr+1, c_o+3, float(Q_mon[rrr]))
#        ws_3.write(rrr+1)

    ws_1.set_column(1,c_o+4,10)
    for rrrr in xrange(len(P_yer)):
        ws_1.write(rrrr+1, c_o+4, D_yer[rrrr], date_format)
        ws_1.write(rrrr+1, c_o+5, float(P_yer[rrrr]))
        ws_1.write(rrrr+1, c_o+6, float(Q_yer[rrrr]))

    ### Write water balance components expressed in mm
    ws_1.write(1, c_o+7,  float(P_mm))
    ws_1.write(1, c_o+8,  float(Qsim_mm))
    ws_1.write(1, c_o+9,  float(ETa_mm))
    ws_1.write(1, c_o+10, float(Ech_mm))
    ws_1.write(1, c_o+11, float(dVs_mm))
    if Dams_on:
        ws_1.write(1, c_o+12, float(D_dam_loss_mm))
        ws_1.write(1, c_o+13, float(P_mm - Qsim_mm - ETa_mm - Ech_mm - dVs_mm -\
                   D_dam_loss_mm))
    else:
        ws_1.write(1, c_o+13, float(P_mm - Qsim_mm - ETa_mm - Ech_mm - dVs_mm))
    ws_1.write(1, c_o+14, float(ar_Qs_mm))
    ws_1.write(1, c_o+15, float(ar_Qo_mm))
    
    ### Write water balance components expressed in % of P
    ws_1.write(5, c_o+7,  float(P_mm/P_mm)*100.)
    ws_1.write(5, c_o+8,  float(Qsim_mm/P_mm)*100.)
    ws_1.write(5, c_o+9,  float(ETa_mm/P_mm)*100.)
    ws_1.write(5, c_o+10, float(Ech_mm/P_mm)*100.)
    ws_1.write(5, c_o+11, float(dVs_mm/P_mm)*100.)
    if Dams_on:
        ws_1.write(5, c_o+12, float((D_dam_loss_mm/P_mm))*100.)
        ws_1.write(5, c_o+13, float((P_mm - Qsim_mm - ETa_mm - Ech_mm - dVs_mm -\
                   D_dam_loss_mm)/P_mm)*100.)
    else:
        ws_1.write(5, c_o+13, float((P_mm - Qsim_mm - ETa_mm - Ech_mm - dVs_mm)\
                                     /P_mm)*100.)
    ws_1.write(5, c_o+14, float(ar_Qs_mm/P_mm)*100.)
    ws_1.write(5, c_o+15, float(ar_Qo_mm/P_mm)*100.)

    percentage_hdr = ['P','Q','ETa','Ech','dVsoil','DAM_loss','Error','Qs','Qo']
    
    for h in xrange(len(percentage_hdr)):
        ws_1.write(4, c_o+7+h,percentage_hdr[h])
    
#    for c in xrange(12):
#        mon_tmp = []
#        for r in xrange(n_months):
#            Q_mon2D[r, c].append(mon_tmp)
    #print Q_mon_p.shape
    #print Q_mon2D
    if Q_mon2D.shape[0]:
        m_names = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct',
                   'Nov','Dec']
        rows    = [r for r in xrange(1,13)]
        i = -1
        for m in [10, 11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9]:
            i += 1
            Q_mon_avg_tmp = np.average(Q_mon2D[:,m-1])
            if np.isnan(Q_mon_avg_tmp):
                Q_mon_avg_tmp = -99999.9
            ws_2.write(rows[i],0, m_names[m-1])
            #ws_3.write(rows[i],0, m_names[m-1])
            for p in xrange(1,11):
                if len(yrs) > 0:
                    if p == 10:
                        Q_mon_p[m-1, p-1] = np.percentile(Q_mon2D[:,m-1], 1)
                    else:
                        #print 'Q_mon2D.shape', Q_mon2D.shape
                        #print Q_mon2D
                        Q_mon_p[m-1, p-1] = np.percentile(Q_mon2D[:,m-1], 100-p*10)
                else:
                    Q_mon_p[m-1, p-1] = 0.0
                ws_2.write(rows[i], p, Q_mon_p[m-1, p-1], st_sc_frmat)
            ws_2.write(rows[i], 13, Q_mon_avg_tmp, st_sc_frmat)
    
    for p in xrange(1,11):
        if p == 10:
            ws_2.write(0, p, str(99)+'%')
        else:
            ws_2.write(0, p, str(p*10)+'%')
    ws_2.write(0, 13, 'Q avg')
    ws_2.write(0, 0, 'Month')
    ws_2.write(0, 14, 'MAP')
    ws_2.write(1, 14, P_MA)
    ws_2.write(0, 15, 'MAQ')
    ws_2.write(1, 15, Q_MA)


#==============================================================================
# Write WWB yearly data on harddrive
#==============================================================================

    # First create the right shape for the WWB - yearly data:
    # Each unit has 13 n's and 13 i's both for mm and m3 make 52 rows.
    # columns = max nb of units = 200:
    #ETa, Vs, Qs_in, Qs_out, Qo_in, Qo_out, Qs_in_r, Qo_in_r, Rain, WBd, SSI, N_yrs, Catchment
    #rt = run type either n or i    
    
    # Do the same for month averages which are averaged sums (m3/month) for each month
    # i.e. all Jan's, Feb's, etc...
    rt = run_name[0]
    
    print 'Run type =', rt
#    print 'n_rec =', n_rec
    
    if n_rec > 500 and (rt == 'n' or rt == 'i') and not run_name[1] == 'P':
        # Cheap fix to get some of the nWWBU into iWWBU as they are unimpacted
        # Disable this for other projects than AMF
        if WWBU_transfer_n2i == True:
            # list of n WWB's which must get into the i database
            WWB_n2i = [76,77,78,80,25]
        cols  = 200
        colsm = cols*10
        rows  = 14*2*2
        rowsm = 23
        spxy  = (rows, cols)
        spxym = (rowsm, colsm)
        
        AR_WWB_y_fn = os.getcwd()+'\\..\\..'+'\\Final_Results\\AR_WWB_y'
        AR_WWB_m_fn = os.getcwd()+'\\..\\..'+'\\Final_Results\\AR_WWB_m'
    
        if os.path.exists(AR_WWB_y_fn):
            AR_WWB_y = np.fromfile(AR_WWB_y_fn,dtype=ftype).reshape(spxy)
        else:
            # Create a 2D array which can hold all the data
            AR_WWB_y = np.ones(spxy,dtype=ftype)*-999.
        
        # same for AR_WWB_m_fn
        if os.path.exists(AR_WWB_m_fn):
            AR_WWB_m = np.fromfile(AR_WWB_m_fn,dtype=ftype).reshape(spxym)
        else:
            # Create a 2D array which can hold all the data
            AR_WWB_m = np.ones(spxym,dtype=ftype)*-999.
        
        # Make a var list to lable vars in left most column
        hdr_AR_WWB_m_n = ['WWB_#','Month','Rain','ETa','Qs_in','Qo_in',
                          'Qs_out','Qo_out','Qs_in_r','Qo_in_r','SSI',
                          'NB_months_avg']
        hdr_AR_WWB_m_i = ['Month','Rain','ETa','Qs_in','Qo_in',
                          'Qs_out','Qo_out','Qs_in_r','Qo_in_r','SSI',
                          'NB_months_avg']

#==============================================================================
#     Sheet 3 = WWB (Wetland water balance)
#==============================================================================

    ws_hdr = ['ETa','Vs','Qs_in','Qs_inII','Qs_out','Qs_outII','Qo_in',
              'Qo_out','Qs_in_r','Qs_in_rII','Qo_in_r','Rain','WB','ERR','SSI']
    um3    = ['(m3/','(m3)' ,'(m3/','(m3/','(m3/','(m3/','(m3/','(m3/',
              '(m3/','(m3/','(m3/','(m3/','(m3/','(m3/','-']
    umm    = ['(mm/','(mm)' ,'(mm/','(mm/','(mm/','(mm/','(mm/','(mm/',
              '(mm/','(mm/','(mm/','(mm/','(mm/','(mm/','-']
    
    nb_var = len(ws_hdr)

    u_m3_d, u_m3_m, u_m3_y, u_mm_d, u_mm_m, u_mm_y = [],[],[],[],[],[]
    cnt = 0
    for u in xrange(len(um3)):
        if cnt != 1 and cnt != nb_var-1:
            u_m3_d.append(um3[u]+'d)')
            u_m3_m.append(um3[u]+'m)')
            u_m3_y.append(um3[u]+'y)')
            u_mm_d.append(umm[u]+'d)')
            u_mm_m.append(umm[u]+'m)')
            u_mm_y.append(umm[u]+'y)')
        else:
            u_m3_d.append(um3[u])
            u_m3_m.append(um3[u])
            u_m3_y.append(um3[u])
            u_mm_d.append(umm[u])
            u_mm_m.append(umm[u])
            u_mm_y.append(umm[u])
        cnt += 1
    
    er_f = '=(INDIRECT(ADDRESS(ROW(),COLUMN()-9))-INDIRECT(ADDRESS(ROW()-1,\
COLUMN()-9)))-INDIRECT(ADDRESS(ROW(),COLUMN()-1))'
    
    # divide the areas wwbu_A_w by 1000 to obtain mm instead of meter
    
    A_2_mm_w = np.zeros(wwbu_A_w.shape)
    for i in xrange(len(wwbu_A_w)):
        A_2_mm_w[i] = wwbu_A_w[i]*1E-3

    ws_3.set_column(1,0,10)
    ws_4.set_column(1,0,10)
    ws_5.set_column(1,0,10)
    ws_51.set_column(1,11,14)
    ws_52.set_column(1,11,14)
    ws_55.set_column(1,0,10)
    ws_56.set_column(1,0,10)
    ws_57.set_column(1,0,10)

    for h in xrange(nb_wwb):
        # Write headers for all WWB sheets
        for i in xrange(nb_var):
            ws_3.write(0, nb_var*h+i+1, ws_hdr[i]+'_'+str(ar_wwb_u[h]))
            ws_3.write(1, nb_var*h+i+1, u_m3_d[i])
            ws_4.write(0, nb_var*h+i+1, ws_hdr[i]+'_'+str(ar_wwb_u[h]))
            ws_4.write(1, nb_var*h+i+1, u_m3_m[i])
            ws_5.write(0, nb_var*h+i+1, ws_hdr[i]+'_'+str(ar_wwb_u[h]))
            ws_5.write(1, nb_var*h+i+1, u_m3_y[i])
            ws_55.write(0, nb_var*h+i+1, ws_hdr[i]+'_'+str(ar_wwb_u[h]))
            ws_55.write(1, nb_var*h+i+1, u_mm_d[i])
            ws_56.write(0, nb_var*h+i+1, ws_hdr[i]+'_'+str(ar_wwb_u[h]))
            ws_56.write(1, nb_var*h+i+1, u_mm_m[i])
            ws_57.write(0, nb_var*h+i+1, ws_hdr[i]+'_'+str(ar_wwb_u[h]))
            ws_57.write(1, nb_var*h+i+1, u_mm_y[i])

        # Write WWB daily data - absolute volumes
        for l in xrange(n_rec):
            wb_tmp_d = Rainfl_w[h,l]-float(ET_out_w[h,l])-Qs_ot_w[h,l]-\
            Qs_ot_w_[h,l]+Qs_in_w[h,l]+Qs_in_w_[h,l]-Qo_ot_w[h,l]+Qo_in_w[h,l]
            if h == 0:
                ws_3.write(l+2,        0,  D_day[l], date_format)
            ws_3.write(l+2, nb_var*h+1+0,  float(ET_out_w[h,l]))
            ws_3.write(l+2, nb_var*h+1+1,  float(Vs_cel_w[h,l]))
            ws_3.write(l+2, nb_var*h+1+2,  Qs_in_w[h,l])
            ws_3.write(l+2, nb_var*h+1+3,  Qs_in_w_[h,l])
            ws_3.write(l+2, nb_var*h+1+4,  Qs_ot_w[h,l])
            ws_3.write(l+2, nb_var*h+1+5,  Qs_ot_w_[h,l])
            ws_3.write(l+2, nb_var*h+1+6,  Qo_in_w[h,l])
            ws_3.write(l+2, nb_var*h+1+7,  Qo_ot_w[h,l])
            ws_3.write(l+2, nb_var*h+1+8,  Qs_in_r[h,l])
            ws_3.write(l+2, nb_var*h+1+9,  Qs_in_r_[h,l])
            ws_3.write(l+2, nb_var*h+1+10, Qo_in_r[h,l])
#            try:
            ws_3.write(l+2, nb_var*h+1+8, float(Rainfl_w[h,l]))
#            except TypeError:
#                print l+2
#                print nb_var*h+1+8
#                print h, l
#                print Rainfl_w[h,l]
#                print type(Rainfl_w[h,l])
#                raise ValueError('here you go.. that is wrong')
            ws_3.write(l+2, nb_var*h+1+9,  wb_tmp_d)
            if l:
                ws_3.write_formula(l+2, nb_var*h+1+10, er_f)
            ws_3.write(l+2, nb_var*h+1+11, SSIndx_w[h,l])
            ws_3.write(l+2, nb_var*h+1+12, SSIndx_w_[h,l])

            # Write WWB daily data - values in mm
            if h == 0:
                ws_55.write(l+2,        0,  D_day[l], date_format)
            ws_55.write(l+2, nb_var*h+1+0,  float(ET_out_w[h,l])/A_2_mm_w[h])
            ws_55.write(l+2, nb_var*h+1+1,  float(Vs_cel_w[h,l])/A_2_mm_w[h])
            ws_55.write(l+2, nb_var*h+1+2,  Qs_in_w[h,l]/A_2_mm_w[h])
            ws_55.write(l+2, nb_var*h+1+3,  Qs_in_w_[h,l]/A_2_mm_w[h])
            ws_55.write(l+2, nb_var*h+1+4,  Qs_ot_w[h,l]/A_2_mm_w[h])
            ws_55.write(l+2, nb_var*h+1+5,  Qs_ot_w_[h,l]/A_2_mm_w[h])
            ws_55.write(l+2, nb_var*h+1+6,  Qo_in_w[h,l]/A_2_mm_w[h])
            ws_55.write(l+2, nb_var*h+1+7,  Qo_ot_w[h,l]/A_2_mm_w[h])
            ws_55.write(l+2, nb_var*h+1+8,  Qs_in_r[h,l]/A_2_mm_w[h])
            ws_55.write(l+2, nb_var*h+1+9,  Qs_in_r_[h,l]/A_2_mm_w[h])
            ws_55.write(l+2, nb_var*h+1+10, Qo_in_r[h,l]/A_2_mm_w[h])
#            try:
            ws_55.write(l+2, nb_var*h+1+8, float(Rainfl_w[h,l])/A_2_mm_w[h])
#            except TypeError:
#                print l+2
#                print nb_var*h+1+8
#                print h, l
#                print Rainfl_w[h,l]
#                print type(Rainfl_w[h,l])
#                raise ValueError('here you go.. that is wrong')
            ws_55.write(l+2, nb_var*h+1+9,  wb_tmp_d/A_2_mm_w[h])
            if l:
                ws_55.write_formula(l+2, nb_var*h+1+10, er_f)
            ws_55.write(l+2, nb_var*h+1+11, SSIndx_w[h,l])
            ws_55.write(l+2, nb_var*h+1+12, SSIndx_w_[h,l])
        
        # Write WWB monthly data - absolute volumes
        for l in xrange(len(D_mon)):
            wb_tmp_m = Pt_m_w[h,l]-float(Ea_m_w[h,l])-Qs_ot_m_w[h,l]+\
            Qs_in_m_w[h,l]-Qo_ot_m_w[h,l]+Qo_in_m_w[h,l]
            if h == 0:
                ws_4.write(l+2,        0,  D_mon[l], date_format)
            ws_4.write(l+2, nb_var*h+1+0,  float(Ea_m_w[h,l]))
            ws_4.write(l+2, nb_var*h+1+1,  float(Vs_m_w[h,l]))
            ws_4.write(l+2, nb_var*h+1+2,  Qs_in_m_w[h,l])
            ws_4.write(l+2, nb_var*h+1+3,  Qs_in_m_w_[h,l])
            ws_4.write(l+2, nb_var*h+1+4,  Qs_ot_m_w[h,l])
            ws_4.write(l+2, nb_var*h+1+5,  Qs_ot_m_w_[h,l])
            ws_4.write(l+2, nb_var*h+1+6,  Qo_in_m_w[h,l])
            ws_4.write(l+2, nb_var*h+1+7,  Qo_ot_m_w[h,l])
            ws_4.write(l+2, nb_var*h+1+8,  Qs_in_m_r[h,l])
            ws_4.write(l+2, nb_var*h+1+9,  Qs_in_m_r_[h,l])
            ws_4.write(l+2, nb_var*h+1+10, Qo_in_m_r[h,l])
            ws_4.write(l+2, nb_var*h+1+11, Pt_m_w[h,l])
            ws_4.write(l+2, nb_var*h+1+12, wb_tmp_m)
            if l:
                ws_4.write_formula(l+2, nb_var*h+1+13, er_f)
            ws_4.write(l+2, nb_var*h+1+14, SSInx_m_w[h,l])
            ws_4.write(l+2, nb_var*h+1+15, SSInx_m_w_[h,l])

            # Write WWB monthly data - values in mm
            if h == 0:
                ws_56.write(l+2,        0,  D_mon[l], date_format)
            ws_56.write(l+2, nb_var*h+1+0,  float(Ea_m_w[h,l])/A_2_mm_w[h])
            ws_56.write(l+2, nb_var*h+1+1,  float(Vs_m_w[h,l])/A_2_mm_w[h])
            ws_56.write(l+2, nb_var*h+1+2,  Qs_in_m_w[h,l]/A_2_mm_w[h])
            ws_56.write(l+2, nb_var*h+1+3,  Qs_in_m_w_[h,l]/A_2_mm_w[h])
            ws_56.write(l+2, nb_var*h+1+4,  Qs_ot_m_w[h,l]/A_2_mm_w[h])
            ws_56.write(l+2, nb_var*h+1+5,  Qs_ot_m_w_[h,l]/A_2_mm_w[h])
            ws_56.write(l+2, nb_var*h+1+6,  Qo_in_m_w[h,l]/A_2_mm_w[h])
            ws_56.write(l+2, nb_var*h+1+7,  Qo_ot_m_w[h,l]/A_2_mm_w[h])
            ws_56.write(l+2, nb_var*h+1+8,  Qs_in_m_r[h,l]/A_2_mm_w[h])
            ws_56.write(l+2, nb_var*h+1+9,  Qs_in_m_r_[h,l]/A_2_mm_w[h])
            ws_56.write(l+2, nb_var*h+1+10, Qo_in_m_r[h,l]/A_2_mm_w[h])
            ws_56.write(l+2, nb_var*h+1+11, Pt_m_w[h,l]/A_2_mm_w[h])
            ws_56.write(l+2, nb_var*h+1+12, wb_tmp_m/A_2_mm_w[h])
            if l:
                ws_56.write_formula(l+2, nb_var*h+1+13, er_f)
            ws_56.write(l+2, nb_var*h+1+14, SSInx_m_w[h,l])
            ws_56.write(l+2, nb_var*h+1+15, SSInx_m_w_[h,l])

        # Write WWB yearly data - absolute volumes
        for l in xrange(len(D_yer)):
            wb_tmp_y = Pt_y_w[h,l]-float(Ea_y_w[h,l])-Qs_ot_y_w[h,l]-\
            Qs_ot_y_w_[h,l]+Qs_in_y_w[h,l]+Qs_in_y_w_[h,l]-Qo_ot_y_w[h,l]+\
            Qo_in_y_w[h,l]
            if h == 0:
                ws_5.write(l+2,        0,  D_yer[l], date_format)
            ws_5.write(l+2, nb_var*h+1+0,  float(Ea_y_w[h,l]))
            ws_5.write(l+2, nb_var*h+1+1,  float(Vs_y_w[h,l]))
            ws_5.write(l+2, nb_var*h+1+2,  float(Vs_y_w_[h,l]))
            ws_5.write(l+2, nb_var*h+1+3,  Qs_in_y_w[h,l])
            ws_5.write(l+2, nb_var*h+1+4,  Qs_in_y_w_[h,l])
            ws_5.write(l+2, nb_var*h+1+5,  Qs_ot_y_w[h,l])
            ws_5.write(l+2, nb_var*h+1+6,  Qs_ot_y_w_[h,l])
            ws_5.write(l+2, nb_var*h+1+7,  Qo_in_y_w[h,l])
            ws_5.write(l+2, nb_var*h+1+8,  Qo_ot_y_w[h,l])
            ws_5.write(l+2, nb_var*h+1+9,  Qs_in_y_r[h,l])
            ws_5.write(l+2, nb_var*h+1+10, Qs_in_y_r_[h,l])
            ws_5.write(l+2, nb_var*h+1+11, Qo_in_y_r[h,l])
            ws_5.write(l+2, nb_var*h+1+12, Pt_y_w[h,l])
            ws_5.write(l+2, nb_var*h+1+13, wb_tmp_y)
            if l:
                ws_5.write_formula(l+2, nb_var*h+1+14, er_f)
            ws_5.write(l+2, nb_var*h+1+15, SSInx_y_w[h,l])

            # Write WWB yearly data - values in mm
            if h == 0:
                ws_57.write(l+2,        0,  D_yer[l], date_format)
            ws_57.write(l+2, nb_var*h+1+0,  float(Ea_y_w[h,l])/A_2_mm_w[h])
            ws_57.write(l+2, nb_var*h+1+1,  float(Vs_y_w[h,l])/A_2_mm_w[h])
            ws_57.write(l+2, nb_var*h+1+2,  float(Vs_y_w_[h,l])/A_2_mm_w[h])
            ws_57.write(l+2, nb_var*h+1+3,  Qs_in_y_w[h,l]/A_2_mm_w[h])
            ws_57.write(l+2, nb_var*h+1+4,  Qs_in_y_w_[h,l]/A_2_mm_w[h])
            ws_57.write(l+2, nb_var*h+1+5,  Qs_ot_y_w[h,l]/A_2_mm_w[h])
            ws_57.write(l+2, nb_var*h+1+6,  Qs_ot_y_w_[h,l]/A_2_mm_w[h])
            ws_57.write(l+2, nb_var*h+1+7,  Qo_in_y_w[h,l]/A_2_mm_w[h])
            ws_57.write(l+2, nb_var*h+1+8,  Qo_ot_y_w[h,l]/A_2_mm_w[h])
            ws_57.write(l+2, nb_var*h+1+9,  Qs_in_y_r[h,l]/A_2_mm_w[h])
            ws_57.write(l+2, nb_var*h+1+10, Qs_in_y_r_[h,l]/A_2_mm_w[h])
            ws_57.write(l+2, nb_var*h+1+11, Qo_in_y_r[h,l]/A_2_mm_w[h])
            ws_57.write(l+2, nb_var*h+1+12, Pt_y_w[h,l]/A_2_mm_w[h])
            ws_57.write(l+2, nb_var*h+1+13, wb_tmp_y/A_2_mm_w[h])
            if l:
                ws_57.write_formula(l+2, nb_var*h+1+14, er_f)
            ws_57.write(l+2, nb_var*h+1+15, SSInx_y_w[h,l])

        # all time averages - absolute volumes
        for i in xrange(nb_var):
            ws_51.write(0, i+1, ws_hdr[i])       # write variable names
            ws_51.write(1, i+1, u_m3_y[i])       # write units
            ws_52.write(0, i+1, ws_hdr[i])       # write variable names
            ws_52.write(1, i+1, u_mm_y[i])       # write units
        #print 'Ea_y_w.size:', Ea_y_w.size
#        print 'cnt_y:', cnt_y
        if cnt_y:
            wwb_y_tmp_m3 = []
            wwb_y_tmp_mm = []
            wb_tmp_y_avg = np.average(Pt_y_w[h,:])-np.average(Ea_y_w[h,:])-\
            np.average(Qs_ot_y_w[h,:])-np.average(Qs_ot_y_w_[h,:])+\
            np.average(Qs_in_y_w[h,:])+np.average(Qs_in_y_w_[h,:])-\
            np.average(Qo_ot_y_w[h,:])+np.average(Qo_in_y_w[h,:])
            ws_51.write(h+2, 0,  ar_wwb_u[h])    # write the wetland number
            ws_51.write(h+2, 1,  np.average(Ea_y_w[h,:]))
            wwb_y_tmp_m3.append( np.average(Ea_y_w[h,:]))
            ws_51.write(h+2, 2,  np.average(Vs_y_w[h,:]))
            wwb_y_tmp_m3.append( np.average(Vs_y_w[h,:]))
            ws_51.write(h+2, 3,  np.average(Vs_y_w_[h,:]))
            wwb_y_tmp_m3.append( np.average(Vs_y_w_[h,:]))
            ws_51.write(h+2, 4,  np.average(Qs_in_y_w[h,:]))
            wwb_y_tmp_m3.append( np.average(Qs_in_y_w[h,:]))
            ws_51.write(h+2, 5,  np.average(Qs_in_y_w_[h,:]))
            wwb_y_tmp_m3.append( np.average(Qs_in_y_w_[h,:]))
            ws_51.write(h+2, 6,  np.average(Qs_ot_y_w[h,:]))
            wwb_y_tmp_m3.append( np.average(Qs_ot_y_w[h,:]))
            ws_51.write(h+2, 7,  np.average(Qs_ot_y_w_[h,:]))
            wwb_y_tmp_m3.append( np.average(Qs_ot_y_w_[h,:]))
            ws_51.write(h+2, 8,  np.average(Qo_in_y_w[h,:]))
            wwb_y_tmp_m3.append( np.average(Qo_in_y_w[h,:]))
            ws_51.write(h+2, 9,  np.average(Qo_ot_y_w[h,:]))
            wwb_y_tmp_m3.append( np.average(Qo_ot_y_w[h,:]))
            ws_51.write(h+2, 10, np.average(Qs_in_y_r[h,:]))
            wwb_y_tmp_m3.append( np.average(Qs_in_y_r[h,:]))
            ws_51.write(h+2, 11, np.average(Qs_in_y_r_[h,:]))
            wwb_y_tmp_m3.append( np.average(Qs_in_y_r_[h,:]))
            ws_51.write(h+2, 12, np.average(Qo_in_y_r[h,:]))
            wwb_y_tmp_m3.append( np.average(Qo_in_y_r[h,:]))
            ws_51.write(h+2, 13, np.average(Pt_y_w[h,:]))
            wwb_y_tmp_m3.append( np.average(Pt_y_w[h,:]))
            ws_51.write(h+2, 14, wb_tmp_y_avg)
            wwb_y_tmp_m3.append( wb_tmp_y_avg)
            if l:
                ws_51.write(h+2, 15, 'NaN')
            ws_51.write(h+2, 16, np.average(SSInx_y_w[h,:]))
            wwb_y_tmp_m3.append( np.average(SSInx_y_w[h,:]))
            ws_51.write(h+2, 17, np.average(SSInx_y_w_[h,:]))
            wwb_y_tmp_m3.append( np.average(SSInx_y_w_[h,:]))
            ws_51.write(h+2, 18, cnt_y)
            wwb_y_tmp_m3.append( cnt_y)
            ws_51.write(h+2, 19, run_name)
            wwb_y_tmp_m3.append(-999.)
            # wwb_y_tmp_m3.append(run_name) only floats allowed

            # all time averages - values in mm
            ws_52.write(h+2, 0,  ar_wwb_u[h])    # write the wetland number
            ws_52.write(h+2, 1,  np.average(Ea_y_w[h,:])/A_2_mm_w[h])
            wwb_y_tmp_mm.append( np.average(Ea_y_w[h,:])/A_2_mm_w[h])
            ws_52.write(h+2, 2,  np.average(Vs_y_w[h,:])/A_2_mm_w[h])
            wwb_y_tmp_mm.append( np.average(Vs_y_w[h,:])/A_2_mm_w[h])
            ws_52.write(h+2, 3,  np.average(Vs_y_w_[h,:])/A_2_mm_w[h])
            wwb_y_tmp_mm.append( np.average(Vs_y_w_[h,:])/A_2_mm_w[h])
            ws_52.write(h+2, 4,  np.average(Qs_in_y_w[h,:])/A_2_mm_w[h])
            wwb_y_tmp_mm.append( np.average(Qs_in_y_w[h,:])/A_2_mm_w[h])
            ws_52.write(h+2, 5,  np.average(Qs_in_y_w_[h,:])/A_2_mm_w[h])
            wwb_y_tmp_mm.append( np.average(Qs_in_y_w_[h,:])/A_2_mm_w[h])
            ws_52.write(h+2, 6,  np.average(Qs_ot_y_w[h,:])/A_2_mm_w[h])
            wwb_y_tmp_mm.append( np.average(Qs_ot_y_w[h,:])/A_2_mm_w[h])
            ws_52.write(h+2, 7,  np.average(Qs_ot_y_w_[h,:])/A_2_mm_w[h])
            wwb_y_tmp_mm.append( np.average(Qs_ot_y_w_[h,:])/A_2_mm_w[h])
            ws_52.write(h+2, 8,  np.average(Qo_in_y_w[h,:])/A_2_mm_w[h])
            wwb_y_tmp_mm.append( np.average(Qo_in_y_w[h,:])/A_2_mm_w[h])
            ws_52.write(h+2, 9,  np.average(Qo_ot_y_w[h,:])/A_2_mm_w[h])
            wwb_y_tmp_mm.append( np.average(Qo_ot_y_w[h,:])/A_2_mm_w[h])
            ws_52.write(h+2, 10, np.average(Qs_in_y_r[h,:])/A_2_mm_w[h])
            wwb_y_tmp_mm.append( np.average(Qs_in_y_r[h,:])/A_2_mm_w[h])
            ws_52.write(h+2, 11, np.average(Qs_in_y_r_[h,:])/A_2_mm_w[h])
            wwb_y_tmp_mm.append( np.average(Qs_in_y_r_[h,:])/A_2_mm_w[h])
            ws_52.write(h+2, 12, np.average(Qo_in_y_r[h,:])/A_2_mm_w[h])
            wwb_y_tmp_mm.append( np.average(Qo_in_y_r[h,:])/A_2_mm_w[h])
            ws_52.write(h+2, 13, np.average(Pt_y_w[h,:])/A_2_mm_w[h])
            wwb_y_tmp_mm.append( np.average(Pt_y_w[h,:])/A_2_mm_w[h])
            ws_52.write(h+2, 14, wb_tmp_y_avg/A_2_mm_w[h])
            wwb_y_tmp_mm.append( wb_tmp_y_avg/A_2_mm_w[h])
            if l:
                ws_52.write(h+2, 15, 'NaN')
            ws_52.write(h+2, 16, np.average(SSInx_y_w[h,:]))
            wwb_y_tmp_mm.append( np.average(SSInx_y_w[h,:]))
            ws_52.write(h+2, 17, np.average(SSInx_y_w_[h,:]))
            wwb_y_tmp_mm.append( np.average(SSInx_y_w_[h,:]))
            ws_52.write(h+2, 18, cnt_y)
            wwb_y_tmp_mm.append( cnt_y)
            ws_52.write(h+2, 19, run_name)
            wwb_y_tmp_mm.append( -999.)
            # wwb_y_tmp_mm.append(run_name) only floats allowed

#            print 'wwb_unit:', ar_wwb_u[h]
#            print 'len(wwb_y_tmp_m3):', len(wwb_y_tmp_m3)
#            print 'len(wwb_y_tmp_mm):', len(wwb_y_tmp_mm)

#            M_WWB_V holds the month averages for from the current run.
#            Shape is (12,nb_wwb,nb_vars)
            wwb = ar_wwb_u[h]
            if rt == 'n' and not run_name[1] == 'P':
                AR_WWB_y[0:13,wwb] = wwb_y_tmp_m3
                for m in xrange(12):
                    AR_WWB_m[0,wwb*12+m] = wwb
                    AR_WWB_m[1,wwb*12+m] = m+1
                    AR_WWB_m[2:11,wwb*12+m] = M_WWB_V[m,h,:]
                    AR_WWB_m[11,wwb*12+m] = M_count[m]
#                print wwb_y_tmp_m3
#                print 'len of AR_WWB_y[ 0:13,wwb]:', len(AR_WWB_y[ 0:13,wwb])
                AR_WWB_y[26:39,wwb] = wwb_y_tmp_mm
                wwb_y_tmp_mm
#                print 'len of AR_WWB_y[26:39,wwb]:', len(AR_WWB_y[26:39,wwb])
                if WWBU_transfer_n2i:
                    if wwb in WWB_n2i:
                        AR_WWB_y[13:26,wwb] = wwb_y_tmp_m3
                        AR_WWB_y[39:52,wwb] = wwb_y_tmp_mm
                        for m in xrange(12):
                            AR_WWB_m[0,wwb*12+m] = wwb
                            AR_WWB_m[1,wwb*12+m] = m+1
                            AR_WWB_m[13:22,wwb*12+m] = M_WWB_V[m,h,:]
                            AR_WWB_m[22,wwb*12+m] = M_count[m]
            elif rt == 'i' and not run_name[1] == 'P':
                AR_WWB_y[13:26,wwb] = wwb_y_tmp_m3
                AR_WWB_y[39:52,wwb] = wwb_y_tmp_mm
                for m in xrange(12):
                    AR_WWB_m[0,wwb*12+m] = wwb
                    AR_WWB_m[1,wwb*12+m] = m+1
                    AR_WWB_m[13:22,wwb*12+m] = M_WWB_V[m,h,:]
                    AR_WWB_m[22,wwb*12+m] = M_count[m]
                    
#                raise ValueError('This is neither a current stat run nor a \
#impact run and therefore funciton is terminating...')
#            print 'shape of AR_WWB_y:', AR_WWB_y.shape
    if (rt == 'n' or rt == 'i') and not run_name[1] == 'P':
#        AR_WWB_m[0:12,0] = hdr_AR_WWB_m_n
#        AR_WWB_m[12:23,0] = hdr_AR_WWB_m_i
        AR_WWB_y.tofile(AR_WWB_y_fn)
        AR_WWB_m.tofile(AR_WWB_m_fn)
        np.savetxt(AR_WWB_y_fn+'_txt.dat', AR_WWB_y)
        np.savetxt(AR_WWB_m_fn+'_txt.dat', AR_WWB_m, delimiter='\t', fmt='%g')
        
        print 'Adding data from %s to the WWB-yearly database' %run_name
        print 'Adding data from %s to the WWB-month database' %run_name
        

    # Write some cell lables
    ws_51.write(0,  0, 'WWB')
    ws_51.write(1,  0, 'Unit#')
    ws_51.write(0, 13, 'N_yrs')
    ws_51.write(1, 13, 'averaged')
    ws_51.write(0, 14, 'Catchment')
    ws_52.write(0,  0, 'WWB')
    ws_52.write(1,  0, 'Unit#')
    ws_52.write(0, 13, 'N_yrs')
    ws_52.write(1, 13, 'averaged')
    ws_52.write(0, 14, 'Catchment')
    
    print n_rec, 'records written to the file', fn_Qout
    
#    plt.clf()
#    fig = plt.figure()
#    ax  = fig.add_subplot(111)
#    for w in xrange(nb_wwb):
#        ax.plot(SSIndx_w[w])
#    plt.savefig('SSI_test.png')

#==============================================================================
# SHEET 6: WRITE COORDINATES OF BWin AND BWout INTO XLSX FILE
#==============================================================================
    coords6  = [ar_bcw_x, ar_bcw_y, bcw_wwb_l]
    b_hdrs6  = ['BWin_X','BWin_Y', 'WWB_unit']

    coords7  = [ar_bwc_x, ar_bwc_y, bwc_wwb_l]
    b_hdrs7  = ['BWout_X', 'BWout_Y', 'WWB_unit']
    
    coords8  = [ar_wwb_x, ar_wwb_y, ar_wwb_l]
    b_hdrs8  = ['WWB_X', 'WWB_Y', 'WWB_unit']
    
    coords9  = [ar_wwb_C_x, ar_wwb_C_y, ar_wwb_C_l]
    b_hdrs9  = ['TRD_X', 'TRD_c_Y', 'Cell_Label']

# this is dublicate of coords7..    
#    coords10 = [ar_wwb_W_x, ar_wwb_W_y, ar_wwb_W_l]
#    b_hdrs10 = ['WWB_X_W' , 'WWB_W_Y' , 'WWBu_W'  ]
#    
#    coords11 = [ar_cdC_x, ar_cdC_y, ar_cdC_l]
#    b_hdrs11 = ['cd_C_X', 'cd_C_Y', 'cd_C'  ]
#    
#    coords12 = [ar_cd_x, ar_cd_y, ar_cd_l]
#    b_hdrs12 = ['cd_X' , 'cd_Y' , 'cd'   ]
    
    coords13 = [ar_cl_x, ar_cl_y, ar_cl_l]
    b_hdrs13 = ['CellLabel_X', 'CellLabel_Y', 'CellLabel']

    for c in xrange(len(b_hdrs6)):
        ws_6.write(0, c, b_hdrs6[c])
        for l in xrange(len(coords6[c])):
            ws_6.write(l+1, c, coords6[c][l])
    for c in xrange(len(b_hdrs7)):
        ws_7.write(0, c, b_hdrs7[c])
        for l in xrange(len(coords7[c])):
            ws_7.write(l+1, c, coords7[c][l])
    for c in xrange(len(b_hdrs8)):
        ws_8.write(0, c, b_hdrs8[c])
        for l in xrange(len(coords8[c])):
            ws_8.write(l+1, c, coords8[c][l])
    for c in xrange(len(b_hdrs9)):
        ws_9.write(0, c, b_hdrs9[c])
        for l in xrange(len(coords9[c])):
            ws_9.write(l+1, c, coords9[c][l])
#    for c in xrange(len(b_hdrs10)):
#        ws_10.write(0, c, b_hdrs10[c])
#        for l in xrange(len(coords10[c])):
#            ws_10.write(l+1, c, coords10[c][l])
#    for c in xrange(len(b_hdrs11)):
#        ws_11.write(0, c, b_hdrs11[c])
#        for l in xrange(len(coords11[c])):
#            ws_11.write(l+1, c, coords11[c][l])
#    for c in xrange(len(b_hdrs12)):
#        ws_12.write(0, c, b_hdrs12[c])
#        for l in xrange(len(coords12[c])):
#            ws_12.write(l+1, c, coords12[c][l])
    for c in xrange(len(b_hdrs13)):
        ws_13.write(0, c, b_hdrs13[c])
        for l in xrange(len(coords13[c])):
            ws_13.write(l+1, c, coords13[c][l])


#==============================================================================
# Write a date column into the WS 14 and 15 sheet (DAMs and PZM)
#==============================================================================
    #print 'len(D_day)', len(D_day)
    if Dams_on:
        for l in xrange(n_rec-2):
            ws_14.write(l+4, 1, D_day[l+2], date_format)
    if Piezos_on:
        for ll in xrange(n_rec):
            ws_15.write(ll+2, 1, D_day[ll], date_format)

    wb.close()
    
    shutil.copyfile(fn_Qout, fn_Qout2)
    

#===============================================================================
#     # RIVER - FLOW DURATION CURVE
#===============================================================================
    if np.sum(ar_Qsim) > 0:
        Q_sim_sorted = np.sort(ar_Qsim)
        Q_sim_descnd = Q_sim_sorted[::-1]
        ranks = np.zeros(n_rec)
        frqcy = np.zeros(n_rec)
        
        for i in xrange(n_rec):
            ranks[i] = i
            frqcy[i] = 100 * (i/(n_rec+1.0))
    
        #print 'rank min, max', np.min(ranks), np.max(ranks)
        #print 'Q_sim_descnd min, max, len:', np.min(Q_sim_descnd), np.max(Q_sim_descnd), len(Q_sim_descnd)
        #plt.plot(frqcy, Q_sim_descnd)
    
        fig_x, fig_y = 8, 3
    
        plt.clf()
        fig     = plt.figure()
        fig_p   = plt.gcf()
        fig_p.set_size_inches(fig_x, fig_y)
        #ax1 = fig.add_axes(rect)
        
        ax1 = fig.add_subplot(111)
        '''
        left = 0.125
            the left side of the subplots of the figure
        right = 0.9
            the right side of the subplots of the figure
        bottom = 0.1
            the bottom of the subplots of the figure
        top = 0.9
            the top of the subplots of the figure
        wspace = 0.2
            the amount of width reserved for blank space between subplots
        hspace = 0.2
            the amount of height reserved for white space between subplots
        '''
        fig.subplots_adjust(left=0.085, right = 0.97, bottom = 0.14, top = 0.97)
        
        #ax1.set_xscale('log')
        #ax1.set_yscale('log')
        ax1.plot(frqcy, Q_sim_descnd, 'gray', linewidth=1)
        
        #xticklabel = ax1.get_xticks()
        
        # solution example from http://old.nabble.com/Bold-Latex-Tick-Labels-td28037900.html
        #tick_locs = range(start, stop, increment)
        #plt.xticks(tick_locs, [r"$\mathbf{%s}$" % x for x in tick_locs]) 
        
        ymax_V = ax1.get_ylim()[1]
        ymin_V = ax1.get_ylim()[0]
        
        xmin_V = ax1.get_xlim()[0]
        xmax_V = ax1.get_xlim()[1]
        
        #ax1.axvlines([10,40,60,90], [x*0.0 for x in xrange(4)], [ymax_V for x in xrange(4)])
        #ax1.vlines([10,40,60,90], ymin_V, ymax_V, color='k', linestyles='solid')
        for v in [10,40,60,90]:
            plt.axvline(v, color='k', linestyle='dotted')
        #print 'xmax_V', xmax_V
        ax1.set_xlim(xmin_V,xmax_V)
        #ax1.boxplot(dataT7C, vert=1, positions = headsT7C, widths=barwidth, whis=1.5)
        #ax1.set_ylim(ymin_V,ymax_V)
        ax1.set_yscale('log')
        
        #print 'Volume Grpah has limits: '
        #print ax1.get_ylim()
        #print '\n\n'
        #print 'Limit from confic are: ' + str(eval(Vol_name))
        #ax1.set_xlim(xlim)
        ax1.set_ylabel(r'$\mathbf{Mean~Monthly~Flow~(m^3~s^{-1})}$')
        ax1.set_xlabel(r'$\mathbf{Exceedance~Probability~}$(\%)', fontweight='bold')
        ax1.set_xticks([x*10 for x in xrange(11)])
        xticklabels, yticklabels = [], []
        for ix in ax1.get_xticks():
            xticklabels.append(r'$\mathbf{%s}$' % int(ix))
        for iy in ax1.get_yticks():
            iy_e = log10(iy)
            yticklabels.append(r'$\mathbf{10^{%s}}$' % int(iy_e))
        ax1.set_xticklabels(xticklabels)
        ax1.set_yticklabels(yticklabels)
    #     lowest 20 % of y
        yu = (ymax_V-ymin_V)*1e-2
        yl = (ymax_V-ymin_V)*.4e-3
        xy_1 = [(0,yu),(10,yu),(40,yu),(60,yu),(90,yu)]
        xy_2 = [(10,yu),(40,yu),(60,yu),(90,yu),(100,yu)]
        xy_3 = [(5,yl),(25,yl),(50,yl),(75,yl),(95,yl)]
        text = ['High\nflows', 'Moist conditions', 'Mid range flows',
                'Dry conditions', 'Low\nflows']
    #    a1=()
    #    a2=()
    #    a3=()
    #    a4=()
    #    a5=()
        for ar in xrange(len(xy_1)):
            ax1.annotate('', xy=xy_1[ar], xytext=xy_2[ar], xycoords='data',
                    arrowprops=dict(arrowstyle='<->'))
            ax1.annotate(text[ar], xy=xy_3[ar],  xycoords='data',
                         horizontalalignment='center', verticalalignment='bottom',
                         fontsize=10)
        #ax1.arrow(5, , 0.5, 0.5, head_width=0.05, head_length=0.1, fc='k', ec='k')
        #ax1.grid(True)
        
        #yticks = ax1.set_yticks([0,       depthA      , 200     , 400     , depthB      , 600     , 800     , ymax])
        #ytickslables = ax1.set_yticklabels([r'$Sur.~0$',  r'$A~ 113$', r'$200$', r'$400$', r'$B~ 454$', r'$600$', r'$800$', r'$1000$'])
        
        #xyax  = mplt.gca()
        #xyax.set_ylim(xyax.get_ylim()[::-1])
        #xticklabel = ax1.get_xticklabels()
        #plt.setp(xticklabel, visible=False)
        
        
        #ax1b   = ax1.twiny()
        #ax1b.set_xlabel(r'$Number~ of~ Measurements$')
        ##ax1b.set_xscale('log')
        #ax1b.set_xticks(headsT7C)
        #ax1b.set_xlim(ax1.get_xlim())
        #ax1b.set_xticklabels(nT7C)
        #ax1b.text(txtX,txtY,r'$T7C$', fontsize=11)
        fig_fn = 'Results\\Flow_Duration_Curve.png'
        #ax1.relim()
        plt.savefig(fig_fn, dpi=150)
        #print 'Ylim: ', ax1.get_ylim()
    
    
    #t = np.arange(0.0, 100, 0.1)
    #s = np.sin(2*np.pi*t)
    #
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    #ax.plot(t, s, '-o')
    #ax.vlines([10,40,60,90], [x*0.0 for x in xrange(4)], [1.0 for x in xrange(4)])


#==============================================================================
#     SATURATION DURATION CURVES FOR WETLANDS
#==============================================================================

    if n_rec > 500 and (rt == 'n' or rt == 'i') and not run_name[1] == 'P':
        shapexy = (n_rec,200)
        
        # Numpy binary file names
        if rt == 'n':
            AR_SSI_WWB_d_fn = os.getcwd() + '\\..\\..' +\
            '\\Saturation_Duration_Curves\\' + 'AR_SSI_WWB_n'
            AR_FQC_WWB_d_fn = os.getcwd() + '\\..\\..' +\
            '\\Saturation_Duration_Curves\\' + 'AR_FQC_WWB_n'
        elif rt == 'i':
            AR_SSI_WWB_d_fn = os.getcwd() + '\\..\\..' +\
            '\\Saturation_Duration_Curves\\' + 'AR_SSI_WWB_i'
            AR_FQC_WWB_d_fn = os.getcwd() + '\\..\\..' +\
            '\\Saturation_Duration_Curves\\' + 'AR_FQC_WWB_i'
        if WWBU_transfer_n2i and rt == 'n':
            AR_SSI_WWB_d_fn_i = os.getcwd() + '\\..\\..' +\
            '\\Saturation_Duration_Curves\\' + 'AR_SSI_WWB_i'
            AR_FQC_WWB_d_fn_i = os.getcwd() + '\\..\\..' +\
            '\\Saturation_Duration_Curves\\' + 'AR_FQC_WWB_i'
            if os.path.exists(AR_SSI_WWB_d_fn_i) and os.path.exists(AR_FQC_WWB_d_fn_i):
                AR_SSI_WWB_i = np.fromfile(AR_SSI_WWB_d_fn_i,dtype=ftype).reshape(shapexy)
                AR_FQC_WWB_i = np.fromfile(AR_FQC_WWB_d_fn_i,dtype=ftype).reshape(shapexy)
            else:
                AR_SSI_WWB_i = np.ones(shapexy,dtype=ftype)*-999.
                AR_FQC_WWB_i = np.ones(shapexy,dtype=ftype)*-999.
        if os.path.exists(AR_SSI_WWB_d_fn) and os.path.exists(AR_FQC_WWB_d_fn):
            AR_SSI_WWB = np.fromfile(AR_SSI_WWB_d_fn,dtype=ftype).reshape(shapexy)
            AR_FQC_WWB = np.fromfile(AR_FQC_WWB_d_fn,dtype=ftype).reshape(shapexy)
        else:
            # Create a 2D array which can hold all the data
            AR_SSI_WWB = np.ones(shapexy,dtype=ftype)*-999.
            AR_FQC_WWB = np.ones(shapexy,dtype=ftype)*-999.
        
        lines   = []
        tab_leg = []
    
        fig_x, fig_y = 8, 8
    
        plt.clf()
        fig     = plt.figure()
        fig_p   = plt.gcf()
        fig_p.set_size_inches(fig_x, fig_y)
        
        ax1 = fig.add_subplot(111)
    
        fig.subplots_adjust(left=0.085, right=0.97, bottom=0.14, top=0.97)
    
        for w in xrange(nb_wwb):
            wu = ar_wwb_u[w]
#            print 'Wetland Unit Number:', wu
            #SSIndx_w[w][SSIndx_w[w]<0] = np.nan
            SSI_sorted = np.sort(SSIndx_w[w])
            SSI_descnd = SSI_sorted[::-1]
            ranks = np.zeros(n_rec)
            frqcy = np.zeros(n_rec)
            
            for i in xrange(n_rec):
                ranks[i] = i
                frqcy[i] = 100 * (i/(n_rec+1.0))
    
            #ra = lambda: random.randint(0,255)
            rgb = '#%02X%02X%02X' % (tuple([random.randint(0,255) for x in xrange(3)]))
        
            lines += ax1.plot(frqcy, SSI_descnd, color=rgb, linewidth=1)
            
            tab_leg.append(r'$\mathbf{Unit~%s}$' %str(ar_wwb_u[w]))
            tab_leg = tab_leg[::-1]
            
            # Following conditions will add data to the SSI_i and SSI_n data
            # bases based on whether or not n-data needs to be written to the
            # SSI_i data base. If so then i_data is prevented to overwrite
            # when i_runs are processed afterwards. Condition is that i_runs
            # have the WWBU_transfer_n2i argument set to true..
            # ToDo: Simplify this somehow!!! Make input file with tranfer units
            
            if rt == 'i' and not run_name[1] == 'P':
                if WWBU_transfer_n2i:
                    if wu not in WWB_n2i:
                        AR_SSI_WWB[:,wu] = SSI_descnd
                        AR_FQC_WWB[:,wu] = frqcy
                else:
                    # wirte all i stuff to i database
                    AR_SSI_WWB[:,wu] = SSI_descnd
                    AR_FQC_WWB[:,wu] = frqcy
            else:
                AR_SSI_WWB[:,wu] = SSI_descnd
                AR_FQC_WWB[:,wu] = frqcy

            if WWBU_transfer_n2i and rt == 'n' and not run_name[1] == 'P':
                if wu in WWB_n2i:
                    AR_SSI_WWB_i[:,wu] = SSI_descnd
                    AR_FQC_WWB_i[:,wu] = frqcy

        ax1.legend(lines, tab_leg, loc='upper right', fancybox=True)
        leg = ax1.get_legend()
        leg.get_frame().set_alpha(0.75)
        
        #xticklabel = ax1.get_xticks()
        
        # solution example from http://old.nabble.com/Bold-Latex-Tick-Labels-td28037900.html
        #tick_locs = range(start, stop, increment)
        #plt.xticks(tick_locs, [r"$\mathbf{%s}$" % x for x in tick_locs]) 
        
#        ymax_V = ax1.get_ylim()[1]
#        ymin_V = ax1.get_ylim()[0]
        ymin_V, ymax_V = 0, 100
        xmin_V = ax1.get_xlim()[0]
        xmax_V = ax1.get_xlim()[1]
        
        for v in [10,40,60,90]:
            plt.axvline(v, color='k', linestyle='dotted')
    
        ax1.set_xlim(xmin_V,xmax_V)
        ax1.set_ylim(ymin_V,ymax_V)    
        #ax1.set_yscale('log')
    
        ax1.set_ylabel(r'$\mathbf{Daily~Saturation~}(\%)$', fontweight='bold')
        ax1.set_xlabel(r'$\mathbf{Exceedance~Probability~}(\%)$', fontweight='bold')
        ax1.set_xticks([x*10 for x in xrange(11)])
        xticklabels, yticklabels = [], []
        for ix in ax1.get_xticks():
            xticklabels.append(r'$\mathbf{%s}$' % int(ix))
        for iy in ax1.get_yticks():
            #iy_e = log10(iy)
            yticklabels.append(r'$\mathbf{%s}$' % int(iy))
        ax1.set_xticklabels(xticklabels)
        ax1.set_yticklabels(yticklabels)
    
        yu = (ymax_V-ymin_V)*1e-8
        yl = (ymax_V-ymin_V)*.4e-7
        xy_1 = [(0,yu),(10,yu),(40,yu),(60,yu),(90,yu)]
        xy_2 = [(10,yu),(40,yu),(60,yu),(90,yu),(100,yu)]
        xy_3 = [(5,yl),(25,yl),(50,yl),(75,yl),(95,yl)]
        text = ['Very\nWet', 'Wet', 'Mid Range',
                'Dry', 'Very\nDry']
    
        for ar in xrange(len(xy_1)):
            ax1.annotate('', xy=xy_1[ar], xytext=xy_2[ar], xycoords='data',
                    arrowprops=dict(arrowstyle='<->'))
            ax1.annotate(text[ar], xy=xy_3[ar],  xycoords='data',\
                         horizontalalignment='center', verticalalignment='bottom',
                         fontsize=10)
        wwb_u_names = ''
        for u in ar_wwb_u:
            wwb_u_names += str(u)+'_'
        fig_fn = os.getcwd()+'\\..\\..'+'\\Saturation_Duration_Curves\\'+\
                 run_name+'_'+wwb_u_names+'Saturation_Duration.png'
        
        if len(fig_fn) > 250:
            first_len = len(fig_fn) - (len(fig_fn) - 255 + 23)
            fig_fn = fig_fn[:first_len] + fig_fn[-23:]
    
        plt.savefig(fig_fn, dpi=150)
        

        print 'Adding data from %s to the SSI-Database' %run_name

        AR_SSI_WWB.tofile(AR_SSI_WWB_d_fn)
        AR_FQC_WWB.tofile(AR_FQC_WWB_d_fn)
        if WWBU_transfer_n2i and rt == 'n' and not run_name[1] == 'P':
            AR_SSI_WWB_i.tofile(AR_SSI_WWB_d_fn_i)
            AR_FQC_WWB_i.tofile(AR_FQC_WWB_d_fn_i)



#    a = raw_input('Press Enter to quit')


#def some_name(date, data, wwb_index):
#    pass
    # 3-D nd.array to hold data