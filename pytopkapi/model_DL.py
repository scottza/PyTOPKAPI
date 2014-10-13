"""Core logic of PyTOPKAPI.

The `run` function in this module contains the logic to run a TOPKAPI
simulation based on the parameters specified in an INI file.

"""

#General module importation
import os
from ConfigParser import SafeConfigParser, NoOptionError

import numpy as np
import tables as h5

#Personnal module importation
import pytopkapi
import utils as ut
import pretreatment as pm
import fluxes as fl
import ode as om
import evap as em
from infiltration import green_ampt_cum_infiltration
import datetime as dtm
import warnings
import exceptions

warnings.filterwarnings('error', message='invalid value encountered in \
double_scalars', category=exceptions.RuntimeWarning)

def run_DL(ini_file='TOPKAPI.ini'):
    """Run the model with the set-up defined by `ini_file`.

    """
    t_ini = dtm.datetime.today()

    ##================================##
    ##  Read the input file (*.ini)   ##
    ##================================##
    config = SafeConfigParser()
    config.read(ini_file)
    print 'Reading\'',ini_file,'\''

    ##~~~~~~ Numerical_options ~~~~~~##
    solve_s = config.getfloat('numerical_options', 'solve_s')
    solve_o = config.getfloat('numerical_options', 'solve_o')
    solve_c = config.getfloat('numerical_options', 'solve_c')
    only_channel_output = config.getboolean('numerical_options',
                                            'only_channel_output')
    try:
        SSI_ET_treshold = config.getfloat('forcing_options', 'SSI_ET_switch')
    except NoOptionError:
        SSI_ET_treshold = 0.75

    ##~~~~~~~~~~~ input files ~~~~~~~~~~~##
    #Param
    file_global_param = config.get('input_files', 'file_global_param')
    file_cell_param   = config.get('input_files', 'file_cell_param')
    #Rain
    file_rain         = config.get('input_files', 'file_rain')
    #ETP
    file_ET           = config.get('input_files', 'file_ET')

    #~~~~~~~~~~~ Group (simulated event) ~~~~~~~~~~~##
    group_name = config.get('groups', 'group_name')

    ##~~~~~~ Calibration ~~~~~~##
    fac_L     = config.getfloat('calib_params', 'fac_L')
    fac_L_    = config.getfloat('calib_params', 'fac_L_II')
    fac_Ks    = config.getfloat('calib_params', 'fac_Ks')
    fac_Ks_   = config.getfloat('calib_params', 'fac_Ks_II')
    fac_n_o   = config.getfloat('calib_params', 'fac_n_o')
    fac_n_c   = config.getfloat('calib_params', 'fac_n_c')
    fac_th_s  = config.getfloat('calib_params', 'fac_th_s')
    fac_th_s_ = config.getfloat('calib_params', 'fac_th_s_II')

    ##~~~~~~ Starting time ~~~~~~~~##
    try:
        start_date = config.get('start_date', 'start_date')
        stmp = []
        for s in start_date.split('/'):
            stmp.append(int(s))
        start_date = dtm.datetime(stmp[2],stmp[1],stmp[0])
    except NoOptionError:
        print 'WARNING: No starting date was provided. This is necessary if \
you want to use the AVC evaporation method. You also need to specify the \
the \'evaporation_options\'.'

    try:
        mon_i = config.getint('evaporation_options', 'month_increase')
        mon_d= config.getint('evaporation_options', 'month_decrease')
    except:
        pass

    ##~~~~~~ Forcing options ~~~~~~~~##
    rain_distributed = config.getboolean('forcing_options', 'rain_dist')
    ET_distributed   = config.getboolean('forcing_options', 'ET_dist')
    try:
        Dams_on      = config.getboolean('forcing_options', 'Dams_on')
        Qseep_dam_mm = config.getfloat('forcing_options', 'Seepage_Dam_mm')
        Dams_ini     = config.getfloat('forcing_options', 'Dams_initial')
        WL_dam_avg   = config.getfloat('forcing_options', 'Average_WL_dam')

    except NoOptionError:
        print 'WARNING: At least one of the parameters for DAMs were not \
provided in the model initiation file. The DAM functionallity was switched \
off.\n\The required parameters for dams are under \'forcing_options\':\n\
Dams_on - Boolean (1=on 0=off)\nQseep_dam_mm - Seepage rate for dam water \
loss (mm/Timestep)\nDams_ini - Initial water level as percentage of FSL\n\
WL_dam_avg - Average depth of water in dam at FSL.'
        Dams_on      = 0
        WL_dam_avg   = 0.
        Dams_ini     = 0.

    if Dams_on:
        print 'Dam module active...'
    else:
        print 'Dam module is off...'

    try:
        ET_module_option = config.get('numerical_options', 'ET_module')
    except NoOptionError:
        print 'WARNING: No option for the evaporation module specified in \
\'%s\'. The options for \'ET_module\' under \'Numerical_options\' can be:\n\
Standard - for the \'evapot_soil_overland\' function\n\
Feddes   - for the \'evapot_soil_overland_feddes\' function in the \'evap\' \
module.\n\
In order to maintain backwards compatability the ET option was set to \
\'Standard\'.' %(ini_file)
        ET_module_option = 'Standard'


    ##~~~~~~ External flows ~~~~~~##
    external_flow = config.getboolean('external_flow', 'external_flow')
    if external_flow:
        file_Qexternal_flow = config.get('external_flow',
                                         'file_Qexternal_flow')
        Xexternal_flow = config.getfloat('external_flow', 'Xexternal_flow')
        Yexternal_flow = config.getfloat('external_flow', 'Yexternal_flow')

    ##~~~~~~~~~~~ output files ~~~~~~~~~~##
    file_out = config.get('output_files', 'file_out')
    ut.check_file_exist(file_out) #create path_out if it doesn't exist
    if os.path.exists(file_out):
        first_run = False
    else:
        first_run = True

    append_output = config.getboolean('output_files', 'append_output')
    if append_output is True:
        fmode = 'a'
    else:
        fmode = 'w'

    ##============================##
    ##   Read the forcing data    ##
    ##============================##
    print 'Reading the forcing data...'

    #~~~~Rainfall
    h5file_in = h5.openFile(file_rain,mode='r')
    group = '/'+group_name+'/'
    node = h5file_in.getNode(group+'rainfall')
    ndar_rain = node.read()
    h5file_in.close()

    n_cols_rain = len(ndar_rain[0,:])
    if not rain_distributed and n_cols_rain > 1:
        raise ValueError('You have chosen non-distributed rainfall, however\
 you provide more than one column of rainfall in the rain forcing file')

    if rain_distributed and n_cols_rain < 2:
        raise ValueError('You have chosen distributed rainfall, however you\
 provide only one cell with rainfall in the rain forcing file')


    #~~~~ETr - Reference crop ET
    h5file_in = h5.openFile(file_ET,mode='r')
    group = '/'+group_name+'/'
    node = h5file_in.getNode(group+'ETr')
    ndar_ETr = node.read()
    h5file_in.close()

    #~~~~ETo - Open water potential evap.
    h5file_in = h5.openFile(file_ET,mode='r')
    group = '/'+group_name+'/'
    node = h5file_in.getNode(group+'ETo')
    ndar_ETo = node.read()
    h5file_in.close()

    #~~~~external_flow flows
    if external_flow:
        ar_Qexternal_flow = np.loadtxt(file_Qexternal_flow)[:, 5]


    ##============================##
    ## Pretreatment of input data ##
    ##============================##
    print 'Processing pretreatment of input data...'

    #~~~~Read Global parameters file
    X, Dt, alpha_s, \
    alpha_o, alpha_c, \
    A_thres, W_min, W_max = pm.read_global_parameters(file_global_param)

    #~~~~Read Cell parameters file
    ar_cell_label, ar_coorx, \
    ar_coory, ar_lambda, \
    ar_Xc, ar_wwb, \
    ar_tan_beta, ar_tan_beta_channel, \
    ar_L0, ar_Ks0, \
    ar_theta_r, ar_theta_s0, \
    ar_n_o0, ar_n_c0, \
    ar_cell_down, ar_pVs_t0, \
    ar_Vo_t0, ar_Qc_t0, \
    ar_kc, psi_b, lamda,\
    ar_L0_, ar_Ks0_, ar_theta_r_, \
    ar_theta_s0_, ar_pVs_t0_, psi_b_,\
    lamda_ = pm.read_cell_parameters_DL(file_cell_param)


    #~~~~Number of cell in the catchment
    nb_cell = len(ar_cell_label)

    #~~~~Computation of cell order
    ar_label_sort = pm.sort_cell(ar_cell_label, ar_cell_down)

    #~~~~Computation of upcells
    li_cell_up = pm.direct_up_cell(ar_cell_label, ar_cell_down, ar_label_sort)

    #~~~~Computation of drained area
    ar_A_drained = pm.drained_area(ar_label_sort, li_cell_up, X)

    #~~~~Apply calibration factors to the parameter values
    ar_L        = ar_L0*fac_L
    ar_L_       = ar_L0_*fac_L_
    ar_Ks       = ar_Ks0*fac_Ks
    ar_Ks_      = ar_Ks0_*fac_Ks_
    ar_n_o      = ar_n_o0*fac_n_o
    ar_n_c      = ar_n_c0*fac_n_c
    ar_theta_s  = ar_theta_s0 * fac_th_s
    ar_theta_s_ = ar_theta_s0_ * fac_th_s_


#    print 'Max L=', max(ar_L)
#    print 'Max Ks=', max(ar_Ks)
#    print 'Max n_o=', max(ar_n_o)
#    print 'Max n_c=', max(ar_n_c)
#    print '\nMin L =', min(ar_L)
#    print 'Min Ks =', min(ar_Ks)
#    print 'Min n_o =', min(ar_n_o)
#    print 'Min n_c =', min(ar_n_c)


    # Check the soil imput files for validity:

    soil_min = min(ar_L)<=0. or min(ar_Ks)<=0. or min(ar_theta_r)<=0. or\
    min(ar_theta_s)<=0. or min(psi_b)<=0. or min(lamda)<=0.

    if soil_min:
        err_msg = 'One of the soil input files contains values below or equal to \
zero. Make sure the soil input rasters have the right extend and \
contain a valid value range.'
        raise ValueError(err_msg)


    #~~~~Computation of model parameters from physical parameters
    ar_Vsm, ar_b_s, ar_b_o, ar_W, ar_b_c, ar_Vsm_, ar_b_s_ \
    = pm.compute_cell_param_DL(X, ar_Xc, Dt, alpha_s, alpha_o, alpha_c,
                               nb_cell, A_thres, W_max, W_min, ar_lambda,
                               ar_tan_beta, ar_tan_beta_channel, ar_L, ar_Ks,
                               ar_theta_r, ar_theta_s, ar_n_o, ar_n_c,
                               ar_A_drained, ar_L_, ar_Ks_, ar_theta_r_,
                               ar_theta_s_)

    #~~~~Look for the cell of external_flow tunnel
    if external_flow:
        cell_external_flow = ut.find_cell_coordinates(ar_cell_label,
                                                      Xexternal_flow,
                                                      Yexternal_flow,
                                                      ar_coorx,
                                                      ar_coory,
                                                      ar_lambda)

        print 'external flows will be taken into account for cell no',\
            cell_external_flow, ' coordinates ('\
            ,Xexternal_flow,',',Yexternal_flow,')'

    #~~~~Number of simulation time steps
    if rain_distributed:
        nb_time_step = len(ndar_rain[:,0])
    else:
        nb_time_step = len(ndar_rain)


    ##=============================##
    ##  Variable array definition  ##
    ##=============================##

    ## Initialisation of the reservoirs
    #Matrix of soil,overland and channel store at the begining of the time step
    if append_output and not first_run:
        print 'Initializing from file...'
        # read from file
        h5file_in = h5.openFile(file_out, mode='r')

        node         = h5file_in.getNode('/Soil/V_s')
        ar_Vs0       = node.read()[-1, :]
        ar_Vs0_      = node.read()[-1, :]
        ar_Vs0_scott = ar_Vs0

        node = h5file_in.getNode('/Overland/V_o')
        ar_Vo0 = node.read()[-1, :]
        ar_Vo0_scott = ar_Vo0

        node = h5file_in.getNode('/Channel/V_c')
        ar_Vc0 = node.read()[-1, :]

        h5file_in.close()
    else:
        print 'Initializing from parms...'
        ar_Vs0  = fl.initial_volume_soil(ar_pVs_t0, ar_Vsm)
        ar_Vs0_ = fl.initial_volume_soil(ar_pVs_t0_, ar_Vsm_)
        ar_Vo0  = ar_Vo_t0
        ar_Vc0  = fl.initial_volume_channel(ar_Qc_t0, ar_W, X, ar_n_c)
        ar_Vd0  = np.ones(nb_cell)*ar_wwb*WL_dam_avg*Dams_ini*1e-2

        ar_Vs0_scott = ar_Vs0
        ar_Vo0_scott = ar_Vo0


    ## Computed variables
    #Matrix of soil,overland and channel store at the end of the time step
    ar_Vs1       = np.ones(nb_cell)*-99.9
    ar_Vs1_      = np.ones(nb_cell)*-99.9
    ar_Vs1_scott = np.ones(nb_cell)*-99.9
    ar_Vo1       = np.ones(nb_cell)*-99.9
    ar_Vo1_scott = np.ones(nb_cell)*-99.9
    ar_Vc1       = np.ones(nb_cell)*-99.9
    ar_Vd1       = np.ones(nb_cell)*-99.9

    #Matrix of outflows between two time steps
    ar_Qs_out  = np.ones(nb_cell)*-99.9
    ar_Qs_out_ = np.ones(nb_cell)*-99.9
    ar_Qo_out  = np.ones(nb_cell)*-99.9
    ar_Qc_out  = np.zeros(nb_cell)
    ar_Qd_inp  = np.ones(nb_cell)*-99.9
    ar_Qd_out  = np.ones(nb_cell)*-99.9

    ## Intermediate variables
    ar_a_s              = np.ones(nb_cell)*-99.9
    ar_a_s_             = np.ones(nb_cell)*-99.9
    ar_a_o              = np.ones(nb_cell)*-99.9
    ar_a_c              = np.ones(nb_cell)*-99.9
    ar_Q_to_next_cell   = np.ones(nb_cell)*-99.9
    ar_Q_to_next_cell_  = np.ones(nb_cell)*-99.9
    ar_Q_to_channel     = np.ones(nb_cell)*-99.9
    ar_Q_to_channel_sub = np.zeros(nb_cell)
    ar_Qc_cell_up       = np.zeros(nb_cell)
    ar_ETa              = np.zeros(nb_cell)
    ar_ETa_scott        = np.zeros(nb_cell)
    ar_ET_channel       = np.zeros(nb_cell)
    ar_ETc              = np.zeros(nb_cell)


    ##=============================##
    ## HDF5 output file definition ##
    ##=============================##
    h5file = h5.openFile(file_out, mode=fmode, title='TOPKAPI_out')

    root = h5file.getNode('/')
    root._v_attrs.pytopkapi_version = pytopkapi.__version__
    root._v_attrs.pytopkapi_git_revision = pytopkapi.__git_revision__

    atom = h5.Float32Atom()
    h5filter = h5.Filters(9)# maximum compression

    # create file structure as necessary
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    grp_name = '/Soil'
    if grp_name not in h5file:
        h5file.createGroup('/', 'Soil', 'Soil arrays')
    if grp_name+'/Qs_out' not in h5file:
        array_Qs_out = h5file.createEArray(grp_name, 'Qs_out',
                                           atom, shape=(0,nb_cell),
                                           title='m3/s', filters=h5filter,
                                           expectedrows=nb_time_step)
    else:
        array_Qs_out = h5file.getNode(grp_name+'/Qs_out')

    if grp_name+'/Qs_outII' not in h5file:
        array_Qs_out_ = h5file.createEArray(grp_name, 'Qs_outII',
                                           atom, shape=(0,nb_cell),
                                           title='m3/s', filters=h5filter,
                                           expectedrows=nb_time_step)
    else:
        array_Qs_out_ = h5file.getNode(grp_name+'/Qs_outII')

    if grp_name+'/V_s' not in h5file:
        array_Vs = h5file.createEArray(grp_name, 'V_s',
                                       atom, shape=(0, nb_cell),
                                       title='m3', filters=h5filter,
                                       expectedrows=nb_time_step+1)
    else:
        array_Vs = h5file.getNode(grp_name+'/V_s')

    if grp_name+'/V_sII' not in h5file:
        array_Vs_ = h5file.createEArray(grp_name, 'V_sII',
                                       atom, shape=(0, nb_cell),
                                       title='m3', filters=h5filter,
                                       expectedrows=nb_time_step+1)
    else:
        array_Vs_ = h5file.getNode(grp_name+'/V_sII')

    if grp_name+'/V_s_scott' not in h5file:
        array_Vs_scott = h5file.createEArray(grp_name, 'V_s_scott',
                                       atom, shape=(0, nb_cell),
                                       title='m3', filters=h5filter,
                                       expectedrows=nb_time_step+1)
    else:
        array_Vs_scott = h5file.getNode(grp_name+'/V_s_scott')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    grp_name = '/Overland'
    if grp_name not in h5file:
        h5file.createGroup('/', 'Overland', 'Overland arrays')
    if grp_name+'/Qo_out' not in h5file:
        array_Qo_out = h5file.createEArray(grp_name, 'Qo_out',
                                           atom, shape=(0,nb_cell),
                                           title='m3/s', filters=h5filter,
                                           expectedrows=nb_time_step)
    else:
        array_Qo_out = h5file.getNode(grp_name+'/Qo_out')

    if grp_name+'/V_o' not in h5file:
        array_Vo = h5file.createEArray(grp_name, 'V_o',
                                       atom, shape=(0,nb_cell),
                                       title='m3', filters=h5filter,
                                       expectedrows=nb_time_step+1)
    else:
        array_Vo = h5file.getNode(grp_name+'/V_o')

    if grp_name+'/V_o_scott' not in h5file:
        array_Vo_scott = h5file.createEArray(grp_name, 'V_o_scott',
                                       atom, shape=(0,nb_cell),
                                       title='m3', filters=h5filter,
                                       expectedrows=nb_time_step+1)
    else:
        array_Vo_scott = h5file.getNode(grp_name+'/V_o_scott')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    grp_name = '/Channel'
    if grp_name not in h5file:
        h5file.createGroup('/', 'Channel', 'Channel arrays')

    if grp_name+'/Qc_out' not in h5file:
        array_Qc_out = h5file.createEArray(grp_name, 'Qc_out',
                                           atom, shape=(0,nb_cell),
                                           title='m3/s', filters=h5filter,
                                           expectedrows=nb_time_step)
    else:
        array_Qc_out = h5file.getNode(grp_name+'/Qc_out')

    if grp_name+'/V_c' not in h5file:
        array_Vc = h5file.createEArray(grp_name, 'V_c',
                                       atom, shape=(0,nb_cell),
                                       title='m3', filters=h5filter,
                                       expectedrows=nb_time_step)
    else:
        array_Vc = h5file.getNode(grp_name+'/V_c')

    if grp_name+'/Ec_out' not in h5file:
        array_Ec_out = h5file.createEArray(grp_name, 'Ec_out',
                                           atom, shape=(0,nb_cell),
                                           title='m3', filters=h5filter,
                                           expectedrows=nb_time_step)
    else:
        array_Ec_out = h5file.getNode(grp_name+'/Ec_out')

    if '/ET_out' not in h5file:
        array_ET_out = h5file.createEArray('/', 'ET_out',
                                           atom, shape=(0,nb_cell),
                                           title='mm', filters=h5filter,
                                           expectedrows=nb_time_step)
    else:
        array_ET_out = h5file.getNode('/ET_out')

    if '/ET_out_scott' not in h5file:
        array_ET_out_scott = h5file.createEArray('/', 'ET_out_scott',
                                           atom, shape=(0,nb_cell),
                                           title='mm', filters=h5filter,
                                           expectedrows=nb_time_step)
    else:
        array_ET_out_scott = h5file.getNode('/ET_out_scott')

    if '/Q_down' not in h5file:
        array_Q_down = h5file.createEArray('/', 'Q_down',
                                           atom, shape=(0,nb_cell),
                                           title='m3', filters=h5filter,
                                           expectedrows=nb_time_step)
    else:
        array_Q_down = h5file.getNode('/Q_down')

    if '/Q_downII' not in h5file:
        array_Q_down_ = h5file.createEArray('/', 'Q_downII',
                                           atom, shape=(0,nb_cell),
                                           title='m3', filters=h5filter,
                                           expectedrows=nb_time_step)
    else:
        array_Q_down_ = h5file.getNode('/Q_down__')

    if '/ETc' not in h5file:
        array_ETc = h5file.createEArray('/', 'ETc',
                                           atom, shape=(0,nb_cell),
                                           title='mm', filters=h5filter,
                                           expectedrows=nb_time_step)
    else:
        array_ETc = h5file.getNode('/ETc')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    grp_name = '/Dam'
    if grp_name not in h5file:
        h5file.createGroup('/', 'Dam', 'Dam arrays')
    if grp_name+'/V_d' not in h5file:
        array_Vd = h5file.createEArray(grp_name, 'V_d',
                                       atom, shape=(0,nb_cell),
                                       title='m3', filters=h5filter,
                                       expectedrows=nb_time_step)
    else:
        array_Vd = h5file.getNode(grp_name+'/V_d')

    if grp_name+'/Qd_in' not in h5file:
        array_Qd_inp = h5file.createEArray(grp_name, 'Qd_in',
                                       atom, shape=(0,nb_cell),
                                       title='m3/s', filters=h5filter,
                                       expectedrows=nb_time_step)
    else:
        array_Qd_inp = h5file.getNode(grp_name+'/Qd_in')

    if grp_name+'/Qd_out' not in h5file:
        array_Qd_out = h5file.createEArray(grp_name, 'Qd_out',
                                       atom, shape=(0,nb_cell),
                                       title='m3/s', filters=h5filter,
                                       expectedrows=nb_time_step)
    else:
        array_Qd_out = h5file.getNode(grp_name+'/Qd_out')

    if append_output is False or first_run is True:
        #Write the initial values into the output file
        array_Vs.append(ar_Vs0.reshape((1,nb_cell)))
        array_Vs_.append(ar_Vs0_.reshape((1,nb_cell)))

        array_Vo.append(ar_Vo0.reshape((1,nb_cell)))

        array_Vs_scott.append(ar_Vs0_scott.reshape((1,nb_cell)))
        array_Vo_scott.append(ar_Vo0_scott.reshape((1,nb_cell)))

        array_Vc.append(ar_Vc0.reshape((1,nb_cell)))
        array_Vd.append(ar_Vd0.reshape((1,nb_cell)))

        array_Qs_out.append(ar_Qs_out.reshape((1,nb_cell)))
        array_Qs_out_.append(ar_Qs_out_.reshape((1,nb_cell)))
        array_Qo_out.append(ar_Qo_out.reshape((1,nb_cell)))
        array_Qc_out.append(ar_Qc_out.reshape((1,nb_cell)))
        array_Qd_inp.append(ar_Qd_inp.reshape((1,nb_cell)))
        array_Qd_out.append(ar_Qd_out.reshape((1,nb_cell)))

        array_Q_down.append( ar_Q_to_next_cell.reshape( (1,nb_cell)))
        array_Q_down_.append(ar_Q_to_next_cell_.reshape((1,nb_cell)))

        array_ET_out.append(ar_ETa.reshape((1,nb_cell)))
        array_ET_out_scott.append(ar_ETa.reshape((1,nb_cell)))

        array_ETc.append(ar_ETc.reshape((1,nb_cell)))

        E_vol = ar_ET_channel*1e-3 * ar_W * ar_Xc
        array_Ec_out.append(E_vol.reshape((1,nb_cell)))

    eff_theta  = ar_theta_s  - ar_theta_r
    eff_theta_ = ar_theta_s_ - ar_theta_r_

    ar_d_t = np.zeros((nb_time_step))
    t_avg_nr = 10

    date_time = start_date

    ##===========================##
    ##     Core of the Model     ##
    ##===========================##
    print '** NB_CELL=',nb_cell
    print '** NB_TIME_STEP=',nb_time_step
    print '--> SIMULATIONS <--'

    ## Loop on time
    for t in range(nb_time_step):
        date_time += dtm.timedelta(seconds=Dt)

        if t<=t_avg_nr:
            print t+1, '/', nb_time_step
        if t:
            d_t = t_1 - t_0
            d_t = d_t.total_seconds()
            s_r = nb_time_step - t
            ar_d_t[t] = d_t
            if t > t_avg_nr:
                t_r = np.average(ar_d_t[0:t+1])*s_r
#                print 'd_t, s_r, t_r', d_t, s_r, t_r
                dys, hrs, mns, scs = ut.seconds_to_dhms(t_r)
                now = dtm.datetime.now()
                ETA = (now + dtm.timedelta(seconds=t_r)).strftime('%H:%M:%S')
                if int(dys) > 0:
                    r_t_str = dys + ' days, '+ hrs+':' + mns + ':' + scs
                else:
                    r_t_str = hrs + ':' + mns + ':' + scs
                print t+1, '/', nb_time_step, r_t_str, ETA

        t_0 = dtm.datetime.now()

        eff_sat  = ar_Vs0 /ar_Vsm
        eff_sat_ = ar_Vs0_/ar_Vsm_

        # estimate soil suction head using Brookes and Corey (1964)
        # psi_b in mm and hence psi in mm too
        psi  = psi_b /np.power(eff_sat , 1.0/lamda)
        psi_ = psi_b_/np.power(eff_sat_, 1.0/lamda_)

        ## Loop on cells
        n=-1
        for cell1 in ar_label_sort:
            cell=np.where(ar_cell_label==cell1)[0][0]
            n=n+1
#            print 'cell:', cell

#==============================================================================
#             Drainage from upper to lower soil layer
#==============================================================================
            # disconnectedness factor = 2*lamda + 2.5
            ar_c = 2*lamda+2.5
            # Pr_prim = preliminary percolation to deeper soil layer controlled
            # by free drainage from the upper soil layer (m/s).
            Pr_prim_rate = fl.perc_prim(ar_Vs0[cell],ar_Ks[cell]*1e-3,
                                        eff_sat[cell], ar_c[cell])

#==============================================================================
#             Flow from uppter to lower soil layer
#==============================================================================
            # The flow is estimated using Green Ampt
            # Convert Pr_prim (m/s) to (mm/s)
#            print 'ar_Vs0[cell] ar_Ks[cell] eff_sat[cell] ar_c[cell]'
#            print ar_Vs0[cell], ar_Ks[cell]*1e-3, eff_sat[cell], ar_c[cell]

            Pr_prim_rate = Pr_prim_rate*1e3
#            print 'Pr_prim_rate (mm/s):', Pr_prim_rate

            infiltration_depth_ = green_ampt_cum_infiltration(Pr_prim_rate,
                                                              psi_[cell],
                                                              eff_theta_[cell],
                                                              eff_sat_[cell],
                                                              ar_Ks_[cell], Dt,
                                                              cell, t)

#            print 'infiltration_depth_', infiltration_depth_
            # Check if infiltration_depth_ is negative
            if infiltration_depth_ < 0:
                print 'WARNING: Infiltration depth from the \
green_ampt_cum_infiltration function has resulted a negative value: %0.2f mm.'\
                %(infiltration_depth_)

            # Infiltration__excess should not be necessary
#            infiltration__excess = Pr_prim_rate*1e-3 - infiltration_depth_/Dt
            # m/s = m/s - m/s

            ## ============================ ##
            ## ===== Lower SOIL STORE ===== ##
            ## ============================ ##

            ar_a_s_[cell] = fl.input_soil(infiltration_depth_, Dt, X,
                                         ar_Q_to_next_cell_, li_cell_up[cell])

            Vs_prim_ = om.solve_storage_eq(cell, ar_a_s_[cell], ar_b_s_[cell],
                                          alpha_s, ar_Vs0_[cell], Dt, solve_s)

            ar_Qs_out_[cell], ar_Vs1_[cell] = fl.output_soil(ar_Vs0_[cell],
                                                             Vs_prim_,
                                                             ar_Vsm_[cell],
                                                             ar_a_s_[cell],
                                                             ar_b_s_[cell],
                                                             alpha_s, Dt)

            return_2_shallow = max(0, ar_a_s_[cell]                   \
                               - ((ar_Vs1_[cell]-ar_Vs0_[cell])/Dt  \
                               + ar_Qs_out_[cell]))
#                               + infiltration__excess*X**2)    should not be necessary

            # Update Vs0 after drainage into second layer
            ar_Vs0[cell] = ar_Vs0[cell] - infiltration_depth_*1e-3*X**2
            ## ======================== ##
            ## ===== INTERCEPTION ===== ##
            ## ======================== ##
            ## No interception for the moment

            ## ======================== ##
            ## ===== INFILTRATION ===== ##
            ## ======================== ##
            if rain_distributed:
                rain_rate = ndar_rain[t, cell]/Dt
            else:
                rain_rate = ndar_rain[t]/Dt

            infiltration_depth = green_ampt_cum_infiltration(rain_rate,
                                                             psi[cell],
                                                             eff_theta[cell],
                                                             eff_sat[cell],
                                                             ar_Ks[cell], Dt,
                                                             cell, t)
            if infiltration_depth < 0:
                print 'WARNING: Infiltration depth from the \
green_ampt_cum_infiltration function has resulted a negative value: %0.2f mm.'\
                %(infiltration_depth)
                print 'rain_rate,psi[cell],eff_theta[cell],eff_sat[cell]\n',\
                rain_rate,',',psi[cell],',',eff_theta[cell],',',\
                eff_sat[cell],',\n',\
            ## ====================== ##
            ## ===== SOIL STORE ===== ##
            ## ====================== ##
            #~~~~ Computation of soil input
            # ar_a_s is the combination of P and O and Vs from upstream cells
            ar_a_s[cell] = fl.input_soil_shallow(infiltration_depth, Dt, X,
                                                 ar_Q_to_next_cell,
                                                 li_cell_up[cell],
                                                 return_2_shallow)

            #~~~~ Resolution of the equation dV/dt=a_s-b_s*V^alpha_s
            # Calculate the volume in the soil store at the end of the
            # current time-step.
#            if (cell == 5938 and t+1 == 188):# or (cell == 2356 and t == 61) or (cell == 2356 and t == 62):
#                print 't,cell,ar_a_s[cell],ar_b_s[cell],alpha_s,\
#                ar_Vsm[cell],ar_Vs0[cell],Dt,cell,solve_s,ar_Ks[cell],\
#                infiltration_depth,li_cell_up[cell],'
#                print t,',',cell,',',ar_a_s[cell],',',ar_b_s[cell],',',\
#                alpha_s,',',ar_Vsm[cell],',',ar_Vs0[cell],',',Dt,',',cell,',',\
#                solve_s,',',ar_Ks[cell],',',infiltration_depth,',',\
#                li_cell_up[cell]
#                h5file.close()

            Vs_prim = om.solve_storage_eq(cell, ar_a_s[cell], ar_b_s[cell],
                                          alpha_s, ar_Vs0[cell], Dt, solve_s)

            #~~~~ Computation of soil outflow and overland input
            ar_Qs_out[cell], ar_Vs1[cell] = fl.output_soil(ar_Vs0[cell],
                                                           Vs_prim,
                                                           ar_Vsm[cell],
                                                           ar_a_s[cell],
                                                           ar_b_s[cell],
                                                           alpha_s, Dt)

            if ar_Qs_out[cell] < 0 and abs(ar_Qs_out[cell]) < 1e-20:
                #ar_Qs_out[cell] = 1e-20
                ar_Qs_out[cell] = 0.0
            elif ar_Qs_out[cell] < 0 and abs(ar_Qs_out[cell]) > 1e-20:
                print 'Problem Soil: output greater than input....'
                print 'n =', n, ' ;  label =', cell
                h5file.close()
                stop

            ## ========================== ##
            ## ===== OVERLAND STORE ===== ##
            ## ========================== ##
            #~~~~ Computation of overland input
            if rain_distributed:
                rain_excess = ndar_rain[t, cell] - infiltration_depth
            else:
                rain_excess = ndar_rain[t] - infiltration_depth
            # convert mm to m^3/s
            rain_excess = max(0, (rain_excess*(10**-3)/Dt)*X**2)

            # The following calculates the saturated overland flow according to:
            # Soil inflos (P,I,O) - Storage - Soil outflow + hortonian overl. flow.
            ar_a_o[cell] = max(0,
                               ar_a_s[cell]                          \
                               - ((ar_Vs1[cell]-ar_Vs0[cell])/Dt     \
                               + ar_Qs_out[cell])                    \
                               + rain_excess)

            #~~~~ Resolution of the equation dV/dt=a_o-b_o*V^alpha_o

            ar_Vo1[cell] = om.solve_storage_eq(cell, ar_a_o[cell],
                                               ar_b_o[cell], alpha_o,
                                               ar_Vo0[cell], Dt, solve_o)

            #~~~~ Computation of overland outflows
            ar_Qo_out[cell] = fl.Qout_computing(ar_Vo0[cell], ar_Vo1[cell],
                                                ar_a_o[cell], Dt)

            if ar_Qo_out[cell] < 0 and abs(ar_Qo_out[cell]) < 1e-20:
                ar_Qo_out[cell] = 0.0
            elif ar_Qo_out[cell] < 0 and abs(ar_Qo_out[cell]) > 1e-20:
                print 'Problem Overland: output greater than input....'
                print 'ar_Qo_out=', ar_Qo_out[cell]
                print 'ar_Vo0=', ar_Vo0[cell], 'ar_Vo1=', ar_Vo1[cell], 'ar_a_o=',ar_a_o[cell]
                print 'n=', n, 'cell Nr.=', cell
                h5file.close()
                stop



            ## ============================= ##
            ## ===== FLOW PARTITIONING ===== ##
            ## ============================= ##
            # ar_Q_to_channel_sub doesn't get used for anything?
            '''
            if len(li_cell_up[cell]) < 1 and ar_lambda[cell] == 1:
                ar_lambda[cell] = 0
                print 'There is no catchment or river cell above the river \
cell %i. To prevent crashing this cell will be treated as terrestrial.' %cell
            '''
            ar_Q_to_next_cell[cell], ar_Q_to_next_cell_[cell],\
            ar_Q_to_channel[cell], ar_Q_to_channel_sub[cell] \
            = fl.flow_partitioning_DL(ar_lambda[cell], ar_Qs_out[cell],
                                      ar_Qs_out[cell], ar_Qo_out[cell],
                                      ar_W[cell], X, ar_Xc[cell])

            if ar_Q_to_channel[cell] < 0.:
                print 'ar_Q_to_channel[cell]', ar_Q_to_channel[cell]
                print 'ar_lambda[cell]', ar_lambda[cell]
                print 'ar_Qo_out[cell]', ar_Qo_out[cell]
                print 'ar_W[cell]', ar_W[cell]
                print 'X,', X
                print 'ar_Xc[cell]', ar_Xc[cell]

            # ET input data needed for dam ET, hence moved up
            if ET_distributed:
                ETr_t = ndar_ETr[t, cell]
                ETo_t = ndar_ETo[t, cell]
            else:
                ETr_t = ndar_ETr[t]
                ETo_t = ndar_ETo[t]

            ## ======================== ##
            ## ===== CHANNEL STORE ==== ##
            ## ======================== ##


            if ar_lambda[cell] == 1:
                if ar_cell_down[cell] >= 0 \
                   and ar_lambda[ar_cell_down[cell]] == 0:
                    print 'Problem: the present cell has a channel but not the\
cell down...'
                    h5file.close()
                    stop

                #~~~~ Computation of channel input
                ar_a_c[cell], \
                ar_Qc_cell_up[cell] = fl.input_channel(ar_Qc_out,
                                                       ar_Q_to_channel[cell],
                                                       li_cell_up[cell])

                if external_flow \
                and cell == np.where(ar_cell_label==cell_external_flow)[0][0]:
                    ar_a_c[cell] = ar_a_c[cell] + ar_Qexternal_flow[t]

                #~~~~ Resolution of the equation dV/dt=a_c-b_c*V^alpha_c

                ar_Vc1[cell] = om.solve_storage_eq(cell, ar_a_c[cell],
                                                   ar_b_c[cell], alpha_c,
                                                   ar_Vc0[cell], Dt, solve_c)

                #~~~~ Computation of channel outflows
                ar_Qc_out[cell] = fl.Qout_computing(ar_Vc0[cell],
                                                    ar_Vc1[cell],
                                                    ar_a_c[cell], Dt)

                if ar_Qc_out[cell] < 0:
                    print 'Problem Channel: output greater than input....'
                    print 'cell, ar_lambda[cell], ar_Qc_out[cell], ar_Vc0[cell]\
, ar_Vc1[cell], ar_a_c[cell]\n', cell, ar_lambda[cell], ar_Qc_out[cell],\
                    ar_Vc0[cell], ar_Vc1[cell], ar_a_c[cell]
                    print 'ar_Q_to_channel[cell] = ', ar_Q_to_channel[cell]
                    x_t,y_t=ut.show_cell_cords(ar_cell_label,cell,ar_coorx,ar_coory)
                    print 'this happens in cell %d with coords x=%0.2f, y=%0.2f'\
                    %(cell, x_t, y_t)
                    h5file.close()
                    stop
                if str(ar_Qc_out[cell]).count('N') > 0:
                    print ar_Qc_out[cell]
                    print 'Problem Channel: Non authorized operand....'
                    h5file.close()
                    stop

                # ~~~~ Computation of storage in Dam if present
                # The label for a dam cell is has to be >= 200 (in ar_wwb)
                ar_Qd_inp[cell] = ar_Qc_out[cell]

                if ar_wwb[cell] >= 200 and Dams_on:
                    ar_Qc_out[cell], \
                    ar_Vd1[cell] = fl.Qout_dam(Dt, WL_dam_avg,
                                               Qseep_dam_mm, ar_Vd0[cell],
                                               ar_Qd_inp[cell], ETo_t,
                                               ar_wwb[cell])
                    ar_Qd_out[cell] = ar_Qc_out[cell]
                else:
                    ar_Vd1[cell] = 0.

            else:
                ar_a_c[cell]     = 0.
                ar_Vc1[cell]     = 0.
                ar_Qc_out[cell]  = 0.


            ## ============================== ##
            ## ===== EVAPOTRANSPIRATION ===== ##
            ## ============================== ##

            #~~~~~ From soil
            ar_ETa_scott[cell], \
            ar_Vs1_scott[cell], \
            ar_Vo1_scott[cell] = em.evapot_soil_overland(ar_Vo1[cell],
                                                   ar_Vs1[cell],
                                                   ar_Vsm[cell],
                                                   ar_kc[cell],
                                                   ETr_t, X)

            # Switch for different ET-Functions:
            if ET_module_option == 'Feddes':
                Vs1_before = ar_Vs1[cell]
                Vo1_before = ar_Vo1[cell]
                ar_ETa[cell], \
                ar_Vs1[cell], \
                ar_Vo1[cell]  = em.evapot_soil_overland_feddes(ar_Vo1[cell],
                                                               ar_Vs1[cell],
                                                               ar_Vsm[cell],
                                                               ar_kc[cell],
                                                               ETr_t, ETo_t, X,
                                                               psi[cell],
                                                               ar_wwb[cell])
            elif ET_module_option == 'Standard':
                ar_ETa[cell], \
                ar_Vs1[cell], \
                ar_Vo1[cell] = em.evapot_soil_overland(ar_Vo1[cell],
                                                       ar_Vs1[cell],
                                                       ar_Vsm[cell],
                                                       ar_kc[cell],
                                                       ETr_t, X)
            elif ET_module_option == 'AVC':
                ar_ETa[cell], \
                ar_Vs1[cell], \
                ar_Vo1[cell]=em.evapot_soil_overland_feddes_AVC(ar_Vo1[cell],
                                                                ar_Vs1[cell],
                                                                ar_Vsm[cell],
                                                                ar_kc[cell],
                                                                ETr_t, ETo_t,
                                                                X, psi[cell],
                                                                ar_wwb[cell],
                                                                eff_sat[cell],
                                                                date_time,
                                                                mon_i, mon_d,
                                                                SSI_ET_treshold)




            #Test whether the soil ET function has produced valid outputs
            if ar_Vs1[cell] > ar_Vsm[cell]:
                print 'Evaporation function overestimates soil volume.'
                print 'ar_ETa_scott[cell],Vs1_before,Vo1_before,ar_Vs1_scott[cell],ar_Vo1_scott[cell]\n',\
                       ar_ETa_scott[cell],',',Vs1_before,',',Vo1_before,',',ar_Vs1_scott[cell],',',ar_Vo1_scott[cell],'\n' # ,',',fsr,'\n'
                print 'ar_ETa[cell],Vs1_before,Vo1_before,ar_Vs1[cell],ar_Vo1[cell]\n',\
                       ar_ETa[cell],',',Vs1_before,',',Vo1_before,',',ar_Vs1[cell],',',ar_Vo1[cell],'\n' #,',',fsr,'\n'
                print 'ar_kc[cell],psi[cell],ETr_t,ETo_t,X,ar_wwb[cell]\n',\
                       ar_kc[cell],',',psi[cell],',',ETr_t,',',ETo_t,',',X,',',ar_wwb[cell]#,',',ETa_s,',',ETa_o

                raise ValueError('Soil volume (%0.4f m3) exceeds porosity or\
 maximum soil volume (%0.4f m3)' %(ar_Vs1[cell], ar_Vsm[cell]))

            #~~~~~ Evaporation from channel
            #ar_ET_channel[cell] in mm of water removed from channel
            if ar_lambda[cell] == 1 and ar_wwb[cell] < 200:
                ar_ET_channel[cell], \
                ar_Vc1[cell] = em.evapor_channel(ar_Vc1[cell],
                                                 ETo_t,
                                                 ar_W[cell], ar_Xc[cell])

            #~~~~~ Potential ET (ETr)
            ar_ETc[cell] = ETr_t * ar_kc[cell]
            t_1 = dtm.datetime.now()

        ####===================================####
        #### Affectation of new vector values  ####
        ####===================================####
        ar_Vs0  = np.array(ar_Vs1)
        ar_Vs0_ = np.array(ar_Vs1_)
        ar_Vo0  = np.array(ar_Vo1)
        ar_Vc0  = np.array(ar_Vc1)
        ar_Vd0  = np.array(ar_Vd1)
        ar_Vs0_scott = np.array(ar_Vs1_scott)
        ar_Vo0_scott = np.array(ar_Vo1_scott)

        ####===================================####
        #### Results writing at each time step ####
        ####===================================####
        array_Vs.append( ar_Vs1.reshape( (1,nb_cell))) # Volume soil store
        array_Vs_.append(ar_Vs1_.reshape((1,nb_cell))) # Volume soil in 2nd layer
        array_Vo.append(ar_Vo1.reshape((1,nb_cell))) # Volume overland store
        array_Vc.append(ar_Vc1.reshape((1,nb_cell))) # Volume in Channel
        array_Vd.append(ar_Vd1.reshape((1,nb_cell))) # Volume in Dam
        array_Vs_scott.append(ar_Vs1_scott.reshape((1,nb_cell))) # Volume soil store
        array_Vo_scott.append(ar_Vo1_scott.reshape((1,nb_cell))) # Volume overland s

        array_Qs_out.append( ar_Qs_out.reshape( (1,nb_cell)))
        array_Qs_out_.append(ar_Qs_out_.reshape((1,nb_cell)))
        array_Qo_out.append(ar_Qo_out.reshape((1,nb_cell)))
        array_Qc_out.append(ar_Qc_out.reshape((1,nb_cell))) # channel outflow m3/s
        array_Qd_inp.append(ar_Qd_inp.reshape((1,nb_cell))) # dam input from channel
        array_Qd_out.append(ar_Qd_out.reshape((1,nb_cell))) # dam out to channel

        array_Q_down.append( ar_Q_to_next_cell.reshape( (1,nb_cell)))
        array_Q_down_.append(ar_Q_to_next_cell_.reshape((1,nb_cell)))

        array_ET_out.append(ar_ETa.reshape((1,nb_cell)))
        array_ET_out_scott.append(ar_ETa_scott.reshape((1,nb_cell)))

        array_ETc.append(ar_ETc.reshape((1,nb_cell)))

        E_vol = ar_ET_channel*1e-3 * ar_W * ar_Xc  # m3/Dt
        array_Ec_out.append(E_vol.reshape((1,nb_cell)))

    h5file.close()

    t_end = dtm.datetime.today()
    t_zro = dtm.datetime(2000,1,1,0,0,0)
    t_dsp = t_zro + dtm.timedelta(seconds=(t_end-t_ini).seconds)
    t_dff =  t_dsp.strftime('%H:%M:%S')

    print 'Model run %s took %s' %(os.getcwd(), t_dff)

    print ' '
    print '***** THE END *****'
    a = raw_input('Enter any value to quit. Enter initializes post processing')
    if a:
        import sys
        sys.exit()
