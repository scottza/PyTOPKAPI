import datetime as dtm
from ConfigParser import SafeConfigParser

import numpy as np
import tables as h5
import matplotlib.pyplot as plt
from matplotlib.dates import date2num
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA

import pytopkapi.utils as ut
import pytopkapi.pretreatment as pm

def run(ini_file='plot_fluxes.ini', start_date = '1/1/1910', TKP_ini = 'TOPKAPI.ini'):
    config = SafeConfigParser()
    config.read(ini_file)
    print 'Read the file ',ini_file

#    file_Qsim=config.get('files','file_Qsim')
    file_results     = config.get('files', 'file_results')
#    file_Qobs        = config.get('files','file_Qobs')
    file_rain        = config.get('files', 'file_rain')
    image_out        = config.get('files', 'image_out')
    parameter_fn     = config.get('files', 'file_cell_param')

    group_name = config.get('groups', 'group_name')
    
    #outlet_ID = config.getint('parameters', 'outlet_ID')
    graph_format = config.get('parameters', 'graph_format')
    
    rain = config.getboolean('flags', 'rain')
    ET   = config.getboolean('flags', 'ET')  # ETc = potential crop ref. ET (ETc = ETr * Kc)
    Vs   = config.getboolean('flags', 'Vs')
    Qs   = config.getboolean('flags', 'Qs')


    config.read(TKP_ini)
    print 'Read the file ',TKP_ini
    
    fn_global_param = config.get('input_files', 'file_global_param')
    
    #~~~~Read Global parameters file
    X, Dt, alpha_s, \
    alpha_o, alpha_c, \
    A_thres, W_min, W_max = pm.read_global_parameters(fn_global_param)

    #~~~~Read Cell parameters file
    ar_cell_label, ar_coorx, \
    ar_coory, ar_lambda, \
    ar_Xc, ar_dam, \
    ar_tan_beta, ar_tan_beta_channel, \
    ar_L0, ar_Ks0, \
    ar_theta_r, ar_theta_s0, \
    ar_n_o0, ar_n_c0, \
    ar_cell_down, ar_pVs_t0, \
    ar_Vo_t0, ar_Qc_t0, \
    ar_kc, psi_b, lamda = pm.read_cell_parameters(parameter_fn)


    # Read the run version number and add it to graph file name
    config.read('run_version.ini')
    ver_n = config.getint('version_number', 'version')
    
    if ver_n < 100 and ver_n >= 10:
        ver_n = '0'+ str(ver_n)
    elif ver_n < 10:
        ver_n = '00' + str(ver_n)
    else: ver_n = str(ver_n)
    
    image_out = image_out + '_' + str(ver_n) + '.' + graph_format

#   Order        ETa      ETc         Vs        Qs
    tab_col=['#FF8000','#FFBF00', '#C56600', '#31007F']
    tab_style=['solid', 'dashed', 'solid', 'dotted']
    tab_width=['1', '1', '1', '1']
    color_P='b'
    transparency_P = 0.5

    #create path_out if it does'nt exist
    ut.check_file_exist(image_out)

    #Rain
    if rain:
        h5file_in = h5.openFile(file_rain, mode='r')
        group = '/'+group_name+'/'
        node = h5file_in.getNode(group + 'rainfall')
        ndar_rain = node.read()
        h5file_in.close()
        #Compute the mean catchment rainfall
        ar_rain = np.average(ndar_rain, axis=1)
        print 'ar_rain.shape', ar_rain.shape
    
    if ET:
        h5file_in = h5.openFile(file_results, mode = 'r')
        #group = '/'+group_name_r+'/'
        node = h5file_in.getNode('/', 'ETc')
        ndar_ETc = node.read()
        np.array(ndar_ETc)
        # Comput the mean ETc in the catchment
        ar_ETc = np.average(ndar_ETc,axis=1)
        ar_ETc = ar_ETc[1:]
        
        # Now load the ETa data
        node = h5file_in.getNode('/', 'ET_out')
        ndar_ETa = node.read()
        h5file_in.close()
        np.array(ndar_ETc)
        # Compute the mean ETa in the catchment
        ar_ETa = np.average(ndar_ETa,axis=1)
        ar_ETa = ar_ETa[1:]
    
    if Vs:
        h5file_in = h5.openFile(file_results, mode = 'r')
        node = h5file_in.getNode('/Soil/', 'V_s')
        ndar_Vs = node.read()
        np.array(ndar_Vs)
        h5file_in.close()
        ar_Vs = np.average(ndar_Vs,axis=1)
        ar_Vs = ar_Vs[1:]
    
    if Vs:
        h5file_in = h5.openFile(file_results, mode = 'r')
        node = h5file_in.getNode('/Soil/', 'Qs_out')
        ndar_Qs = node.read()
        np.array(ndar_Qs)
        h5file_in.close()
        ar_Qs = np.average(ndar_Qs,axis=1)
        ar_Qs = ar_Qs[1:]

    #Read the simulated data ETa-data
#    file_h5 = file_results
#    ndar_ETa_out = ut.read_one_array_hdf(file_h5,'ET_out')
#    ar_ETa_out = np.ones((ndar_ETa_out.shape[0])) * -99.9
#    for l in xrange(ndar_ETa_out.shape[0]):
#        ar_ETa_out[l] = np.average(ndar_ETa_out[l, :])
    #ndar_Qc_out = ut.read_one_array_hdf(file_h5,'/Channel/','Qc_out')
    #ar_Qsim = ndar_Qc_out[1:,outlet_ID]
    #print len(ndar_Qc_out[1:,outlet_ID])

    ar_date = []
    n_rec = len(ar_ETa)
    
    ar_sd = start_date.split('/')
    
    if len(ar_sd[2]) == 4 and len(ar_sd[1]) < 3 and len(ar_sd[0]) < 3:
        pass
    else:
        raise ValueError('The date format for the start date has to be \
entered in the format \'DD/MM/YYYY\'.\
 ')

    stmp = []
    for s in ar_sd:
        stmp.append(int(s))        
    start_date = dtm.datetime(stmp[2],stmp[1],stmp[0])
    
    print 'start date:', start_date
    
    curr_dt = start_date
    for d in xrange(n_rec):
        curr_dt += dtm.timedelta(seconds=Dt)
        ar_date.append(curr_dt)

    #Read the obs
    #Qobs
# Used to load the date from the Q_obs file. However if now flow data is
# present then this doesn't work. Hence date will be generated with starting date.
#    ar_date, ar_Qobs = read_observed_flow(file_Qobs)

    delta = date2num(ar_date[1]) - date2num(ar_date[0])

#===============================================================================
#     Graph
#===============================================================================

    #fig, ax = plt.subplots()
    plt.clf()
    ax = host_subplot(111, axes_class=AA.Axes)
    #fig = plt.gcf()

    lines = []
    tab_leg = []
    if ET:
        lines += ax.plot(ar_date, ar_ETa,
                         color=tab_col[0],
                         linestyle=tab_style[0], linewidth=tab_width[0])
        tab_leg.append('ETa')
        tab_leg = tab_leg[::-1]

        lines += ax.plot(ar_date, ar_ETc,
                         color=tab_col[1],
                         linestyle=tab_style[1], linewidth=tab_width[1])
        tab_leg.append('ETc')

    '''
    if nash:
        nash_value = ut.Nash(ar_Qsim,ar_Qobs)
        lines += ax.plot(ar_date[0:1], ar_Qsim[0:1], 'w:')
        tab_leg.append(('Eff = '+str(nash_value)[0:5]))
    '''
    
    #ax.set_xlim(ar_date[0], ar_date[-1])
    ytitle=r'$Evapotranspiration \ (mm)$'
    ax.set_ylabel(ytitle, fontsize=25)
    ax.set_title('Version ' + ver_n)

    ax2 = ax.twinx()
    ax2.set_ylabel(r'$Rainfall \ (mm)$', fontsize=25, color=color_P)
    ax2.bar(ar_date, ar_rain, width=delta,
            facecolor='blue', edgecolor='blue', alpha=transparency_P)
    ax2.set_ylim(max(ar_rain)*2, min(ar_rain))
    
    if Vs:
        ax3 = ax.twinx()
        offset = 50
        new_fixed_axis3 = ax3.get_grid_helper().new_fixed_axis
        ax3.axis['right'] = new_fixed_axis3(loc='right', axes=ax3, offset=(offset, 0))
        ax3.axis['right'].toggle(all=True)
        lines += ax3.plot(ar_date, ar_Vs, color=tab_col[2], linestyle=tab_style[2], linewidth=tab_width[2])
        tab_leg.append('Vs')
        ax3.set_xlim(ar_date[0], ar_date[-1])
        ytitle=r'$Volume \ soil \ (m^3)$'
        ax3.set_ylabel(ytitle, fontsize=25)

    if Qs:
        ax4 = ax.twinx()
        new_fixed_axis4 = ax4.get_grid_helper().new_fixed_axis
        ax4.axis['right'] = new_fixed_axis4(loc='right', axes=ax4, offset=(100, 0))
        ax4.axis['right'].toggle(all=True)
        lines += ax4.plot(ar_date, ar_Qs,
                         color=tab_col[3],
                         linestyle=tab_style[3], linewidth=tab_width[3])
        tab_leg.append(r'$Qs$')
        ax4.set_ylabel(r'$Soil~flow~(m^3/s)$')
        #print max(ar_Qs), min(ar_Qs)

    ax2.legend(lines, tab_leg, loc='center right', fancybox=True)
    leg = ax2.get_legend()
    leg.get_frame().set_alpha(0.75)

    # rotate and align the tick labels so they look better,
    # unfortunately autofmt_xdate doesn't work with twinx due to a bug
    # in matplotlib <= 1.0.0 so we do it manually
    ## fig.autofmt_xdate()

    fig_x = 14
    fig_y = 11

    cf = plt.gcf()
    cf.set_size_inches(fig_x, fig_y)

    bottom   = 0.08
    top      = 0.97
    right    = 0.85
    left     = 0.082
    rotation = 30
    ha='right'

    for ax in cf.get_axes():
        if hasattr(ax, 'is_last_row') and ax.is_last_row():
            for label in ax.get_xticklabels():
                label.set_ha(ha)
                label.set_rotation(rotation)
        else:
            for label in ax.get_xticklabels():
                label.set_visible(False)
            ax.set_xlabel('')

    cf.subplots_adjust(left=left, bottom=bottom, right=right, top=top)

    plt.savefig(image_out)
    #plt.show()

def read_observed_flow(file_name):
    """Read the observed flow from a data file.

    """
    date = np.loadtxt(file_name, dtype=np.int, usecols=(0, 1, 2, 3, 4))
    dates = [dtm.datetime(yr, mon, dy, hr, mn) for yr, mon, dy, hr, mn in date]

    Q = np.loadtxt(file_name, usecols=(5,))

    return dates, Q