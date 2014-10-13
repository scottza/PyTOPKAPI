import datetime as dt
from ConfigParser import SafeConfigParser

import numpy as np
import tables as h5
import matplotlib.pyplot as plt
from matplotlib.dates import date2num, DateFormatter
import datetime as dtm
from matplotlib import rc, font_manager

import pytopkapi.utils as ut

def run(ini_file = 'plot-flow-precip.ini'):
    rc('text', usetex=False)
    ar_Qsim, ar_Qobs = flux_plot(ini_file)
    write_validation_results2file(ar_Qsim, ar_Qobs)
    flux_plot_calibration(ini_file)


def flux_plot(ini_file):
    rc('text', usetex=False)
    config = SafeConfigParser()
    config.read(ini_file)
    print 'Read the file ', ini_file

    file_Qsim = config.get('files','file_Qsim')
    file_Qobs = config.get('files','file_Qobs')
    file_rain = config.get('files','file_rain')
    image_out = config.get('files','image_out')

    group_name=config.get('groups','group_name')

    outlet_ID = config.getint('parameters', 'outlet_ID')
    graph_format = config.get('parameters', 'graph_format')

    Qobs       = config.getboolean('flags','Qobs')
    Pobs       = config.getboolean('flags','Pobs')
    nash       = config.getboolean('flags','nash')
    R2         = config.getboolean('flags','R2')
    RMSE       = config.getboolean('flags','RMSE')
    RMSE_norm  = config.getboolean('flags','RMSE_norm')
    Diff_cumul = config.getboolean('flags','Diff_cumul')
    Bias_cumul = config.getboolean('flags','Bias_cumul')
    Err_cumul  = config.getboolean('flags','Err_cumul')
    Abs_cumul  = config.getboolean('flags','Abs_cumul')

    # Read the run version number and add it to graph file name
    config.read('run_version.ini')
    ver_n = config.getint('version_number', 'version')

    if ver_n < 100 and ver_n >= 10:
        ver_n = '0'+ str(ver_n)
    elif ver_n < 10:
        ver_n = '00' + str(ver_n)
    else: ver_n = str(ver_n)

    image_out = image_out + '_' + str(ver_n) + '.' + graph_format

    tab_col=['k','r']
    tab_style=['-','-']
    tab_width=['1','1']
    color_P='b'
    transparency_P = 0.5#(0 for invisible)

    #create path_out if it does'nt exist
    ut.check_file_exist(image_out)

    #Read the obs
    #Qobs
    ar_date, ar_Qobs = read_observed_flow(file_Qobs)

    delta = date2num(ar_date[1]) - date2num(ar_date[0])

    #Rain
    if Pobs:
        h5file_in=h5.openFile(file_rain,mode='r')
        group='/'+group_name+'/'
        node = h5file_in.getNode(group+'rainfall')
        ndar_rain=node.read()
        h5file_in.close()
        #Compute the mean catchment rainfall
        ar_rain=np.average(ndar_rain,axis=1)

    #Read the simulated data Q
    file_h5=file_Qsim
    ndar_Qc_out=ut.read_one_array_hdf(file_h5,'/Channel/','Qc_out')
    ar_Qsim=ndar_Qc_out[1:,outlet_ID]
    #print len(ndar_Qc_out[1:,outlet_ID])

    ##Graph
    fig, ax = plt.subplots()

    fig_x = 14
    fig_y = 11

    lines = []
    tab_leg = []

    if Qobs:
        lines += ax.plot(ar_date, ar_Qobs,
                         color=tab_col[-1],
                         linestyle=tab_style[-1], linewidth=tab_width[-1])
        tab_leg.append(('Observation'))
        tab_leg = tab_leg[::-1]

    lines += ax.plot(ar_date, ar_Qsim,
                     color=tab_col[0],
                     linestyle=tab_style[0], linewidth=tab_width[0])
    tab_leg.append('Model')

    if nash:
        nash_value = ut.Nash(ar_Qsim,ar_Qobs)
        lines += ax.plot(ar_date[0:1], ar_Qsim[0:1], 'w:')
        tab_leg.append(('Eff       = '+str(nash_value)[0:5]))

    if R2:
        r2_value = ut.R2(ar_Qsim,ar_Qobs)
        lines += ax.plot(ar_date[0:1], ar_Qsim[0:1], 'w:')
        tab_leg.append('R^2       = '+str(r2_value)[0:5])

    if RMSE:
        RMSE_value = ut.RMSE(ar_Qsim,ar_Qobs)
        lines += ax.plot(ar_date[0:1], ar_Qsim[0:1], 'w:')
        tab_leg.append('RMSE      = '+str(RMSE_value)[0:5])

    if RMSE_norm:
        RMSE_norm_value = ut.RMSE_norm(ar_Qsim,ar_Qobs)
        lines += ax.plot(ar_date[0:1], ar_Qsim[0:1], 'w:')
        tab_leg.append('RMSE_norm = '+str(RMSE_norm_value)[0:5])

    if Diff_cumul:
        Diff_cumul_value = ut.Diff_cumul(ar_Qsim,ar_Qobs)
        lines += ax.plot(ar_date[0:1], ar_Qsim[0:1], 'w:')
        tab_leg.append('Diff_cum  = '+str(Diff_cumul_value)[0:5])

    if Bias_cumul:
        Bias_cumul_value = ut.Bias_cumul(ar_Qsim,ar_Qobs)
        lines += ax.plot(ar_date[0:1], ar_Qsim[0:1], 'w:')
        tab_leg.append('Bias_cum  = '+str(Bias_cumul_value)[0:5])

    if Err_cumul:
        Err_cumul_value = ut.Err_cumul(ar_Qsim,ar_Qobs)
        lines += ax.plot(ar_date[0:1], ar_Qsim[0:1], 'w:')
        tab_leg.append('Err_cum   = '+str(Err_cumul_value)[0:5])

    if Abs_cumul:
        Abs_cumul_value = ut.Abs_cumul(ar_Qsim,ar_Qobs)
        lines += ax.plot(ar_date[0:1], ar_Qsim[0:1], 'w:')
        tab_leg.append('Abs_cum   = '+str(Abs_cumul_value)[0:5])

    ax.set_xlim(ar_date[0], ar_date[-1])
    ytitle=r'$Q \  (m^3/s)$'
    ax.set_ylabel(ytitle, fontsize=18)
    ax.set_title('Version ' + ver_n)

    ax2 = ax.twinx()

    ax2.set_ylabel(r'$Rainfall~(mm)$', fontsize=18, color=color_P)
    ax2.bar(ar_date, ar_rain, width=delta,
            facecolor='blue', edgecolor='blue', alpha=transparency_P)
    ax2.set_ylim(max(ar_rain)*2, min(ar_rain))

    ax2.legend(lines, tab_leg, loc='center right', fancybox=True)
    leg = ax2.get_legend()
    leg.get_frame().set_alpha(0.75)

    # rotate and align the tick labels so they look better,
    # unfortunately autofmt_xdate doesn't work with twinx due to a bug
    # in matplotlib <= 1.0.0 so we do it manually
    ## fig.autofmt_xdate()

#    cf = plt.gcf()
    fig.set_size_inches(fig_x, fig_y)

    bottom = 0.08
    top    = 0.97
    right  = 0.96
    left   = 0.062
    rotation=30
    ha='right'

    for ax in fig.get_axes():
        if hasattr(ax, 'is_last_row') and ax.is_last_row():
            for label in ax.get_xticklabels():
                label.set_ha(ha)
                label.set_rotation(rotation)
        else:
            for label in ax.get_xticklabels():
                label.set_visible(False)
            ax.set_xlabel('')

    fig.subplots_adjust(left=left, bottom=bottom, right=right, top=top)

    plt.savefig(image_out)
    #plt.show()
    return ar_Qsim, ar_Qobs


def read_observed_flow(file_name):
    """Read the observed flow from a data file.

    """
    date = np.loadtxt(file_name, dtype=np.int, usecols=(0, 1, 2, 3, 4))
    dates = [dt.datetime(yr, mon, dy, hr, mn) for yr, mon, dy, hr, mn in date]

    Q = np.loadtxt(file_name, usecols=(5,))

    return dates, Q

def write_validation_results2file(ar_Qsim, ar_Qobs, topkapi_ini_file = 'TOPKAPI.ini'):


    config = SafeConfigParser()
    config.read(topkapi_ini_file)
    fn_change_log_out = config.get('output_files', 'file_change_log_out')

    f_out = open(fn_change_log_out, 'a')

    values = []
    values.append(ut.Nash(ar_Qsim,ar_Qobs))
    values.append(ut.R2(ar_Qsim, ar_Qobs))
    values.append(ut.RMSE(ar_Qsim, ar_Qobs))
    values.append(ut.RMSE_norm(ar_Qsim, ar_Qobs))
    values.append(ut.Diff_cumul(ar_Qsim, ar_Qobs))
    values.append(ut.Bias_cumul(ar_Qsim, ar_Qobs))
    values.append(ut.Err_cumul(ar_Qsim, ar_Qobs))
    values.append(ut.Abs_cumul(ar_Qsim, ar_Qobs))

    names = ['nash_value', 'r2_value', 'RMSE_value', 'RMSE_norm_value', 'Diff_cumul_value', 'Bias_cumul_value', 'Err_cumul_value', 'Abs_cumul_value']

    for v in values:
        f_out.write(names[np.where(values==v)[0][0]] + ' = ' + str(round(v,3)) + '\n')

    f_out.write('______________________________\n')
    f_out.close()


def flux_plot_calibration(ini_file):
    #rc('text.latex', preamble=r'\usepackage{cmbright}')
    #sizeOfFont = 11
    #fontProperties = {'family':'sans-serif','sans-serif':['Arial'], 'weight':'normal', 'size':sizeOfFont}
    #ticks_font = font_manager.FontProperties(family='Arial', style='normal', size=sizeOfFont, weight='normal', stretch='normal')
    rc('text', usetex=True, fontsize=11)
    #rc('font', **fontProperties)

    config = SafeConfigParser()
    config.read(ini_file)
    print 'Read the file ',ini_file

    file_Qsim = config.get('files','file_Qsim')
    file_Qobs = config.get('files','file_Qobs')
    file_rain = config.get('files','file_rain')
    image_out = config.get('files','image_out')

    group_name=config.get('groups','group_name')

    outlet_ID    = config.getint('parameters', 'outlet_ID')
    graph_format = config.get('parameters', 'graph_format')
    start_cal    = config.get('parameters', 'start_calibration')

    Qobs       = config.getboolean('flags','Qobs')
    Pobs       = config.getboolean('flags','Pobs')
    nash       = config.getboolean('flags','nash')
    R2         = config.getboolean('flags','R2')
    RMSE       = config.getboolean('flags','RMSE')
    RMSE_norm  = config.getboolean('flags','RMSE_norm')
    Diff_cumul = config.getboolean('flags','Diff_cumul')
    Bias_cumul = config.getboolean('flags','Bias_cumul')
    Err_cumul  = config.getboolean('flags','Err_cumul')
    Abs_cumul  = config.getboolean('flags','Abs_cumul')

    # Read the run version number and add it to graph file name
    config.read('run_version.ini')
    ver_n = config.getint('version_number', 'version')

    if ver_n < 100 and ver_n >= 10:
        ver_n = '0'+ str(ver_n)
    elif ver_n < 10:
        ver_n = '00' + str(ver_n)
    else: ver_n = str(ver_n)

    image_out = image_out[0:-3] + '/Calibration/PPQ_' + str(ver_n) + '_CAL.' + graph_format

    tab_col=['k','r']
    tab_style=['-','-']
    tab_width=['1','1']
    color_P='b'
    transparency_P = 0.5#(0 for invisible)

    #create path_out if it does'nt exist
    ut.check_file_exist(image_out)

    #Read the obs
    #Qobs
    ar_date, ar_Qobs = read_observed_flow(file_Qobs)

    start_d = []
    for s in start_cal.split('/'):
        start_d.append(int(s))

    start_date = dtm.datetime(start_d[2],start_d[1],start_d[0])
    print 'Calibration period commences on %s' %(start_date.strftime('%d %b \'%y'))

    cc = 0    # commencing calibration

    for d in ar_date:
        #print start_date, '==', d
        #print 'cc =', cc
        cc += 1
        if start_date == d:
            break

#    print cc

    #del ar_date[:] does not work for arrays
    delta = date2num(ar_date[1]) - date2num(ar_date[0])

    #Rain
    if Pobs:
        h5file_in=h5.openFile(file_rain,mode='r')
        group='/'+group_name+'/'
        node = h5file_in.getNode(group+'rainfall')
        ndar_rain=node.read()
        h5file_in.close()
        #Compute the mean catchment rainfall
#        print ndar_rain.shape
        ar_rain=np.average(ndar_rain[cc:],axis=1)

    #Read the simulated data Q
    file_h5=file_Qsim
    ndar_Qc_out=ut.read_one_array_hdf(file_h5,'/Channel/','Qc_out')
    ar_Qsim=ndar_Qc_out[1:,outlet_ID]
    ar_Qsim = ar_Qsim[cc:]
    #print len(ndar_Qc_out[1:,outlet_ID])

    ar_Qobs = ar_Qobs[cc:]

    ##Graph
    plt.clf()
    fig, ax = plt.subplots()

    fig_x = 9
    fig_y = 4

    lines   = []
    tab_leg = []

#    print len(ar_date[cc:])
#    print len(ar_Qobs)
#    print len(ar_Qsim)
#    print len(ar_rain)

    if Qobs:
        lines += ax.plot(ar_date[cc:], ar_Qobs,
                         color=tab_col[-1],
                         linestyle=tab_style[-1], linewidth=tab_width[-1])
        tab_leg.append(('Observation'))
        tab_leg = tab_leg[::-1]

    lines += ax.plot(ar_date[cc:], ar_Qsim,
                     color=tab_col[0],
                     linestyle=tab_style[0], linewidth=tab_width[0])
    tab_leg.append('Model')

#    if nash:
#        nash_value = ut.Nash(ar_Qsim,ar_Qobs)
#        lines += ax.plot(ar_date[cc:cc+1], ar_Qsim[0:1], 'w:')
#        tab_leg.append(('Eff       = '+str(nash_value)[0:5]))

    if R2:
        r2_value = ut.R2(ar_Qsim,ar_Qobs)
        lines += ax.plot(ar_date[cc:cc+1], ar_Qsim[0:1], 'w:')
        tab_leg.append((r'$\mathbf{R^2 = %s}$' % str(r2_value)[0:4]))
#
#    if RMSE:
#        RMSE_value = ut.RMSE(ar_Qsim,ar_Qobs)
#        lines += ax.plot(ar_date[cc:cc+1], ar_Qsim[0:1], 'w:')
#        tab_leg.append(('RMSE      = '+str(RMSE_value)[0:5]))
#
#    if RMSE_norm:
#        RMSE_norm_value = ut.RMSE_norm(ar_Qsim,ar_Qobs)
#        lines += ax.plot(ar_date[cc:cc+1], ar_Qsim[0:1], 'w:')
#        tab_leg.append(('RMSE_norm = '+str(RMSE_norm_value)[0:5]))

    if Diff_cumul:
        Diff_cumul_value = ut.Diff_cumul(ar_Qsim,ar_Qobs)
        lines += ax.plot(ar_date[cc:cc+1], ar_Qsim[0:1], 'w:')
        tab_leg.append((r'$\mathbf{Cumulative~Difference = %s (m^3/s)}$' %str(Diff_cumul_value)[0:4]))

#    if Bias_cumul:
#        Bias_cumul_value = ut.Bias_cumul(ar_Qsim,ar_Qobs)
#        lines += ax.plot(ar_date[0:1], ar_Qsim[0:1], 'w:')
#        tab_leg.append(('Bias_cum  = '+str(Bias_cumul_value)[0:5]))

    if Err_cumul:
        Err_cumul_value = ut.Err_cumul(ar_Qsim,ar_Qobs)
        lines += ax.plot(ar_date[cc:cc+1], ar_Qsim[0:1], 'w:')
        tab_leg.append(r'$\mathbf{Cumulative~Error = %s (%s)}$'\
        %(str(Err_cumul_value)[0:4], r'\%'))

#    if Abs_cumul:
#        Abs_cumul_value = ut.Abs_cumul(ar_Qsim,ar_Qobs)
#        lines += ax.plot(ar_date[0:1], ar_Qsim[0:1], 'w:')
#        tab_leg.append(('Abs_cum   = '+str(Abs_cumul_value)[0:5]))


    ax.set_xlim(ar_date[cc], ar_date[-1])
    ytitle=r'$\mathbf{Q~(m^3/s)}$'
    ax.set_ylabel(ytitle, fontsize=13)
    #ax.set_title('Version ' + ver_n)

    ax2 = ax.twinx()

    ax2.set_ylabel(r'$\mathbf{Rainfall~(mm)}$', fontsize=13, color=color_P)
    ax2.bar(ar_date[cc:], ar_rain, width=delta,
            facecolor='blue', edgecolor='blue', alpha=transparency_P)
    ax2.set_ylim(max(ar_rain)*2, min(ar_rain))

    ax2.legend(lines, tab_leg, loc='center right', fontsize = 12, fancybox=True)
    leg = ax2.get_legend()
    leg.get_frame().set_alpha(0.75)

    # rotate and align the tick labels so they look better,
    # unfortunately autofmt_xdate doesn't work with twinx due to a bug
    # in matplotlib <= 1.0.0 so we do it manually
    ## fig.autofmt_xdate()

#    cf = plt.gcf()
    fig.set_size_inches(fig_x, fig_y)

    bottom   = 0.13
    top      = 0.98
    right    = 0.95
    left     = 0.07
    rotation = 30
    ha = 'right'

    formatter = DateFormatter('%d %b \'%y')
    plt.gcf().axes[0].xaxis.set_major_formatter(formatter)

    fig.canvas.draw()

    tmp = []
#    for l in ax2.get_yticks():
#        print l
    for label in ax2.get_yticklabels():
#        print label.get_text().replace('$','')
        tmp.append(r'$\mathbf{%s}$' %label.get_text().replace('$',''))
    ax2.set_yticklabels(tmp)
        #label.set_fontproperties(ticks_font)

    tmp2 = []
    for l in ax.get_yticks():
        tmp2.append(r'$\mathbf{%s}$' %l)
        #label.set_fontproperties(ticks_font)
    ax.set_yticklabels(tmp2)
#
#    for label in ax2.get_yticklabels():
        #label.set_fontproperties(ticks_font)

    for ax in fig.get_axes():
        if hasattr(ax, 'is_last_row') and ax.is_last_row():
            for label in ax.get_xticklabels():
                #label.strftime('%d %b \'%y')
                label.set_ha(ha)
                label.set_rotation(rotation)
        else:
            for label in ax.get_xticklabels():
                label.set_visible(False)
            ax.set_xlabel('')

    fig.subplots_adjust(left=left, bottom=bottom, right=right, top=top)

    plt.savefig(image_out, dpi=180)
    #plt.show()
