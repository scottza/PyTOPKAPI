import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import random
#import pytopkapi
from pytopkapi.results_analysis import Extract_flow_at_outlet, plot_Qsim_Qobs_Rain_MOD
import pytopkapi.utils as ut
import gc

#reload(pytopkapi.results_analysis)

#if __name__ == '__main__':

# For NL:
#desired =  ["NL_C03_10_24", "NL_C04_12_12", "NL_C06_12_12", "NL_C07_12_12",\
#            "NL_C10_12_12", "NL_C1112_11_03", "NL_C1314_10_24",\
#            "NL_C17_10_24", "NL_C89_10_24", "NL_iC06_2001",\
#            "NL_iC07_1993_e", "NL_iC10_1993_e", "NL_iC1112_1993_e",\
#            "NL_iC1314_1993_e", "NL_iC17_1993", "NL_iC89_1993_e",\
#            "NL_nC03_10_24", "NL_nC04_12_12", "NL_nC06_1993",\
#            "NL_nC07_1993", "NL_nC10_1993_e", "NL_nC1112_1993_e",\
#            "NL_nC1314_1993_e", "NL_nC17_1993", "NL_nC89_1993_e",\
#            "NL_iC03_1993", "NL_iC04_1993", "NL_nC04_1993","NL_nC03_1993"]

# for AMF:

desired = ['AMF_iC01_1970_daily_e', 'AMF_iC02_1970_daily_e', 'AMF_nC01_1970_daily_e',
           'AMF_nC02_1970_daily_e']

#project = 'NL'
#project = 'AMF'


def batch_Extract_flow_at_outlet(group, basepath, desired=desired,
                                 project='', WWBU_transfer_n2i=True,
                                 SD_n='', SD_i=''):
    '''
    Options for execution:
    'a' = all starting with 'Project_'
    'n' = current state (natural) starting with 'Project_n'
    'i' = impacted starting with 'Project_i'
    'c' = calibration runs starting with 'Project_C'
    '''
    os.chdir(basepath)
    dirs_a = os.listdir(os.getcwd())
    d_NL_a = []
    d_NL_n = []
    d_NL_i = []
    d_NL_C = []
    sday_a, sday_n, sday_i, sday_C = [], [], [], []
    
    if group == 'n':
        year_desired = int(SD_n.split('/')[-1])
    elif group == 'i':
        year_desired = int(SD_i.split('/')[-1])
    
    def find_year_in_string(string_in):
        iy = ''
        year = 0
        for i in string_in:
            try:
                ii =  int(i)
                iy += i
                if len(iy) == 4:
                    year = int(iy)
            except ValueError:
                iy = ''
        return year


    for d in dirs_a:
        # find the year in the dir
        year = find_year_in_string(d)
        d_splt = d.split('_')
        end5 = (d[-5:-1]+d[-1]).split('_')
        #print end5
        if d in desired:
#            print d        
            #print 'current dir is', d, 'last 5 split:', end5
            if d[0:3] == project+'_' and len(end5) >= 2:
                d_NL_a.append(d)
                if d[0:4] == project+'_C' and len(d) > 6:
                    sd_tmp = end5
                    if len(sd_tmp) < 2:
                        sday_a.append('1/7/1993')
                    else:
                        sday_a.append(sd_tmp[1] + '/' + sd_tmp[0] + '/2012')
                else:
                    sday_a.append('1/7/1993')
            if d_splt[0]+d_splt[1][0] == project+'n' and year==year_desired:
                d_NL_n.append(d)
                sday_n.append(SD_n)
            if d_splt[0]+d_splt[1][0] == project+'i' and year==year_desired:
                d_NL_i.append(d)
                sday_i.append(SD_i)
            if d_splt[0]+d_splt[1][0] == project+'C' and len(d)>6:
                d_NL_C.append(d)
                sd_tmp = end5
                if len(sd_tmp) < 2:
                        sday_C.append('1/7/1993')
                else:
                    sday_C.append(sd_tmp[1] + '/' + sd_tmp[0] + '/2012')

    if group == 'a':
        execdirs = d_NL_a
        sd       = sday_a
    if group == 'n':
        execdirs = d_NL_n
        sd       = sday_n
    if group == 'i':
        execdirs = d_NL_i
        sd       = sday_i
    if group == 'c':
        execdirs = d_NL_C
        sd       = sday_C

    for dd in xrange(len(execdirs)):
        print 'basepath: ', basepath
        print 'execdirs: ', execdirs[dd]
        os.chdir(os.path.join(basepath,execdirs[dd]+'\\run_the_model\\'))
        Extract_flow_at_outlet('plot-flow-precip.ini', start_stamp=sd[dd],
                               split_month='10', topkapi_ini='TOPKAPI.ini',
                               WWBU_transfer_n2i=WWBU_transfer_n2i)
        gc.collect()



def cal_graphs(group, basepath, desired=desired):
    
    os.chdir(basepath)
    dirs_a = os.listdir(os.getcwd())
    d_NL_a = []
    d_NL_n = []
    d_NL_i = []
    d_NL_C = []
    sday_a, sday_n, sday_i, sday_C = [], [], [], []


    for d in dirs_a:
        end5 = (d[-5:-1]+d[-1]).split('_')
        #print end5
        if d in desired:
            print d            
            #print 'current dir is', d, 'last 5 split:', end5
            if d[0:3] == 'NL_' and len(end5) >= 2:
                d_NL_a.append(d)
            if d[0:4] == 'NL_n':
                d_NL_n.append(d)
            if d[0:4] == 'NL_i':
                d_NL_i.append(d)
            if d[0:4] == 'NL_C' and len(d)>6:
                d_NL_C.append(d)

    if group == 'a':
        execdirs = d_NL_a
    if group == 'n':
        execdirs = d_NL_n
    if group == 'i':
        execdirs = d_NL_i
    if group == 'c':
        execdirs = d_NL_C

    for dd in xrange(len(execdirs)):
        print 'basepath: ', basepath
        print 'execdirs: ', execdirs[dd]
        os.chdir(os.path.join(basepath,execdirs[dd]+'\\run_the_model\\'))
        plot_Qsim_Qobs_Rain_MOD.run('plot-flow-precip.ini')


def SSI_dur_curves_plotter(path, scenario, nb_lines=12, compare=1,
                           forcing_fn=''):
    if nb_lines > 12:
        raise ValueError('Too many lines per graph. Only a maximum of 12 lines\
                          per figure are allowed...')
    plt.rc('text', usetex=True)
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'], 'size':'10'})
    
    ftype           = np.float32
#    n_rec = 7217
    
    if compare:
        AR_SSI_WWB_d_fn = path + 'AR_SSI_WWB_n'
        AR_FQC_WWB_d_fn = path + 'AR_FQC_WWB_n'
        # For NL only
#        n_rec = ut.file_len('C:\TOPKAPI\\NL_iC89_1993_e\\run_the_model\\forcing_variables\\Create_ETo\\ETo_1993.dat')
        n_rec = ut.file_len(forcing_fn)
        
        shapexy         = (n_rec, 200)
        AR_SSI_WWB_n    = np.fromfile(AR_SSI_WWB_d_fn,dtype=ftype).reshape(shapexy)
        AR_FQC_WWB_n    = np.fromfile(AR_FQC_WWB_d_fn,dtype=ftype).reshape(shapexy)
        AR_SSI_WWB_d_fn = path + 'AR_SSI_WWB_i'
        AR_FQC_WWB_d_fn = path + 'AR_FQC_WWB_i'
        AR_SSI_WWB_i    = np.fromfile(AR_SSI_WWB_d_fn,dtype=ftype).reshape(shapexy)
        AR_FQC_WWB_i    = np.fromfile(AR_FQC_WWB_d_fn,dtype=ftype).reshape(shapexy)
    elif scenario == 'n':
        AR_SSI_WWB_d_fn = path + 'AR_SSI_WWB_n'
        AR_FQC_WWB_d_fn = path + 'AR_FQC_WWB_n'
        n_rec = ut.file_len('C:\TOPKAPI\\NL_nC89_1993_e\\run_the_model\\forcing_variables\\Create_ETo\\ETo_1993.dat')
        shapexy         = (n_rec, 200)
        AR_SSI_WWB    = np.fromfile(AR_SSI_WWB_d_fn,dtype=ftype).reshape(shapexy)
        AR_FQC_WWB    = np.fromfile(AR_FQC_WWB_d_fn,dtype=ftype).reshape(shapexy)
    elif scenario == 'i':
        AR_SSI_WWB_d_fn = path + 'AR_SSI_WWB_i'
        AR_FQC_WWB_d_fn = path + 'AR_FQC_WWB_i'
        n_rec = ut.file_len('C:\TOPKAPI\\NL_iC89_1993_e\\run_the_model\\forcing_variables\\Create_ETo\\ETo_1993.dat')
        shapexy         = (n_rec, 200)
        AR_SSI_WWB    = np.fromfile(AR_SSI_WWB_d_fn,dtype=ftype).reshape(shapexy)
        AR_FQC_WWB    = np.fromfile(AR_FQC_WWB_d_fn,dtype=ftype).reshape(shapexy)

    lines   = []
    tab_leg = []
    fig_x, fig_y = 8, 5
    plt.clf()
    fig     = plt.figure()
    fig_p   = plt.gcf()
    fig_p.set_size_inches(fig_x, fig_y)
    ax1 = fig.add_subplot(111)
    ll = 0.08
    rr = 0.97
    bb = 0.08
    tt = 0.97
    fig.subplots_adjust(left=ll, right=rr, bottom=bb, top=tt)

    #Totally random colors:
    #cl =   ['648AFF','000000','02E12D','DF3C1C','FFFF00','7E0112','9FFFFF',\
#            'FFA200','268E00','DAFB14','A623B7','CA1F74','0000CD']
    
    #Distinct different:
#    cl=["FF0000", "00FF00", "0000FF", "FFFF00", "FF00FF", "00FFFF", "000000", 
#        "800000", "008000", "000080", "808000", "800080", "008080", "808080", 
#        "C00000", "00C000", "0000C0", "C0C000", "C000C0", "00C0C0", "C0C0C0", 
#        "400000", "004000", "000040", "404000", "400040", "004040", "404040", 
#        "200000", "002000", "000020", "202000", "200020", "002020", "202020", 
#        "600000", "006000", "000060", "606000", "600060", "006060", "606060", 
#        "A00000", "00A000", "0000A0", "A0A000", "A000A0", "00A0A0", "A0A0A0", 
#        "E00000", "00E000", "0000E0", "E0E000", "E000E0", "00E0E0", "E0E0E0"]
    
    #Color groups
#    cl=['FFC500','A68000','FFFA00','A6A300','FF8B00',
#        'A65A00','2419B2','352F85','100873','766FD5','',
#        '','','','','','','','']
    
#    cl=['FF9000','A65E00','FFC273','FFC900','A68200','FFE173','FF2300',
#        'A61700','FF8673','0095FF','0061A6','73C5FF']
    
    cl     = ['FF9500','661400','FFCF00','5C7C00','FF2800',
              'BF123A','085EAB','01294C','006561','33CDC7']
#    cl_n   = ['1921B1','009999','00BD39','95EC00']
#    cl_i   = ['FFBC00','FF7400','FF2300','D2006B']B4F200
    cl_n   = ['A60000','25D500','024A68','FF7400','2F3485','A66F00']
    cl_i   = ['FF4040','59EA3A','225E79','FF9640','7277D8','BF8F30']

    #ls   = [':',':','-.','-.','-','-','--','--']
    ls   = ['-','--','-','--','-','--','-','--']
    ls_n = '-'
    ls_i = '--'

    # NL
#    exclusions = [4,5,6,11,24,30,31,33,34,35,37,46,51,52,53,54,55,56,57,58,59,60,61,
#                  62,63,68,72,74,75,76,78,86,87,88,89,93,94,95,97]
#    w_splitter = [2,3,10,14,20,23,27,32,40,45,50,67,73,81,85,96]

    
    WWBU_mined_fn_e = path+'\\..\\_Forcing_Data\\WWBU\\WWBU_mined.dat'
    WWBU_mined_fn_s = path+'\\..\\_Forcing_Data\\WWBU\\SSI_splitter.dat'
    exclusions = np.loadtxt(WWBU_mined_fn_e, dtype=np.int)
    w_splitter = np.loadtxt(WWBU_mined_fn_s, dtype=np.int)
    w_splitter = list(w_splitter)

    if compare:
        nb_lines = nb_lines/2
        AR_SSI_WWB = AR_SSI_WWB_n

    cnt = -1
    wwb_u_names = ''
    nb_cols = AR_SSI_WWB.shape[1]
    ccnt = -1
    gr_cnt = 0

    for i in xrange(nb_cols):
        if i+1 == nb_cols:
#            cnt = nb_lines+1
            print 'last column...', 'i=',i,' Forcing to create last graph'
            w_splitter.append(i)
        if AR_SSI_WWB[0,i] <= 0:
           #print AR_SSI_WWB[0,i]
            pass
        
        if AR_SSI_WWB[0,i] > 0 and i not in exclusions:
            print 'adding WWB unit lines', i, 'to figure.'
            ccnt += 1
            cnt  += 1
            wwb_u_names += str(i)+'_'
            #rgb         = '#%02X%02X%02X' % (tuple([random.randint(0,255) for x in xrange(3)]))
            if compare:
                SSI_descnd_n  = AR_SSI_WWB_n[:,i]
                frqcy_n       = AR_FQC_WWB_n[:,i]
                SSI_descnd_i  = AR_SSI_WWB_i[:,i]
                frqcy_i       = AR_FQC_WWB_i[:,i]
                SSI_descnd_n[SSI_descnd_n<0] = np.nan
                frqcy_n[frqcy_n<0]           = np.nan
                SSI_descnd_i[SSI_descnd_i<0] = np.nan
                frqcy_i[frqcy_i<0]           = np.nan
                tab_leg.append(r'$\mathbf{Unit~%s~Current~State}$' %str(i))
                if np.sum(np.isnan(SSI_descnd_i)) == len(SSI_descnd_i):
                    print 'unit %d lost and set to zero..' %i
                    SSI_descnd_i = np.ones(len(SSI_descnd_i))*5e-1
                    frqcy_i = np.linspace(0,100, num=len(SSI_descnd_i))
                    tab_leg.append(r'$\mathbf{Unit~%s~Lost}$' %str(i))
                else:
                    tab_leg.append(r'$\mathbf{Unit~%s~Impacted}$' %str(i))
                print ccnt
                lines += ax1.plot(frqcy_n, SSI_descnd_n, ls=ls_n, color='#'+cl_n[ccnt], linewidth=1.5)
                lines += ax1.plot(frqcy_i, SSI_descnd_i, ls=ls_i, color='#'+cl_i[ccnt], linewidth=1.5)
                
            else:
                SSI_descnd  = AR_SSI_WWB[:,i]
                frqcy       = AR_FQC_WWB[:,i]
                SSI_descnd[SSI_descnd<0] = np.nan
                frqcy[frqcy<0]           = np.nan
                #print 'plotting wwb unit', i, 'ccnt is', ccnt
                lines      += ax1.plot(frqcy, SSI_descnd, ls=ls[ccnt], color='#'+cl[ccnt], linewidth=1)
                if scenario == 'n':
                    l = ' Current State'
                elif scenario == 'i':
                    l = ' Impacted'
                tab_leg.append(r'$\mathbf{Unit~%s}$' %str(i)+l)

        if AR_SSI_WWB[0,i] > 0:
#            if cnt < nb_lines:    # this case skipps the save graph things
            current_graph = i
            if current_graph in w_splitter or ccnt+1 == nb_lines:
                print 'making graph. Count:',cnt, 'i:',i
                gr_cnt += 1
                if gr_cnt > 3:
#                    break
                    pass
                ax1.legend(lines, tab_leg, loc='upper right', fancybox=True,
                           fontsize=11)
                leg     = ax1.get_legend()
                leg.get_frame().set_alpha(0.75)
#                ymax_V  = ax1.get_ylim()[1]
#                ymin_V  = ax1.get_ylim()[0]
                ymin_V, ymax_V = 0, 100
                xmin_V  = ax1.get_xlim()[0]
                xmax_V  = ax1.get_xlim()[1]
                for v in [10,40,60,90]:
                    plt.axvline(v, color='k', ls=':')
                ax1.set_xlim(xmin_V,xmax_V)
                ax1.set_ylim(ymin_V,ymax_V)
                #ax1.set_yscale('log')
                ax1.set_ylabel(r'$\mathbf{Soil~Saturation~Index~}(\%)$', fontweight='bold')
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
                yu = (ymax_V-ymin_V)*1e-2
                yl = (ymax_V-ymin_V)*2e-2
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
                ut.check_folder_exist(path+'\\10\\')
                if compare:
                    fig_fn = path+'\\10\\N_vs_I_'+wwb_u_names+'SSI_Duration.png'
                else:
                    fig_fn = path+'\\10\\'+scenario+'_'+wwb_u_names+'SSI_Duration.png'
            
                plt.savefig(fig_fn, dpi=200)
                print 'figure', fig_fn, 'saved...'
                
                #___Now clear current figure and create new one..___#
                cnt = -1
                ccnt = -1
#                f_nb = i+1
                lines   = []
                tab_leg = []
                plt.clf()
                fig     = plt.figure()
                fig_p   = plt.gcf()
                fig_p.set_size_inches(fig_x, fig_y)
                ax1     = fig.add_subplot(111)
                fig.subplots_adjust(left=ll, right=rr, bottom=bb, top=tt)
                wwb_u_names = ''

    a = raw_input('Press Enter to Quit')
    print a


def plot_bars(path):
    '''
    Plots the water balace (WB) components of all wetland units for the
    current state as well as for the impact scenario.
    This function reads columns for all WWB units from one large np array
    saved on disk. The array is arranged as such that each WWB unit is
    stored in one seperate column. Each row contains the various water balance
    component. The different scenarios and units are specified below.
    '''
    
    sn_m3 = 0     # Start row index for natural runs in m3.
    en_m3 = 13    # End  row index for natural runs in m3.
    sn_mm = 26    # Start row index for natural runs in mm.
    en_mm = 39    # End row index for natural runs in mm.
    si_m3 = 13    # Start row index for impact runs in m3.
    ei_m3 = 26    # End row index for impact runs in m3.
    si_mm = 39    # Start row index for impact runs in mm.
    ei_mm = 52    # End row index for impact runs in  mm.
    
    plt.rc('text', usetex=True)
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'], 'size':'8'})
    
    AR_WWB_y_fn = path + 'AR_WWB_y'
    print '..Opening file:', AR_WWB_y_fn
    
    ftype = np.float32
    cols  = 200
    rows  = 14*2*2
    xy    = (rows, cols)
    
    AR_WWB_y    = np.fromfile(AR_WWB_y_fn, dtype=ftype).reshape(xy)
    print 'shape of AR_WWB_y:', AR_WWB_y.shape
    
#    coords = np.loadtxt(path+'Coords_Barplots.txt',delimiter=',',\
#                        skiprows=1, usecols=(2,3))
#    #print coords
#    res    = 3.5

    path_out = path + '\\WWB_barplots\\'
        
    ut.check_folder_exist(path_out)
    
    WWBU_mined_fn_e = path+'\\..\\_Forcing_Data\\WWBU\\WWBU_mined.dat'
    exclusions = np.loadtxt(WWBU_mined_fn_e)
    
    width  = 0.35
    
    # HGM for New Largo:
#    hgm    = ['','CVB','HS','CVB','HS','HS','HS','HS','CVB','HS','CVB','HS',\
#              'CVB','HS','HS','CVB','CVB','CVB','HS','HS','HS','CVB','HS','HS',\
#              'HS','CVB','HS','HS','HS','CVB','CVB','HS','HS','HS','CVB','HS',\
#              'HS','CVB','HS','HS','CVB','HS','HS','CVB','HS','HS','CVB','HS',\
#              'HS','CVB','HS','HS','HS','HS','PAN','PAN','PAN','HS','PAN','PAN',\
#              'DAM','DAM','CVB','HS','HS','HS','HS','HS','HS','CVB','HS','HS',\
#              'HS','HS','CVB','HS','HS','CVB','HS','HS','HS','HS','HS','HS',\
#              'HS','HS','HS','CVB','HS','HS','CVB','HS','HS','HS','PAN','HS',\
#              'CVB','CVB']
    
    # HGM input changed to read it in as text file
    
    HGM_fn = path+'\\..\\_Forcing_Data\\WWBU\\HGM.dat'
    with open(HGM_fn, 'r') as f:
        hgm = f.read().splitlines()
    
    
    co     = ['#0090FF','#D4471E','#BF9A30','#F2D88F','#40BD18','#94E77A','#004DFF','#00174D']
    #co     = ['#3515B0','#FF8B00','#FF2800','#B4F200','#086CA2','#FFD100','#BE008A','#00BB3F']
    #xl     = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i','j','k','l','m']
    xl     = [r'$\mathbf{Rain}$',r'$\mathbf{ET_a}$',r'$\mathbf{Qs_{in}}$',\
              r'$\mathbf{Qs_{out}}$',r'$\mathbf{Qo_{in}}$',r'$\mathbf{Qo_{out}}$',\
              r'$\mathbf{Riv_{Qs}}$',r'$\mathbf{Riv_{Qo}}$',r'$\mathbf{Rain}$']
    x0     = np.arange(13)
    #x1 = [x for x in xrange(1,14)]
    w_idx = -1

    for w in xrange(200):
        if w > 200:
            break
        if AR_WWB_y[0,w] > 0 and w not in exclusions:
            w_idx += 1
            print '..plotting WWB unit', w
            s = 9   # slicer
            n_m3 = np.delete(AR_WWB_y[sn_m3:en_m3,w][0:s],1)
            i_m3 = np.delete(AR_WWB_y[si_m3:ei_m3,w][0:s],1)
            n_mm = np.delete(AR_WWB_y[sn_mm:en_mm,w][0:s],1)
            i_mm = np.delete(AR_WWB_y[si_mm:ei_mm,w][0:s],1)
            #print '1st:', n_m3

            #now the data should look like:
            #ETa,Qs_in,Qs_out,Qo_in,Qo_out,Qs_in_r,Qo_in_r,Rain
            #Insert rain into first column, then delete last value
            ri    = 7    # rainfall index
            n_m3 = np.insert(n_m3,0,n_m3[ri])[0:s-1]
            i_m3 = np.insert(i_m3,0,i_m3[ri])[0:s-1]
            n_mm = np.insert(n_mm,0,n_mm[ri])[0:s-1]
            i_mm = np.insert(i_mm,0,i_mm[ri])[0:s-1]
            #print '2st:', n_m3

            ss   = len(n_m3)
            co   = co[0:ss]
            xl   = xl[0:ss]
            x0   = x0[0:ss]
            #print xl

            n_m3[n_m3<0] = 0.
            i_m3[i_m3<0] = 0.
            n_mm[n_mm<0] = 0.
            i_mm[i_mm<0] = 0.
            
            n_m3[n_m3==0.] = 0.001
            i_m3[i_m3==0.] = 0.001
            n_mm[n_mm==0.] = 0.001
            i_mm[i_mm==0.] = 0.001

            #print len(n_m3), len(i_m3), len(n_mm), len(i_mm)

            #barplot_mm(x0, n_mm, i_mm, width, xl, co, w, path_out, hgm[w])

            barplot_mm_small(x0, n_mm, i_mm, width, xl, co, w, path_out, hgm[w])
            
#            fn_wdf = path_out+'\\TIF\\sWWB_unit_'+str(w)+'.tfw'

#            create_world_file(coords[w_idx][0], coords[w_idx][1], res, fn_wdf)




def barplot_mm(x0,a,b,width,xl,color,w,path_out,hgm):
    
    fig, ax = plt.subplots()
    
    fig_x, fig_y = 2.8, 2.8
    fig_p   = plt.gcf()
    fig_p.set_size_inches(fig_x, fig_y)
    
    bar_n_mm = ax.bar(x0,       a, width, color=color)
    bar_i_mm = ax.bar(x0+width, b, width, color=color, hatch='/')
    bar_f2   = ax.bar(0, 0, 0, color='w', hatch='/')
    bar_f3   = ax.bar(0, 0, 0, color='w', edgecolor='w', lw=0)
    
    ax.set_ylabel(r'$\mathbf{Wetland~Water~Balance~(mm)}$')
    ax.set_title(r'$\mathbf{Unit~%s}$' %str(w))
    ax.set_xticks(x0+width)
    ax.set_xticklabels(xl)
    
    #ax.legend((bar_n_mm[0], bar_i_mm[0]), ('a', 'b'))
    
    fig_fn = path_out+'WWB_unit_'+str(w)+'.png'
    s = plt.savefig(fig_fn, dpi=150)




def barplot_mm_small(x0, a, b, width, xl, color, w, path_out, hgm):
    
    fig, ax = plt.subplots()
    
    fig_x, fig_y = 4, 2.7
    fig_p   = plt.gcf()
    fig_p.set_size_inches(fig_x, fig_y)
    
    #fig.subplots_adjust(left=0.085, right = 0.97, bottom = 0.14, top = 0.97)
    bar_n_mm = ax.bar(x0,       a, width, color=color)
    bar_i_mm = ax.bar(x0+width, b, width, color=color, hatch='/')
    bar_f1   = ax.bar(0, 0, 0, color='w')
    bar_f2   = ax.bar(0, 0, 0, color='w', hatch='/')
    bar_f3   = ax.bar(0, 0, 0, color='w', edgecolor='w', lw=0)
    
    ax.set_ylabel(r'$\mathbf{Annual~WWB~Components~(mm)}$')
#    ax.set_title(r'$\mathbf{Unit~%s}$' %str(w))
    
    # New largo:
    # 8000 for VB's
    # 4000 for HS's
    
    if hgm == 'HS':        
        ax.set_ylim(0,2000)
    elif hgm == 'CVB' or hgm == 'UVB':
        ax.set_ylim(0,5000)
    elif hgm == 'PAN':
        ax.set_ylim(0,2000)
    ax.set_xticks(x0+width, )
    ax.set_xticklabels(xl)
    yticklabels = []
    for iy in ax.get_yticks():
        yticklabels.append(r'$\mathbf{%s}$' % int(iy))
    ax.set_yticklabels(yticklabels)
    ax.tick_params()
    
#    ax.text(.5,0.93,r'$\mathbf{Unit~%s}$' %str(w), horizontalalignment='center',
#            transform=ax.transAxes, fontsize=8,
#            bbox=dict(facecolor='w',alpha=0.5))
    
    ll,rr,bb,tt = 0.14, 0.97, 0.07, 0.97
    fig.subplots_adjust(left=ll, right=rr, bottom=bb, top=tt)
    
#    axx = plt.gca()
#    wxt = axx.get_window_extent()
    #print wxt.width
    #wwt = wxt.width
    #wht = wxt.height
    #print wxt
    ccc = b[2:len(b)-1]
    if np.average(ccc) > 0.008:
        i_label = r'$\mathbf{Impacted}$'
    else:
        i_label = r'$\mathbf{Lost}$'
    n_label = r'$\mathbf{Unimpacted}$'
    title   = r'$\mathbf{%s~Unit~%s}$' %(hgm, str(w))
    
    if hgm == 'HS':
        hgm_c = '#3EE400'
    elif hgm == 'CVB' or hgm == 'UVB':
        hgm_c = '#03600F'
    elif hgm == 'PAN':
        hgm_c = '#00C7FF'
    else:
        hgm_c = 'k'
    
    legend = ax.legend((bar_f1[0], bar_f2[0]),(n_label, i_label),fancybox=True,
                       title=title, bbox_to_anchor=(0.47,0.97))
    # bbox can be a BboxBase instance, a tuple of [left, bottom, width, height]
    # in the given transform (normalized axes coordinate if None), or a tuple of
    # [left, bottom] where the width and height will be assumed to be zero
    
    plt.setp(legend.get_title(),fontsize='10',color=hgm_c)
    leg = ax.get_legend()
    leg.get_frame().set_alpha(0.75)
    
    # for AMF where one graph is peking 12K mm of Qo
    def autolabel(rects):
    # attach some text labels
        for rect in rects:
#            print rect
            height = rect.get_height()
            y_max  = ax.get_ylim()[1]
            if height > y_max:
                print w
                print rect
                ax.text(rect.get_x()+rect.get_width()/2., int(y_max*0.82),
                        r'$\mathbf{%d}$'%int(height),ha='center',
                        va='bottom',rotation=90,size=8,
                        bbox=dict(boxstyle='rarrow,pad=0.1',lw=0.5,
                                  fc='w',alpha=0.7,ec='gray'))
    autolabel(bar_n_mm)
    
    autoAxis = ax.axis()
    #print 'autoAxis', autoAxis
    wwt_m = autoAxis[1]-autoAxis[0]
    wht_m = autoAxis[3]-autoAxis[2]
    wwt  = wwt_m/(1-ll-(1-rr))
    wht  = wht_m/(1-bb-(1-tt))
    #print 'window width, height:', wwt, wht
    #rec = plt.Rectangle((autoAxis[0]-0.7,autoAxis[2]-0.2),(autoAxis[1]-autoAxis[0])+1,(autoAxis[3]-autoAxis[2])+0.4,fill=False,lw=1)
    rec = mpl.patches.Rectangle((0-ll*wwt,0-bb*wht),width=wwt,height=wht,
                                fill=False,lw=1.5)
    rec = ax.add_patch(rec)
    rec.set_clip_on(False)
    
    # to find out which file formats for savefig are supported:
    #print plt.gcf().canvas.get_supported_filetypes_grouped()
    
    fig_fn = path_out+'sWWB_unit_'+str(w)+'.png'
    s = plt.savefig(fig_fn, dpi=300)

#    fig_fn_tif = path_out+'\\TIF\\sWWB_unit_'+str(w)+'.tif'
#    plt.savefig(fig_fn_tif, dpi=300)



def create_world_file(x, y, res, fn_wdf):
    
    "Write a world file for each image for automated georeferencing of the\
    bar plots."
    '''
    Line 1: A: pixel size in the x-direction in map units/pixel
    Line 2: D: rotation about y-axis
    Line 3: B: rotation about x-axis
    Line 4: E: pixel size in the y-direction in map units, almost always negative[3]
    Line 5: C: x-coordinate of the center of the upper left pixel
    Line 6: F: y-coordinate of the center of the upper left pixel
    '''    
    
#    a_out = np.zeros((1,6))
    lines = [res,0,-0,res*-1,int(x),int(y)]
    #frmts = ['%0.2f %i %i %0.2f %i %i']
    frmts = '%0.2f'
    print type(lines)
    print lines
    
    print type(frmts)
    print frmts

#    for l in xrange(len(lines)):
#        a_out[0,l] = lines[l]

    np.savetxt(fn_wdf,lines,fmt=frmts)









