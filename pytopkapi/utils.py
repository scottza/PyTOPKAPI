import os
import os.path
from subprocess import Popen, PIPE

import tables as h5
import numpy as np
import traceback
import h5py
import pytopkapi
import pytopkapi.pretreatment as pm
import numpy as np
import datetime as dtm
#from pytopkapi.parameter_utils.create_file_04_2013_Kc import read_raster

from ConfigParser import SafeConfigParser, NoOptionError

# System utility functions
def exec_command(cmd_args):
    """Execute a shell command in a subprocess

    Convenience wrapper around subprocess to execute a shell command
    and pass back stdout, stderr, and the return code. This function
    waits for the subprocess to complete, before returning.

    Usage example:
    >>> stdout, stderr, retcode = exec_command(['ls', '-lhot'])

    Parameters
    ----------
    cmd_args : list of strings
        The args to pass to subprocess. The first arg is the program
        name.

    Returns
    -------
    stdout : string
        The contents of stdout produced by the shell command
    stderr : string
        The contents of stderr produced by the shell command
    retcode : int
        The return code produced by the shell command

    """
    proc = Popen(cmd_args, stdout=PIPE, stderr=PIPE)
    stdout, stderr = proc.communicate()
    proc.wait()

    return stdout, stderr, proc.returncode

########################
##   For graphics     ##
########################
def CRange(ar_x):
    '''
    Returns the range of an array
    '''
    lim=np.array([min(ar_x),max(ar_x)])
    return lim

def f_axe(p,xc):
    '''
    Returns the value in the array xc
    associated to a relative value inside [0,1]
    '''
    xc=np.sort(xc)
    pos=xc[0]+p*(xc[-1]-xc[0])
    return pos

def string(integer,len_str_out):
    """
    From a given integer, return an string of length len_str_out completed by zero
    Example:
    ut.string(1,3)-->'001'
    """
    str_zero='0'
    str_int=str(integer)
    len_int=len(str_int)
    if len_str_out-len_int<0:
        print '****ERROR: length of string too short'
        str_out=''
    else:
        str_out=(len_str_out-len_int)*str_zero+str_int
    return str_out

def from_float_array_to_string_array(ar_float,unique=False):
    if unique:
        a=str(np.unique(ar_float)).split()[1:]
    else:
        a=str(ar_float).split()[1:]
    a[-1]=a[-1].replace(']','')
    ar_string=a
    return ar_string

##############################
##   For file management    ##
##############################

def check_file_exist(filename):
    path_name, file_name = os.path.split(filename)
    if not os.path.exists(path_name) and path_name != '':
        print path_name, 'has been created'
        os.makedirs(path_name)

def check_folder_exist(folder_name):
    if not os.path.exists(folder_name):
        print folder_name, 'has been created'
        os.mkdir(folder_name)

def read_one_array_hdf(file_h5,group,name):
    h5file_in=h5.openFile(file_h5,mode='r')
    node = h5file_in.getNode(group+name)
    array=node.read()
    h5file_in.close()
    return array

##############################
##        Statistics        ##
##############################

def mov_avg(ar_float,period):
    '''
    period is a multiple of 2
    '''

    nb_ind=len(ar_float)-period
    ar_out=np.zeros(nb_ind)
    for i in range(nb_ind):
        n=period/2
        ind_mid=i+n
        ar_out[i]=np.average(ar_float[ind_mid-n:ind_mid+n])

    return ar_out

##~~~   Comparison of 2 vectors   ~~~~##
# Functions defining useful criteria comparing two vectors
## REFERENCE is Y
def R(ar_x,ar_y):
    R=np.corrcoef(ar_x,ar_y)
    return R[0,1]

def R2(ar_x,ar_y):
    R=np.corrcoef(ar_x,ar_y)
    return R[0,1]**2

def Nash(ar_x,ar_y):
    eff=1-sum((ar_y-ar_x)**2)/sum((ar_y-np.mean(ar_y))**2)
    return eff

def RMSE(ar_x,ar_y):
    rmserr=(np.mean((ar_y-ar_x)**2))**0.5
    return rmserr

def RMSE_norm(ar_x,ar_y):
    rmserr=(np.mean((ar_y-ar_x)**2))**0.5
    rmsenorm=rmserr/np.mean(ar_y)
    return rmsenorm

def Bias_cumul(ar_x,ar_y):
    b=sum(ar_x)/sum(ar_y)
    return b

def Diff_cumul(ar_x,ar_y):
    diff=sum(ar_x)-sum(ar_y)
    return diff

def Abs_cumul(ar_x,ar_y):
    abs_diff=abs(sum(ar_x)-sum(ar_y))
    return abs_diff

def Err_cumul(ar_x,ar_y):
    err_rel=abs(sum(ar_x)-sum(ar_y))/sum(ar_y)
    return err_rel


##############################
##        HDF5 files        ##
##############################

####HOW to remove a group
#h5file.removeNode('/', groupname)

##############################
##      Works on vectors    ##
##############################

def find_dist_max(ar_coorx,ar_coory):
    """
    Compute the maximum distance between several points defined by their coordinates ar_coorx and ar_coory
    """
    nb_cell=len(ar_coorx)
    max_dist=0.
    for i in range(nb_cell):
        for j in range(nb_cell):
            max_dist=max(max_dist,distance(ar_coorx[i],ar_coory[i],ar_coorx[j],ar_coory[j]))
    return max_dist

def distance(x1,y1,x2,y2):
    """
    Compute the distance between two points
    """
    dist=((x1-x2)**2+(y1-y2)**2)**0.5
    return dist

def find_cell_coordinates(ar_cell_label, Xoutlet, Youtlet,
                          ar_coorx, ar_coory, ar_lambda, channel=True):
    """Find the label of the cell closest to (Xoutlet, Youtlet).

    Find the label of the model cell containing the specified location. The
    co-ordinates of the location must be given in the same co-ordinate system
    as that specifying the model catchment.

    Parameters
    ----------
    ar_cell_label : (N,) int array
        Numbers labelling each cell.
    Xoutlet : float
        The x co-ordinate of a point. This is the Longitude expressed in
        metres using the same projection as `ar_coorx`.
    Youtlet : float
        The y co-ordinate of a point. This is the Longitude expressed in
        metres using the same projection as `ar_coory`.
    ar_coorx : (N,) float array
        The x co-ordinate of the centre of each cell (m). This is the Longitude
        expressed in metres using a Transverse Mercator projection, but any
        appropriate projection can be used.
    ar_coory : (N,) float array
        The y co-ordinate of the centre of each cell (m). This is the Latitude
        expressed in metres using a Transverse Mercator projection, but any
        appropriate projection can be used.
    ar_lambda : (N,) int array
        Switch indicating whether the current cell contains a channel. A value
        of `1` indicates a channel cell, `0` indicates no channel.
    channel : boolean (default=True)
        Allows cells with or without channels to be chosen.

    Returns
    -------
    cell_outlet : int
        The label for the cell closest to the defined location.

    """
    tab_x=np.unique(ar_coorx);X=abs(tab_x[0]-tab_x[1])
    dist_max=3*X
    dist_min=dist_max
    nb_cell=len(ar_cell_label)
    cell_outlet=-999.9
    for i in range(nb_cell):
        dist=distance(Xoutlet,Youtlet,ar_coorx[i],ar_coory[i])
        if channel:
            if dist < dist_min and ar_lambda[i]==1.:
                dist_min=dist
                cell_outlet=ar_cell_label[i]
        else:
            if dist<dist_min:
                dist_min=dist
                cell_outlet=ar_cell_label[i]


    if cell_outlet<0:
        print "Wrong coordinates"
        stop
    return cell_outlet

def show_cell_cords(ar_cell_label, cell, ar_coorx, ar_coory):
    cor_x = []
    cor_y = []
    cor_x.append(ar_coorx[np.where(ar_cell_label == cell)])
    cor_y.append(ar_coory[np.where(ar_cell_label == cell)])
    return cor_x[0], cor_y[0]

def pv(var):
    (filename,line_number,function_name,text)=traceback.extract_stack()[-2]
    print('%s: %s'%(text[text.find('(')+1:-1],var))

def CFF(fn_in_r, fn_in_e, len_cells, rain_distributed=0, TKP_ini='TOPKAPI.ini',
        CRF_ini='..\\create_the_parameter_files\\create_file.ini',
        NB_periods = 1, splits_i=(0)):
    """ CFF stands for: Create Forcing Files Checks if forcing files (rain and ET) exist, if not it creates them.

     Parameters
     ----------
     fn_in:
     filename of input rainfall file (string)

     len_cells:
     number of cells of the catchment (int)

     rain_distributed:
     bool to switch between same rainfall in the catchment
     or different rainfall for each cell. The rain cell parameter is in the
     GIS_bin_files folder
     
     Returns
     -------
     Nothing
    """
    
    # Test wether the rain or ET input files are windows shortcuts or real files
    
#    for fn in os.listdir(os.path.join(os.getcwd(),os.path.split(fn_in_r)[0])):
#        if fn.endswith('.lnk'):
#            if os.path.split(fn_in_r)[1] + '.lnk' == fn:
#                print 'Rainfall input file is a shortcut: ', fn
#                from win32com.client import Dispatch
#                shell = Dispatch('WScript.Shell')
#                shortcut = shell.CreateShortCut(fn)
#                fn_in_r = shortcut.Targetpath
#                print fn_in_r
#                print 'The target %s of the shortcut will be used' %(fn_in_r)

#    for fn in os.listdir(os.path.join(os.getcwd(),os.path.split(fn_in_e)[0])):
#        if fn.endswith('.lnk'):
#            if os.path.split(fn_in_e)[1] + '.lnk' == fn:
#                print 'Rainfall input file is a shortcut: ', fn
#                from win32com.client import Dispatch
#                shell = Dispatch('WScript.Shell')
#                shortcut = shell.CreateShortCut(fn)
#                fn_in_e = shortcut.Targetpath
#                print 'The target %s of the shortcut will be used' %(fn_in_e)
    

    print 'current directory: ', os.getcwd()
    # check if the forcing files already exist
    fn_h5_r = 'forcing_variables/rainfields.h5'
    fn_h5_e = 'forcing_variables/ET.h5'
    force_d = 'forcing_variables/'
    
    if rain_distributed:
        config = SafeConfigParser()
        config.read(CRF_ini)
        #print config.has_section('raster_files')
        try:
            if NB_periods == 2:
                rain_fname1 = config.get('raster_files', 'rain_fname1')
                rain_fname2 = config.get('raster_files', 'rain_fname2')
            else:
                rain_fname = config.get('raster_files', 'rain_fname')
            mask_fname = config.get('raster_files', 'mask_fname')
        except NoOptionError:
            print 'No raster file name was specified in the create file ini-file \
\'%s\'. Therefore distributed rainfall was switched off.'%(CRF_ini)
            rain_distributed = False
    if rain_distributed:
        if NB_periods == 2:
            rain_fname1 = '..\\'+rain_fname1
            rain_fname2 = '..\\'+rain_fname2
        else:
            rain_fname = '..\\'+rain_fname
        mask_fname = '..\\'+mask_fname
        
        config.read(TKP_ini)
        try:
            param_fname = config.get('input_files', 'file_cell_param')
        except NoOptionError:
            raise ValueError('The cell parameter file does not have an entry \
for \'file_cell_param\'.Therefore distributed rainfall was switched off.')
            rain_distributed = False
    if rain_distributed:
        print os.getcwd()
        print rain_fname1,', ',rain_fname2
        if NB_periods == 2:
            try:
                with open(rain_fname1): pass
            except IOError:
               raise ValueError('Rain distributed but raster file \'%s\' \
does not exist.'%(rain_fname1.split('/')[-1]))
            try:
                with open(rain_fname2): pass
            except IOError:
               raise ValueError('Rain distributed but raster file \'%s\' \
does not exist.'%(rain_fname2.split('/')[-1]))
        else:
            try:
                with open(rain_fname): pass
            except IOError:
               raise ValueError('Rain distributed but raster file \'%s\' \
does not exist.'%(rain_fname.split('/')[-1]))
        '''
        ar_cell_label, ar_coorx, \
        ar_coory, ar_lambda, \
        ar_Xc, ar_wwb, \
        ar_tan_beta, ar_tan_beta_channel, \
        ar_L0, ar_Ks0, \
        ar_theta_r, ar_theta_s0, \
        ar_n_o0, ar_n_c0, \
        ar_cell_down, ar_pVs_t0, \
        ar_Vo_t0, ar_Qc_t0, \
        ar_kc, psi_b, lamda = pm.read_cell_parameters(param_fname)
        '''
        from pytopkapi.parameter_utils.create_file_04_2013_Kc import read_raster
        ar_mask = read_raster(mask_fname)
        
        if NB_periods == 2:
            ar_rain1 = read_raster(rain_fname1)
            ar_rain2 = read_raster(rain_fname2)
            ar_rain1 = ar_rain1[ar_mask==1]
            ar_rain2 = ar_rain2[ar_mask==1]
        else:
            ar_rain = read_raster(rain_fname)
            ar_rain = ar_rain[ar_mask==1]
#        ar_rain = pytopkapi.parameter_utils.create_file_04_2013_Kc.read_raster(rain_fname)
#        ar_mask = pytopkapi.parameter_utils.create_file_04_2013_Kc.read_raster(mask_fname)

    if fn_h5_r.split('/')[-1] not in os.listdir(os.path.join(os.getcwd(),force_d)):
        if rain_distributed:
#            rain_forcing_dist(fn_in_r, ar_rain, ar_cell_label,splits_i)
            if NB_periods == 2:
                rain_forcing_dist(fn_in_r,[ar_rain1,ar_rain2],NB_periods,splits_i)
            else:
                rain_forcing_dist(fn_in_r, ar_rain,splits_i)
        else:
            rain_forcing(fn_in_r, len_cells)
        print 'Rain forcing file created:\n%s' %fn_h5_r
    else:
        print 'Rain forcing file already exists..'
    if fn_h5_e.split('/')[-1] not in os.listdir(os.path.join(os.getcwd(),force_d)):
        et_forcing(fn_in_e, len_cells)
        print 'ET forcing file created:\n%s' %fn_h5_e
    else:
        print 'ET forcing file already exists..'

def rain_forcing(fn_in, len_cells):
    # Writes the rainfal forcing data into a hd5 file
    # Input:
    # fn_in: filename of input rainfall file (string)
    # len_cells: number of cells of the catchment (int)
    
    #### Output file
    fn_h5 = 'forcing_variables/rainfields.h5'
    with open(fn_in, 'r') as f_in:
        rain = []
        for l in f_in:
            rain.append(float(l))

    raindata = []
    for i in range(len(rain)):
        line_tmp = []
        for c in range(len_cells):
            line_tmp.append(rain[i])
        raindata.append(line_tmp)
    rain = np.array(raindata)

    ### Create the H5 file ###
    f_h5 = h5py.File(fn_h5, 'w')
    sub_group = f_h5.create_group('sample_event')    
    sub_group.create_dataset('rainfall', data=rain)    
    f_h5.close()


#def rain_forcing_dist(fn_in_r, ar_rain, ar_cell_label):
def rain_forcing_dist(fn_in_r, ar_rain, NB_periods=1, splits_i=(0)):
    print 'rain data file is being split at index value %i' %(splits_i)
    # Writes the rainfal forcing data into a hd5 file
    # Input:
    # fn_in_r: filename of input rainfall file (string)
    # ar_rain: array of the rain raster file
    # ar_cell_label: array of the naming of all cells
    
#    if NB_periods > 1:
        # specify the split lines for the various periods
        #splits_i = 634
        #158
#        cols1    = [0,1]
#        cols2    = [1,2,3,4]
    
    P_rec = np.loadtxt(fn_in_r,dtype=np.float32,delimiter='\t',skiprows=1, ndmin=2)
    print 'Shape of rain input file:', P_rec.shape
    
    if NB_periods > 1:
        Len_cells = len(ar_rain[0])
    else:
        Len_cells = len(ar_rain)
    Len_P_rec = len(P_rec[:,0])
    print 'Len_cells:', Len_cells
    print 'Len_P_rec:', Len_P_rec
    
    ndar_rain = np.zeros((Len_P_rec,Len_cells),dtype=np.float32)
    
    #for _ in xrange(NB_periods):
    '''
    The column indices of the rainfall records (P_rec) are assigned to every
    cell in the ar_rain raster files. The splits_i will switch between the
    tow (in this case) different periods: ar_rain[0] and ar_rain[1]. Each rain
    value will be then written to ndar_rain[y,x] which is the final ndarray
    holding all reanfall records (y) for every cell (x).
    '''
    for y in xrange(splits_i):
        for x in xrange(Len_cells):
            ndar_rain[y,x] = P_rec[y, int(ar_rain[0][x]-1)]
            P_tmp1 = P_rec[y, int(ar_rain[0][x]-1)]
            if P_tmp1 == -999:
                print '-999 discovered. y, x, ar_rain[0][x]-1:', y, x, ar_rain[0][x]
                stop
    for y in xrange(splits_i, Len_P_rec):
        for x in xrange(Len_cells):
            ndar_rain[y,x] = P_rec[y, int(ar_rain[1][x]-1)]
            P_tmp2 = P_rec[y, int(ar_rain[1][x]-1)]
            if P_tmp2 == -999:
                print '-999 discovered. y, x, ar_rain[0][x]-1:', y, x, ar_rain[1][x]
                stop

    #### Output file
    fn_h5 = 'forcing_variables/rainfields.h5'

    ### Create the H5 file ###
    f_h5 = h5py.File(fn_h5, 'w')
    sub_group = f_h5.create_group('sample_event')
    sub_group.create_dataset('rainfall', data=ndar_rain)    
    f_h5.close()

def et_forcing(fn_in, len_cells):
    #### Output file
    fn_h5 = 'forcing_variables/ET.h5'
    with open(fn_in, 'r') as f_in:
        et = []
        for l in f_in:
            et.append(float(l))

    etdata_o = []
    etdata_r = []    
    for i in range(len(et)):
        line_tmp = []
        for c in range(len_cells):
            line_tmp.append(et[i])
        etdata_r.append(line_tmp)

    et_r = np.array(etdata_r)
    for i in range(len(et)):
        line_tmp = []
        for c in range(len_cells):
            line_tmp.append(et[i]/1.2)
        etdata_o.append(line_tmp)
    
    et_o = np.array(etdata_o)
    
    ### Create the H5 file ###
    f_h5 = h5py.File(fn_h5, 'w')
    sub_group = f_h5.create_group('sample_event')
    sub_group.create_dataset('ETr', data=et_r)
    sub_group.create_dataset('ETo', data=et_o)
    f_h5.close()


def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def cell_debugger(TPKPI_ini, cell):
    
    config = SafeConfigParser()
    config.read(TPKPI_ini)
    file_cell_param   = config.get('input_files', 'file_cell_param')
    
    ar_cell_label, ar_coorx, \
    ar_coory, ar_lambda, \
    ar_Xc, ar_wwb, \
    ar_tan_beta, ar_tan_beta_channel, \
    ar_L0, ar_Ks0, \
    ar_theta_r, ar_theta_s0, \
    ar_n_o0, ar_n_c0, \
    ar_cell_down, ar_pVs_t0, \
    ar_Vo_t0, ar_Qc_t0, \
    ar_kc, psi_b, lamda = pm.read_cell_parameters(file_cell_param)
    
    
    x_t, y_t = show_cell_cords(ar_cell_label,cell,ar_coorx,ar_coory)
    print 'The affected cell %i has the coordinates x=%0.2f, y=%0.2f'%(cell, x_t, y_t)


def seconds_to_dhms(seconds):
    days = str(int(seconds // (3600 * 24)))
    hours = str(int((seconds // 3600) % 24))
    minutes = str(int((seconds // 60) % 60))
    seconds = str(int(round(seconds % 60)))
    if len(days) < 2:
        days = '0' + days
    if len(hours) < 2:
        hours = '0' + hours
    if len(minutes) < 2:
        minutes = '0' + minutes
    if len(seconds) < 2:
        seconds = '0' + seconds
    return days, hours, minutes, seconds


def last_day_of_month(date):
    if date.month == 12:
        return date.replace(day=31)
    return (date.replace(month=date.month+1, day=1) - dtm.timedelta(days=1)).day