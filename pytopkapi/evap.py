import numpy as np
import datetime as dtm
import pytopkapi.utils as ut
"""Evaporation and evapotranspiration routines.

A selection of routines that can be used to compute the
evapotranspiration losses in PyTOPKAPI.

"""

def evapot_soil_Liu_and_Todini_ETc(Vs0,Vsm,kc,ETr,X):
    """
    The evapotranspiration taken up from the soil:
      - at the reference crop rate ETc (potential evapotranspiration rate)
      - without overland runoff infiltration
    """
    #ETr is in mm
    #From Reference crop evapotranspiration
    #to crop evapotranspiration
    ETc=kc*ETr
    #to actual evapotranspiration
    ETa=ETc
    #From mm to m
    ETa=ETa*1e-3
    #Compute the new soil volume
    Vs1=max(0.,Vs0-ETa*(X**2))
    #From m to mm
    ETa=ETa*1e3
    return ETa, Vs1

def evapot_soil_Liu_and_Todini(Vs0,Vsm,kc,ETr,X):
    """
    The evapotranspiration taken up from the soil:
      - at the rate ks*ETc ks being the soil saturation rate
      - without overland runoff infiltration
    """
    #ETr is in mm
    #From Reference crop evapotranspiration
    #to crop evapotranspiration
    ETc=kc*ETr
    #to actual evapotranspiration
    ks=Vs0/Vsm
    ETa=ks*ETc
    #From mm to m
    ETa=ETa*1e-3
    #Compute the new soil volume
    Vs1=max(0.,Vs0-ETa*(X**2))
    #From m to mm
    ETa=ETa*1e3
    return ETa, Vs1

def evapot_soil_overland(Vo0, Vs0, Vsm, kc, ETr, X):
    """Compute the evapotranspiration from a model cell.

    The evapotranspiration loss is first taken from the overland
    store, if the storage in the overland store cannot satisfy the ET
    demand then the water is extracted from the soil store. In cases
    where the ET demand is greater than the available water, both the
    soil and overland store are totally drained.

    Parameters
    ----------
    Vo0 : scalar
        Volume in the overland store before ET removal (:math:`m^3`)
    Vs0 : scalar
        Volume in the soil store before ET removal (:math:`m^3`)
    Vsm : scalar
        Volume of soil moisture for the store under saturated
        conditions (:math:`m^3`)
    kc : scalar
        Dimensionless crop co-efficient for the model cell
    ETr : scalar
        Reference crop ET for the cell during the current time-step
        (:math:`mm`)
    X : scalar
        The lateral dimension of the grid-cell (:math:`m`)

    Returns
    -------
    ETa : scalar
        The actual ET removed from the soil store, calculated as
        :math:`K_c K_s ET_r` with :math:`K_s` the relative saturation
        of the soil store before ET is removed (:math:`mm`)
    Vs1 : scalar
        Volume in the soil store after ET removal (:math:`m^3`)
    Vo1 : scalar
        Volume in the overland store after ET removal (:math:`m^3`)

    """
    #ETr is in mm
    #From Reference crop evapotranspiration
    #to crop evapotranspiration
    ETc = kc*ETr
    #to actual evapotranspiration
    ks = Vs0/Vsm
    ETa = ks*ETc
    #From mm to m
    ETa = ETa*1e-3
    if Vo0 > 0:
        if Vo0-ETa*(X**2) >= 0:
            Vo1 = Vo0-ETa*(X**2)
            Vs1 = Vs0
        else:
            Vo1 = 0.
            Vs1 = max(0., Vs0-(ETa*(X**2)-Vo0))
    else:
        Vo1 = 0.
        Vs1 = max(0., Vs0-ETa*(X**2))
    #From m to mm
    ETa = ETa*1e3

    return ETa, Vs1, Vo1

def evapor_channel(Vc0, ETo, W, X):
    """Compute the evaporation from a channel store.

    Calculate the evaporation loss from a channel store. The function
    attempts to remove the demand from the open water potential
    evaporation from the channel store, if there isn't sufficient
    water in the store then the total volume in the channel is
    removed.

    Parameters
    ----------
    Vc0 : scalar
        Volume in the channel store before evaporation removal
        (:math:`m^3`)
    ETo : scalar
        The potential evaporation from an open water surface
        (:math:`mm`)
    W : scalar
        Width of the channel (:math:`m`)
    X : scalar
        The length of the channel cell, this can be different from the
        cell dimension if the channel runs along the cell diagonal
        (:math:`m`)

    Returns
    -------
    ET_channel : scalar
        The actual depth of water removed from the channel store
        (:math:`mm`)
    Vc1 : scalar
        Volume in the channel store after evaporation removal
        (:math:`m^3`)

    """
    ETo = ETo*1e-3
    if Vc0-ETo*W*X > 0:
        Vc1 = Vc0-ETo*W*X
        Vevap_c = ETo*W*X
    else:
        Vc1 = 0
        Vevap_c = Vc0
    ET_channel = Vevap_c*1e3/(W*X)

    return ET_channel, Vc1

def intercept_rain_ET(P,ETr,kc):
    """
    The evapotranspiration is taken from the precipitation:
      - at the reference crop rate ETc
    """
    ETc=kc*ETr
    if P-ETc>=0:
        Pn=P-ETc
        ETcn=0.
        ETa=ETc
    else:
        Pn=0.
        ETcn=(ETc-P)
        ETa=P
    ETrn=ETcn/kc
    return Pn,ETrn,ETa


def evapot_soil_overland_feddes(Vo0, Vs0, Vsm, kc, ETr, ETo, X, psi, WWB):
    """Compute the evapotranspiration from a model cell.
    The evapotranspiration loss is first taken from the overland
    store, if the storage in the overland store cannot satisfy the ET
    demand then the water is extracted from the soil store. Water from the
    overland store is evaporated at the rate of open water ET (ETo).
    In cases where the ET demand is greater than the available water, both the
    soil and overland store are drained up to the wilting point.
    The potential transpiration rate is reduced according to the Feddes
    water stress function which is dependent on the matric suction (psi).

    Parameters
    ----------
    Vo0 : scalar
        Volume in the overland store before ET removal (:math:`m^3`)
    Vs0 : scalar
        Volume in the soil store before ET removal (:math:`m^3`)
    Vsm : scalar
        Volume of soil moisture for the store under saturated
        conditions (:math:`m^3`)
    kc : scalar
        Dimensionless crop co-efficient for the model cell
    ETr : scalar
        Reference ET for the cell during the current time-step
        (:math:`mm`)
    ETo : scalar
        Open water evaporation, here used for overland flow ET
        (:math:`mm`)
    X : scalar
        The lateral dimension of the grid-cell (:math:`m`)
    psi : scalar
        The soil suction value before ET removal estimated from psi_b and the
        Lambda parameter from  Brooks (:math:`mm`)
    WWB : scalar
        Contains the wetland unit number or dam value respectively. It will be
        used to switch between wetland and terrestrial soil. Usually wetland
        plants are able to transpire at potential rate despite water logging.

    Returns
    -------
    ETa : scalar
        The actual ET removed from the soil store, calculated using
         psi and the stress coefficient from Feddes (:math:`mm`)
    Vs1 : scalar
        Volume in the soil store after ET removal (:math:`m^3`)
    Vo1 : scalar
        Volume in the overland store after ET removal (:math:`m^3`)

    """
    # From Reference crop ET to crop ET in m
    ETc = kc*ETr*1e-3

    # Calculate the Feddes stress response parameter
    fsr = Feddes_stress_response(psi, WWB)

    # Convert open water et (ETo) to m
    ETo = ETo*1e-3

    if Vo0 > 0:
        # case where ET demand (ETo) is met from overland store (Vo0)
        if Vo0-ETo*(X**2) >= 0:
            ETa_o = ETo*(X**2)
            ETa_s = 0
            Vo1 = Vo0-ETo*(X**2)
            Vs1 = Vs0
        else:
            # case where ET demand (ETo) is partly met by Vo0.
            # Once Vo0 store is depleted then the remainder of ET demand is
            # calculated as the remainder of the ratio of ETo to Vo0.
            # The remaining "percentage" is then taken from the soil store.
            R_ETo_Vo0 = (ETo*X**2-Vo0)/(ETo*X**2)
            ETa_o = Vo0
            ETa_s = ETc*fsr*(X**2)*R_ETo_Vo0
            Vo1 = 0.
            Vs1 = max(0., Vs0-ETa_s)
    else:
        # case where overland store is 0 and hence ET demand is taken from Vs0
        Vo1 = 0.
        ETa_s = ETc*fsr*(X**2)
        Vs1 = max(0., Vs0-ETa_s)
        ETa_o = 0

    # Lump ET from overland and soil store together and convert from m3 to mm
    ETa = ((ETa_o+ETa_s)/X**2)*1e3

    return ETa, Vs1, Vo1


def Feddes_stress_response(psi, WWB):
    '''Compute the Feddes sink term which is the dimensionless stress parameter
    dependent on the matric pressure and plant properties i.e. plant response
    to water and oxygen stress.

    Parameters
    ----------
    psi : scalar
        The soil suction value before ET removal estimated from the Brooks
        Lambda parameter and the relative saturation (mm)
    WWB : scalar
        Contains the wetland unit number or dam value respectively. It will be
        used to switch between wetland and no terrestrial plants. Usually wetland
        plants are able to transpire at potential rate despite water logging
        and therefore for wetland areas oxygen stress is neglected.

    Returns
    -------
    fsr : scalar
        Stress response coefficient (0-1)
    '''

    #convert psi to be negative and in m
    psi = psi/-1000.

    # Feddes parameters (plant dependent redcution due to water or oxygen stress)
    # Idealy there should be a cataloug of differnt parameter sets for differnt
    # land use classes. E.g. pasture, millies, trees, wetlands, ect.
    # However this would entail to have an additional distributed parameter.
    # The parameters for non-wetland areas are from
    # All tresholds are in m
    if WWB > 0 and WWB < 200:
        h00 = 1.
        h0p = 0.5
    else:
        h00 = -0.1
        h0p = -0.25
    h02 = -2.
    h03 = -8.
    h04 = -80.

    y = 0
    if psi <= h00 and psi >= h0p:
#        print 'case 1'
        m   = 1./(h0p - h00)
        h   = psi
        b   = m*h00* - 1
        y   = m*h + b

    if psi < h0p and psi >= h02:
#        print 'case 2'
        y = 1

    if psi < h02 and psi >= h03:
#        print 'case 3'
        m = -0.2/(h03 - h02)
        h   = psi
        b   = -m*h02 + 1
        y   = m*h + b

    if psi < h03 and psi >= h04:
#        print 'case 4'
        m = -0.8/(h04 - h03)
        h   = psi
        b   = -m*h03 + 0.8
        y   = m*h + b

    if y > 1:
        raise ValueError('Feddes stress function does not perform correctly. \
The stress response coefficient (0-1) has yielded %0.2f' %y)
    return y

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    h_min, h_max = -10, 100000
    y = np.zeros((h_max-h_min))
    x = np.zeros((h_max-h_min))
    cnt = 0
    for i in xrange(h_min, h_max):
#        i = i/10.
        y[cnt] = Feddes_stress_response(i, 0)
        x[cnt] = i/-1000.
#        print 'x, y:', i/-1000., y[cnt]
        cnt += 1
    plt.clf()

    plt.plot(x,y)
#    a = plt.gca()
#    a.set_xscale('log')
    plt.show()


def evapot_soil_overland_feddes_AVC(Vo0, Vs0, Vsm, kc, ETr, ETo, X, psi, WWB,
                                    SSI, dat, m_i, m_d, SSI_ET_treshold):

    '''This function reduces the evapotranspiration during certain months
    of the year when vegetation is inactive. For the Highfeld this period is
    approx. between March and October where vegetation is inactive due to
    cold conditions, not necessarily due to water stress. However if the soil
    is very wet i.e. SSI > 85 % then direct soil evaporation does occur.
    In this case evaporation takes place at the ETo - open water rate.

        Parameters
    ----------
    Vo0 : scalar
        Volume in the overland store before ET removal (:math:`m^3`)
    Vs0 : scalar
        Volume in the soil store before ET removal (:math:`m^3`)
    Vsm : scalar
        Volume of soil moisture for the store under saturated
        conditions (:math:`m^3`)
    kc : scalar
        Dimensionless crop co-efficient for the model cell
    ETr : scalar
        Reference ET for the cell during the current time-step
        (:math:`mm`)
    ETo : scalar
        Open water evaporation, here used for overland flow ET
        (:math:`mm`)
    X : scalar
        The lateral dimension of the grid-cell (:math:`m`)
    psi : scalar
        The soil suction value before ET removal estimated from psi_b and the
        Lambda parameter from  Brooks (:math:`mm`)
    WWB : scalar
        Contains the wetland unit number or dam value respectively. It will be
        used to switch between wetland and terrestrial soil. Usually wetland
        plants are able to transpire at potential rate despite water logging.
    SSI : scalar
        Soil saturation index (0-1), dimensionless.
    dat : datetime.datetime object
        date_time stamp of the current time step.
    m_i : scalar
        Month when vegetation starts to increase, i.e. spring.
    m_d : scalar
        Month when vegetation starts to decrease, i.e. autum.

    '''

    K_owa = 0.8

    # Conversion from reference crop ET to crop ET in m
    ETc = kc*ETr*1e-3

    # Calculate the Feddes stress response parameter
    fsr = Feddes_stress_response(psi, WWB)

    # Convert open water et (ETo) to m
    ETo  = ETo * 1e-3
    fac_t, fac_e, s = f_v_AVC(SSI, fsr, dat, m_i, m_d, SSI_ET_treshold)

    if Vo0 > 0:
        # case where ET demand (ETo) is met from overland store (Vo0)
        if Vo0-ETo*(X**2) >= 0:
            ETa_o = ETo*(X**2)
            ETa_s = 0
            Vo1 = Vo0-ETo*(X**2)
            Vs1 = Vs0
        else:
            # case where ET demand (ETo) is partly met by Vo0.
            # Once Vo0 store is depleted then the remainder of ET demand is
            # calculated as the remainder of the ratio of ETo to Vo0.
            # The remaining "percentage" is then taken from the soil store.
            R_ETo_Vo0 = (ETo*X**2-Vo0)/(ETo*X**2)
            ETa_o     = Vo0
            Ea        = ETo * fac_e * K_owa * (X**2) * R_ETo_Vo0
            Ta        = ETc * fac_t * (X**2) * R_ETo_Vo0
            ETa_s     = Ea + Ta

            Vo1 = 0.
            Vs1 = max(0., Vs0-ETa_s)
    else:
        # case where overland store is 0 and hence ET demand is taken from Vs0
        Vo1   = 0.
        Ea    = ETo * fac_e * K_owa * (X**2)
        Ta    = ETc * fac_t * (X**2)
        ETa_s = Ea + Ta
        Vs1   = max(0., Vs0-ETa_s)
        ETa_o = 0

    # Lump ET from overland and soil store together and convert from m3 to mm
    ETa = ((ETa_o+ETa_s)/X**2)*1e3
    if __name__ == "__main__":
        return ETa, Vs1, Vo1, s
    else:
        return ETa, Vs1, Vo1


def f_v_AVC(SSI, fsr, dat, m_i, m_d, SSI_ET_treshold):
    '''See documentation of evapot_soil_overland_feddes_AVC function

    Returns
    f_Es_Tv : scalar
    factor Evaporation soil Transpiration vegetation - A dimensionless factor
    to reduce the transpiration if vegetation if it is inactive during winter.
    '''
    # Define boundaries i.e. max and min factor
    hi    = 1.
    lo    = 0.1
    yer   = dat.year
    mon   = dat.month
    day   = dat.day
    t_ssi = SSI_ET_treshold
    seasn = 0

    # get day of year from the dat
    dat_doy = dat.timetuple().tm_yday
#    m_i_doy = dtm.date(dat.year, m_i, 1).timetuple().tm_yday
#    m_d_doy = dtm.date(dat.year, m_d, 1).timetuple().tm_yday
    doy_su_s = dtm.datetime(yer,m_i+1,1).timetuple().tm_yday
    doy_su_e = \
    dtm.datetime(yer,m_d-1,ut.last_day_of_month(dtm.datetime(yer,m_d-1,1)))\
    .timetuple().tm_yday
    doy_wi_s = dtm.datetime(yer,m_d+1,1).timetuple().tm_yday
    doy_wi_e = \
    dtm.datetime(yer,m_i-1,ut.last_day_of_month(dtm.datetime(yer,m_i-1,1)))\
    .timetuple().tm_yday

    summer1 = range(doy_su_s, 367)
    summer2 = range(1, doy_su_e+1)
    winter1 = range(doy_wi_s, doy_wi_e+1)

#    print 'doy_su_s', doy_su_s
#    print 'doy_su_e', doy_su_e+1
#    print 'doy_wi_s', doy_wi_s
#    print 'doy_wi_e', doy_wi_e+1
#    print 'mon', mon

    if SSI <= t_ssi:
        # case where SSI is below the soil evaporation treshold
        if mon == m_i:
            seasn = 1
            # case where the factor is increasing
            nb_days = ut.last_day_of_month(dat)
            t_fac   = (hi-lo) * (float(day)/nb_days) + lo
        elif mon == m_d:
            seasn = -1
            # case where the factor is decreasing
            nb_days = ut.last_day_of_month(dat)
            t_fac   = (hi-lo) * (1-float(day)/nb_days) + lo
        elif dat_doy in summer1 or dat_doy in summer2:
            seasn = 8
            # case where veg is fully active
            t_fac   = hi
        elif dat_doy in winter1:
            seasn = 4
            # case where veg is inactive
            t_fac   = lo
        else:
            raise ValueError('Date fell through ET algorithem. Date not \
covered here is %s.' %dat.strftime("%d %b %Y"))
        t_fac = t_fac * fsr
        e_fac = 0.
    else:
        # Soil evaporation takes place at a rate of open water ET
        # when SSI is larger than a treshold
        if mon == m_i:
            seasn = 1
            # case where the factor is increasing
            nb_days = ut.last_day_of_month(dat)
            t_fac   = (hi-lo) * (float(day)/nb_days) + lo
            e_fac   = (1-t_fac)
        elif mon == m_d:
            seasn = -1
            # case where the factor is decreasing
            nb_days = ut.last_day_of_month(dat)
            t_fac   = (hi-lo) * (1-float(day)/nb_days) + lo
            e_fac   = (1-t_fac)
        elif dat_doy in summer1 or dat_doy in summer2:
            seasn = 8
            # case where veg is fully active
            t_fac   = hi
            e_fac   = (1-t_fac)
        elif dat_doy in winter1:
            seasn = 4
            # case where veg is inactive
            t_fac   = lo
            e_fac   = (1-t_fac)
        else:
            raise ValueError('Date fell through ET algorithem. Date not \
covered here is %s.' %dat.strftime("%d %b %Y"))
        t_fac = t_fac * fsr
    return t_fac, e_fac, seasn

if __name__ == "__main__":

#    dat = dtm.datetime(2012,2,29)
#    print 'doy of dat', dat.timetuple().tm_yday
#    f_v_AVC(0.6,1.,dat,10,3)




#    stop
    import time

    # define start date
    start_date = dtm.date.today()
    dat = start_date
    # constants:
    kc  = 1.
    X   = 10.
    WWB = 0
    m_i = 10
    m_d = 3
    nts = 400 # number of time steps
    Vsm = 20.
    lda = 0.25
    psb = 280.8 # psi_b
    Dt  = 86400
    thd = 0.5


    # generate random data
    ar_Vo0 = np.random.random(size=nts)*0.
    ar_Vs0 = np.random.random(size=nts)*Vsm
    ar_ssi = ar_Vs0/Vsm
    ar_psi = psb/np.power(ar_ssi, 1.0/lda)

    ar_ETr = np.random.random(size=nts)*6
    ar_ETo = ar_ETr * 0.8

    ar_ETa = np.ones((nts))*-999
    ar_Vs1 = np.ones((nts))*-999
    ar_Vo1 = np.ones((nts))*-999
    ar_dat = []
    ar_t_fac = np.ones((nts))*-999
    ar_e_fac = np.ones((nts))*-999
    ar_seasn = np.ones((nts))*-999

    t0 = time.time()
    # loop through time
    for t in xrange(nts):
        dat += dtm.timedelta(seconds=Dt)
        ar_dat.append(dat)

        fsr = Feddes_stress_response(ar_psi[t], WWB)
        ar_t_fac[t], ar_e_fac[t], s0 = f_v_AVC(ar_ssi[t],fsr,dat,m_i,m_d,thd)
        ar_ETa[t], ar_Vs1[t], ar_Vo1[t], s =\
        evapot_soil_overland_feddes_AVC(ar_Vo0[t], ar_Vs0[t], Vsm, kc,
                                        ar_ETr[t], ar_ETo[t], X, ar_psi[t],
                                        WWB, ar_ssi[t], dat, m_i, m_d, thd)
        if ar_t_fac[t] + ar_e_fac[t] > 1.0:
            raise ValueError('The f_v_AVC function procuces invalid results \
for t_fac (%s) and e_fac (%s) respectively. Their sum cannot be > 1.' \
        %(ar_ETa[t], ar_Vs1[t]))
#        print dat, s
        print ar_psi[t]/1000., fsr


    t1 = time.time()
    dt = t1-t0
    print 'AVC execution of %i runs has taken %0.4f seconds.' %(nts, dt)

    fig   = plt.figure()

    ax1   = fig.add_subplot(111)
    ax2   = ax1.twinx()
    lines = []
    t_leg = []
    color = ['orange', 'black', 'red', 'brown', 'blue', 'green', 'brown']

    lines += ax1.plot(ar_dat, ar_ETr,   color[0], marker='o', linewidth=1)
    lines += ax1.plot(ar_dat, ar_ETo,   color[1], marker='o', linewidth=1)
    lines += ax1.plot(ar_dat, ar_ETa,   color[2], marker='o', linewidth=1)
    lines += ax2.plot(ar_dat, ar_Vo0,   color[3], marker='o', linewidth=1)
    lines += ax2.plot(ar_dat, ar_ssi,   color[4], marker='o', linewidth=1)
    lines += ax2.plot(ar_dat, ar_t_fac, color[5], marker='o', linewidth=0)
    lines += ax2.plot(ar_dat, ar_e_fac, color[6], marker='o', linewidth=0)

    t_leg.append(r'$ET_r$')
    t_leg.append(r'$ET_o$')
    t_leg.append(r'$ET_a$')
    t_leg.append(r'$Vo0$')
    t_leg.append(r'$SSI$')
    t_leg.append(r'$T-fac$')
    t_leg.append(r'$E-fac$')

    #ax_A.set_xlim(0,4)
    for v in [0.75]:
        ax2.axhline(v, color='k', linestyle='dotted')
    ax1.legend(lines, t_leg, loc='best', fancybox=True)
    leg = ax1.get_legend()
    leg.get_frame().set_alpha(0.75)

    fig.savefig('P:\\Mooifontein\\PyTopkapi\\_Forcing_Data\\ET\\test_AVC', dpi=200)

#==============================================================================
#     SECOND TEST WITH CONSTANT INPUT VALUES
#==============================================================================

    ssi    = 0.76
    psi    = psb/np.power(ssi, 1.0/lda)

    dat    = start_date
    ar_dat = []
    fsr    = Feddes_stress_response(psi, WWB)

    print 'psi (mm):', psi
    print 'fsr:', fsr

    test_len = 365

    ar_t_fac = np.ones((test_len))*-999.
    ar_e_fac = np.ones((test_len))*-999.
    ar_seasn = np.ones((test_len))*-999.

    for i in xrange(test_len):
        dat += dtm.timedelta(seconds=Dt)
        ar_dat.append(dat)

        ar_t_fac[i],ar_e_fac[i],ar_seasn[i] = f_v_AVC(ssi,fsr,dat,m_i,m_d,thd)


    fig2  = plt.figure()
    ax1   = fig2.add_subplot(111)
    ax2   = ax1.twinx()
#    ax2   = ax1.twinx()
    lines = []
    t_leg = []
    color = ['orange', 'black', 'red', 'brown', 'blue', 'green', 'brown']

    lines += ax1.plot(ar_dat, ar_t_fac, color[5], marker='.', linewidth=0)
    lines += ax1.plot(ar_dat, ar_e_fac, color[6], marker='.', linewidth=0)
    lines += ax2.plot(ar_dat, ar_seasn, color[1], marker='.', linewidth=0)

    t_leg.append('Tfac')
    t_leg.append('Tfac')
    t_leg.append('Season')

    ax2.set_ylim(-2,10)
    ax1.set_ylim(0,1)
    ax1.axhline(.75, color='k', linestyle='dotted')
    ax1.axhline(ssi, color='b', linestyle='solid')
    ax1.axhline(fsr, color='g', linestyle='solid')
    ax1.legend(lines, t_leg, loc='best', fancybox=True)
    leg = ax1.get_legend()
    leg.get_frame().set_alpha(0.75)

    fig2.savefig('P:\\Mooifontein\\PyTopkapi\\_Forcing_Data\\ET\\test_AVC2', dpi=200)

'''
    dat = dtm.date(2014,9,28)
    hi = 1.
    lo = 0.1

    for i in xrange(10):
        dat += dtm.timedelta(seconds=Dt)
        d = float(dat.day)
        nb_days = ut.last_day_of_month(dat)
        t_fac   = (hi-lo) * (d/nb_days) + lo
        print 'day',d
        print 'nb_days', nb_days
        print 't_fac', t_fac
'''
