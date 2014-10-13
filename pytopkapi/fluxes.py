"""
Functions required by the TOPKAPI model for the management of input
and output of cells.

.. note::

  The subroutines solving the differential equation are not in this
  module (see :mod:`~TOPKAPI.ode` for more information)

"""

import numpy as np

##        ROUTINES FOR SOIL STORE
#```````````````````````````````````````````
def initial_volume_soil(ar_pVs_t0, ar_Vsm):
    """ initialize_volume_soil
    Compute the intial content of water (volume of water) from the
    initial saturated values given in the parameter file.
    """
    ar_Vs_t0=ar_pVs_t0/100.*ar_Vsm
    return ar_Vs_t0

def input_soil(P, Dt, X, ar_Q_to_next_cell, ar_cell_up):
    """Compute the total input to a soil store.

    Calculate the total rate of input to a single soil store. This
    comprises the sum of rainfall input, subsurface contribution from
    upstream cells and the overland contribution from upstream cells.

    Parameters
    ----------
    P : scalar
        Precipitation input to the cell during the current time-step
        (:math:`mm`)
    Dt : scalar
        The length of the current time-step in seconds
    X : scalar
        The lateral dimension of the grid-cell (:math:`m`)
    ar_Q_to_next_cell : (N,) Numpy array
        The total contribution from each cell to it's downstream
        neighbour as a result of subsurface and overland fluxes
        calculated during the previous timestep (:math:`m^3/s`)
    ar_cell_up : list of `int`
        List of integer indices into `ar_Q_to_next_cell`. The indices
        point to the cells upstream of the current cell

    Returns
    -------
    a_s : scalar
        The input flux to the soil store during the current time-step
        in :math:`m^3/s`

    """
    #Transformation of P in mm to P_flux in m^3/s
    P_flux=(P*(10**-3)/Dt)*X**2
    #Case 1: cell without up)
    ind=ar_cell_up[ar_cell_up>-90.]
    ar_sel=ar_Q_to_next_cell[ind]
    if ar_sel[ar_sel<0].size!=0:
        print ''
        print 'STOP-ERROR: By computing -->Upstream from cell n. ',ind[ar_sel<0],' missing'
        print ''
        a_s='Error on upstream'
        return a_s
    else:
        a_s=P_flux+ar_sel.sum()
        if a_s < 0:
            print 'WARNING: a_s is negative'
            print 'ind, ar_cell_up, ar_sel, ar_sel.sum(), P_flux,',\
            ind,',',ar_cell_up,',',ar_sel,',',ar_sel.sum(),',',P_flux
        return a_s

def input_soil_shallow(P, Dt, X, ar_Q_to_next_cell, ar_cell_up,
                       return_2_shallow):
    """Compute the total input to a soil store.

    Calculate the total rate of input to a single soil store. This
    comprises the sum of rainfall input, subsurface contribution from
    upstream cells and the overland contribution from upstream cells.

    Parameters
    ----------
    P : scalar
        Precipitation input to the cell during the current time-step
        (:math:`mm`)
    Dt : scalar
        The length of the current time-step in seconds
    X : scalar
        The lateral dimension of the grid-cell (:math:`m`)
    ar_Q_to_next_cell : (N,) Numpy array
        The total contribution from each cell to it's downstream
        neighbour as a result of subsurface and overland fluxes
        calculated during the previous timestep (:math:`m^3/s`)
    ar_cell_up : list of `int`
        List of integer indices into `ar_Q_to_next_cell`. The indices
        point to the cells upstream of the current cell
    return_2_shallow : scalar
        Return flow from the deep to the shallow layer (m^3/s)

    Returns
    -------
    a_s : scalar
        The input flux to the soil store during the current time-step
        in :math:`m^3/s`

    """
    #Transformation of P in mm to P_flux in m^3/s
    P_flux=(P*(10**-3)/Dt)*X**2
    #Case 1: cell without up)
    ind=ar_cell_up[ar_cell_up>-90.]
    ar_sel=ar_Q_to_next_cell[ind]
    if ar_sel[ar_sel<0].size!=0:
        print ''
        print 'STOP-ERROR: By computing -->Upstream from cell n. ',ind[ar_sel<0],' missing'
        print ''
        a_s='Error on upstream'
        return a_s
    else:
        a_s=P_flux+ar_sel.sum()+return_2_shallow
        if a_s < 0:
            print 'WARNING: a_s is negative'
            print 'ind, ar_cell_up, ar_sel, ar_sel.sum(), P_flux, ,return_lower',\
            ind,',',ar_cell_up,',',ar_sel,',',ar_sel.sum(),',',P_flux,',',return_2_shallow
        return a_s



def output_soil(Vs_t0, Vs_t1_prim, Vsm, a_s, b_s, alpha_s, Dt):
    """Compute the outflow from and volume in a soil store.

    Calculate the outflow rate and final volume in a soil store after
    transferring volume greater than the saturated soil moisture
    content to the overland store of the same model cell.

    Parameters
    ----------
    Vs_t0 : scalar
        Volume of water in the soil store at the beginning of the
        current time- step (:math:`m^3`)
    Vs_t1_prim : scalar
        Volume of water in the soil store at the end of the current
        time-step as a result of combined input to and drainage from
        the non-linear reservoir (:math:`m^3`)
    Vsm : scalar
        Volume of soil moisture for the store under saturated
        conditions (:math:`m^3`)
    a_s : scalar
        The input flux to the soil store during the current time-step
        (:math:`m^3/s`)
    b_s : scalar
        The constant term of the non-linear differential equation
        :math:`dV_s/dt = a_s - b_sV_s^{\\alpha_s}`
    alpha_s : scalar
        The dimensionless pore-size distribution parameter for the
        soil store
    Dt : scalar
        The length of the current time-step (:math:`s`)

    Returns
    -------
    Qs_out : scalar
        Rate of flow out of the soil store during the time-step
        (:math:`m^3/s`)
    Vs_out : scalar
        Volume remaining in the soil store at the end of the time-step
        (:math:`m^3`)

    """
    if Vs_t1_prim > Vsm:
        Q_max = b_s*Vsm**alpha_s

        Qs_out=Q_max
        Vs_out=Vsm
    else:
        Qs_out = Qout_computing(Vs_t0, Vs_t1_prim, a_s, Dt)
        Vs_out = Vs_t1_prim

    if Qs_out < 0:
        print 'a=', a_s, 'Vs_t1_prim=', Vs_t1_prim, 'Vs_t0=', Vs_t0, 'Vsm=', Vsm
        print 'Qs=', Qs_out

    return Qs_out, Vs_out


def output_soil2(Vs_t0_2, Vs_t1_prim2, Vsm2, a_s2, b_s2, alpha_s2, Dt):
    """Compute the outflow from and volume in a soil store.

    Calculate the outflow rate of the deep soil layer and final volume of the
    soil store after transferring volume greater than the saturated soil to
    the shallower soil layer.

    Parameters
    ----------
    Vs_t0_2 : scalar
        Volume of water in the soil store at the beginning of the
        current time- step (:math:`m^3`)
    Vs_t1_prim2 : scalar
        Volume of water in the soil store at the end of the current
        time-step as a result of combined input to and drainage from
        the non-linear reservoir (:math:`m^3`)
    Vsm2 : scalar
        Volume of soil moisture for the store under saturated
        conditions (:math:`m^3`)
    a_s2 : scalar
        The input flux to the soil store during the current time-step
        (:math:`m^3/s`)
    b_s2 : scalar
        The constant term of the non-linear differential equation
        :math:`dV_s/dt = a_s - b_sV_s^{\\alpha_s}`
    alpha_s2 : scalar
        The dimensionless pore-size distribution parameter for the
        soil store
    Dt : scalar
        The length of the current time-step (:math:`s`)

    Returns
    -------
    Qs_out2 : scalar
        Rate of flow out of the soil store during the time-step
        (:math:`m^3/s`)
    Vs_out2 : scalar
        Volume remaining in the soil store at the end of the time-step
        (:math:`m^3`)

    """
    if Vs_t1_prim2 > Vsm2:
        Q_max = b_s2*Vsm2**alpha_s2

        Qs_out2=Q_max
        Vs_out2=Vsm2
    else:
        Qs_out2 = Qout_computing(Vs_t0_2, Vs_t1_prim2, a_s2, Dt)
        Vs_out2 = Vs_t1_prim2

    if Qs_out2 < 0:
        print 'a_s2=', a_s2, 'Vs_t1_prim2=', Vs_t1_prim2, 'Vs_t0_2=', Vs_t0_2, 'Vsm2=', Vsm2
        print 'Qs2=', Qs_out2

    return Qs_out2, Vs_out2


def output_soil_parak(Vs_t0, Vs_t1_prim, Vsm, b_s, alpha_s):
    if Vs_t1_prim> Vsm:
        Vs_out=Vsm
    else:
        Vs_out=Vs_t1_prim
    Qs_out=Qout_computing2(Vs_t0,Vs_out, b_s, alpha_s)

    return Qs_out

def input_overland(Vs_t0,Vs_prim,Vsm,a_s,b_s,alpha_s,Dt):
    Q_prim=Qout_computing2(Vs_t0,Vs_prim,a_s,Dt)
    Q_max=b_s*Vsm**alpha_s
    a_o=max(0,Q_prim-Q_max)

    return a_o, Q_max

def Qout_computing(V_t0, V_t1_prim, a, Dt):
    """Compute the outflow `Qout` from a generic water store.

    Calculate the mean outflow rate during the current time-step from
    a generic water store. The outflow is calculated as the difference
    between the inflow and rate of change in storage during the
    time-step.

    Parameters
    ----------
    V_t0 : scalar
        Volume of water at the start of the current time-step
        (:math:`m^3`)
    V_t1_prim : scalar
        Volume of water at the end of the current time-step
        (:math:`m^3`)
    a : scalar
        The inflow rate to the water store during the time-step
        (:math:`m^3/s`)
    Dt : scalar
        The length of the current time-step (:math:`s`)

    Returns
    -------
    Qout : scalar
        The outflow rate from the water store during the time-step
        (:math:`m^3/s`)

    """
    Qout = a - (V_t1_prim - V_t0)/Dt
    return Qout

def Qout_computing2(V_t0,V_t1,b,alpha):
    """ Compute the output flows Qout from the computed water volumes:

    Returns
    -------
    Qout : scalar
        The outflow rate from the water store during the time-step
        (:math:`m^3/s`)

    """
    Qout=b/2.*(V_t1**alpha+V_t0**alpha)

    return Qout

def perc_prim(Vs0,Ks,eff_sat,c):
    """Computes the drainage from the upper to lower soil layer.

    This function calculates the free drainage flux from the upper soil layer
    using the Brooks and Corey function Pr(\theta) = Ksat Se^c

    Parameters
    ----------
    Vs0 : scalar
        Initial soil water volume (m^3)
    Ks : scalar
        Saturated hydraulic conductivity (m/s)
    eff_sat : scalar
        Effective saturation
    c : scalar
        Disconnectedness index, exponent of the Brooks and Corey funciton

    Returns
    -------
    Pr_prim : scalar
        Preliminary percolation rate to the deeper soil layer (m/s)
    """
    return Ks*eff_sat**c

def Qout_dam(Dt, WL_dam_avg, Qseep_mm, Vd0, Qc_in, ETo_t, A_dam):
    """
    Dam routing function.
    Computes the losses (Evaporation (E) and seepage) and the volume leaving
    the dam in case it is full to the downstream river cell.
    E and seepage losses are handled in a hirachial manner. If enough water
    is available then E and seepage are at potential rate times the area of the
    dam. If not then seepage losses take prevalence and E will be the remainder
    of water left in the dam until it is dry. Water is only passed on to the
    downstream river cell if the dam reaches FSL.

    Parameters
    ----------
    Dt : scalar
        Time step (seconds)
    WL_dam_avg : scalar
        Average water height in dam at full supply level (FSL)
    Qseep_mm : scalar
        Seepage rate losses for dam (mm/TMSP)
    Vd0 : scalar
        Volume in dam at the beginning of the time step
    Qc_in : scalar
        Inflow rate into dam (m3/TMSP)
    ETo_t : scalar
        Open water evaporation (mm/TMSP)
    A_dam : scalar
        Area of the dam (m2)

    Returns
    -------
    Qd_out : scalar
        Outflow rate from dam to downstream channel cell (m3/TMSP)
    Vd1 : scalar
        Volume of water in the dam at the end of the time step.
    """

    ETo_d_pot  = ETo_t * 1e-3 * A_dam        # [m3]
    Q_seep_pot = Qseep_mm * 1e-3 * A_dam     # [m3]
    Q_in       = Qc_in * Dt                  # [m3]
    Vd_max     = A_dam * WL_dam_avg          # [m3]

    if Vd0 + Q_in > Q_seep_pot:
        Q_seep = Q_seep_pot
        if Vd0 + Q_in - Q_seep > ETo_d_pot:
            ETo_dam = ETo_d_pot
        else:
            ETo_dam = Vd0 + Q_in - Q_seep
    else:
        ETo_dam = 0.
        Q_seep = Vd0 + Q_in

    V_bal = Q_in + Vd0 - Q_seep - ETo_dam
    if V_bal > 0:
        if V_bal > Vd_max:
            Qd_out = V_bal - Vd_max
            Vd1    = Vd_max
        else:
            Qd_out = 0.
            Vd1    = V_bal
    else:
        Qd_out = 0.
        Vd1    = 0.

    Qd_out = Qd_out/Dt

    return Qd_out, Vd1


def flow_partitioning(Lambda, Qs_out, Qo_out, W, X, Xc):
    """Partition the outflows from soil and overland stores.

    Calculates the correct partitioning of outflows from the soil and
    overland stores into contributions going to the downstream cell
    and/or the channel store of the current cell (if this exists).

    Parameters
    ----------
    Lambda : int
        Switch indicating whether the current cell contains a
        channel. A value of `1` indicates a channel cell, `0`
        indicates no channel
    Qs_out : scalar
        Outflow from the soil store during the current time-step
        (:math:`m^3/s`)
    Qo_out : scalar
        Outflow from the overland store during the current time-step
        (:math:`m^3/s`)
    W : scalar
        Width of the channel (:math:`m`)
    X : scalar
        The lateral dimension of the grid-cell (:math:`m`)
    Xc : scalar
        The length of the channel cell, this can be different from `X`
        if the channel runs along the cell diagonal (:math:`m`)

    Returns
    -------
    Q_to_next_cell : scalar
        Combined outflow from soil and overland to the downstream cell
        (:math:`m^3/s`)
    Q_to_channel : scalar
        Combined outflow from soil and overland to the channel store
        (:math:`m^3/s`)
    Q_to_channel_sub : scalar
        Outflow from soil store to the channel store (:math:`m^3/s`)

    """
    if Lambda != 0: #ToDo: check for invalid values of Lambda.
        if Lambda != 1:
            print 'Warning: The channel switch indicator \'Lamda\' is not 1\
but %d. Fix you river input file according to the specs.' %Lambda
        Q_to_next_cell   = (1-Lambda*W*X/(X**2))*(Qs_out+Qo_out)
        Q_to_channel     = (Lambda*W*X/(X**2))*(Qs_out+Qo_out)
        Q_to_channel_sub = (Lambda*W*X/(X**2))*(Qs_out)
    else:
        Q_to_next_cell = (Qs_out+Qo_out)
        Q_to_channel = 0.
        Q_to_channel_sub = 0.

    return Q_to_next_cell, Q_to_channel, Q_to_channel_sub


def flow_partitioning_DL(Lambda, Qs_out, Qs_out_, Qo_out, W, X, Xc):
    """Partition the outflows from soil and overland stores.

    Calculates the correct partitioning of outflows from the soil and
    overland stores into contributions going to the downstream cell
    and/or the channel store of the current cell (if this exists).

    Parameters
    ----------
    Lambda : int
        Switch indicating whether the current cell contains a
        channel. A value of `1` indicates a channel cell, `0`
        indicates no channel
    Qs_out : scalar
        Outflow from the soil store during the current time-step
        (:math:`m^3/s`)
    Qs_out_ : scalar
        Outflow from the soil store of the deep soil layer
        (:math:`m^3/s`)
    Qo_out : scalar
        Outflow from the overland store during the current time-step
        (:math:`m^3/s`)
    W : scalar
        Width of the channel (:math:`m`)
    X : scalar
        The lateral dimension of the grid-cell (:math:`m`)
    Xc : scalar
        The length of the channel cell, this can be different from `X`
        if the channel runs along the cell diagonal (:math:`m`)

    Returns
    -------
    Q_to_next_cell : scalar
        Combined outflow from soil and overland to the downstream cell
        (:math:`m^3/s`)
    Q_to_channel : scalar
        Combined outflow from soil and overland to the channel store
        (:math:`m^3/s`)
    Q_to_channel_sub : scalar
        Outflow from soil store to the channel store (:math:`m^3/s`)

    """
    if Lambda != 0: #ToDo: check for invalid values of Lambda.
        if Lambda != 1:
            print 'Warning: The channel switch indicator \'Lamda\' is not 1\
but %d. Fix you river input file according to the specs.' %Lambda
        Q_to_next_cell   = (1-Lambda*W*X/(X**2))*(Qs_out +Qo_out)
        Q_to_next_cell_  = (1-Lambda*W*X/(X**2))*(Qs_out_)
        Q_to_channel     = (Lambda*W*X/(X**2))*(Qs_out+Qs_out_+Qo_out)
        Q_to_channel_sub = (Lambda*W*X/(X**2))*(Qs_out+Qs_out_)
    else:
        Q_to_next_cell   = (Qs_out +Qo_out)
        Q_to_next_cell_  = (Qs_out_+Qo_out)
        Q_to_channel     = 0.
        Q_to_channel_sub = 0.

    return Q_to_next_cell, Q_to_next_cell_, Q_to_channel, Q_to_channel_sub

def input_channel(ar_Qc_out, Q_to_channel, ar_cell_up):
    """Compute the total inflow to the channel of a channel cell.

    Calculate the total inflow to the channel as the sum of inflows
    from upstream cells and inflows from the soil and overland stores
    in the current cell.

    Parameters
    ----------
    ar_Qc_out : (N,) Numpy array
        Array of channel outflows from each cell in the catchment for
        the current time-step (:math:`m^3/s`)
    Q_to_channel : scalar
        Combined outflow from soil and overland to the channel store
        (:math:`m^3/s`)
    ar_cell_up : list of `int`
        List of integer indices into `ar_Qc_out`. The indices point to
        the cells upstream of the current cell.

    Returns
    -------
    a_c : scalar
        The total inflow to the channel store during the current
        time-step (:math:`m^3/s`)
    Qc_cell_up : scalar
        The contribution to the total inflow to the channel store from
        upstream cells during the current time-step (:math:`m^3/s`)

    """
    ind=ar_cell_up[ar_cell_up>-90.]
    if len(ind)>0:
        ar_Qc_cell_up=ar_Qc_out[ind]
        a_c=Q_to_channel+ar_Qc_cell_up.sum()
    else:
        ar_Qc_cell_up=np.array([0.])
        if Q_to_channel < 0:
            a_c = 0.
            print 'set a_c = 0. It would have been %f' %Q_to_channel
            print 'ar_cell_up = ', ar_cell_up
        else:
            a_c=Q_to_channel

    Qc_cell_up = ar_Qc_cell_up.sum()

    return a_c, Qc_cell_up

def manning_depth(Q,W,X,n):
    """Compute Manning depth for flow `Q`.

    Compute with the manning equation (high width hypothesis) the
    water depth h corresponding to the flow Q for the channel
    characteristics: width W, manning n. The volume is then returned
    V=hWX

    """
    h=(n*Q/(W**(3./2.)))**(2./5.)

    return h

def initial_volume_channel(ar_Q,ar_W,X,ar_n_c):
    """Compute initial channel volume.

    The initial volume ar_Vc_t0 is computed for all the channel cells,
    from a known initial flow.

    """
    nb_cell=len(ar_W)
    ar_Vc_t0=np.zeros(nb_cell)
    for i in  range(nb_cell):
        if ar_W[i]>0.:
            h=manning_depth(ar_Q[i],ar_W[i],X,ar_n_c[i])
            ar_Vc_t0[i]=h*ar_W[i]*X

    return ar_Vc_t0
