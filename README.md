# Multiscale_Ulva
Code and data of multiscale Ulva model

## Introduction

    This is a multi-reactor, algae farm, simulation base function that will be 
    solved in time. The multiple reactors are built along the streamwise (x) direction and numbered
    in the order of appearance, i.e. x=0 is the first reactor that meets the original nutrient stream. 
    Every following reactor will receive a different nutrient concentration, Nenv(x) due to the uptake 
    by the previous reactor and dilution, in case of a diluting environment.

    The main construction here is that each reactor is 4 ODEs: 
    1. dNenv/dt (Nenv = Nsea)
    2. dNext/dt
    3. dNint/dt, and 
    4. dm/dt

    we develop a list of such using a loop and in the loop we create the coupling by propagating the Nenv[i] 
    from the first reactor to the following ones.

## List of variables and parameters

    State variables
    Nsea (Nenv):  nutrient content in the stream of water, [umol N/l]
    Nsea_upstream : nutrient content upstream the reactor, [umol N/l]
    Next : nutrient content in the reactor [umol N/l]
    Nint : nutrient content in the algae biomass [% g N/g DW]
    m : biomass density in reactor [g DW/l]
        
    Parameters (given and fitted)*:
    miu:       Maximum specific growth rate, [1/h]
    losses20:  Specific rate of biomass losses, [1/h]
    teta:      An empiric factor of biomass losses, [-] 
    Nintcrit:  Threshold Nint level below which the growth rate slows down (f(N_int )<1), [% g N/g DW]
    Nintmax:   Maximum internal nitrogen concentration, [% g N/g DW]
    Nintmin:   Minimum internal nitrogen concentration, [% g N/g DW]
    Vmax:      Maximum N uptake rate, [μmol N/g DW/ h]
    Ks:        Nitrogen half saturation uptake constant, [μmol N/l]    
    S:         Salinity level in water [PSU]
    T:         Water temperature [C]
    I0:         Incident photon irradiance at water surface, [μmol photons/m2/s]
    Iaverage:   Average photon irradiance in reactor, [μmol photons/m2/s]
    KI:        Light half saturation constant, [μmol photons/m2/s] 
    K0:        Water light extinction coefficient, [1/m]
    Ka:        Ulva light extinction coefficient, [m2/gDW] 
    Tmin:      Minimal temperature allowing Ulva growth, [℃]
    Topt:      Optimal temperature for Ulva growth, [℃]
    Tmax:      Maximal temperature allowing Ulva growth, [℃]
    Smin:      Minimal salinity allowing Ulva growth, [PSU]
    Sopt:      Optimal salinity for Ulva growth, [PSU]
    Smax:      Maximal salinity allowing Ulva growth, [PSU]
    Qsea:       Stream flow through an area equivalent to the reactor narrow-side cross-section, [l/h]
    Qp:        flow rate of the airlift pump [l/h], transfers Sea to the reactor to increase Next, due to conservation, 
               the same flow rate overflows from the reactor back to the Qsea, reducing the concentration in the Qsea stream,
               thus affecting the following reactors
    d:         Dilution ratio between every two reactors, [%]
    Vcage:     Reactor volume, [m3]
    Z:        Maximum water depth is the reactor, [m]
  
    *additional arguments are imported into args
    
## Functions and notebooks

    All functions are writen in 'myfunctions_multi_scale.py'
    The basic function is <multi_N_f_un> that performs the numeric solution in time (t) and space (x). This function is used in the
    notebooks: 'Sensitivity_analysis', Ýear-round_productivities' and others.
    The functions <Reading_val> and <Reading_val_IMS> were used for calibration process using two different types of irradiance data.
    
    Calibration was performed and evaluated in the following notebooks:
    1. 'Calibration_plots'
    2. 'RMSRE1_Error_plots'
    3. 'RMSRE2_Error_plots'
    4. 'RMSRE 1'
    5. 'RMSRE 2'
    
    Some functions were built especially for specific figures:
    1. <plot_result_un> looks at all 4 state variables and plots a line in a different color for every X's cage (i.e every 10th cage). 
       This function is used in the 'Dilution' notebook
    2. <plot_result_Nsea> adds cages till the level of Nsea decreases to a set value. When it finds this cage, it plots the lines for
       dynamics of the 4 state variables in the first and the last cages. This function is used in the '10uM_threshold' notebook.
    3. <plot_result_seasons> plots <plot_result_Nsea> for the three relevant seasons (Autumn, Winter and Spring) to find farm size
       for each season, and plots last-cage dynamics for each season, following a farm-size calculated by winter dynamics. This function
       is used in the '10uM_threshold' notebook.
    4. <plot_result_seasons_days> plots the dynamics of the last cage according to the farm-size calculated by winter dynamics. The
       difference compared to previous plots, is that here, cultivation cycles are adjusted in each season to achieve a concstant Nint
       content rather than the Nsea threshold. This function is used in the '10uM_threshold' notebook.
    5. <plot_result_Qp_2cages> plots the dynamics of the first and last cages in a farm with spefic Qp, enabling to examine the effect of
       different Qp values on spatial distribution in the farm. This function is used in the 'Qp' notebook. 

## Data
    
    The data is distributed along six excel files:
    1. input
       'Parameters_multi-scale' - has all parameter/arguments/initial conditions of model and simulations
       'Parameters_Reading' - has all parameters used for calibration
       'T_multi-scale' - has water temperature data (before interpulation) for simulations
    2. ims_data_2017_PAR, ims_data_2017_umol_photons, ims_T_data_2017 and HOBO - have temperature and light intensity data for calibration
    3. ims_data_2014_umol_photons - has light intensity data for simulations


## How to install 

   We recommend using Anaconda Python Distribution and `conda` environments. After installing Anaconda or Miniconda, please use the following command to reproduce the environment:

         conda env create --file environment.yml
         conda activate multiscale_ulva


## How to cite this work

      Please cite our publication ... 
      Please cite our code using DOI: 
