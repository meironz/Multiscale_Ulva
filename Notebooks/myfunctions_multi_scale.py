import matplotlib.pyplot as plt
import numpy as np

# Functions of numeric solutions: 

# Function used to follow dynamics along cages - here used to follow spatial effects of dilution
def multi_N_f_un(t, y, Nintcrit,Nintmax,Nintmin,Vmax,Ks,KN,miu,S,Z,KI,K0,Ka,Topt,\
            Tmin,Tmax,losses20,teta,Sopt,Smin,Smax, Qp, Qsea, Nsea_upstream,f1,f0,dilution,n,umol_to_percent_DW):

    dydt = []
    
    # let's figure out the number of reactors, every 4 state variables are a single reactor
    n_reactors = len(y)//4
    
    Temp = f1(t)
    I0 = f0(t)
    
    for i_reactor in range(n_reactors): # loop that constructs the coupled ODEs
        
        # state variables: 
        
        i = 4*i_reactor # temporary counter of the state variables
        
        Nsea = y[i]   
        Next = y[i+1]  
        Nint = y[i+2]  
        m = y[i+3]
        
        # Nutrient consumption:
        Neff = (Nintmax - Nint)/(Nintmax - Nintmin)  # units: [ ]
        if Next <= 0:
            uN = 0
        else: 
            uN = Vmax * Next / (Ks + Next)  # units: [umol N/g DW/h]
        if Nint >= Nintcrit:
            fN = 1
        else:
            fN = ((Nint - Nintmin)/Nint) / ((Nintcrit - Nintmin)/Nintcrit)  # units: [ ]
        fP = 1  #(N:P < 12)


        # density - light penetration effects:
        SD = m * 1.785 / 2 * 1000 # Stocking Density. units: [g DW/ m^2]
        I_average = (I0 / (K0 * Z + Ka * SD)) * (1 - np.exp(-(K0 * Z + Ka * SD))) # units: [umul photons/(m^2 s)]
        fI = I_average / (I_average + KI) # units: [-]

        
        # Temperature effects:         
        if Temp <= Topt:
            Tx = Tmin
        else:
            Tx = Tmax

        fT = np.exp(-2.3 * ((Temp - Topt) / (Tx - Topt))**n) # Temp from temperature data

        # S (salinity) effects
        if S < Sopt:
            Sx = Smin
            b = 2.5
            if S < 5:
                fS = ((S - Smin)/(Sopt - Sx))
            elif S >= 5:
                fS = 1 - ((S - Sopt)/(Sx - Sopt)) ** b           
        elif S >= Sopt:
            Sx = Smax
            b = 4.4
            fS = 1 - ((S - Sopt)/(Sx - Sopt)) ** b

        # empirically defined losses
        losses = losses20 * teta ** (Temp - 20)

        # limiting factors:
        g = min(fN,fI,fP) * fT * fS
        
        # Qsea to reactor and back, per reactor, input is from the last one:
        if i_reactor > 0:
            Nsea_upstream = y[i-4] # y is a vector state, 4 state variables backwards is previous Nsea
            
        # Reactor Nsea -> Next -> Nint -> m and feedback
        dNsea = (Qsea * (Nsea_upstream * (1 - dilution) - Nsea) - Qp * (Nsea - Next)) / (1.785 * 1000)
        dNext = Qp / (1.785 * 1000) * (Nsea - Next) - Neff * uN * m # [umol N/l/h]
        dNint = umol_to_percent_DW * Neff * uN - Nint * miu * g #units: [%g N/g DW/h]

        if fI == 0:
            losses = 0
        dm = (miu * g - losses) * m #units: [g DW/l/h]
        
        
        dydt.extend([dNsea,dNext,dNint,dm])
        
    return dydt


# growth function with in-situ (HOBO) data and Reading parameters
def Reading_val(t, y, Nintcrit,Nintmax,Nintmin,Vmax,Ks,KN,miu,S,Z,KI,K0,Ka,Topt,\
            Tmin,Tmax,losses20,teta,Sopt,Smin,Smax, Qp, Qsea, Nsea_upstream,f1,f0,dilution,n,umol_to_percent_DW):

    import matplotlib.pyplot as plt
    import numpy as np
    
    dydt = []
    
    # every 4 state variables are a single reactor
    n_reactors = 1
      
    Temp = f1(t)
    
    I_average = f0(t)
    
    for i_reactor in range(n_reactors): # loop that constructs the coupled ODEs

        
        # state variables: 
        
        i = 4*i_reactor # temporary counter of the state variables
        
        Nsea = y[i]   
        Next = y[i+1]  
        Nint = y[i+2]  
        m = y[i+3]
        
        # Nutrient consumption:
        Neff = (Nintmax - Nint)/(Nintmax - Nintmin)  # units: [ ]
        if Next <= 0:
            uN = 0
        else: 
            uN = Vmax * Next / (Ks + Next)  # units: [umol N/g DW/h]
        if Nint >= Nintcrit:
            fN = 1
        else:
            fN = ((Nint - Nintmin)/Nint) / ((Nintcrit - Nintmin)/Nintcrit)  # units: [ ]
        fP = 1  #(N:P < 12)


        # density - light penetration effects:
        SD = m * 1.785 / 2 * 1000 # Stocking Density. units: [g DW/ m^2]
        fI = I_average / (I_average + KI) # units: [-]

        
        # Temperature effects: 
        
        if Temp <= Topt:
            Tx = Tmin
        else:
            Tx = Tmax

        fT = np.exp(-2.3 * ((Temp - Topt) / (Tx - Topt))**n) # Temp from temperature data

        # S (salinity) effects

        if S < Sopt:
            Sx = Smin
            b = 2.5
            if S < 5:
                fS = ((S - Smin)/(Sopt - Sx))
            elif S >= 5:
                fS = 1 - ((S - Sopt)/(Sx - Sopt)) ** b           
        elif S >= Sopt:
            Sx = Smax
            b = 4.4 # found by solver in fs file
            fS = 1 - ((S - Sopt)/(Sx - Sopt)) ** b

        # empirically defined losses
        losses = losses20 * teta ** (Temp - 20)


        # limiting factors:
        g = min(fN,fI,fP) * fT * fS

        
        # Qsea to reactor and back, per reactor, input is from the last one:

        if i_reactor > 0:
            Nsea_upstream = y[i-4] # y is a vector state, 4 state variables backwards is previous Nsea
            
        dNsea = (Qsea * (Nsea_upstream * (1 - dilution) - Nsea) - Qp * (Nsea - Next)) / (1.785 * 1000)
        
        # Reactor Next -> Nint -> m and feedback
        dNext = Qp / (1.785 * 1000) * (Nsea - Next) - Neff * uN * m # [umol N/l/h] 
        dNint = umol_to_percent_DW * Neff * uN - Nint * miu * g #units: [%g N/g DW/h]

        if fI == 0:
            losses = 0
        dm = (miu * g - losses) * m #units: [g DW/l/h]
        
        
        dydt.extend([dNsea,dNext,dNint,dm])
        
    return dydt


# Growth function with ex-situ (IMS) data and Reading parameters
def Reading_val_IMS(t, y, Nintcrit,Nintmax,Nintmin,Vmax,Ks,KN,miu,S,Z,KI,K0,Ka,Topt,\
            Tmin,Tmax,losses20,teta,Sopt,Smin,Smax, Qp, Qsea, Nsea_upstream,f1,f0,dilution,n,umol_to_percent_DW):
    """
    This is a second version of Reading_val. The difference is that this function works with IMS data and not HOBO data
    """
    
    import matplotlib.pyplot as plt
    import numpy as np
    
    dydt = []
    
    # every 4 state variables are a single reactor
    n_reactors = 1
    
    Temp = f1(t)
    #Temp = 20 #https://isramar.ocean.org.il/isramar_data/TimeSeries.aspx
    I0 = f0(t)
    
    for i_reactor in range(n_reactors): # loop that constructs the coupled ODEs
        
        # state variables: 
        
        i = 4*i_reactor # temporary counter of the state variables
        
        Nsea = y[i]   
        Next = y[i+1]  
        Nint = y[i+2]  
        m = y[i+3]
        
        # Nutrient consumption:
        Neff = (Nintmax - Nint)/(Nintmax - Nintmin)  # units: [ ]
        if Next <= 0:
            uN = 0
        else: 
            uN = Vmax * Next / (Ks + Next)  # units: [umol N/g DW/h]
        if Nint >= Nintcrit:
            fN = 1
        else:
            fN = ((Nint - Nintmin)/Nint) / ((Nintcrit - Nintmin)/Nintcrit)  # units: [ ]
        fP = 1  #(N:P < 12)

        # density - light penetration effects:
        SD = m * 1.785 / 2 * 1000 # Stocking Density. units: [g DW/ m^2]
        I_average = (I0 / (K0 * Z + Ka * SD)) * (1 - np.exp(-(K0 * Z + Ka * SD))) # units: [umul photons/(m^2 s)]
        fI = I_average / (I_average + KI) # units: [-]

        
        # Temperature effects: 
       
        if Temp <= Topt:
            Tx = Tmin
        else:
            Tx = Tmax

        fT = np.exp(-2.3 * ((Temp - Topt) / (Tx - Topt))**n) # Temp from temperature data

        # S (salinity) effects
        if S < Sopt:
            Sx = Smin
            b = 2.5
            if S < 5:
                fS = ((S - Smin)/(Sopt - Sx))
            elif S >= 5:
                fS = 1 - ((S - Sopt)/(Sx - Sopt)) ** b           
        elif S >= Sopt:
            Sx = Smax
            b = 4.4
            fS = 1 - ((S - Sopt)/(Sx - Sopt)) ** b

        # empirically defined losses
        losses = losses20 * teta ** (Temp - 20)


        # limiting factors:
        g = min(fN,fI,fP) * fT * fS
        
        
        # Qsea to reactor and back, per reactor, input is from the last one:
        if i_reactor > 0:
            Nsea_upstream = y[i-4] # y is a vector state, 4 state variables backwards is previous Nsea
            
        dNsea = (Qsea * (Nsea_upstream * (1 - dilution) - Nsea) - Qp * (Nsea - Next)) / (1.785 * 1000)
        
        # Reactor Next -> Nint -> m and feedback
        dNext = Qp / (1.785 * 1000) * (Nsea - Next) - Neff * uN * m # [umol N/l/h]
        dNint = umol_to_percent_DW * Neff * uN - Nint * miu * g #units: [%g N/g DW/h]

        if fI == 0:
            losses = 0
        dm = (miu * g - losses) * m #units: [g DW/l/h]
        
        
        dydt.extend([dNsea,dNext,dNint,dm])
        
    return dydt

# Functions of plots

# plotting NSEA, NEXT, NINT and M vs time - each reactor (or 10th reactor) gets a color
def plot_result_un(T,NSEA,NEXT,NINT,M,resolution,n_reactors,Nsea,fig,ax,season):
    
    import matplotlib.pyplot as plt
    import numpy as np
    import seaborn as sns
    
    sns.set_palette(sns.cubehelix_palette(20, start=2, rot=0, dark=0.2, light=.7, reverse=True))
    fig,ax = plt.subplots(3,1,sharex=True,figsize=(8,8))
        
    xlabels = ['day 0', '', 'day 2', '', 'day 4', '', 'day 6', '', 'day 8', '','day 10', '','day 12','', 'day 14']
    t = T[0]
    I = []
    
    for i in range(0,n_reactors,resolution):
        ax[0].plot(t,NSEA[i],'-',markersize=3)
        ax[0].set_ylim([0,Nsea*1.2])
        #ax[1].plot(t,NEXT[i],'.', markersize=3) # removed Next plot as it is redundant for the figure
        #ax[1].set_ylim([0,Nsea*1.2])
        ax[1].plot(t,NINT[i],'-', markersize=3)
        ax[1].set_ylim([0.5,4.5])
        ax[2].plot(t,M[i],'-', markersize=3)
        
        I.append(i)
        
    ax[0].set_ylabel('Nenv \n [µmol $l^{-1}$]',fontsize=13, weight="bold")
    #ax[1].set_ylabel('Next \n [umol / l]',fontsize=10, weight="bold") # removed Next plot as it is redundant for the figure
    ax[1].set_ylabel('Nint \n [% g N $g DW^{-1}$]',fontsize=13, weight="bold")
    ax[2].set_ylabel('m \n [g DW $l^{-1}$]',fontsize=13, weight="bold")
    
    
    ax[0].set_xticklabels([])
    ax[1].set_xticklabels([])
    ax[2].set_xticklabels([])
    #ax[3].set_xticklabels([])
    ax[0].set_xticks(np.linspace(T[0][0],T[0][0]+14*24,15))
    ax[1].set_xticks(np.linspace(T[0][0],T[0][0]+14*24,15))
    ax[2].set_xticks(np.linspace(T[0][0],T[0][0]+14*24,15))
    #ax[3].set_xticks(np.linspace(T[0][0],T[0][0]+14*24,15))
    ax[2].set_xticklabels([str(i) for i in xlabels], rotation=45,fontsize=12, weight="bold")

# plotting NSEA, NEXT, NINT and M vs time - first cage and cage N with Nsea < Nsea_lim 
def plot_result_Nsea(T,NSEA,NEXT,NINT,M,Nsea,i_reactor):
    
    import matplotlib.pyplot as plt
    import numpy as np
    
    fig,ax = plt.subplots(3,1,sharex=True,figsize=(12,10))
    
    xlabels = ['day 0', 'day 1', 'day 2', 'day 3', 'day 4', 'day 5', 'day 6', 'day 7', 'day 8', 'day 9','day 10', 'day 11',
               'day 12','day 13', 'day 14']
    t = T[0]
    I = []
    colors = ['black', 'red']
    reactors = [0,i_reactor] 
    for i in reactors:
        ax[0].plot(t,NSEA[i],'.',color=colors[reactors.index(i)],markersize=2)
        ax[0].set_ylim([0,Nsea*1.2])
        #ax[1].plot(t,NEXT[i],'.',color='black', markersize=2) # removed Next plot as it is redundant for the figure
        #ax[1].set_ylim([0,Nsea*1.2])
        ax[1].plot(t,NINT[i],'.',color=colors[reactors.index(i)], markersize=2)
        ax[1].set_ylim([0.5,4.5])
        ax[2].plot(t,M[i],'.',color=colors[reactors.index(i)], markersize=2)
        I.append(i)
        
    ax[0].set_ylabel('Nsea \n [µmol $l^{-1}$]',fontsize=12, weight="bold")
    #ax[1].set_ylabel('Next \n [umol / l]',fontsize=10, weight="bold")
    ax[1].set_ylabel('Nint \n [% g N $g DW^{-1}$]',fontsize=12, weight="bold")
    ax[2].set_ylabel('m \n [g DW $l^{-1}$]',fontsize=12, weight="bold")

    
    ax[0].set_xticklabels([])
    ax[1].set_xticklabels([])
    ax[2].set_xticklabels([])
    #ax[3].set_xticklabels([])
    ax[0].set_xticks(np.linspace(T[0][0],T[0][0]+14*24,15))
    ax[1].set_xticks(np.linspace(T[0][0],T[0][0]+14*24,15))
    ax[2].set_xticks(np.linspace(T[0][0],T[0][0]+14*24,15))
    #ax[3].set_xticks(np.linspace(T[0][0],T[0][0]+14*24,15))
    #ax[3].set_xlabel('Time',fontsize=14, weight="bold")
    ax[2].set_xticklabels([str(i) for i in xlabels], rotation=45,fontsize=10, weight="bold")
    
    #ax[0].legend(I)
    ax[0].legend(['first','last'],fontsize=14)


# plotting NSEA, NINT and M vs time - cage N with Nsea < Nsea_lim in winter, for 3 seasons
def plot_result_seasons(T,NSEA_seasons,NINT_seasons,M_seasons,Nsea,NSEA_seasons_f,lengths):
    
    import matplotlib.pyplot as plt
    import numpy as np
    
    fig = plt.figure()
    fig,ax = plt.subplots(2,2,figsize = (16,9))
    plt.rcParams.update({'font.size': 11.5})

    xlabels = ['day 0', '', 'day 2', '', 'day 4', '', 'day 6', '', 'day 8', '','day 10', '','day 12','', 'day 14']
    t = T[0]
    I = []
    colors = ['blue', 'red', 'orange']
    seasons = ['Winter', 'Spring', 'Autumn'] 
    for i in seasons:
        ax.flat[2].plot(t,NSEA_seasons[seasons.index(i)],'-',color=colors[seasons.index(i)],markersize=2)
        ax.flat[2].set_ylim([0,Nsea*1.2])
        ax.flat[1].plot(t,NINT_seasons[seasons.index(i)],'-',color=colors[seasons.index(i)], markersize=2)
        ax.flat[1].set_ylim([1.5,4])
        ax.flat[3].plot(t,M_seasons[seasons.index(i)],'-',color=colors[seasons.index(i)], markersize=2)
        I.append(i)
        
    ax.flat[2].set_ylabel('Nenv \n [µmol $l^{-1}$]',fontsize=14, weight="bold")
    ax.flat[1].set_ylabel('Nint \n [% g N $g DW^{-1}$]',fontsize=14, weight="bold")
    ax.flat[3].set_ylabel('m \n [g DW $l^{-1}$]',fontsize=14, weight="bold")
    
    ax.flat[2].set_xticklabels([])
    ax.flat[3].set_xticklabels([])
    ax.flat[1].set_xticklabels([])
    ax.flat[2].set_xticks(np.linspace(T[0][0],T[0][0]+14*24,15))
    ax.flat[3].set_xticks(np.linspace(T[0][0],T[0][0]+14*24,15))
    ax.flat[1].set_xticks(np.linspace(T[0][0],T[0][0]+14*24,15))
    ax.flat[2].set_xticklabels([str(i) for i in xlabels], rotation=45,fontsize=12)
    ax.flat[3].set_xticklabels([str(i) for i in xlabels], rotation=45,fontsize=12)
    ax.flat[1].set_xticklabels([str(i) for i in xlabels], rotation=45,fontsize=12)
    ax.flat[0].set_xlabel('Number of reactors in farm',fontsize=14, weight="bold")
    ax.flat[0].set_ylabel('Final Nenv \n [µmol N $l^{-1}$]',fontsize=14, weight="bold")
    colors = ['blue', 'red', 'orange']
    seasons = ['Winter', 'Spring', 'Autumn']
    for i in range (3):
        x = range(len(NSEA_seasons_f[i]))
        print(len(NSEA_seasons_f[i]))
        lengths.append(len(x))
        y = NSEA_seasons_f[i]
        ax.flat[0].plot(x,y, '--',color=colors[i],markersize=1.5)
        ax.flat[0].set_xlim([0,1.1*max(lengths)])
    ax.flat[0].legend(seasons,loc='upper right',fontsize=14)  
    return lengths


# plotting NSEA, NINT and M vs time - cage N with Nsea < Nsea_lim in winter, for 3 seasons
def plot_result_seasons_days(T1,T2,T3,NSEA_seasons,NINT_seasons,M_seasons,Nsea):
    """ Plot time series of the results, it's organized as a list of solution structures, sol.t,sol.y 
    
    We will plot 3 axes each for the separate quantity, only the last reactors get a line
    """
    
    import matplotlib.pyplot as plt
    import numpy as np
    
    fig,ax = plt.subplots(3,1,sharex=True,figsize=(12,10))
    
    xlabels = ['day 0', 'day 1', 'day 2', 'day 3', 'day 4', 'day 5', 'day 6', 'day 7', 'day 8', 'day 9','day 10', 'day 11',
               'day 12','day 13', 'day 14']
    
    I = []
    colors = ['blue', 'red', 'orange']
    seasons = ['Winter', 'Spring', 'Autumn'] 
    T_all = [T1, T2, T3]
    for i in seasons:
        T = T_all[seasons.index(i)]
        t = range(len(T[0]))
        ax[0].plot(t,NSEA_seasons[seasons.index(i)],'.',color=colors[seasons.index(i)],markersize=2)
        ax[0].set_ylim([0,Nsea*1.2])
        ax[1].plot(t,NINT_seasons[seasons.index(i)],'.',color=colors[seasons.index(i)], markersize=2)
        ax[1].set_ylim([0.5,4.5])
        ax[2].plot(t,M_seasons[seasons.index(i)],'.',color=colors[seasons.index(i)], markersize=2)
        I.append(i)
        
    ax[0].set_ylabel('Nenv \n [µmol $l^{-1}$]',fontsize=10, weight="bold")
    #ax[1].set_ylabel('Next \n [umol / l]',fontsize=10, weight="bold") # removed Next plot as it is redundant for the figure
    ax[1].set_ylabel('Nint \n [% g N $g DW^{-1}$]',fontsize=10, weight="bold")
    ax[2].set_ylabel('m \n [g DW $l^{-1}$]',fontsize=10, weight="bold")
    
    ax[0].set_xticklabels([])
    ax[1].set_xticklabels([])
    ax[2].set_xticklabels([])
    #ax[3].set_xticklabels([])
    ax[0].set_xticks(np.linspace(T1[0][0],T1[0][0]+14*24,15))
    ax[1].set_xticks(np.linspace(T1[0][0],T1[0][0]+14*24,15))
    ax[2].set_xticks(np.linspace(T1[0][0],T1[0][0]+14*24,15))
    #ax[3].set_xticks(np.linspace(T[0][0],T[0][0]+14*24,15))
    ax[2].set_xticklabels([str(i) for i in xlabels], rotation=45,fontsize=10, weight="bold")
    
    ax[0].legend(seasons)

    
# plotting NSEA, NEXT, NINT and M vs time - cage N with Nsea < Nsea_lim in winter, for 3 different Qps - first and last cage
def plot_result_Qp_2cages(T,NSEA_Qp1,NSEA_Qp2,NEXT_Qp1,NEXT_Qp2,NINT_Qp1,NINT_Qp2,M_Qp1,M_Qp2,Nsea,Qps,season):
    
    import matplotlib.pyplot as plt
    import numpy as np
    import seaborn as sns
    
    fig,ax = plt.subplots(2,2,figsize=(15,8))

    xlabels = ['day 0', '', 'day 2', '', 'day 4', '', 'day 6', '', 'day 8', '','day 10', '','day 12','', 'day 14']
    
    t = T[0]
    I = []
    colors = ['turquoise', 'dodgerblue','blue'] 
    for i in Qps:
        ax.flat[0].plot(t,NSEA_Qp1[Qps.index(i)],'-',markersize=2.5,color=colors[Qps.index(i)],linewidth=2)
        ax.flat[0].plot(t,NSEA_Qp2[Qps.index(i)],':',markersize=3,color=colors[Qps.index(i)],linewidth=2.5)
        ax.flat[0].set_ylim([0,Nsea*1.2])
        ax.flat[2].plot(t,NEXT_Qp1[Qps.index(i)],'-',markersize=3,color=colors[Qps.index(i)])
        ax.flat[2].plot(t,NEXT_Qp2[Qps.index(i)],':',markersize=3,color=colors[Qps.index(i)],linewidth=2.5)
        ax.flat[2].set_ylim([0,Nsea*1.2])
        ax.flat[1].plot(t,NINT_Qp1[Qps.index(i)],'-',markersize=3,color=colors[Qps.index(i)])
        ax.flat[1].plot(t,NINT_Qp2[Qps.index(i)],':',markersize=3,color=colors[Qps.index(i)],linewidth=2.5)
        ax.flat[1].set_ylim([0.5,4.5])
        ax.flat[3].plot(t,M_Qp1[Qps.index(i)],'-',markersize=3,color=colors[Qps.index(i)])
        ax.flat[3].plot(t,M_Qp2[Qps.index(i)],':',markersize=3,color=colors[Qps.index(i)],linewidth=2.5)
        I.append(i)
        
    ax.flat[0].set_ylabel('Nenv \n [µmol $l^{-1}$]',fontsize=14, weight="bold")
    ax.flat[2].set_ylabel('Next \n [µmol $l^{-1}$]',fontsize=14, weight="bold")
    ax.flat[1].set_ylabel('Nint \n [% g N $g DW^{-1}$]',fontsize=14, weight="bold")
    ax.flat[3].set_ylabel('m \n [g DW $l^{-1}$]',fontsize=14, weight="bold")
    
    ax.flat[0].set_xticklabels([])
    ax.flat[1].set_xticklabels([])
    ax.flat[2].set_xticklabels([])
    ax.flat[3].set_xticklabels([])
    ax.flat[0].set_xticks(np.linspace(T[0][0],T[0][0]+14*24,15))
    ax.flat[1].set_xticks(np.linspace(T[0][0],T[0][0]+14*24,15))
    ax.flat[2].set_xticks(np.linspace(T[0][0],T[0][0]+14*24,15))
    ax.flat[3].set_xticks(np.linspace(T[0][0],T[0][0]+14*24,15))
    ax.flat[2].set_xticklabels([str(i) for i in xlabels], rotation=45,fontsize=12, weight="bold")
    ax.flat[3].set_xticklabels([str(i) for i in xlabels], rotation=45,fontsize=12, weight="bold")
    
    legend = ax.flat[3].legend(['Qp = 460 l/h, first reactor', 'Qp = 460 l/h, last reactor','Qp = 15 l/h, first reactor', 'Qp = 15 l/h, last reactor','Qp = 0 l/h, all reactors'],ncol=3,fontsize='x-large',framealpha=0.8,bbox_to_anchor=(0.88, -0.25))