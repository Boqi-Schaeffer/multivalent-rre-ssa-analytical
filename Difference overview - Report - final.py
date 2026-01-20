import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.differentiate
import cmasher as cmr


"""Type of data you want to analyze"""
n_sets = 4
set_types = np.array([3,6]) #arms
n_set_types = set_types.size
include_0_bound = False
r_range_name = "full"

"""Functions"""
def change_color_value(color, amount=0.5): # >1 darkens the color
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = np.array(colorsys.rgb_to_hls(*mc.to_rgb(c)))
    return colorsys.hls_to_rgb(c[0],1-amount * (1-c[1]),c[2])

def idx(i_set_type, i_set): # I want a different number assigned to each dataset, for coloring
    return i_set + i_set_type*n_sets

def color(i_set_type, i_set):
    color = plt.cm.Set2(i_set)
    if set_types[i_set_type] == 3:
        return change_color_value(color, 0.5)
    elif set_types[i_set_type] == 6:
        return change_color_value(color, 1.1)
    
def label(i_set_type, i_set):
    if i_set < 2:
        peak_loc = 4
    else:
        peak_loc = 5

    if i_set == 0 or i_set == 2:
        peak_height = 'high'
    else:
        peak_height = 'low'

    label = f'{peak_height} peak at $10^{peak_loc} - {set_types[i_set_type]}$ arms'
    return label

N_TIMES = 2
mean_markersize = 70
scatter_markersize = 30
scatter_markersize = 4
mean_abs_marker = '+'
line_datapoints = '-'
# line_excluded = (0, (5,5))
line_excluded = '-'
marker_excluded = 'o'
markersize_excluded = 4
color_excluded = change_color_value('crimson', 1.2)
zorder_mean = 2.5

figheight = 8.5
figwidth = 16

plt.rcParams.update({'font.size': 16})
plt.rc('legend', fontsize=15) 


"""Initialization"""
statistics = pd.DataFrame


for i_set_type in range(n_set_types):
    for i_set in range(n_sets):
        N_ARMS = set_types[i_set_type]
        # applying the above scaling factor to high rate 1 arm makes the selectivity difference almost nonexistent, just like the percent difference in bound fraction
        R_range_sf = 0
        automate_avg_time = True
        show_species = True


        """Data file info"""
        parent_folder = f'BEP runs\\Data files'
        folder_name = f'{parent_folder}\\{r_range_name} receptor range\\{N_ARMS} arms\\set {i_set}'


        """CONSTANTS"""
        '''taken from stochastic model (hard-coded)'''
        kons = np.zeros((2, 4))

        '''define from predetermined group of 4 sets of kon_sol and kon_sur, see table in BEP report'''
        if N_ARMS == 3:
            #high peak at 10^4, low receptor depletion
            kons[0,0] = 1e-12
            kons[1,0] = 5e-3

            #low peak at 10^4, high receptor depletion
            kons[0,1] = 5e-8
            kons[1,1] = 1e-4

            #high peak at 10^5, low receptor depletion
            kons[0,2] = 1e-11
            kons[1,2] = 3e-4

            # low peak at 10^5, high receptor depletion
            kons[0,3] = 1e-7
            kons[1,3] = 1e-5

        elif N_ARMS == 6:
            #high peak at 10^4, low receptor depletion
            kons[0,0] = 1e-15
            kons[1,0] = 5e-4

            #low peak at 10^4, high receptor depletion (already done)
            kons[0,1] = 1e-8
            kons[1,1] = 5e-5

            #high peak at 10^5, low receptor depletion
            kons[0,2] = 1e-13
            kons[1,2] = 4e-5

            #low peak at 10^5, high receptor depletion
            kons[0,3] = 1e-8
            kons[1,3] = 5e-6


        LOWPOWER = 2 #TODO: take from stochastic data
        if r_range_name == 'full':
            HIGHPOWER = 7
            cells = 25
        elif r_range_name == 'low': 
            HIGHPOWER = 4.5
            cells = 13
        kon_sol = kons[0, i_set]
        kon_sur = kons[1, i_set]
        koff = 1
        AREA = 1
        VOLUME = 1*AREA

        # kon_sol = 1e-8
        # kon_sur = 1e-4
        # koff = 1
        # AREA = 1
        # VOLUME = 1*AREA

        # kon_sol = 1e-6
        # kon_sur = 1e-2
        # koff = 100
        # AREA = 1
        # VOLUME = 1*AREA


        '''Settings for data extraction'''
        AVERAGING_TIME = 10                  # seconds you want to average over (use instead of averaging fraction)
        CUTOFF_TIME = False


        """EXTRACTING AND AVERAGING REACTANT NUMBERS AT DIFFERENT TIME POINTS"""
        '''Opening receptor data'''
        # file_path = f'.\\{parent_folder}\\{folder_name}\\Testrun_3_arms_concentration_of_species_receptA.dat'
        file_path = f'.\\{folder_name}\\Testrun_{N_ARMS}_arms_concentration_of_species_receptA.dat'
        dfR = pd.read_csv(file_path, sep=r'\s{3,}', skiprows = 1, engine = 'python')


        if dfR['Time (t)'].is_monotonic_increasing:
            print('Time is increasing')
        else:
            print('Time goes up and down')
            for i_time in range(dfR['Time (t)'].size-1):
                if dfR['Time (t)'].iloc[i_time]>dfR['Time (t)'].iloc[i_time+1]:
                    print("difference is at step ", dfR['# Step (n)'].iloc[i_time], '\nIf this is empty, that means that the row is Null values.')
            exit()

        '''Extracting general data from data file'''
        N_cells = dfR.shape[1]-2
        N_PARTICLES = 1e6
        max_time = dfR.iat[-1, 1]
        if CUTOFF_TIME and CUTOFF_TIME <= max_time:
            max_time = CUTOFF_TIME
        if automate_avg_time:
            AVERAGING_TIME = max_time/N_TIMES/2
        max_step = dfR.shape[0]
        N_receptors = np.logspace(LOWPOWER, HIGHPOWER, N_cells)
        step_size = int(max_step/N_TIMES)

        '''Check that my simulation times will not go into the negatives'''
        while max_time/N_TIMES - AVERAGING_TIME < 0:
            print(f"N_TIMES (= {N_TIMES}) or AVERAGING_TIME (= {AVERAGING_TIME}) must be decreased")
            N_TIMES = int(input('N_TIMES = '))
            AVERAGING_TIME = float(input('AVERAGING_TIME = '))

        averages = np.zeros((N_ARMS+2, N_cells, N_TIMES)) # so (type, N_receptors, timepoint) like in my numerical model
        sim_times = np.zeros(N_TIMES)

        '''Receptors'''
        for i in range(N_TIMES):
            low_time = max_time/N_TIMES*(i+1) - AVERAGING_TIME #TODO is <0 for small max_time but large N_TIMES
            high_time = max_time/N_TIMES*(i+1)
            sim_times[i] = (high_time + low_time)/2

            temp = dfR['Time (t)'].between(low_time, high_time)
            times = dfR[dfR['Time (t)'].between(low_time, high_time)] # 2D df with (time, cell)
            averages[-1, :, i] = times.iloc[:, 2:].mean(0) # take average over the time axis

        '''Particles'''
        for j in range(0, N_ARMS+1):
            # file_path = f'.\\{parent_folder}\\{folder_name}\\Testrun_3_arms_concentration_of_species_n{j}0.dat'
            file_path = f'.\\{folder_name}\\Testrun_{N_ARMS}_arms_concentration_of_species_n{j}0.dat'
            df = pd.read_csv(file_path, sep=r'\s{3,}', skiprows = 1, engine = 'python')

            for i in range(N_TIMES):
                low_time = max_time/N_TIMES*(i+1) - AVERAGING_TIME
                high_time = max_time/N_TIMES*(i+1)
                times = df[df['Time (t)'].between(low_time, high_time)] # 2D df with (time, cell)
                averages[j, :, i] = times.iloc[:, 2:].mean(0) # take average over the time axis

        """CONVERTING AVERAGES"""
        '''To bound fraction'''
        bound_fraction_sto = np.sum(averages[1:N_ARMS+1], 0)/N_PARTICLES # output: array(cell, time)
        #Bookmark, I'm trying to transpose the bound fraction
        '''To selectivity'''
        selectivity_sto = np.zeros((N_receptors.size, N_TIMES))
        for i in range(N_TIMES):
            selectivity_sto[:, i] = np.gradient(np.log(bound_fraction_sto[:,i]), np.log(N_receptors))


        """         """
        """NUMERICAL"""
        """         """ 

        from scipy.integrate import solve_ivp
        import matplotlib.pyplot as plt
        import numpy as np

        import time

        """SYSTEM"""
        initial_values = np.zeros(N_ARMS + 2)
        initial_values[0] = N_PARTICLES

        """Rate constants k"""
        k = np.zeros((N_ARMS+1, N_ARMS+1))

        for i in range(k.shape[0]):
            for j in range(k.shape[1]):
                k[i,j] = AREA
                if i == 0 and j == 1:
                    k[i,j] *= kon_sol
                    k[i,j] /= AREA*VOLUME # reaction of a surface molecule with a solution molecule
                elif j-i == 1: # if it's a binding reaction
                    k[i,j] *= kon_sur
                    k[i,j] /= AREA*AREA # reaction of a surface molecule with a surface molecule
                elif j-i == -1: # if it's an unbinding reaction
                    k[i,j] *= koff
                    k[i,j] /= AREA # reaction of a surface molecule

        '''setting rate constants corrected for number of arms'''
        c = np.zeros(k.shape)

        for i in range(k.shape[0]):
            for j in range(k.shape[1]):
                if i<j:
                    c[i,j] = k[i,j]*(N_ARMS - i)  # to bind an arm, multiply k by number of unbound arms
                elif i>j:
                    c[i,j] = k[i,j]*i # to release an arm, multiply k by number of bound arms

        """PLOTTING CONSTANTS"""
        tmin = 0 # When does time start running, NOT the lowest simulated time I graph my selectivity for.
        tmax = max_time
        nsteps = 100
        solver = 'LSODA'

        """CALCULATIONS"""
        """System of ordinary differential equations for all molecules and their configurations"""
        # f[0...n] is state with n arms bound, f[-1] is receptor concentration
        def rhs(t, f):
            ODE = np.zeros(N_ARMS+2)
            for i in range(N_ARMS+1): # we make the ODE system components by going through each state and adding the change to all states that a reaction from that state induces
                if i<N_ARMS: #binding reaction
                    change = c[i,i+1]*f[-1]*f[i]
                    ODE[i] -= change #subtract from i_arms state
                    ODE[i+1] += change #add to i+1_arms state
                    ODE[-1] -= change #subtract from number of receptors

                if i>0: #unbinding reaction
                    change = c[i,i-1]*f[i]
                    ODE[i] -= change #subtract from i_arms state
                    ODE[i-1] += change #add to i-1_arms state
                    ODE[-1] += change #add to number of receptors
            # ODE[-1]=0
            return ODE


        """Solve system of ODEs"""
        #maybe TODO: parallelize. You may need to use Numba or some other integrator than solve_ivp https://stackoverflow.com/questions/57706940/solving-ode-with-large-set-of-initial-conditions-in-parallel-through-python

        bound_fraction_num = np.zeros((N_receptors.size, sim_times.size)) # bound_fraction[N, t] so N on the y-axis and t on the x-axis
        results = np.zeros((initial_values.size, N_receptors.size, sim_times.size)) # results [type, N, t] with type being the same as the initial values (N0, N1, N2, N3, R)

        for i in range (N_receptors.size):
            initial_values[-1] = N_receptors[i]
            res = solve_ivp(rhs, (tmin, tmax), initial_values, t_eval = sim_times, method=solver)
            results[:, i] = res.y
            bound_fraction_num[i] = np.sum(res.y[1:-1], 0)/N_PARTICLES

        '''converting to selectivity'''
        selectivity_num = np.zeros((N_receptors.size, N_TIMES))

        for i in range(sim_times.size):
            selectivity_num[:, i] = np.gradient(np.log(bound_fraction_num[:,i]), np.log(N_receptors))

        """Calculating and saving differences"""
        if i_set_type == 0 and i_set == 0:
            selectivity_differences = np.zeros((n_set_types, n_sets, N_receptors.size, N_TIMES))
            percent_differences = np.zeros((n_set_types, n_sets, N_receptors.size, N_TIMES))

        '''selectivity_difference over receptors'''
        selectivity_difference = abs(-selectivity_sto+selectivity_num)/(selectivity_sto+selectivity_num)*2*100
        # selectivity_difference = abs(selectivity_sto/selectivity_num)
        selectivity_differences[i_set_type, i_set] = selectivity_difference

        '''Bound fraction over receptors difference'''
        percent_difference = 100*abs(-bound_fraction_sto+bound_fraction_num)/((bound_fraction_num+bound_fraction_sto)/2)
        # percent_difference = abs(bound_fraction_sto/bound_fraction_num) 
        percent_differences[i_set_type, i_set] = percent_difference

        
        print(f"[set type = {i_set_type} - set = {i_set} - range = {r_range_name}] complete")

"""PLOTTING"""
f, axs = plt.subplots(2, 2, sharex=True, sharey=True)
f.set_figheight(figheight)
f.set_figwidth(figwidth)

axs[0,0].set_xscale('log')
axs[1,0].set_xscale('log')
axs[0,1].set_xscale('log')
axs[1,1].set_xscale('log')

'''Selectivity difference'''
for i_set in range(n_sets):
    for i_set_type in range(n_set_types):
        axs[0,0].semilogx(N_receptors, 
                       selectivity_differences[i_set_type, i_set, :, -1], 
                       color = color(i_set_type, i_set),
                       linestyle = line_datapoints,
                       label = label(i_set_type, i_set))

mask = np.isnan(selectivity_differences)+np.isinf(selectivity_differences)
masked_selectivity_differences = np.ma.masked_array(selectivity_differences, mask) #ignore infinity or NaN selectivities
mean = np.ma.mean(masked_selectivity_differences[:, :, :, -1], axis=(0,1))
print(f'Selectivity mean {r_range_name} range\n{mean}')
axs[0,0].scatter(N_receptors, 
                 mean, marker = '_', 
                 color = 'k', 
                 s = mean_markersize, 
                 label = "Mean",
                 zorder=zorder_mean)

"""datapoints left out"""
twin00 = axs[0,0].twinx()
twin00.semilogx(N_receptors, 
              np.sum(mask[:,:,:,-1], (0,1)), 
              color = color_excluded,
              marker = marker_excluded,
              markersize = markersize_excluded)
twin00.set(ylim = (0,8), ylabel = "Datapoints excluded\nfrom mean")
twin00.yaxis.label.set_color(color_excluded)
twin00.tick_params(axis = 'y', colors = color_excluded)


'''percent differences of bound fraction'''
for i_set in range(n_sets):
    for i_set_type in range(n_set_types):
        axs[0,1].semilogx(N_receptors, 
                       percent_differences[i_set_type, i_set, :, -1], 
                       color = color(i_set_type, i_set),
                       linestyle = line_datapoints,
                       label = label(i_set_type, i_set))

if include_0_bound == True:
    masked_percent_differences = np.ma.masked_array(percent_differences, np.isnan(percent_differences)+np.isinf(percent_differences)) #ignore infinity or NaN selectivities
else:
    mask = np.isnan(percent_differences)+np.isinf(percent_differences)+percent_differences > 199.999999999
    masked_percent_differences = np.ma.masked_array(percent_differences, 
                                                       mask=mask) #ignore difference if stochastic bound fraction is 0
mean = np.ma.mean(masked_percent_differences[:, :, :, -1], axis=(0,1))
print(f'Bound fraction mean {r_range_name} range\n{mean}')
axs[0,1].scatter(N_receptors, mean, marker = '_', color = 'k', s = mean_markersize, label = "Mean", zorder = zorder_mean)

"""datapoints left out"""
twin10 = axs[0,1].twinx()
twin10.semilogx(N_receptors, 
              np.sum(mask[:,:,:,-1], (0,1)), 
              color = color_excluded,
              marker = marker_excluded,
              markersize = markersize_excluded)
twin10.set(ylim = (0,8), ylabel = "Datapoints excluded\nfrom mean")
twin10.yaxis.label.set_color(color_excluded)
twin10.tick_params(axis = 'y', colors = color_excluded)
























"""     """
"""AGAIN"""
"""     """

r_range_name = "low"

for i_set in range(n_sets):
    for i_set_type in range(n_set_types):
        N_ARMS = set_types[i_set_type]
        # applying the above scaling factor to high rate 1 arm makes the selectivity difference almost nonexistent, just like the percent difference in bound fraction
        R_range_sf = 0
        automate_avg_time = True
        show_species = True


        """Data file info"""
        parent_folder = f'BEP runs\\Data files'
        if set_types[i_set_type] == 6 and i_set == 2:
            folder_name = f'{parent_folder}\\low receptor range\\{N_ARMS} arms\\set {i_set} break_iter 1'
            time_weighted = True
        else:
            folder_name = f'{parent_folder}\\{r_range_name} receptor range\\{N_ARMS} arms\\set {i_set}'
            time_weighted = False


        """CONSTANTS"""
        '''taken from stochastic model (hard-coded)'''
        kons = np.zeros((2, 4))

        '''define from predetermined group of 4 sets of kon_sol and kon_sur, see table in BEP report'''
        if N_ARMS == 3:
            #high peak at 10^4, low receptor depletion
            kons[0,0] = 1e-12
            kons[1,0] = 5e-3

            #low peak at 10^4, high receptor depletion
            kons[0,1] = 5e-8
            kons[1,1] = 1e-4

            #high peak at 10^5, low receptor depletion
            kons[0,2] = 1e-11
            kons[1,2] = 3e-4

            # low peak at 10^5, high receptor depletion
            kons[0,3] = 1e-7
            kons[1,3] = 1e-5

        elif N_ARMS == 6:
            #high peak at 10^4, low receptor depletion
            kons[0,0] = 1e-15
            kons[1,0] = 5e-4

            #low peak at 10^4, high receptor depletion (already done)
            kons[0,1] = 1e-8
            kons[1,1] = 5e-5

            #high peak at 10^5, low receptor depletion
            kons[0,2] = 1e-13
            kons[1,2] = 4e-5

            #low peak at 10^5, high receptor depletion
            kons[0,3] = 1e-8
            kons[1,3] = 5e-6


        LOWPOWER = 2 #TODO: take from stochastic data
        if r_range_name == 'full':
            HIGHPOWER = 9
            cells = 25
        elif r_range_name == 'low': 
            HIGHPOWER = 4.5
            cells = 13
        kon_sol = kons[0, i_set]
        kon_sur = kons[1, i_set]
        koff = 1
        AREA = 1
        VOLUME = 1*AREA


        '''Settings for data extraction'''
        AVERAGING_TIME = 10                  # seconds you want to average over (use instead of averaging fraction)
        CUTOFF_TIME = False


        """EXTRACTING AND AVERAGING REACTANT NUMBERS AT DIFFERENT TIME POINTS"""
        '''Opening receptor data'''
        # file_path = f'.\\{parent_folder}\\{folder_name}\\Testrun_3_arms_concentration_of_species_receptA.dat'
        file_path = f'.\\{folder_name}\\Testrun_{N_ARMS}_arms_concentration_of_species_receptA.dat'
        dfR = pd.read_csv(file_path, sep=r'\s{3,}', skiprows = 1, engine = 'python')


        if dfR['Time (t)'].is_monotonic_increasing:
            print('Time is increasing')
        else:
            print('Time goes up and down')
            for i_time in range(dfR['Time (t)'].size-1):
                if dfR['Time (t)'].iloc[i_time]>dfR['Time (t)'].iloc[i_time+1]:
                    print("difference is at step ", dfR['# Step (n)'].iloc[i_time], '\nIf this is empty, that means that the row is Null values.')
            exit()

        '''Extracting general data from data file'''
        N_cells = dfR.shape[1]-2
        N_PARTICLES = 1e6
        max_time = dfR.iat[-1, 1]
        if CUTOFF_TIME and CUTOFF_TIME <= max_time:
            max_time = CUTOFF_TIME
        if automate_avg_time:
            AVERAGING_TIME = max_time/N_TIMES/2
        max_step = dfR.shape[0]
        N_receptors = np.logspace(LOWPOWER, HIGHPOWER, N_cells)
        step_size = int(max_step/N_TIMES)

        '''Check that my simulation times will not go into the negatives'''
        while max_time/N_TIMES - AVERAGING_TIME < 0:
            print(f"N_TIMES (= {N_TIMES}) or AVERAGING_TIME (= {AVERAGING_TIME}) must be decreased")
            N_TIMES = int(input('N_TIMES = '))
            AVERAGING_TIME = float(input('AVERAGING_TIME = '))

        averages = np.zeros((N_ARMS+2, N_cells, N_TIMES)) # so (type, N_receptors, timepoint) like in my numerical model
        sim_times = np.zeros(N_TIMES)

        '''Receptors'''
        for i in range(N_TIMES):
            low_time = max_time/N_TIMES*(i+1) - AVERAGING_TIME #TODO is <0 for small max_time but large N_TIMES
            high_time = max_time/N_TIMES*(i+1)
            sim_times[i] = (high_time + low_time)/2

            if time_weighted:
                slice = dfR[dfR['Time (t)'].between(low_time, high_time)] # 2D df with (time, cell)
                dt = slice['Time (t)'].iloc[1:].to_numpy() - slice['Time (t)'].iloc[:-1].to_numpy() # the weight is the time before the next reaction, aka the weight is the dwell time of the system.
                slice = slice.iloc[1:, 2:].to_numpy()
                averages[-1, :, i] = np.average(slice, axis = 0, weights = dt)
            else:
                slice = dfR[dfR['Time (t)'].between(low_time, high_time)] # 2D df with (time, cell)
                averages[-1, :, i] = slice.iloc[:, 2:].mean(0) # take average over the time axis

        '''Particles'''
        for j in range(0, N_ARMS+1):
            # file_path = f'.\\{parent_folder}\\{folder_name}\\Testrun_3_arms_concentration_of_species_n{j}0.dat'
            file_path = f'.\\{folder_name}\\Testrun_{N_ARMS}_arms_concentration_of_species_n{j}0.dat'
            df = pd.read_csv(file_path, sep=r'\s{3,}', skiprows = 1, engine = 'python')

            for i in range(N_TIMES):
                low_time = max_time/N_TIMES*(i+1) - AVERAGING_TIME
                high_time = max_time/N_TIMES*(i+1)
                if time_weighted:
                    slice = df[df['Time (t)'].between(low_time, high_time)] # 2D df with (time, cell)
                    dt = slice['Time (t)'].iloc[1:].to_numpy() - slice['Time (t)'].iloc[:-1].to_numpy() # the weight is the time before the next reaction, aka the weight is the dwell time of the system.
                    slice = slice.iloc[1:, 2:].to_numpy()    
                    averages[j, :, i] = np.average(a=slice, axis = 0, weights = dt)
                else:
                    slice = df[df['Time (t)'].between(low_time, high_time)] # 2D df with (time, cell)
                    averages[j, :, i] = slice.iloc[:, 2:].mean(0) # take average over the time axis

        """CONVERTING AVERAGES"""
        '''To bound fraction'''
        bound_fraction_sto = np.sum(averages[1:N_ARMS+1], 0)/N_PARTICLES # output: array(cell, time)
        #Bookmark, I'm trying to transpose the bound fraction
        '''To selectivity'''
        selectivity_sto = np.zeros((N_receptors.size, N_TIMES))
        for i in range(N_TIMES):
            selectivity_sto[:, i] = np.gradient(np.log(bound_fraction_sto[:,i]), np.log(N_receptors))


        """         """
        """NUMERICAL"""
        """         """ 

        from scipy.integrate import solve_ivp
        import matplotlib.pyplot as plt
        import numpy as np

        import time

        """SYSTEM"""
        initial_values = np.zeros(N_ARMS + 2)
        initial_values[0] = N_PARTICLES

        """Rate constants k"""
        k = np.zeros((N_ARMS+1, N_ARMS+1))

        for i in range(k.shape[0]):
            for j in range(k.shape[1]):
                k[i,j] = AREA
                if i == 0 and j == 1:
                    k[i,j] *= kon_sol
                    k[i,j] /= AREA*VOLUME # reaction of a surface molecule with a solution molecule
                elif j-i == 1: # if it's a binding reaction
                    k[i,j] *= kon_sur
                    k[i,j] /= AREA*AREA # reaction of a surface molecule with a surface molecule
                elif j-i == -1: # if it's an unbinding reaction
                    k[i,j] *= koff
                    k[i,j] /= AREA # reaction of a surface molecule

        '''setting rate constants corrected for number of arms'''
        c = np.zeros(k.shape)

        for i in range(k.shape[0]):
            for j in range(k.shape[1]):
                if i<j:
                    c[i,j] = k[i,j]*(N_ARMS - i)  # to bind an arm, multiply k by number of unbound arms
                elif i>j:
                    c[i,j] = k[i,j]*i # to release an arm, multiply k by number of bound arms

        """PLOTTING CONSTANTS"""
        tmin = 0 # When does time start running, NOT the lowest simulated time I graph my selectivity for.
        tmax = max_time
        nsteps = 100
        solver = 'LSODA'

        """CALCULATIONS"""
        """System of ordinary differential equations for all molecules and their configurations"""
        # f[0...n] is state with n arms bound, f[-1] is receptor concentration
        def rhs(t, f):
            ODE = np.zeros(N_ARMS+2)
            for i in range(N_ARMS+1): # we make the ODE system components by going through each state and adding the change to all states that a reaction from that state induces
                if i<N_ARMS: #binding reaction
                    change = c[i,i+1]*f[-1]*f[i]
                    ODE[i] -= change #subtract from i_arms state
                    ODE[i+1] += change #add to i+1_arms state
                    ODE[-1] -= change #subtract from number of receptors

                if i>0: #unbinding reaction
                    change = c[i,i-1]*f[i]
                    ODE[i] -= change #subtract from i_arms state
                    ODE[i-1] += change #add to i-1_arms state
                    ODE[-1] += change #add to number of receptors
            # ODE[-1]=0
            return ODE


        """Solve system of ODEs"""
        #maybe TODO: parallelize. You may need to use Numba or some other integrator than solve_ivp https://stackoverflow.com/questions/57706940/solving-ode-with-large-set-of-initial-conditions-in-parallel-through-python

        bound_fraction_num = np.zeros((N_receptors.size, sim_times.size)) # bound_fraction[N, t] so N on the y-axis and t on the x-axis
        results = np.zeros((initial_values.size, N_receptors.size, sim_times.size)) # results [type, N, t] with type being the same as the initial values (N0, N1, N2, N3, R)

        for i in range (N_receptors.size):
            initial_values[-1] = N_receptors[i]
            res = solve_ivp(rhs, (tmin, tmax), initial_values, t_eval = sim_times, method=solver)
            results[:, i] = res.y
            bound_fraction_num[i] = np.sum(res.y[1:-1], 0)/N_PARTICLES

        '''converting to selectivity'''
        selectivity_num = np.zeros((N_receptors.size, N_TIMES))

        for i in range(sim_times.size):
            selectivity_num[:, i] = np.gradient(np.log(bound_fraction_num[:,i]), np.log(N_receptors))

        """Calculating and saving differences"""
        if i_set_type == 0 and i_set == 0:
            selectivity_differences = np.zeros((n_set_types, n_sets, N_receptors.size, N_TIMES))
            percent_differences = np.zeros((n_set_types, n_sets, N_receptors.size, N_TIMES))

        '''selectivity_difference over receptors'''
        selectivity_difference = abs(-selectivity_sto+selectivity_num)/(selectivity_sto+selectivity_num)*2*100
        # selectivity_difference = abs(selectivity_sto/selectivity_num)
        selectivity_differences[i_set_type, i_set] = selectivity_difference

        '''Bound fraction over receptors difference'''
        percent_difference = 100*abs(-bound_fraction_sto+bound_fraction_num)/((bound_fraction_num+bound_fraction_sto)/2)
        # percent_difference = abs(bound_fraction_sto/bound_fraction_num)       
        percent_differences[i_set_type, i_set] = percent_difference

        
        print(f"[set type = {i_set_type} - set = {i_set} - range = {r_range_name}] complete")

"""PLOTTING"""

'''Selectivity difference'''
for i_set in range(n_sets):
    for i_set_type in range(n_set_types):
        axs[1,0].semilogx(N_receptors, 
                       selectivity_differences[i_set_type, i_set, :, -1], 
                       color = color(i_set_type, i_set),
                       linestyle = line_datapoints,
                       label = label(i_set_type, i_set))

mask = np.isnan(selectivity_differences)+np.isinf(selectivity_differences)
masked_selectivity_differences = np.ma.masked_array(selectivity_differences, mask=mask) #ignore infinity or NaN selectivities
mean = np.ma.mean(masked_selectivity_differences[:, :, :, -1], axis=(0,1))
print(f'Selectivity mean {r_range_name} range\n{mean}')

axs[1,0].scatter(N_receptors, mean, marker = '_', color = 'k', s = mean_markersize, label = "Mean", zorder = zorder_mean)
axs[1,0].title.set_text(r"Selectivity difference low $N_R$ range")

"""datapoints left out"""
twin10 = axs[1,0].twinx()
twin10.semilogx(N_receptors, 
              np.sum(mask[:,:,:,-1],(0,1)), 
              color = color_excluded,
              marker = marker_excluded,
              markersize = markersize_excluded)
twin10.set(ylim = (0,8), ylabel = "Datapoints excluded\nfrom mean")
twin10.yaxis.label.set_color(color_excluded)
twin10.tick_params(axis = 'y', colors = color_excluded)


'''percent differences of bound fraction'''
for i_set in range(n_sets):
    for i_set_type in range(n_set_types):
        axs[1,1].semilogx(N_receptors, 
                       percent_differences[i_set_type, i_set, :, -1], 
                       color = color(i_set_type, i_set),
                       linestyle = line_datapoints,
                       label = label(i_set_type, i_set))

if include_0_bound == True:
    masked_percent_differences = np.ma.masked_array(percent_differences, np.isnan(percent_differences)+np.isinf(percent_differences)) #ignore infinity or NaN selectivities
else:
    mask = np.isnan(percent_differences)+np.isinf(percent_differences)+percent_differences > 199.999999999
    masked_percent_differences = np.ma.masked_array(percent_differences, mask=mask) #ignore difference if stochastic bound fraction is 0
mean = np.ma.mean(masked_percent_differences[:, :, :, -1], axis=(0,1))
print(f'Bound fraction mean {r_range_name} range\n{mean}')

axs[1,1].scatter(N_receptors, mean, marker = '_', color = 'k', s = mean_markersize, label = "Mean", zorder = zorder_mean)

axs[1,1].title.set_text(r"Bound fraction difference low $N_R$ range") #bookmark

"""datapoints left out"""
twin11 = axs[1,1].twinx()
twin11.semilogx(N_receptors, 
              np.sum(mask[:,:,:,-1],(0,1)), 
              color = color_excluded,
              marker = marker_excluded,
              markersize = markersize_excluded)
twin11.set(ylim = (0,8), ylabel = "Datapoints excluded\nfrom mean")
twin11.yaxis.label.set_color(color_excluded)
twin11.tick_params(axis = 'y', colors = color_excluded)

"""Legend"""
# axs[0,1].legend()
# axs[1,1].legend()
handles, labels = axs[0,0].get_legend_handles_labels()
f.legend(handles,
         labels, 
         loc='lower center', 
         bbox_to_anchor=(0.5, 0.02),
         fancybox=False, 
         shadow=False, 
         ncols=n_sets+1)

axs[0,0].title.set_text(r"Selectivity difference full $N_R$ range")
axs[0,1].title.set_text(r"Bound fraction difference full $N_R$ range")

"""Axis labels"""
for column in range(2):
    axs[0,column].set_xlabel(r"Number of receptors $N_R$")
    axs[1,column].set_xlabel(r"Number of receptors $N_R$")
    axs[0,column].set_ylabel(r"Percent difference (%)")
    axs[1,column].set_ylabel(r"Percent difference (%)")

    # axs[0,column].set_ylim(bottom = -0.1, top = 5)
    # axs[1,column].set_ylim(bottom = -0.1, top = 5)

    axs[0,column].grid()
    axs[1,column].grid()

plt.tight_layout()

plt.subplots_adjust(left=0.07, right=0.92, wspace=0.31, bottom=0.2)

plt.show()