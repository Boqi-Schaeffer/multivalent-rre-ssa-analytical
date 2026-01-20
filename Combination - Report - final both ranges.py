import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.differentiate
import cmasher as cmr

"""Type of data you want to analyze"""
N_ARMS = 6
i_set = 3
R_range_sf = 0
automate_avg_time = True
show_species = True
r_range_name = "full"


"""Data file info"""
parent_folder = f'BEP runs\\Data files'
folder_name = f'{parent_folder}\\{r_range_name} receptor range\\{N_ARMS} arms\\set {i_set}'
# folder_name = f'Data2'
# folder_name = f'{parent_folder}\\{r_range_name} receptor range\\{N_ARMS} arms\\set {i_set} 600 particles'

if N_ARMS == 6 and i_set == 2:
    folder_name2 = f'{parent_folder}\\low receptor range\\{N_ARMS} arms\\set {i_set} break_iter 1'
else:
    folder_name2 = f'{parent_folder}\\low receptor range\\{N_ARMS} arms\\set {i_set}'



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
N_TIMES = 2
AVERAGING_TIME = 10                  # seconds you want to average over (use instead of averaging fraction)
CUTOFF_TIME = False

'''Colors for plotting'''
color = plt.cm.berlin(np.linspace(0, np.min([1, N_TIMES*0.1]), N_TIMES))
color_num = plt.cm.Blues(np.linspace(0.3, 0.85, N_TIMES))
# color = plt.cm.hot(1-np.linspace(0.2, 1, N_TIMES))
color_sto = plt.cm.Oranges(np.linspace(0.3, 0.85, N_TIMES))
# color_num = cmr.ocean(1-np.linspace(0.3, 0.8, N_TIMES))
# color_sto = cmr.amber(1-np.linspace(0.1, 0.6, N_TIMES))
# color_num = color
# color_sto = color
if N_TIMES == 1:
    color_num = plt.cm.Blues(np.linspace(0.85, 0.85, N_TIMES))
    color_sto = plt.cm.Oranges(np.linspace(0.85, 0.85, N_TIMES))

'''markers for plotting'''
marker_sto= 'x'
marker_num = 'o'
marker_size_sto = 10
marker_size_num = 7

line_ana = (0, (5,5))
color_ana = 'yellowgreen'
line_receptors = 'dashdot'
color_receptors = 'purple'

'''other plotting params'''
figheight = 7.5
figwidth = 16

plt.rcParams.update({'font.size': 16})
plt.rc('legend', fontsize=15) 



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
N_receptors = np.logspace(LOWPOWER, HIGHPOWER, N_cells)
powers = np.log10(N_receptors)


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
bound_fraction_sto = np.sum(averages[1:N_ARMS+1], 0)/N_PARTICLES # output: array(time, cell)

'''To selectivity'''
selectivity_sto = np.zeros((N_TIMES, N_receptors.size))
for i in range(N_TIMES):
    logfraction = np.log(bound_fraction_sto[:,i])
    logreceptors = np.log(N_receptors) 
    selectivity_sto[i] = np.gradient(np.log(bound_fraction_sto[:,i]), np.log(N_receptors))


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
    print('i is ', i)
    results[:, i] = res.y
    bound_fraction_num[i] = np.sum(res.y[1:-1], 0)/N_PARTICLES

'''converting to selectivity'''
selectivity_num = np.zeros((N_TIMES, N_receptors.size))

for i in range(sim_times.size):
    selectivity_num[i] = np.gradient(np.log(bound_fraction_num[:,i]), np.log(N_receptors))

"""PLOTTING"""
f, axs = plt.subplots(2, 2, sharex=True, sharey='col')
f.set_figheight(figheight)
f.set_figwidth(figwidth)

'''selectivity over receptors semilog for many simulation times'''
'''Legend is methods together'''
if N_TIMES == 1:
    axs[0,0].semilogx(N_receptors, selectivity_num[i], 
                      marker = marker_num, 
                      label = f"Numerical", 
                      c = color_num[i],
                      markersize = marker_size_num)
    axs[0,0].semilogx(N_receptors, 
                      selectivity_sto[i], 
                      marker = marker_sto, 
                      label = f"Stochastic", 
                      c = color_sto[i],
                      markersize = marker_size_sto)
else:
    for i in range(sim_times.size):
        time = format(sim_times[i], '.2f')
        axs[0,0].semilogx(N_receptors, 
                          selectivity_num[i], 
                          marker = marker_num, 
                          label = f"RRE - {time}s", 
                          c = color_num[i],
                          markersize = marker_size_num)
    for i in range(sim_times.size):
        time = format(sim_times[i], '.2f')
        axs[0,0].semilogx(N_receptors, 
                          selectivity_sto[i], 
                          marker = marker_sto, 
                          label = f"SSA - {time}s", 
                          c = color_sto[i],
                          markersize = marker_size_sto)
        

# plt.figtext(0.2, 0.2, f"\
#     kon_sol = {kon_sol} \n\
#     kon_sur = {kon_sur} \n\
#     koff = {koff} \n\
#     Number of arms = {N_ARMS} \n\
#     Number of particles = {initial_values[0]} \n\n\
#     ODE solver = {solver} \n\n\
#     Averaged over {AVERAGING_TIME} seconds")
axs[0,0].set_xlabel(r"Number of receptors $N_R$")
axs[0,0].set_ylabel(r"Selectivity $\alpha$")
axs[0,0].axhline(0, ls=':', c='grey')
axs[0,0].axvline(0, ls=':', c='grey')
axs[0,0].grid()

# '''Difference at endpoint'''
# difference = selectivity_sto[-1]-selectivity_num[-1]
# axs[0,1].semilogx(N_receptors, difference, marker = marker_num, color = 'k')
# axs[0,1].set_xlabel("Number of receptors")
# axs[0,1].set_ylabel(f"Selectivity")
# axs[0,1].axhline(0, ls=':', c='grey')
# axs[0,1].axvline(0, ls=':', c='grey')
# axs[0,1].legend()
# axs[0,1].grid()

# '''Difference over receptors'''
# for i in range(N_TIMES):
#     difference = -selectivity_sto[i, round(cells*R_range_sf):]+selectivity_num[i, round(cells*R_range_sf):]
#     axs[0,1].semilogx(N_receptors[round(cells*R_range_sf):], 
#                       difference, 
#                       label = f'{format(sim_times[i], '.2f')}s',
#                       marker = marker_num, 
#                       color = color[i],
#                       markersize = marker_size_num)
# axs[0,1].set_xlabel("Number of receptors")
# axs[0,1].set_ylabel(f"Selectivity difference N-S")
# axs[0,1].axhline(0, ls=':', c='grey')
# axs[0,1].axvline(0, ls=':', c='grey')
# axs[0,1].legend()
# axs[0,1].grid()

'''Bound fraction over receptors'''
# if show_species:
    # for i_species in range(1, N_ARMS+1):
    #     axs[0,1].loglog(N_receptors, 
    #                     averages[i_species, :, -1]/N_PARTICLES, 
    #                     label = f"{i_species} arms bound - Stochastic",
    #                     color = color_sto[i_species])
    #     axs[0,1].loglog(N_receptors, 
    #                     results[i_species, :, -1]/N_PARTICLES, 
    #                     label = f"{i_species} arms bound - Numerical",
    #                     color = color_num[i_species])
    # axs[0,1].legend()

    

if N_TIMES == 1:
    axs[0,1].loglog(N_receptors, 
                    bound_fraction_num[:,-1], 
                    marker = marker_num, 
                    label = f"Numerical", 
                    color = color_num[0],
                    markersize = marker_size_num)
    axs[0,1].loglog(N_receptors, 
                    bound_fraction_sto[:,-1], 
                    marker = marker_sto, 
                    label = f"Stochastic", 
                    color = color_sto[0],
                    markersize = marker_size_sto)
else:
    # axs[0,1].loglog(N_receptors, bound_fraction_num[:,-1], marker = marker_num, label = f"Numerical", color = 'purple')
    # axs[0,1].loglog(N_receptors, bound_fraction_sto[:,-1], marker = marker_sto, label = f"Stochastic", color = 'green')
    for i in range(N_TIMES):
        axs[0,1].loglog(N_receptors, 
                        bound_fraction_num[:,i], 
                        marker = marker_num, 
                        label = f"Numerical", 
                        color = color_num[i],
                        markersize = marker_size_num)
    for i in range(N_TIMES):
        axs[0,1].loglog(N_receptors, 
                        bound_fraction_sto[:,i], 
                        marker = marker_sto, 
                        label = f"Stochastic", 
                        color = color_sto[i],
                        markersize = marker_size_sto)

#Free receptors
twin = axs[0,1].twinx()
twin.semilogx(N_receptors, 
              results[-1, :, -1]/N_receptors, 
              color = 'purple',
              linestyle = line_receptors)
twin.set(ylim = (0,1.1), ylabel = "Free receptor fraction")
twin.yaxis.label.set_color('purple')
twin.tick_params(axis = 'y', colors = 'purple')

# fully_bound_fraction_s = averages[N_ARMS, :, -1]/np.sum(averages[1:N_ARMS+1, :, -1], axis=0)
# fully_bound_fraction_n = results[N_ARMS, :, -1]/np.sum(results[1:N_ARMS+1, :, -1], axis=0)

# axs2 = axs[0,1].twinx()
# axs2.semilogx(N_receptors, 
#                 fully_bound_fraction_s, 
#                 label = f"fully bound fraction - Stochastic",
#                 marker = marker_sto,
#                 color = color_sto[-1]) 
# axs2.semilogx(N_receptors, 
#                 fully_bound_fraction_n, 
#                 label = f"fully bound fraction - Numerical",
#                 marker = marker_num,
#                 color = color_num[-1]) 
# axs2.legend()

axs[0,1].set_xlabel(r"Number of receptors $N_R$")
axs[0,1].set_ylabel(r"Bound fraction $\theta$")
axs[0,1].axhline(0, ls=':', c='grey')
axs[0,1].axvline(0, ls=':', c='grey')
axs[0,1].grid()



# '''Bound fraction over receptors difference'''
# for i in range(sim_times.size):
#     percentage_difference = 100*(-bound_fraction_sto[:,i]+bound_fraction_num[:,i])/((bound_fraction_num[:,i]+bound_fraction_sto[:,i])/2)
#     time = format(sim_times[i], '.2f')
#     axs[1,1].semilogx(N_receptors, 
#                       percentage_difference, 
#                       marker = marker_num, 
#                       color = color[i], 
#                       label = f'{time}s')

# percentage_difference = 100*(-bound_fraction_sto[:,-1]+bound_fraction_num[:,-1])/((bound_fraction_num[:,-1]+bound_fraction_sto[:,-1])/2)
# axs[1,1].semilogx(N_receptors, percentage_difference, marker = marker_num, color = 'k')
# axs[1,1].set_xlabel("Number of receptors")
# axs[1,1].set_ylabel(f"Percentage difference (%) N-S")
# axs[1,1].axhline(0, ls=':', c='grey')
# axs[1,1].axvline(0, ls=':', c='grey')
# axs[1,1].legend()
# axs[1,1].grid()
# plt.show()





# """GRAPH OVER TIME"""
# sim_times = df['Time (t)'] # array of all times
# cell = 19
# print(N_receptors[cell])
# concentrations_sto = np.zeros((N_ARMS+2, sim_times.size))

# '''Particles'''
# for bound_arms in range(0, N_ARMS+1):

#     file_path = f'.\\{parent_folder}\\{folder_name}\\Testrun_{N_ARMS}_arms_concentration_of_species_n{bound_arms}0.dat'
#     df = pd.read_csv(file_path, sep=r'\s{3,}', skiprows = 1, engine = 'python')
    
#     concentration = df.iloc[:, cell+2]
#     concentrations_sto[bound_arms] = concentration # df(time, cell) take the values of all times at a specific cell.

# '''To bound fraction'''
# bound_fraction_sto = np.sum(concentrations_sto[1:-1], 0)/N_PARTICLES # output: array(time, cell)

# """Solve system of ODEs"""
# bound_fraction_num = np.zeros((N_receptors.size, sim_times.size)) # bound_fraction[N, t] so N on the y-axis and t on the x-axis
# results = np.zeros((initial_values.size, N_receptors.size, sim_times.size)) # results [type, N, t] with the order of the types being the same as the initial values (N0, N1, N2, N3, R)

# tmax = tmax*2

# for i in range (N_receptors.size):
#     initial_values[-1] = N_receptors[i]
#     res = solve_ivp(rhs, (tmin, tmax), initial_values, t_eval = sim_times, method=solver)
#     print('i is ', i)
#     results[:, i] = res.y
#     bound_fraction_num[i] = np.sum(res.y[1:-1], 0)/N_PARTICLES

# '''plot bound fraction over time'''
# print('N receptors is ', N_receptors)
# plt.figure(figsize=(5,4))
# plt.plot(sim_times, bound_fraction_num[cell], label = 'Numerical')
# plt.plot(sim_times, bound_fraction_sto, label = 'Stochastic')
# plt.xlabel("Time (s)")
# plt.ylabel("Bound fraction")
# plt.axhline(0, ls=':', c='grey')
# plt.axvline(0, ls=':', c='grey')
# plt.legend()
# plt.grid()
# plt.show()





"""              """
"""  ANALYTICAL  """
"""              """
# A = 22e6
rho_0 = 1e-8
rho_0 = 1e100
k = N_ARMS

rho_0 = N_PARTICLES/VOLUME
rho_n = rho_0 #? Is this correct?

K_intra = 1/AREA*kon_sur/koff
K_A = 1/rho_0*kon_sol/koff

K_intra = kon_sur/koff*(1/AREA)
K_A = 1/rho_0*kon_sol/koff*(1/VOLUME)

'''KA and Kintra should be unrelated to things like concentration and surface, no?'''
# rho_n = 1e-8
# A = 22e6

def calc_bound_fraction_analytical(sigma_R):
    Kav_A = K_A/K_intra*( (1+sigma_R*AREA*K_intra)**N_ARMS - 1)
    bound_fraction = rho_n*Kav_A/(1+rho_n*Kav_A)
    return bound_fraction

sigma_R = N_receptors/AREA
bound_fraction_analytical = calc_bound_fraction_analytical(sigma_R)

axs[0,1].loglog(N_receptors, bound_fraction_analytical, label = f"Analytical - steady state", color = color_ana, linestyle = line_ana)




# temp = scipy.differentiate.derivative(calc_bound_fraction_analytical, N_receptors)
# selectivity_analytical = temp.df
selectivity_analytical = np.gradient(np.log(bound_fraction_analytical), np.log(N_receptors))
axs[0,0].semilogx(N_receptors, selectivity_analytical, label = f"Analytical - steady state", color = color_ana, linestyle = line_ana)

"""Legend for both graphs"""
handles, labels = axs[0,0].get_legend_handles_labels()
f.legend(handles,
         labels, 
         loc='lower right', 
         bbox_to_anchor=(0.99, 0.63),
         fancybox=False, 
         shadow=False, 
         ncols=1)




















"""     """
"""LOW RANGE"""
"""     """

r_range_name = "low"

folder_name = folder_name2
# folder_name = f'{parent_folder}\\{r_range_name} receptor range\\{N_ARMS} arms\\set {i_set} retry 2gb'
# folder_name = f'{parent_folder}\\{r_range_name} receptor range\\{N_ARMS} arms\\set {i_set} retry 2gb'


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
    # kons[1,2] = 4e-5*1.5 ##trying to figure out why set 2 is different

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
max_time = dfR.iat[-1, 1]
if CUTOFF_TIME and CUTOFF_TIME <= max_time:
    max_time = CUTOFF_TIME
if automate_avg_time:
    AVERAGING_TIME = max_time/N_TIMES/2

N_receptors = np.logspace(LOWPOWER, HIGHPOWER, N_cells)
powers = np.log10(N_receptors)


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

    if 'break_iter 1' in folder_name2:
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
        if 'break_iter 1' in folder_name2:
            slice = df[df['Time (t)'].between(low_time, high_time)] # 2D df with (time, cell)
            dt = slice['Time (t)'].iloc[1:].to_numpy() - slice['Time (t)'].iloc[:-1].to_numpy() # the weight is the time before the next reaction, aka the weight is the dwell time of the system.
            slice = slice.iloc[1:, 2:].to_numpy()    
            averages[j, :, i] = np.average(a=slice, axis = 0, weights = dt)
        else:
            slice = df[df['Time (t)'].between(low_time, high_time)] # 2D df with (time, cell)
            averages[j, :, i] = slice.iloc[:, 2:].mean(0) # take average over the time axis

        

"""CONVERTING AVERAGES"""
'''To bound fraction'''
bound_fraction_sto = np.sum(averages[1:N_ARMS+1], 0)/N_PARTICLES # output: array(time, cell)

'''To selectivity'''
selectivity_sto = np.zeros((N_TIMES, N_receptors.size))
for i in range(N_TIMES):
    logfraction = np.log(bound_fraction_sto[:,i])
    logreceptors = np.log(N_receptors) 
    selectivity_sto[i] = np.gradient(np.log(bound_fraction_sto[:,i]), np.log(N_receptors))


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
    print('i is ', i)
    results[:, i] = res.y
    bound_fraction_num[i] = np.sum(res.y[1:-1], 0)/N_PARTICLES

'''converting to selectivity'''
selectivity_num = np.zeros((N_TIMES, N_receptors.size))

for i in range(sim_times.size):
    selectivity_num[i] = np.gradient(np.log(bound_fraction_num[:,i]), np.log(N_receptors))

"""PLOTTING"""
'''selectivity over receptors semilog for many simulation times'''
'''Legend is methods together'''
if N_TIMES == 1:
    axs[1,0].semilogx(N_receptors, selectivity_num[i], 
                      marker = marker_num, 
                      label = f"Numerical", 
                      c = color_num[i],
                      markersize = marker_size_num)
    axs[1,0].semilogx(N_receptors, 
                      selectivity_sto[i], 
                      marker = marker_sto, 
                      label = f"Stochastic", 
                      c = color_sto[i],
                      markersize = marker_size_sto)
else:
    for i in range(sim_times.size):
        time = format(sim_times[i], '.2f')
        axs[1,0].semilogx(N_receptors, 
                          selectivity_num[i], 
                          marker = marker_num, 
                          label = f"RRE - {time}s", 
                          c = color_num[i],
                          markersize = marker_size_num)
    for i in range(sim_times.size):
        time = format(sim_times[i], '.2f')
        axs[1,0].semilogx(N_receptors, 
                          selectivity_sto[i], 
                          marker = marker_sto, 
                          label = f"SSA - {time}s", 
                          c = color_sto[i],
                          markersize = marker_size_sto)
        

# plt.figtext(0.2, 0.2, f"\
#     kon_sol = {kon_sol} \n\
#     kon_sur = {kon_sur} \n\
#     koff = {koff} \n\
#     Number of arms = {N_ARMS} \n\
#     Number of particles = {initial_values[0]} \n\n\
#     ODE solver = {solver} \n\n\
#     Averaged over {AVERAGING_TIME} seconds")
axs[1,0].set_xlabel(r"Number of receptors $N_R$")
axs[1,0].set_ylabel(r"Selectivity $\alpha$")
axs[1,0].axhline(0, ls=':', c='grey')
axs[1,0].axvline(0, ls=':', c='grey')
axs[1,0].grid()

# '''Difference at endpoint'''
# difference = selectivity_sto[-1]-selectivity_num[-1]
# axs[0,1].semilogx(N_receptors, difference, marker = marker_num, color = 'k')
# axs[0,1].set_xlabel("Number of receptors")
# axs[0,1].set_ylabel(f"Selectivity")
# axs[0,1].axhline(0, ls=':', c='grey')
# axs[0,1].axvline(0, ls=':', c='grey')
# axs[0,1].legend()
# axs[0,1].grid()

# '''Difference over receptors'''
# for i in range(N_TIMES):
#     difference = -selectivity_sto[i, round(cells*R_range_sf):]+selectivity_num[i, round(cells*R_range_sf):]
#     axs[0,1].semilogx(N_receptors[round(cells*R_range_sf):], 
#                       difference, 
#                       label = f'{format(sim_times[i], '.2f')}s',
#                       marker = marker_num, 
#                       color = color[i],
#                       markersize = marker_size_num)
# axs[0,1].set_xlabel("Number of receptors")
# axs[0,1].set_ylabel(f"Selectivity difference N-S")
# axs[0,1].axhline(0, ls=':', c='grey')
# axs[0,1].axvline(0, ls=':', c='grey')
# axs[0,1].legend()
# axs[0,1].grid()

'''Bound fraction over receptors'''
# if show_species:
    # for i_species in range(1, N_ARMS+1):
    #     axs[0,1].loglog(N_receptors, 
    #                     averages[i_species, :, -1]/N_PARTICLES, 
    #                     label = f"{i_species} arms bound - SSA",
    #                     color = color_sto[i_species])
    #     axs[0,1].loglog(N_receptors, 
    #                     results[i_species, :, -1]/N_PARTICLES, 
    #                     label = f"{i_species} arms bound - Numerical",
    #                     color = color_num[i_species])
    # axs[0,1].legend()

    

if N_TIMES == 1:
    axs[1,1].loglog(N_receptors, 
                    bound_fraction_num[:,-1], 
                    marker = marker_num, 
                    label = f"Numerical", 
                    color = color_num[0],
                    markersize = marker_size_num)
    axs[1,1].loglog(N_receptors, 
                    bound_fraction_sto[:,-1], 
                    marker = marker_sto, 
                    label = f"Stochastic", 
                    color = color_sto[0],
                    markersize = marker_size_sto)
else:
    # axs[0,1].loglog(N_receptors, bound_fraction_num[:,-1], marker = marker_num, label = f"Numerical", color = 'purple')
    # axs[0,1].loglog(N_receptors, bound_fraction_sto[:,-1], marker = marker_sto, label = f"Stochastic", color = 'green')
    for i in range(N_TIMES):
        axs[1,1].loglog(N_receptors, 
                        bound_fraction_num[:,i], 
                        marker = marker_num, 
                        label = f"Numerical", 
                        color = color_num[i],
                        markersize = marker_size_num)
    for i in range(N_TIMES):
        axs[1,1].loglog(N_receptors, 
                        bound_fraction_sto[:,i], 
                        marker = marker_sto, 
                        label = f"Stochastic", 
                        color = color_sto[i],
                        markersize = marker_size_sto)

#Free receptors
twin = axs[1,1].twinx()
twin.semilogx(N_receptors, 
              results[-1, :, -1]/N_receptors, 
              color = 'purple',
              linestyle = line_receptors)
twin.set(ylim = (0,1.1), ylabel = "Free receptor fraction")
twin.yaxis.label.set_color('purple')
twin.tick_params(axis = 'y', colors = 'purple')

# fully_bound_fraction_s = averages[N_ARMS, :, -1]/np.sum(averages[1:N_ARMS+1, :, -1], axis=0)
# fully_bound_fraction_n = results[N_ARMS, :, -1]/np.sum(results[1:N_ARMS+1, :, -1], axis=0)

# axs2 = axs[1,1].twinx()
# axs2.semilogx(N_receptors, 
#                 fully_bound_fraction_s, 
#                 label = f"fully bound fraction - Stochastic",
#                 marker = marker_sto,
#                 color = color_sto[-1]) 
# axs2.semilogx(N_receptors, 
#                 fully_bound_fraction_n, 
#                 label = f"fully bound fraction - Numerical",
#                 marker = marker_num,
#                 color = color_num[-1]) 
# axs2.legend()

axs[1,1].set_xlabel(r"Number of receptors $N_R$")
axs[1,1].set_ylabel(r"Bound fraction $\theta$")
axs[1,1].axhline(0, ls=':', c='grey')
axs[1,1].axvline(0, ls=':', c='grey')
axs[1,1].grid()

# '''Bound fraction over receptors difference'''
# for i in range(sim_times.size):
#     percentage_difference = 100*(-bound_fraction_sto[:,i]+bound_fraction_num[:,i])/((bound_fraction_num[:,i]+bound_fraction_sto[:,i])/2)
#     time = format(sim_times[i], '.2f')
#     axs[1,1].semilogx(N_receptors, 
#                       percentage_difference, 
#                       marker = marker_num, 
#                       color = color[i], 
#                       label = f'{time}s')

# # percentage_difference = 100*(-bound_fraction_sto[:,-1]+bound_fraction_num[:,-1])/((bound_fraction_num[:,-1]+bound_fraction_sto[:,-1])/2)
# # axs[1,1].semilogx(N_receptors, percentage_difference, marker = marker_num, color = 'k')
# axs[1,1].set_xlabel("Number of receptors")
# axs[1,1].set_ylabel(f"Percentage difference (%) N-S")
# axs[1,1].axhline(0, ls=':', c='grey')
# axs[1,1].axvline(0, ls=':', c='grey')
# axs[1,1].legend()
# axs[1,1].grid()
# # plt.show()





# """GRAPH OVER TIME"""
# sim_times = df['Time (t)'] # array of all times
# cell = 19
# print(N_receptors[cell])
# concentrations_sto = np.zeros((N_ARMS+2, sim_times.size))

# '''Particles'''
# for bound_arms in range(0, N_ARMS+1):

#     file_path = f'.\\{parent_folder}\\{folder_name}\\Testrun_{N_ARMS}_arms_concentration_of_species_n{bound_arms}0.dat'
#     df = pd.read_csv(file_path, sep=r'\s{3,}', skiprows = 1, engine = 'python')
    
#     concentration = df.iloc[:, cell+2]
#     concentrations_sto[bound_arms] = concentration # df(time, cell) take the values of all times at a specific cell.

# '''To bound fraction'''
# bound_fraction_sto = np.sum(concentrations_sto[1:-1], 0)/N_PARTICLES # output: array(time, cell)

# """Solve system of ODEs"""
# bound_fraction_num = np.zeros((N_receptors.size, sim_times.size)) # bound_fraction[N, t] so N on the y-axis and t on the x-axis
# results = np.zeros((initial_values.size, N_receptors.size, sim_times.size)) # results [type, N, t] with the order of the types being the same as the initial values (N0, N1, N2, N3, R)

# tmax = tmax*2

# for i in range (N_receptors.size):
#     initial_values[-1] = N_receptors[i]
#     res = solve_ivp(rhs, (tmin, tmax), initial_values, t_eval = sim_times, method=solver)
#     print('i is ', i)
#     results[:, i] = res.y
#     bound_fraction_num[i] = np.sum(res.y[1:-1], 0)/N_PARTICLES

# '''plot bound fraction over time'''
# print('N receptors is ', N_receptors)
# plt.figure(figsize=(5,4))
# plt.plot(sim_times, bound_fraction_num[cell], label = 'Numerical')
# plt.plot(sim_times, bound_fraction_sto, label = 'Stochastic')
# plt.xlabel("Time (s)")
# plt.ylabel("Bound fraction")
# plt.axhline(0, ls=':', c='grey')
# plt.axvline(0, ls=':', c='grey')
# plt.legend()
# plt.grid()
# plt.show()





"""              """
"""  ANALYTICAL  """
"""              """
# A = 22e6
rho_0 = 1e-8
rho_0 = 1e100
k = N_ARMS

rho_0 = N_PARTICLES/VOLUME
rho_n = rho_0 #? Is this correct?

K_intra = 1/AREA*kon_sur/koff
K_A = 1/rho_0*kon_sol/koff

K_intra = kon_sur/koff*(1/AREA)
K_A = 1/rho_0*kon_sol/koff*(1/VOLUME)

'''KA and Kintra should be unrelated to things like concentration and surface, no?'''
# rho_n = 1e-8
# A = 22e6

def calc_bound_fraction_analytical(sigma_R):
    Kav_A = K_A/K_intra*( (1+sigma_R*AREA*K_intra)**N_ARMS - 1)
    bound_fraction = rho_n*Kav_A/(1+rho_n*Kav_A)
    return bound_fraction

sigma_R = N_receptors/AREA
bound_fraction_analytical = calc_bound_fraction_analytical(sigma_R)

axs[1,1].loglog(N_receptors, bound_fraction_analytical, label = f"Analytical - steady state", color = color_ana, linestyle = line_ana)




# temp = scipy.differentiate.derivative(calc_bound_fraction_analytical, N_receptors)
# selectivity_analytical = temp.df
selectivity_analytical = np.gradient(np.log(bound_fraction_analytical), np.log(N_receptors))
axs[1,0].semilogx(N_receptors, selectivity_analytical, label = f"Analytical - steady state", color = color_ana, linestyle = line_ana)

"""Subplot settings"""
axs[0,0].title.set_text(r"Selectivity comparison full $N_R$ range")
axs[0,1].title.set_text(r"Bound fraction comparison full $N_R$ range")
if 'break_iter 1' in folder_name2:
    axs[1,0].title.set_text(r"Selectivity comparison low $N_R$ range")
    axs[1,1].title.set_text(r"Bound fraction comparison low $N_R$ range")
else:
    axs[1,0].title.set_text(r"Selectivity comparison low $N_R$ range")
    axs[1,1].title.set_text(r"Bound fraction comparison low $N_R$ range")

'''legend for both graphs'''
handles, labels = axs[1,0].get_legend_handles_labels()
f.legend(handles,
         labels, 
         loc='upper right', 
         bbox_to_anchor=(0.99, 0.4),
         fancybox=False, 
         shadow=False, 
         ncols=1)


plt.subplots_adjust(left=0.048, right=0.8, top = 0.95, wspace=0.22)

plt.tight_layout()
plt.subplots_adjust(left=0.048, right=0.71, wspace=0.31)


plt.show()