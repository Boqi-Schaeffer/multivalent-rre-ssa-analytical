from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import numpy as np
from decimal import Decimal

"""SYSTEM"""
N_ARMS = 3
N_PARTICLES = 600
i_set = 0
R_range_sf = 0
r_range_name = "full"

"""Rate constants k"""
k = np.zeros((N_ARMS+1, N_ARMS+1))

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

for i in range(k.shape[0]):
    for j in range(k.shape[1]):
        k[i,j] = 1
        if i == 0 and j == 1:
            k[i,j] *= kon_sol # reaction of a surface molecule with a solution molecule
        elif j-i == 1: # if it's a binding reaction
            k[i,j] *= kon_sur # reaction of a surface molecule with a surface molecule
        elif j-i == -1: # if it's an unbinding reaction
            k[i,j] *= koff

'''setting rate constants corrected for number of arms'''
c = np.zeros(k.shape)

for i in range(k.shape[0]):
    for j in range(k.shape[1]):
        if i<j:
            c[i,j] = k[i,j]*(N_ARMS - i)  # to bind an arm, multiply k by number of unbound arms
        elif i>j:
            c[i,j] = k[i,j]*i # to release an arm, multiply k by number of bound arms

"""PLOTTING CONSTANTS"""
'''time axis'''
tmin = 0 # When does time start running, NOT the lowest simulated time I graph my selectivity for.
tmax = 4000
# tmax = 57000
# tmax = 2110000
n_times = 6
sim_times = np.arange(tmax/n_times, tmax + tmax/n_times, tmax/n_times)
sim_times = np.logspace(np.log10(tmax), np.log10(tmax*10**(n_times-1)), n_times)
tmax = sim_times[-1]

def format_e(n):
    a = '%E' % n
    return a.split('E')[0].rstrip('0').rstrip('.') + 'E' + a.split('E')[1]

'''receptor axis'''
steps_pp_10 = 16 # steps per power of 10
p_N_receptors_min = 2
p_N_receptors_max = 7 # power of maximum number of receptors
data_points = steps_pp_10*(p_N_receptors_max-p_N_receptors_min)+1

N_receptors = np.logspace(p_N_receptors_min, p_N_receptors_max, data_points)

sigma_R = N_receptors/AREA

'''other'''
color = plt.cm.Blues(np.linspace(0.3, 0.85, n_times))
line_ana = (0, (5,5))
color_ana = 'yellowgreen'
line_receptors = 'dashdot'
color_receptors = 'purple'
linewidth = 2.5
marker_num = 'o'

figheight = 7.5
figwidth = 16

plt.rcParams.update({'font.size': 16})
plt.rc('legend', fontsize=15) 
plt.rcParams['axes.titlesize'] = 18 

"""INITIALIZING ARRAYS"""
initial_values = np.zeros(N_ARMS + 2)
initial_values[0] = N_PARTICLES
bound_fraction = np.zeros((N_receptors.size, sim_times.size)) # bound_fraction[N, t] so N on the y-axis and t on the x-axis
results = np.zeros((initial_values.size, N_receptors.size, sim_times.size)) # results [type, N, t] with type being the same as the initial values (N0, N1, N2, N3, R)

"""THE SYSTEM OF ODEs"""
def rhs(t, f):
    ODE = np.zeros(N_ARMS+2)
    for i in range(N_ARMS+1): # we make the ODE system components by going through each state and adding the change to all states that a reaction from that state induces
        if i<N_ARMS:
            #binding reaction
            change = c[i,i+1]*f[-1]*f[i]
            ODE[i] -= change #subtract from i_arms state
            ODE[i+1] += change #add to i+1_arms state
            ODE[-1] -= change #subtract from number of receptors

        if i>0:
            #unbinding reaction
            change = c[i,i-1]*f[i]
            ODE[i] -= change #subtract from i_arms state
            ODE[i-1] += change #add to i-1_arms state
            ODE[-1] += change #add to number of receptors
    return ODE

"""SOLVING THE SYSTEM"""
for i in range (N_receptors.size):
    initial_values[-1] = N_receptors[i]
    res = solve_ivp(rhs, (tmin, tmax), initial_values, t_eval = sim_times, method='LSODA')
    print('i is ', i)
    results[:, i] = res.y
    bound_fraction[i] = np.sum(res.y[1:-1], 0)/initial_values[0]

"""PLOTTING"""
f, axs = plt.subplots(2, 2)
f.set_figheight(figheight)
f.set_figwidth(figwidth)

'''bound fraction'''
# plt.figure(figsize=(5,4))
for i in range(sim_times.size):
    time = format_e(sim_times[i])
    axs[1,0].semilogx(N_receptors, bound_fraction[:, i], label = f'RRE - {time}s', c = color[i])

axs[1,0].set_xlabel("Number of receptors")
axs[1,0].set_ylabel("Bound fraction")
axs[1,0].axhline(0, ls=':', c='grey')
axs[1,0].axvline(0, ls=':', c='grey')
axs[1,0].legend()
axs[1,0].grid()

'''loglog bound'''
# plt.figure(figsize=(5,4))
for i in range(sim_times.size):
    time = format_e(sim_times[i])
    axs[0,1].loglog(N_receptors, 
                    bound_fraction[:, i], 
                    c = color[i],
                    label = f'RRE - {time}s',
                    lw = linewidth)

axs[0,1].set_xlabel(r"Number of receptors $N_R$")
axs[0,1].set_ylabel(r"Bound fraction $\theta$")
axs[0,1].axhline(0, ls=':', c='grey')
axs[0,1].axvline(0, ls=':', c='grey')
axs[0,1].grid()



'''selectivity over receptors semilog for many simulation times'''
width = 15
height = 7
# plt.figure(figsize=(width,height))
for i in range(sim_times.size):
    time = format_e(sim_times[i])
    selectivity = np.gradient(np.log(bound_fraction[:,i]), np.log(N_receptors))
    axs[0,0].semilogx(N_receptors, selectivity, label = f'RRE - {time}s', c = color[i],
                      lw = linewidth)
    # axs[0,0].semilogx(N_receptors, selectivity, label = f'{N_ARMS} arms', c = color[i])

plt.figtext(0.2, 0.2, f"\
    kon_sol = {kon_sol} \n\
    kon_sur = {kon_sur} \n\
    koff = {koff} \n\
    number of particles = {initial_values[0]} \n\
    ODE solver = {'LSODA'}")
axs[0,0].set_xlabel(r"Number of receptors $N_R$")
axs[0,0].set_ylabel(r"Selectivity $\alpha$")
axs[0,0].axhline(0, ls=':', c='grey')
axs[0,0].axvline(0, ls=':', c='grey')
axs[0,0].grid()

"""Free receptors"""
twin = axs[0,1].twinx()
twin.semilogx(sigma_R, 
              results[-1, :, -1]/AREA/sigma_R, 
              color = 'purple',
              linestyle = line_receptors)
twin.set(ylim = (0,1.1), ylabel = "Free receptor fraction")
twin.yaxis.label.set_color('purple')
twin.tick_params(axis = 'y', colors = 'purple')



"""              """
"""  ANALYTICAL  """
"""              """

'''defining variables as given in (Linne, 2022)'''
rho_0 = N_PARTICLES/VOLUME
k = N_ARMS

## Only works for A = V = 1
K_intra = 1/AREA*kon_sur/koff
K_A = 1/rho_0*kon_sol/koff

## Works for all A and V
K_intra = kon_sur/koff
K_A = 1/rho_0*kon_sol/koff

def calc_bound_fraction_analytical(sigma_R):
    Kav_A = K_A/K_intra*( (1+sigma_R*AREA*K_intra)**N_ARMS - 1)
    bound_fraction = rho_0*Kav_A/(1+rho_0*Kav_A)
    return bound_fraction

sigma_R = N_receptors/AREA
bound_fraction_analytical = calc_bound_fraction_analytical(sigma_R)

axs[0,1].loglog(N_receptors, 
                bound_fraction_analytical, 
                label = f"Analytical", 
                color = color_ana, 
                linestyle = (0,(5,5)),
                lw = linewidth)

selectivity_analytical = np.gradient(np.log(bound_fraction_analytical), np.log(sigma_R))
axs[0,0].semilogx(N_receptors, 
                  selectivity_analytical, 
                  label = f"Analytical - steady state", 
                  color = color_ana, 
                  linestyle = (0,(5,5)),
                  lw = linewidth)

"""Legend for both graphs"""
handles, labels = axs[0,0].get_legend_handles_labels()
f.legend(handles,
         labels, 
         loc='lower right', 
         bbox_to_anchor=(0.99, 0.57),
         fancybox=False, 
         shadow=False, 
         ncols=1)
plt.subplots_adjust(left=0.060, right=0.77, wspace=0.3, hspace=0.3)

axs[0,0].set_title("Selectivity comparison - long equilibration")
axs[0,1].set_title("Bound fraction comparison - long equilibration")


plt.subplots_adjust(left=0.048, right=0.71, wspace=0.31)


plt.show()