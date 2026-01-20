from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import numpy as np

"""SYSTEM"""
N_ARMS = 3
N_A = 6.022e23 #avogadros number

"""TRYING TO MATCH THE PAPER"""
# To match 3 arm system from paper
AREA = 22e6 #micrometer squared
VOLUME = 0.1e3*AREA #parafilm is about 100 micrometers thick
rho_n = 1e-8
N_PARTICLES = rho_n*N_A*VOLUME*1e-15 # mol/L => particles/L => particles

# N_ARMS = 3
# K_intra = 5e-12
# K_A = 1.2e-5

N_ARMS = 6
K_intra = 1.2e-12
K_A = 1.2e-5

def format_e(n):
    a = '%E' % n
    return a.split('E')[0].rstrip('0').rstrip('.') + 'E' + a.split('E')[1]

"""Rate constants k"""
k = np.zeros((N_ARMS+1, N_ARMS+1))

'''based on behavior (on/off, surface/solution)'''
rho_0 = N_PARTICLES/VOLUME
rho_0 = rho_n


koff = 1

# from paper but does not work
kon_sol = K_A*koff*rho_0 # this scaling factor makes them match # moves bound fraction up and down
kon_sur = K_intra*AREA*koff # moves selevtivity peak left and right

# works
kon_sol = K_A*koff*rho_0 # this scaling factor makes them match # moves bound fraction up and down
kon_sur = K_intra*koff # moves selevtivity peak left and right


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
tmax = 1e-1
n_times = 5
sim_times = np.arange(tmax/n_times, tmax + tmax/n_times, tmax/n_times)
sim_times = np.logspace(np.log10(tmax), np.log10(tmax*10**(n_times-1)), n_times)
tmax = sim_times[-1]

'''receptor axis'''
steps_pp_10 = 16 # steps per power of 10
p_N_receptors_min = 3
p_N_receptors_max = 6 # power of maximum number of receptors
data_points = steps_pp_10*(p_N_receptors_max-p_N_receptors_min)+1

sigma_R = np.logspace(p_N_receptors_min, p_N_receptors_max, data_points)

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
bound_fraction = np.zeros((sigma_R.size, sim_times.size)) # bound_fraction[N, t] so N on the y-axis and t on the x-axis
results = np.zeros((initial_values.size, sigma_R.size, sim_times.size)) # results [type, N, t] with type being the same as the initial values (N0, N1, N2, N3, R)

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
    # ODE[-1] = 0
    return ODE

"""SOLVING THE SYSTEM"""
for i in range (sigma_R.size):
    initial_values[-1] = sigma_R[i]*AREA
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
    # time = format(sim_times[i], '.2f')
    time = format_e(sim_times[i])
    axs[1,0].semilogx(sigma_R, bound_fraction[:, i], label = f'RRE - {time}s', c = color[i])

axs[1,0].set_xlabel(r"Receptor density $\sigma_\mathrm{R}$ ($\mu m^{-2}$)")
axs[1,0].set_ylabel(r"Bound fraction $\theta$")
axs[1,0].axhline(0, ls=':', c='grey')
axs[1,0].axvline(0, ls=':', c='grey')
axs[1,0].grid()

'''loglog bound'''
# plt.figure(figsize=(5,4))
for i in range(sim_times.size):
    axs[0,1].loglog(sigma_R, bound_fraction[:, i], c = color[i], lw = linewidth)

axs[0,1].set_xlabel(r"Receptor density $\sigma_\mathrm{R}$ ($\mu m^{-2}$)")
axs[0,1].set_ylabel(r"Bound fraction $\theta$")
axs[0,1].axhline(0, ls=':', c='grey')
axs[0,1].axvline(0, ls=':', c='grey')
axs[0,1].grid()



'''selectivity over receptors semilog for many simulation times'''
width = 15
height = 7
# plt.figure(figsize=(width,height))
for i in range(sim_times.size):
    # time = format(sim_times[i], '.2f')
    time = format_e(sim_times[i])
    selectivity = np.gradient(np.log(bound_fraction[:,i]), np.log(sigma_R))
    axs[0,0].semilogx(sigma_R, selectivity, label = f"RRE - {time}s", c = color[i],
                      linewidth = linewidth)
    # axs[0,0].semilogx(N_receptors, selectivity, label = f'{N_ARMS} arms', c = color[i])

plt.figtext(0.2, 0.2, f"\
    kon_sol = {kon_sol} \n\
    kon_sur = {kon_sur} \n\
    koff = {koff} \n\
    number of particles = {N_PARTICLES} \n\
    ODE solver = {'LSODA'}")
axs[0,0].set_xlabel(r"Receptor density $\sigma_\mathrm{R}$ ($\mu m^{-2}$)")
axs[0,0].set_ylabel(r"Selectivity $\alpha$")
axs[0,0].axhline(0, ls=':', c='grey')
axs[0,0].axvline(0, ls=':', c='grey')
axs[0,0].grid()

for i in range(sim_times.size):
    axs[1,1].semilogx(sigma_R, results[-1, :, i]/AREA/sigma_R, c = color[i])
    # axs[1,1].loglog(N_receptors, results[-1, :, i]/N_receptors, c = color[i])

axs[1,1].set_xlabel(r"Receptor density $\sigma_\mathrm{R}$ ($\mu m^{-2}$)")
axs[1,1].set_ylabel("fraction of free receptors")
axs[1,1].axhline(0, ls=':', c='grey')
axs[1,1].axvline(0, ls=':', c='grey')
axs[1,1].grid()

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
    
def calc_bound_fraction_analytical(sigma_R):
    Kav_A = K_A/K_intra*( (1+sigma_R*AREA*K_intra)**N_ARMS - 1)
    bound_fraction = rho_n*Kav_A/(1+rho_n*Kav_A)
    return bound_fraction

bound_fraction_analytical = calc_bound_fraction_analytical(sigma_R)

axs[0,1].loglog(sigma_R, bound_fraction_analytical, 
                label = f"Analytical", 
                color = color_ana, 
                linestyle = line_ana, 
                linewidth = linewidth)


# temp = scipy.differentiate.derivative(calc_bound_fraction_analytical, N_receptors)
# selectivity_analytical = temp.df
selectivity_analytical = np.gradient(np.log(bound_fraction_analytical), np.log(sigma_R))
axs[0,0].semilogx(sigma_R, selectivity_analytical, label = f"Analytical - steady state", color = color_ana, linestyle = line_ana,
                  linewidth = linewidth)

"""Legend for both graphs"""
handles, labels = axs[0,0].get_legend_handles_labels()
f.legend(handles,
         labels, 
         loc='lower right', 
         bbox_to_anchor=(0.99, 0.57),
         fancybox=False, 
         shadow=False, 
         ncols=1)
plt.subplots_adjust(left=0.060, right=0.75, wspace=0.27, hspace=0.3)

axs[0,0].set_title("Selectivity comparison - in vitro params")
axs[0,1].set_title("Bound fraction comparison - in vitro params")

plt.subplots_adjust(left=0.048, right=0.71, wspace=0.31)

plt.show()