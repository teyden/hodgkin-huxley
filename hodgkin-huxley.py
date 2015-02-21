import numpy as np
import matplotlib.pyplot as plt

# Define Constants of HH model
## Reversal Potentials
E_Na = 115 #mV
E_K = -12 #mV
E_L = 10.6 #mV
## Channel Max Conductances
g_Na = 120 #mS/cm^2
g_K = 36 #mS/cm^2
g_L = 0.3 #mS/cm^2

# Define Functions of HH model
## Alpha Function for each n, m, and h
def alpha_n(u):
    return (0.1-0.01*u)/(np.exp(1-0.1*u)-1);
def alpha_m(u):
    return (2.5-0.1*u)/(np.exp(2.5-0.1*u)-1);
def alpha_h(u):
    return 0.07*np.exp(-u/20);
    
## Beta Function for each n, m, and h
def beta_n(u):
    return 0.125*np.exp(-u/80);
def beta_m(u):
    return 4*np.exp(-u/18);
def beta_h(u):
    return 1/(np.exp(3-0.1*u)+1);
    
## DEs of m, n, and h
def dn_dt(n,u):
    return alpha_n(u)*(1-n)-beta_n(u)*n;
def dm_dt(m,u):
    return alpha_m(u)*(1-m)-beta_m(u)*m;
def dh_dt(h,u):
    return alpha_h(u)*(1-h)-beta_h(u)*h;
    
## HH Model Equation for determining sum of total ion currents
def sigma_I_k(u,m,h,n):
    return g_Na*(m**3)*h*(u-E_Na)+g_K*(n**4)*(u-E_K)+g_L*(u-E_L);
## DE for u(t)
def du_dt(I_k,I_in,C):
    return (-I_k+I_in)/C;

# Initialize Variables
u = 0 # initial voltage across the membrane (mV)
t = 0 # initial time
delta_t = 0.001 # time step
runtime = 100 # total time to run experiment
I_k = 0 # sum of ionic currents
C = 1 # capacitance

# Initialize Gating Variables
n = alpha_n(u)/(alpha_n(u)+beta_n(u)) # controls probability of open K channels
m = alpha_m(u)/(alpha_m(u)+beta_m(u)) # controls probability of open Na channels
h = alpha_h(u)/(alpha_h(u)+beta_h(u)) # controls probability of open Na channels

# Define Input Current
def InputCurrent(I_type,t,I_init):
    # define sinusoidal input current
    if I_type == "sinusoidal":
        return 5*np.sin(t)+5
    # define 1ms input current spike at t=I_init
    if I_type == "single":
        if t < I_init:
            return 0
        else:
            if t > I_init + 1:
                return 0
            else:
                return 10
    # define constant input current of current I_init
    if I_type == "constant":
        return I_init

# Initialize lists for plotting
y_u = [0] # y-axis for membrane potential
y_I = [0] # y-axis for current
x = [0] # x-axis for time

# Use the Forward Euler Method to determine new values
while t < runtime:
    I_in = InputCurrent("constant",t,20) # what kind of input current are we plotting?
    n = n + delta_t*dn_dt(n,u) # new n, based on old n
    m = m + delta_t*dm_dt(m,u) # new m, based on old m
    h = h + delta_t*dh_dt(h,u) # new h, based on old h
    I_k = sigma_I_k(u,m,h,n) # new I_k based on new n, m, h
    u = u + delta_t*du_dt(I_k,I_in,C) # new u, based on old 
    t = t + delta_t # increase time by time step
    y_u.append(u) # save new value of u
    y_I.append(I_in)
    x.append(t) # save new value of t


counts = []
spikesV = []
for i in y_u :
    if y_u(i) > 0 & y_u(i-1) < 0:
        spikesV = spikesV + [y_u(i)]
        counts = counts + [1]

        
        
        


# Plot results
fig, ax = plt.subplots()

# Set up double y-axes
axes = [ax, ax.twinx()]
axes[1].spines['right'].set_position(('axes', 1))
axes[1].set_frame_on(True)
axes[1].patch.set_visible(False)

# Plot current results
axes[0].plot(x,y_u,x,[0]*len(x),'k--')
axes[0].set_ylabel('Membrane Potential', color='Blue')
axes[0].axis([0,runtime,-20,120])
axes[0].set_xlabel('Time')

axes[1].plot(x,y_I,'g')
axes[1].set_ylabel('Current', color='Green')
axes[1].axis([0,runtime,-20,120])
    
plt.show()


