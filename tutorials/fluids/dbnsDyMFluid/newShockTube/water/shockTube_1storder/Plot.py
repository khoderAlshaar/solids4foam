import numpy as np
import matplotlib.pyplot as plt
from sympy import *
from pathlib import Path

class state:
    def __init__(self, T, u, P, gm, P0, Cv):
        self.T = T
        self.u = u
        self.P = P
        self.gm = gm
        self.P0 = P0
        self.Cv = Cv
        self.rho = (self.P + self.P0)/(self.T*self.Cv*(self.gm - 1))
        self.a = sqrt(self.gm*(self.P + self.P0)/self.rho)

    @classmethod
    def rho_P_u(cls, rho, u, P, gm, P0, Cv):
        return cls((P + P0)/(rho*Cv*(gm - 1)), u, P, gm, P0, Cv)

def F_k(P_star, state_k):
    p = symbols('p')
    if (P_star >= state_k.P):
      A = 2/((state_k.gm + 1)*state_k.rho)
      B = ((state_k.gm -1)/(state_k.gm + 1))*(state_k.P + state_k.P0)
      F = (p - state_k.P)*sqrt(A/(p + state_k.P0 + B))
    if (P_star < state_k.P):
      F = 2*state_k.a/(state_k.gm -1)*( pow((p + state_k.P0)/(state_k.P + state_k.P0), (state_k.gm-1)/(2*state_k.gm)) - 1)
    return F

def finding_P_star(left, right, epsilon_min):
    p = symbols('p')
    P_star = right.P
    while(1):
        fL = F_k(P_star, left)
        fR = F_k(P_star, right)
        dfL = diff(fL, p)
        dfR = diff(fR, p)
        P_star_new = P_star - ((fL.subs(p, P_star) + fR.subs(p, P_star) + right.u - left.u)/(dfL.subs(p, P_star) + dfR.subs(p, P_star)))
        epsilon = (P_star_new - P_star)/P_star
        P_star =  P_star_new
        if(abs(epsilon)<epsilon_min):
            break
    return P_star

def finding_u_star(left, right, P_star):
    p = symbols('p')
    fL = F_k(P_star, left)
    fR = F_k(P_star, right)
    u_star = 0.5*(left.u + right.u) + 0.5*(fR.subs(p, P_star) - fL.subs(p, P_star))
    return u_star

def finding_rho_star_k(state_k, P_star):
    if(P_star>state_k.P):
      return state_k.rho*(((P_star + state_k.P0)/(state_k.P + state_k.P0 ) ) + (state_k.gm-1)/(state_k.gm+1))/(((state_k.gm-1)/(state_k.gm+1))*((P_star + state_k.P0) / (state_k.P + state_k.P0))+1);
    if(P_star<state_k.P):
      return state_k.rho*pow((P_star + state_k.P0)/(state_k.P + state_k.P0 ), 1/state_k.gm)  

def wave_speed_k(P_star, state_k, sign):
    if(P_star<=state_k.P):
        q_k = 1
    if(P_star>state_k.P):
        q_k = sqrt(1 + ((state_k.gm + 1)/(2*state_k.gm))*((P_star + state_k.P0)/(state_k.P + state_k.P0) - 1))
    S_k = state_k.u  + sign*state_k.a*q_k
    return S_k

def sampling_solution(left_star, right_star, left, right, X, x_m, t):
    P_star = 0.5*(left_star.P +  right_star.P)
    u_star = 0.5*(left_star.u +  right_star.u)
    S_L = wave_speed_k(P_star, left, -1)
    S_R = wave_speed_k(P_star, right, 1)
    T, P, rho, u = np.zeros_like(X), np.zeros_like(X), np.zeros_like(X), np.zeros_like(X)
    i = 0
    for x in X:
        speed_rel = (x - x_m)/t
        if(speed_rel<u_star):
             if(speed_rel<S_L):
                T[i], P[i], rho[i], u[i] = left.T, left.P, left.rho, left.u
             else:
                if(P_star>left.P):
                   T[i], P[i], rho[i], u[i] = left_star.T, P_star, left_star.rho, u_star
                if(P_star<left.P):
                   STL = left_star.u - left_star.a
                   if(speed_rel<STL):     
                        P[i] = (left.P + left.P0)*pow(2/(left.gm+1)+((left.gm-1)/((left.gm+1)*left.a))*(left.u - speed_rel),(2*left.gm)/(left.gm-1))  - left.P0
                        rho[i] = left.rho*pow(2/(left.gm+1)+((left.gm-1)/((left.gm+1)*left.a))*(left.u- speed_rel),2/(left.gm-1))
                        u[i] = (2/(left.gm+1))*(left.a + (left.gm-1)*left.u/2 + speed_rel)
                        T[i] =   (P[i] + left.P0)/(rho[i]*left.Cv*(left.gm - 1)) 
                   else:
                        T[i], P[i], rho[i], u[i] = left_star.T, P_star, left_star.rho, u_star
        else:
             if(speed_rel>S_R):
                T[i], P[i], rho[i], u[i] = right.T, right.P, right.rho, right.u
             else:
                if(P_star>right.P):
                   T[i], P[i], rho[i], u[i] = right_star.T, P_star, right_star.rho, u_star
                if(P_star<right.P):
                   STL = right_star.u + right_star.a
                   if(speed_rel>STL):        
                        P[i] = (right.P + right.P0)*pow(2/(right.gm+1)+((right.gm-1)/((right.gm+1)*right.a))*(right.u - speed_rel),(2*right.gm)/(right.gm-1)) - right.P0
                        rho[i] = right.rho*pow(2/(right.gm+1)+((right.gm-1)/((right.gm+1)*right.a))*(right.u- speed_rel),2/(right.gm-1))
                        u[i] = (2/(right.gm+1))*(-right.a + (right.gm-1)*right.u/2 + speed_rel)
                        T[i] =   (P[i] + right.P0)/(rho[i]*right.Cv*(right.gm - 1)) 
                   else:
                        T[i], P[i], rho[i], u[i] =  right_star.T, P_star, right_star.rho, u_star
        i +=1
    return rho, P, u, T

def Find_solution(left, right, X, x_0, time):
    P_star = finding_P_star(left, right, 1e-12)
    u_star = finding_u_star(left, right, P_star)
    rho_star_L = finding_rho_star_k(left, P_star)
    rho_star_R = finding_rho_star_k(right, P_star)
    left_star = state.rho_P_u(rho_star_L, u_star, P_star, left.gm, left.P0, left.Cv)
    right_star = state.rho_P_u(rho_star_R, u_star, P_star, right.gm, right.P0, right.Cv)
    rho, P, u, T = sampling_solution(left_star, right_star, left, right, X, x_0, time)
    return X, rho, P, u, T

def add_plots(ax, rho, P, u, T, X, color_, lineStyle_, label_):
    ax[0, 0].plot(X, rho,  color=color_, linestyle= lineStyle_,label=label_)
    ax[0, 1].plot(X, P/1e06, color=color_, linestyle= lineStyle_, label = label_)
    ax[1, 0].plot(X, u, color=color_, linestyle= lineStyle_, label = label_)
    ax[1, 1].plot(X, T, color=color_, linestyle= lineStyle_, label = label_)
def set_labels(ax):
    ax[0, 0].set_ylabel('Density (kg/m3)')
    ax[0, 1].set_ylabel('Pressure (MPa)')
    ax[1, 0].set_ylabel('Velocity (m/s)')
    ax[1, 1].set_ylabel('Temperature (K)')
    ax[0, 0].set_xlabel('Position (m)')
    ax[0, 1].set_xlabel('Position (m)')
    ax[1, 0].set_xlabel('Position (m)')
    ax[1, 1].set_xlabel('Position (m)')
    
    ax[1, 1].yaxis.get_major_formatter().set_scientific(False)
    ax[1, 1].yaxis.get_major_formatter().set_useOffset(False)

    ax[0, 0].legend(frameon=True, facecolor='white', edgecolor='black')
    ax[0, 1].legend(frameon=True, facecolor='white', edgecolor='black')
    ax[1, 0].legend(frameon=True, facecolor='white', edgecolor='black')
    ax[1, 1].legend(frameon=True, facecolor='white', edgecolor='black')
def read_data(p_rho_T, u_file):
    data = np.loadtxt(p_rho_T, delimiter=' ', skiprows=0)
    x, P, rho, T = data[:, 0], data[:, 1], data[:, 2], data[:, 3]
    data_u = np.loadtxt(u_file, delimiter=None, skiprows=0, unpack=True)
    u = data_u[1]
    return x, P, rho, T, u

####### Initial left and right conditions #####
# left = state(288.15, 0, 0.5e6, 2.782917482129208, 773221168.8164893, 1505.254836659763) 
# right = state(288.15, 0, 0.1e6, 2.782917482129208,773221168.8164893, 1505.254836659763) 

left = state(300, 0, 1e06, 2.8, 8.5e08, 1495) 
right = state(300, 0, 0.1e06, 2.8, 8.5e08, 1495) 

### coordinates for shock tube 
X = np.linspace(0, 3, 301)

### Finding the final solution time from case
postprocess_folder = "./postProcessing/sets"
base_dir = Path(postprocess_folder )
folder = next(d for d in base_dir.iterdir() if d.is_dir())

time = float(folder.name)
result_folder = postprocess_folder  + "/" +folder.name

### Finding the exact solution
x, rho, p, u, T = Find_solution(left, right, X, 2.0, time)

### Ploting numerical solution with exact solution
fig, ax = plt.subplots(2, 2, figsize=(10,8))
add_plots(ax, rho, p, u, T, X, 'red', '--', 'Exact solution')
X_data, P_data, rho_data, T_data, u_data = read_data(result_folder + '/uniform_p_rho_T.xy', result_folder + '/uniform_U.xy')
add_plots(ax, rho_data, P_data, u_data, T_data, X_data, 'blue', 'solid' , 'Numerical solution')
set_labels(ax)
fig.tight_layout()
fig.savefig("./Plot_comparison.png", dpi = 200)

