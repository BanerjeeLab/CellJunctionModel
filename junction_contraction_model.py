#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 17 13:53:18 2021

@author: staddon
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['font.sans-serif'] = ['arial']
matplotlib.rcParams['font.size'] = 16
matplotlib.rcParams['axes.linewidth'] = 2.25
matplotlib.rcParams['xtick.major.width'] = 2.25
matplotlib.rcParams['ytick.major.width'] = 2.25
matplotlib.rcParams['figure.figsize'] = [4, 4]


def solve(t, nx=21, tension=0.4, mu=2.5, E=1, k_top=1, k_bot=1):
    
    y0 = np.linspace(0, 1, nx)
    
    dy = y0[1] - y0[0]
    
    # Equation of motion
    def dydt(t, y):
            
        f = np.zeros(y.shape)
        stress = np.ones(y.shape) * tension * (t >= 0) * (t <= 5)
        friction = np.ones(y.shape) * mu
        
        # Junction strain
        f[1:-1] = E * (y[:-2] - 2 * y[1:-1] + y[2:]) / dy ** 2
        
        # Active contraction
        f[1:-1] += -(stress[:-2] - stress[2:]) / 2 / dy
        
        
        f[0] = -k_bot * y[0] + stress[0] + E * (y[1] - y[0] - dy) / dy
        f[-1] = -k_top * (y[-1] - 1) - stress[-1] - E * (y[-1] - y[-2] - dy) / dy
    
        return f / friction
    
    sol = solve_ivp(dydt, [t.min(), t.max()], y0, max_step=0.1, t_eval=t)
    
    return sol.t, sol.y


if __name__ == '__main__':
    
    
    colors = [plt.cm.magma(0.25), plt.cm.magma(0.75)]
    markers = ['s', 'o']
    
    nx = 21
    
    t = np.linspace(-5, 10)
    
    
    """
    No feedback
    """
    
    # Rough experimental values
    lam = np.linspace(1, 2, nx)
    mu = 40.1
    k = 4.52
    
    t, y = solve(t, nx=nx, tension=lam, mu=mu, E=1, k_bot=k, k_top=k)
    
    # Chymograph
    for i in range(y.shape[0]):
        color = plt.cm.magma(i / y.shape[0] * 0.75 + 0.25)
        plt.plot(t, y[i, :], color=color, lw=3)
        
    plt.axvline(0, ls=':', color=(0.75, 0.75, 0.75), lw=3)
    plt.axvline(5, ls=':', color=(0.75, 0.75, 0.75), lw=3)
    plt.ylabel('Normalized Position')
    plt.xlabel('Time (mins)')
    plt.ylim(0, 1)
    plt.xlim(-5, 10)
    plt.title('No Feedback - Asymmetric')
    plt.show()
    
    
    
    """
    Tension feedback model
    """

    # Feedback parameters
    lam_m = 1
    lam_lm = 2
    
    mu_m = 20.3
    mu_lm = 77.3
    
    k_m = 0.996
    k_lm = 17.5
    
    lam = np.linspace(lam_m, lam_lm, nx)
    mu = mu_m * (mu_lm / mu_m) ** (lam - 1)

    t, y = solve(t, nx=nx, tension=lam, mu=mu, E=1, k_bot=k_m, k_top=k_lm)

    # Chymograph
    for i in range(y.shape[0]):
        color = plt.cm.magma(i / y.shape[0] * 0.75 + 0.25)
        plt.plot(t, y[i, :], color=color, lw=3)
        
    plt.axvline(0, ls=':', color=(0.75, 0.75, 0.75), lw=3)
    plt.axvline(5, ls=':', color=(0.75, 0.75, 0.75), lw=3)
    plt.ylabel('Normalized Position')
    plt.xlabel('Time (mins)')
    plt.ylim(0, 1)
    plt.xlim(-5, 10)
    plt.title('Feedback model')
    plt.show()
    
    
    """
    Half activation
    """
    
    # Feedback parameters
    lam_m = 1
    lam_lm = 2
    
    mu_m = 20.3
    mu_lm = 77.3
    
    k_m = 0.996
    k_lm = 17.5

    lam = np.linspace(0, 0, nx)
    lam[:11] = 1
    
    mu = mu_m * (mu_lm / mu_m) ** (lam - 1)
    
    k_top = k_m * (k_lm / k_m) ** (0 - lam_m)
    k_bot = k_m * (k_lm / k_m) ** (1 - lam_m)

    t, y = solve(t, nx=nx, tension=lam, mu=mu, E=1, k_top=k_top, k_bot=k_bot)
    
    # Chymograph
    for i in range(y.shape[0]):
        color = plt.cm.magma(lam[i] * 0.5 + 0.25)
        plt.plot(t, y[i, :], color=color, lw=3)
        
    plt.axvline(0, ls=':', color=(0.75, 0.75, 0.75), lw=3)
    plt.axvline(5, ls=':', color=(0.75, 0.75, 0.75), lw=3)
    plt.ylabel('Normalized Position')
    plt.xlabel('Time (mins)')
    plt.ylim(0, 1)
    plt.xlim(-5, 10)
    plt.title('Feedback model - half')
    plt.show()
