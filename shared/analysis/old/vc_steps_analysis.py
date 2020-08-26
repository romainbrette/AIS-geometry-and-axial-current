"""
Functions to analyze VC steps.
"""

from brian2 import * 
from scipy import special
from math import factorial
from scipy.optimize import minimize

def cumulative_gaussian_distribution(x, mu, sigma):
    """
    Cumulative gaussian distribution with mean mu and std sigma.
    x: command potentials obtained with the staircase method
    """
    cgb = 0.5*(1 + special.erf((x-mu)/sqrt(2*sigma**2)))
    return cgb

def log_lik(params, x):
    """
    A function that computes the log-likelihood function for our data (3D function)
    params: [mu, sigma, delta]
    x[0]: command potentials obtained with the staircase method
    x[1]: spike or not (spike: 1, no spike: 0)
    """
    llik = 0
    for i in range(len(x[0])):        
        llik += log(1./(factorial(x[1][i]) * factorial(1-x[1][i]))) \
                + x[1][i] * cumulative_gaussian_distribution(x[0][i], params[0], params[1]) \
                + (1-x[1][i]) * (1-cumulative_gaussian_distribution(x[0][i], params[0], params[1])) 
    return -llik