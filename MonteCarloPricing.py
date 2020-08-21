import numpy as np
from math import log, exp, sqrt
from scipy.stats import sem
from scipy.stats import norm
from statistics import stdev
from statistics import mean

# Brownian-Bridge Monte-Carlo 

### - Determines the payoff of a single sample path by simulating the final and minumum values
def FLBP_BSBrownianBridge_Simulation(S_t, M_t, beta, r, q, sigma, tau):

    S = S_t*exp( (r - q - sigma*sigma/2 )*tau + sigma*np.random.normal(0, sqrt(tau)) )
    u = np.random.uniform()
    
    w_t = log(S_t)
    w = log(S)

    w_max = 0.5*(w_t + w + sqrt( (w_t - w)*(w_t - w) - 2*sigma*sigma*tau*log(1-u) ))
    S_max = exp(w_max)

    running_max = max(S_max, M_t)
    payoff = max(beta*running_max - S, 0)
    discounting_factor = exp(-r*tau)

    return discounting_factor*payoff

### - Determine the "brownian-bridge" payoff under several simulation trials and calculates the empirical mean 
def FLBP_BSPricing_MonteCarloBB(n_simulations, S_t, M_t, beta, r, q, sigma, tau):

    results = []

    for i in range(0, n_simulations):
        results.append( FLBP_BSBrownianBridge_Simulation(S_t, M_t, beta, r, q, sigma, tau) )
    
    print('Mean obtained with Brownian-Bridge Monte-Carlo (' + str(n_simulations) + ' simulations): ' + str(mean(results)))
    print('Standard Deviation obtained with Brownian-Bridge Monte-Carlo: (' + str(n_simulations) + ' simulations): ' + str(stdev(results)))
    print('Standard Error obtained with Brownian-Bridge Monte-Carlo: (' + str(n_simulations) + ' simulations): ' + str(sem(results)))