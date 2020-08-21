from scipy.stats import norm
from math import sqrt, exp, log, pow

def h1_plus (r, q, sigma, S_t, M_t, beta, tau):

    f_ext = 1/(sigma*sqrt(tau))
    s1 = log(S_t/(beta*M_t))
    s2 = (r - q + 0.5*sigma*sigma)*tau

    return f_ext*(s1 + s2)

def h1_minus (r, q, sigma, S_t, M_t, beta, tau):

    f_ext = 1/(sigma*sqrt(tau))
    s1 = log(S_t/(beta*M_t))
    s2 = (r - q - 0.5*sigma*sigma)*tau

    return f_ext*(s1 + s2)

def h2_plus (r, q, sigma, S_t, M_t, beta, tau):

    f_ext = 1/(sigma*sqrt(tau))
    s1 = log(M_t/(beta*S_t))
    s2 = (r - q - 0.5*sigma*sigma)*tau

    return f_ext*(s1 + s2)

def h2_minus (r, q, sigma, S_t, M_t, beta, tau):

    f_ext = 1/(sigma*sqrt(tau))
    s1 = log(M_t/(beta*S_t))
    s2 = (r - q + 0.5*sigma*sigma)*tau

    return f_ext*(s1 - s2)

def FLB_Put_Pricing_Formula (r, q, sigma, S_t, M_t, beta, tau):

    gamma = 2*(r-q)/(sigma*sigma)
    h1_plus_ = h1_plus(r, q, sigma, S_t, M_t, beta, tau)
    h1_minus_ = h1_minus(r, q, sigma, S_t, M_t, beta, tau)
    h2_plus_ = h2_plus(r, q, sigma, S_t, M_t, beta, tau)
    h2_minus_ = h2_minus(r, q, sigma, S_t, M_t, beta, tau)

    s1 = beta*M_t*exp(-r*tau)*norm.cdf(-h1_minus_)
    s2 = S_t*exp(-q*tau)*norm.cdf(-h1_plus_)

    first_part = s1 - s2

    if r == q:
        f_ext = beta*S_t*exp(-r*tau)*sigma*sqrt(tau)
        s3 = norm.cdf(-h2_plus_)
        s4 = h2_plus*norm.cdf(-h2_plus_) 
    else:
        f_ext = beta*S_t/gamma
        s3 = exp(-q*tau)*pow(beta,gamma)*norm.cdf(-h2_minus_)
        s4 = exp(-r*tau)*pow(M_t/S_t, gamma)*norm.cdf( -h2_plus_)
    
    second_part = f_ext*(s3 - s4)

    return first_part + second_part