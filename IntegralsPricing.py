from math import sqrt, exp, log, erf, erfc

# I_1
def I_1_1 (a, sigma, p , l, tau):

    f1 = exp( (2*a + sigma)*(2*p + sigma*tau)/2 )
    f2 = 1 + erf( (p - l + tau*(a + sigma))/sqrt(2*tau) )

    return f1*f2

def I_1_2 (a, sigma, p , l, tau):

    f_ext = exp( -(2*a + sigma)*(a*tau -p) )/(2 + sigma/a)
    s1 = exp( tau*(2*a + sigma)*(2*a + sigma)/2)*(1 + erf( (-l + p + tau*(a + sigma))/sqrt(2*tau)) )
    s2 = exp( (2*a + sigma)*(l - p + a*tau) )*erfc( (l - p + a*tau)/sqrt(2*tau) )

    return f_ext*(s1 - s2)

def I_1 (a, sigma, p , l, tau):
    
    return I_1_1(a, sigma, p , l, tau) - I_1_2(a, sigma, p , l, tau)

# I_2
def I_2_1 (a, p , l, tau):

    return 0.5*(1 - erf((a*tau - l - p)/sqrt(2*tau)))

def I_2_2 (a, p , l, tau):

    f_ext = exp(2*a*l)
    
    return 0.5*f_ext*(1 - erf((a*tau + l - p)/sqrt(2*tau)))

def I_2 (a, p , l, tau):
    
    return I_2_1(a, p, l, tau) - I_2_2(a, p , l, tau)

#I_3
def I_3_1 (a, sigma, p , l, tau):

    f_ext = exp(sigma*tau*(2*a + sigma)/2)

    return 0.5*f_ext*(1 - erf((tau*(a + sigma) - l - p)/sqrt(2*tau)))

def I_3_2 (a, sigma, p , l, tau):

    f_ext = exp(sigma*tau*(2*a + sigma)/2 + 2*p*(a+sigma))

    return 0.5*f_ext*(1 + erf((tau*(a + sigma) - l + p)/sqrt(2*tau)))

def I_3 (a, sigma, p , l, tau):

    return I_3_1(a, sigma, p , l, tau) + I_3_2(a, sigma, p , l, tau)

# put price
def FLB_Put_Pricing_Integrals (r, q, sigma, S_t, M_t, beta, tau):
    
    a = (r - q)/sigma - sigma/2
    l = log(M_t/S_t)/sigma

    if beta < 1:
        p = log(beta)/sigma
    else:
        p = 0

    put_price = exp(-r*tau)*S_t*(beta*I_1(a, sigma, p, l, tau) + beta*(M_t/S_t)*I_2(a, p, l, tau) - I_3(a,sigma, p, l, tau))

    return put_price