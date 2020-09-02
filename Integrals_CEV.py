from numpy import array, linspace
from scipy import integrate
from CEV_Max import joint_density

def Integral_x(S_t, M_t, y, tau, b, a, beta, limsup, beta_frac, p):

    if y > p*M_t:
        liminf = y/p
    else:
        liminf = y
    
    def integrand(x):
        
        if x > M_t:
           return (beta_frac*x - y)*joint_density(x, y, S_t, tau, b, a, beta)
        
        return (beta_frac*M_t - y)*joint_density(x, y, S_t, tau, b, a, beta)

    val = integrate.quad(lambda x: integrand(x), liminf, limsup)
    return val[0]

def Double_Int_Out_y (S_t, M_t, tau, b, a, beta, beta_frac, p, limsup, number_y_points):

    y_points = linspace(5, limsup*p, number_y_points)
    x_points = []

    for y in y_points:
        x_points.append(Integral_x(S_t, M_t, y, tau, b, a, beta, limsup, beta_frac, p))
    
    result = integrate.trapz( x_points, x=y_points)
    return result

def Double_Int_Out_y1 (S_t, M_t, tau, b, a, beta, beta_frac, p, limsup):

    return integrate.quad(lambda y: Integral_x(S_t, M_t, y, tau, b, a, beta, limsup, beta_frac, p), 1, limsup*p)[0]