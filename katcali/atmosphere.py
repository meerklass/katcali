import numpy as np


def opacity(T, RH, P, h, f):
    """
    Calculates zenith opacity according to ITU-R P.676-9. This is for 
    elevations > 10 deg. Use as "Tsky * (1 - exp(-opacity/sin(el)))" for 
    elevation dependence. (N.B. Same as used for KAT-7.)
    
    Parameters
    ----------
    temp : float
        Temperature in deg C.
        
    rel_hum : float
        Relative humidity, 0 < RH < 1.
    
    pressure : float
        Dry air pressure in hPa (equiv. mbar).
        
    height : float
        Height above sea level in km.
    
    freq : float
        Frequency in GHz (must be < 55 GHz).
    
    Returns
    -------
    opac : float
        Approximate atmospheric opacity at zenith [Nepers].
    """
    # [hPa] from A. L. Buck research manual 1996
    es = 6.1121*np.exp((18.678-T/234.5) * T / (257.14+T))
    
    # [g/m^3] from A. L. Buck research manual 1996
    # (ITU-R omitted the factor "RH" - a mistake)
    rho = RH*es*216.7/(T+273.15)

    # The following is taken directly from ITU-R P.676-9
    p_tot = P + es # from eq 3

    rho = rho*np.exp(h/2) # Adjust to sea level as per eq 32

    # eq 22
    r_t = 288./(273.+T)
    r_p = p_tot/1013.
    phi = lambda a, b, c, d: r_p**a*r_t**b*np.exp(c*(1-r_p)+d*(1-r_t))
    E_1 = phi(0.0717,-1.8132,0.0156,-1.6515)
    E_2 = phi(0.5146,-4.6368,-0.1921,-5.7416)
    E_3 = phi(0.3414,-6.5851,0.2130,-8.5854)
    
    # Following is valid only for f <= 54 GHz
    yo = ( 7.2*r_t**2.8 / (f**2+0.34*r_p**2*r_t**1.6) + 0.62*E_3 / ((54-f)**(1.16*E_1)+0.83*E_2) ) * f**2 * r_p**2 *1e-3
    
    # eq 23
    n_1 = 0.955*r_p*r_t**0.68 + 0.006*rho
    n_2 = 0.735*r_p*r_t**0.5 + 0.0353*r_t**4*rho
    g = lambda f, f_i: 1+(f-f_i)**2/(f+f_i)**2
    yw = (  3.98*n_1*np.exp(2.23*(1-r_t))/((f-22.235)**2+9.42*n_1**2)*g(f,22) + 11.96*n_1*np.exp(0.7*(1-r_t))/((f-183.31)**2+11.14*n_1**2)
          + 0.081*n_1*np.exp(6.44*(1-r_t))/((f-321.226)**2+6.29*n_1**2) + 3.66*n_1*np.exp(1.6*(1-r_t))/((f-325.153)**2+9.22*n_1**2)
          + 25.37*n_1*np.exp(1.09*(1-r_t))/(f-380)**2 + 17.4*n_1*np.exp(1.46*(1-r_t))/(f-448)**2
          + 844.6*n_1*np.exp(0.17*(1-r_t))/(f-557)**2*g(f,557) + 290*n_1*np.exp(0.41*(1-r_t))/(f-752)**2*g(f,752)
          + 8.3328e4*n_2*np.exp(0.99*(1-r_t))/(f-1780)**2*g(f,1780)
          ) * f**2*r_t**2.5*rho*1e-4

    # eq 25
    t_1 = 4.64/(1+0.066*r_p**-2.3) * np.exp(-((f-59.7)/(2.87+12.4*np.exp(-7.9*r_p)))**2)
    t_2 = 0.14*np.exp(2.12*r_p) / ((f-118.75)**2+0.031*np.exp(2.2*r_p))
    t_3 = 0.0114/(1+0.14*r_p**-2.6) * f * (-0.0247+0.0001*f+1.61e-6*f**2) / (1-0.0169*f+4.1e-5*f**2+3.2e-7*f**3)
    ho = 6.1/(1+0.17*r_p**-1.1)*(1+t_1+t_2+t_3)

    # eq 26
    sigma_w = 1.013/(1+np.exp(-8.6*(r_p-0.57)))
    hw = 1.66*( 1 + 1.39*sigma_w/((f-22.235)**2+2.56*sigma_w) + 3.37*sigma_w/((f-183.31)**2+4.69*sigma_w) + 1.58*sigma_w/((f-325.1)**2+2.89*sigma_w) )

    # Attenuation from dry & wet atmosphere relative to a point outside of the 
    # atmosphere
    A = yo*ho*np.exp(-h/ho) + yw*hw*np.exp(-h/hw) # [dB] from equations 27, 30 & 31

    return A*np.log(10)/10.0 # Convert dB to Nepers
    
