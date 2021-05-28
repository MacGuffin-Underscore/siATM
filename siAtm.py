# siAtm.py
"""
Atmospheric Model in SI Units

Rhett A. Smith 19.05.21

Uses scientific units ONLY, input the altitude in (m) and gamma
"""
import math

def siAtm(hg, gam):
    re = 6.356766e6     # radius of earth (m)
    g0 = 9.8066         # gravitational accelaration (m/s^2)
    R = 287             # gas constant
    h = re/(re+hg)*hg   # geo-potential altitude (m)
    if gam == 0:
        gam = 1.4
        
    # Known Altitudes (m)
    h_n = [0, 11000, 25000, 47000, 53000, 79000, 90000, 105000]

    # Known Tempuratures (K)
    t_n = [288.16, 216.66, 216.66, 282.66, 282.66, 165.66, 165.66, 225.66]
    
    # Lapse Rates
    a01 = (t_n[1]-t_n[0])/(h_n[1]-h_n[0])
    a23 = (t_n[3]-t_n[2])/(h_n[3]-h_n[2])
    a45 = (t_n[5]-t_n[4])/(h_n[5]-h_n[4])
    a67 = (t_n[7]-t_n[6])/(h_n[7]-h_n[6])
    
    # Pressure Calculations
    p0 = 101325;
    p1 = p0*(t_n[1]/t_n[0])**(-g0/(a01*R))
    p2 = p1*math.exp(-(g0/(R*t_n[1]))*(h_n[2]-h_n[1]))
    p3 = p2*(t_n[3]/t_n[2])**(-g0/(a23*R))
    p4 = p3*math.exp(-(g0/(R*t_n[3]))*(h_n[4]-h_n[3]))
    p5 = p4*(t_n[5]/t_n[4])**(-g0/(a45*R));
    p6 = p5*math.exp(-(g0/(R*t_n[5]))*(h_n[6]-h_n[5]))
    p7 = p6*(t_n[7]/t_n[6])**(-g0/(a67*R))

    # Conditionals
    #   Finds the 1st conditions that are true
    if h < h_n[1]:
        t = t_n[0]+a01*h
        p = p0*(t/t_n[0])**(-g0/(a01*R))
    elif h < h_n[2]:
        t = t_n[1]
        p = p1*math.exp(-(g0/(R*t))*(h-h_n[1]))
    elif h < h_n[3]:
        t = t_n[2]+a23*(h-h_n[2]);
        p = p2*(t/t_n[2])**(-g0/(a23*R));
    elif h < h_n[4]:
       t = t_n[3]
       p = p3*math.exp(-(g0/(R*t))*(h-h_n[3]))
    elif h < h_n[5]:
        t = t_n[4]+a45*(h-h_n[4])
        p = p4*(t/t_n[4])**(-g0/(a45*R))
    elif h < h_n[6]:
        t = t_n[5]
        p = p5*math.exp(-(g0/(R*t))*(h-h_n[5]));
    elif h < h_n[7]:
        t = t_n[6]+a67*(h-h_n[6]);
        p = p6*(t/t_n[6])**(-g0/(a67*R));
    else:
        t = t_n[7];
        p = p7*math.exp(-(g0/(R*t))*(h-h_n[7]));     
        
    # rho Calcultion
    rho = p/(R*t)
    
    # Speed of Sound (m/s)
    a = math.sqrt(gam*R*t)
    
    return [t, p, rho, a]
