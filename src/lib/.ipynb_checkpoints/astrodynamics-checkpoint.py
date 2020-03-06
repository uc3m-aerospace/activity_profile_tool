#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 16:56:56 2019

@author: davidmorante
"""
    
import numpy  as np
import pandas as pd 
import datetime as datetime

    
#
# Set global variable
#
RE    = +6378.1363;          # km (Earth Radius)
mu    = +398600.440;         # km^3 s^2
J2    = +0.0010826267;       # (Vallado 2007): in units of RE
J3    = -0.0000025327;       # (Vallado 2007): in units of RE
J4    = -0.0000016196;       # (Vallado 2007): in units of RE
RAD   = 57.2957795131;       # Conversion factor from Rad to Degrees
WE    = 7.29211514670698e-5; # Earth's angular velocity rad/s (GMAT reference Manual)
pi    = np.pi;
    
def csv2coe(csv):
    """
    ## This function is based on Algortihm 4.2 from the book 
    ## Orbital Mechanics for Engineering Students by H.D. Curtis, 
    ## Third Edition (2014). 
    ## The algorithm has been generalised to accept multiple objects.
    """  
    ## Create position and velocity vector
    #
    R = np.array([csv.RX,
                  csv.RY,
                  csv.RZ]).T;
        ##
    V = np.array([csv.VX,
                  csv.VY,
                  csv.VZ]).T;
    
    ## Module of R and V
    r  = np.sum(np.abs(R)**2,axis=-1)**(1./2);
    v  = np.sum(np.abs(V)**2,axis=-1)**(1./2);
    
    ## Radial velocity
    vr = np.sum(R*V, axis=1)/r

    ## Angular momentum
    H  = np.cross(R,V);
    h  = np.sum(np.abs(H)**2,axis=-1)**(0.5);
    kk = np.divide(H.T,h);
    kk = kk.T;
    
    ## Equation 4.7: orbit's inclination
    incl = np.arccos( np.divide(H[:,2],h));
    
    ## Equation 4.10: Eccentricity
    E = 1.0/mu*((v**2 - mu/r)*R.T - r*vr*V.T);
    E = E.T;
    e = np.sum(np.abs(E)**2,axis=-1)**(1./2);
    
    ## Equation 4.62 (a < 0 for a hyperbola): Semimajor axis
    a = (h**2/mu)/(1 - e**2);
    
    ## Equation 4.9: Ascending Node
    ## RA = 0 if n = 0
    
    z  = np.array([[0., 0., 1.]])
    zq = np.ones((1,np.size(R,0)));
    KN = np.matmul(z.T,zq);
    N  = np.cross(KN.T,H);
    nn = np.sum(np.abs(N)**2,axis=-1)**(1./2);
    #
    cond1 = (nn != 0);
    RA = np.arccos(N[:,0]/nn);
    #
    #if any(cond1 == 0):
    #    RA = 0;
    #    nn  = h;
    
    ## Modify only those objects with N(2) < 0
    cond2     = (N[:,1] < 0);
    RA[cond2] = 2*np.pi - RA[cond2];    
        
    ## Equation 4.12: Argument of perigee
    ## w = 0 if n = 0 or e < eps
    eps   = 1.e-12;
    cond3 = (nn != 0)*(e > eps);
    w     = np.arccos(np.sum(N*E, axis=1)/nn/e)*cond3;
        
    ## Modify only those objects with
    cond4    = ((nn != 0)*(e > eps)*(E[:,2] < 0));
    w[cond4] = 2*np.pi - w[cond4];
    w[w == 2*np.pi] = 0;
        
    ## Equation 4.13a (incorporating the case e = 0): True anomaly
    cond5     = (e > eps);
    TA        = np.arccos(np.sum(E*R, axis=1)/e/r)*cond5;
    cond6     = (e > eps)*(vr < 0);
    TA[cond6] = 2*np.pi - TA[cond6];
    
    ## Modify only those objects with e < eps
    cp    = np.cross(N,R);
    cp    = cp.T;
    cond7 = cp[2,:] < 0;
    
    if any(~cond5):
        TA[~cond5]         = np.arccos(np.sum(N[~cond5,:]*R[~cond5,:], axis=1)/nn[~cond5]/r[~cond5]);
        TA[(~cond5)*cond7] = 2*np.pi - TA[(~cond5)*cond7];
        
    ## Convert to a data Frame
    coe = {
    "a" : a,
    "e" : e,
    "incl" : incl,
    "w" : w,
    "RA": RA,
    "TA": TA,
    "h" : h,
    }
    coe = pd.DataFrame(coe) 
    
    ## Output Orbital elements
    return coe
    

def coe2csv(coe):
    """
    ## Rotation matrices
    """
    sz     = np.size(coe.a);
    Rnode  = np.zeros((3,3,sz));
    Ri     = np.zeros((3,3,sz));
    Romega = np.zeros((3,3,sz));
    
    ## Rotation matrix for inclination %
    Ri[0,0,:]     = np.ones(sz);
    Ri[1,1,:]     = np.cos(coe.incl);
    Ri[2,2,:]     = Ri[1,1,:];
    Ri[1,2,:]     = np.sin(coe.incl);
    Ri[2,1,:]     = -Ri[1,2,:];
    
    ## Rotation matrix for ascending node %
    Rnode[2,2,:]  = np.ones(sz);
    Rnode[0,0,:]  = np.cos(coe.RA);
    Rnode[1,1,:]  = Rnode[0,0,:];
    Rnode[0,1,:]  = np.sin(coe.RA);
    Rnode[1,0,:]  = -Rnode[0,1,:];
    
    ## Rotation matrix for argument perigee %
    Romega[2,2,:] = np.ones(sz);
    Romega[0,0,:] = np.cos(coe.w);
    Romega[1,1,:] = Romega[0,0,:];
    Romega[0,1,:] = np.sin(coe.w);
    Romega[1,0,:] = -Romega[0,1,:];
    
    ## Determine semiparameter
    p = coe.a*(1-coe.e**2);
    h = np.sqrt(mu*p)
    ## Define some constants for computational optimization
    sqrMuP = mu / h;
    sinNu  = np.sin(coe.TA);
    cosNu  = np.cos(coe.TA);
    
    ## Perifocal Frame
    rPQW = np.array([p * cosNu / ( 1 + coe.e * cosNu ),
            p * sinNu / ( 1 + coe.e * cosNu ),
            np.zeros(sz)])
            
    vPQW = np.array([-sqrMuP * sinNu,      
            sqrMuP * ( coe.e + cosNu ),
            np.zeros(sz)])
    
    ## Rotate to Earth Frame
    R = np.zeros((sz,3));
    V = np.zeros((sz,3));
    
    for i in range(0,sz):
      rot = np.matmul(np.matmul(Romega[:,:,i],Ri[:,:,i]),Rnode[:,:,i]);
      R[i,:] = np.matmul(rot.T,rPQW[:,i]);
      V[i,:] = np.matmul(rot.T,vPQW[:,i]);

    ## Convert to a data Frame
    csv = {
    "RX" : R[:,0],
    "RY" : R[:,1],
    "RZ" : R[:,2],
    "VX" : V[:,0],
    "VY" : V[:,1],
    "VZ" : V[:,2],
    }
    csv = pd.DataFrame(csv) 
    
    ## Output Cartesian State Vector
    return csv
    
def rk4(dfdt,x0,h):
    #
    # Runge Kutta (4) Integrator
    #
    k1 = dfdt(x0);
    k2 = dfdt(x0 + 0.5*h*k1);
    k3 = dfdt(x0 + 0.5*h*k2);
    k4 = dfdt(x0 + h*k3);
    
    x = x0 + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4)*h;
    
    return x
    
    
def TBP_csv(x):
    """
    # Propagate Orbit with Cowell method including Earth Oblateness effects J2,J3,J4
    """
    ## Create position and velocity vector
    R = np.array([x[0],
                  x[1],
                  x[2]]);
    ##
    V = np.array([x[3],
                  x[4],
                  x[5]]);
    ## Compute radia distance
    r = np.sqrt(R[0]**2 + R[1]**2 + R[2]**2);
    
    ## Potential acceleration
    a0 = -mu*R/r**3.0;
    
    ## Acceleration due to the Earth oblateness
    aJ2cte = (3.0/2.0)*mu*J2*RE**2 / r**4;
    z2r2 = R[2]**2 / r**2;
    aJ2x = R[0]/r*(5*z2r2 - 1);
    aJ2y = R[1]/r*(5*z2r2 - 1);
    aJ2z = R[2]/r*(5*z2r2 - 3);
    aJ2 = aJ2cte*np.array([aJ2x, aJ2y, aJ2z]);
    
    aJ3cte = 0.5*mu*J3*RE**3 / r**5;
    zr   = R[2] / r;
    aJ3x = 5*(R[0]/r) * ( 7*zr**3 - 3*zr);
    aJ3y = 5*(R[1]/r) * ( 7*zr**3 - 3*zr);
    aJ3z = 3*((35.0/3.0)*zr**4 - 10*zr**2 + 1);
    aJ3  = aJ3cte*np.array([aJ3x, aJ3y, aJ3z]);
    
    aJ4cte = (15.0/8.0)*mu*J4*RE**4 / r**7;
    aJ4x = R[0]*(1 - 14*zr**2 + 21*zr**4);
    aJ4y = R[1]*(1 - 14*zr**2 + 21*zr**4);
    aJ4z = R[2]*(5 - (70.0/3.0)*zr**2 + 21*zr**4);
    aJ4 = aJ4cte*np.array([aJ4x, aJ4y, aJ4z]);
    
    ## Join acceleration
    a = a0 + aJ2 + aJ3 + aJ4;
    
    ## Compute derivatives
    dfdt = np.zeros(6);
    dfdt[0:3] = V;
    dfdt[3:6] = a;
    
    ## Ouput derivatives
    return dfdt


def orbit(csv0,t0,tf,h):
    """
    # Integrate an orbit between the initial and final time t0 and tf
    """
    n    = np.floor((tf-t0)/h);
    n    = n.astype(int)
    x    = np.zeros((6,n+1));
    time = np.zeros(n+1);
    #
    # Convert to position and velocity vector
    #
    x0    = np.zeros(6);
    x0[0] = csv0.RX;
    x0[1] = csv0.RY;
    x0[2] = csv0.RZ;
    x0[3] = csv0.VX;
    x0[4] = csv0.VY;
    x0[5] = csv0.VZ;
    #
    # First Step
    #
    x[:,0]  = x0;
    time[0] = t0;
    #
    for i in range(1,n):  
        
        x[:,i]  = rk4(TBP_csv,x[:,i-1],h);
        time[i] = time[i-1] + h;
    #
    # Adjust final step to stop at tf
    #
    hr      = tf-time[i];
    x[:,i+1]  = rk4(TBP_csv,x[:,i],hr)
    time[i+1] = time[i] + hr;
    #
    # Convert to Panda's Data Frame
    #
    csv = {
    "RX" : x[0,:],
    "RY" : x[1,:],
    "RZ" : x[2,:],
    "VX" : x[3,:],
    "VY" : x[4,:],
    "VZ" : x[5,:],
    }
    csv = pd.DataFrame(csv) 
    ## Ouput derivatives
    return time, csv


def getSolarPosition (JD):
    """
    # Compute Solar Vector in ECI given the Julian Date (JD)
    """
    n = JD - 2451545.0;
    n = JD - 2451545.0;
    L = (280.460 + 0.9856474 * n)/RAD;
    g = (357.528 + 0.9856003 * n)/RAD;
     #
    epsilon= (23.439 - 4e-7 * n) / RAD;
    lamb= L + (1.915 * np.sin(g) + 0.020 * np.sin(2*g)) / RAD;

    RS = np.array([np.cos(lamb), np.cos(epsilon)*np.sin(lamb), np.sin(epsilon)*np.sin(lamb)]).T;

    return RS
    
def getElipseCondition(time, csv):  
    """
    # Function to get the Earth Eclipse condition Given a trajectory
    # Time input is given in Julian Days
    # return +1 when the spacecraft is in eclipse
    # return -1 when the spaceraft is in sunlight
    """   
    # Create position and velocity vector
    #
    R = np.array([csv.RX,
                  csv.RY,
                  csv.RZ]).T;
    #
    shadow = np.zeros(np.size(R,0)) 
    #
    SunPos = getSolarPosition (time)
    #
    prj = np.sum(SunPos*R, axis=1);
    # 
    shadow = -1.0 * (prj > 0);
    #
    # Modulus of R
    #
    r  = np.sum(np.abs(R)**2,axis=-1)**(1./2);
    #
    delta  = prj + np.sqrt(r**2 - RE**2);
    #
    # Outside cylinder condition
    #
    shadow[delta < 0] = +1; 
    shadow[delta > 0] = +0;
    #
    return shadow
    
def getJulianDate(date):
    """
    Convert a datetime object into julian float.
    Args:
        date: datetime-object of date in question

    Returns: float - Julian calculated datetime.
    Raises: 
        TypeError : Incorrect parameter type
        ValueError: Date out of range of equation
    """

    # Ensure correct format
    if not isinstance(date, datetime.datetime):
        raise TypeError('Invalid type for parameter "date" - expecting datetime')
    elif date.year < 1801 or date.year > 2099:
        raise ValueError('Datetime must be between year 1801 and 2099')

    # Perform the calculation
    julian_datetime = 367 * date.year - int((7 * (date.year + int((date.month + 9) / 
        12.0))) / 4.0) + int(
        (275 * date.month) / 9.0) + date.day + 1721013.5 + (
        date.hour + date.minute / 60.0 + date.second /(60**2))/ 24.0 - 0.5 * np.copysign(
        1, 100 * date.year + date.month - 190002.5) + 0.5

    return julian_datetime    


def getVisivility(lat, lon, time, csv, el_min):
    """
    # Function to compute the visibility from a given point on Earth
    # return +1 when the ground station is visible
    # return -1 when the ground station is NOT visible
    """  
    # Create position and velocity vector
    #
    R = np.array([csv.RX,
                  csv.RY,
                  csv.RZ]).T;
    #
    r = np.sum(np.abs(R)**2,axis=-1)**(1./2);
    #
    #
    # Convert lat, lon to cartesian coodinates in ECEF Referenfe Frame
    #
    R_ECEF = np.array([RE*np.cos(lon)*np.cos(lat), RE*np.sin(lon)*np.cos(lat), RE*np.sin(lat)]);
    #
    # Convert cartesian ECEF to cartesian ECI
    #
    R_ECI  = ecef2eci(R_ECEF, time)
    #
    # Distance form station to spacecraft
    #
    drel  = R_ECI.T-R;
    d     = np.sum(np.abs(drel)**2,axis=-1)**(1./2);
    el    = np.arccos((r**2 - RE**2 - d**2)/(-2*RE*d)) - np.pi/2
    drel  = R_ECI.T-R

#
    #
    visibility = el > el_min;

    return visibility


def ecef2eci(R_ECEF,time):
    """
    # Function to compute rotation  matrix from ECEF to ECI (simple model)
    # Formulas taken from the US Naval Observatory
    """      
    #
    # T is the Julian Date in julian centuries
    #
    d    = time - 2451545.0;
    T    = d/ 36525;
    #
    # Compute Greenwich Mean sidereal Time (in hours)
    #
    GMST  = 2*np.pi*(0.7790572732640 + 1.00273781191125448*d)
    # 
    # Compute Rotation Matrix
    #
    R_ECI = np.zeros((3,np.size(T)))
    #
    for i in range(0,np.size(T)):
             RT = np.array([[+np.cos(GMST[i]), -np.sin(GMST[i]), 0], 
                  [ +np.sin(GMST[i]),+np.cos(GMST[i]), 0],  
                  [0          , 0           , 1  ]]);
    #
             R_ECI[:,i] = np.matmul(RT,R_ECEF.T)
    #
        
    return R_ECI








