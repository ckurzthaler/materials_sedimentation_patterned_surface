## packages
import numpy as np
import math
import scipy as sp
import matplotlib.pyplot as plt
from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show
#import mpmath
from scipy import optimize
from numpy import random
from scipy import special
from scipy.stats import norm
import sys
from scipy.special import legendre
from fractions import Fraction
from decimal import Decimal
from scipy.special import jv
from scipy.special import j1

# mathematical functions
## bessel functions
def Pn(n,x):
    pn = legendre(n)
    return(pn(x))

def Pndiff(n,k,x):
    pnd = np.polyder(legendre(n),k)
    return(pnd(x))

## gegenbauer polynomials
def gegenbauer_p(x,n):
    if(n==0):
        return(1.)
    elif(n==1):
        return(-x)
    else:
        return((Pn(n-2,x)-Pn(n,x))/(2.*n-1.))

def gegenbauer_p_diff(n,k,x):
    if(n==0):
        return(0.)
    elif(n==1):
        if(k==1):
            return(-1.)
        else:
            return(0.)
    else:
        return((Pndiff(n-2,k,x)-Pndiff(n,k,x))/(2.*n-1.))

# bispherical coordinates
def rr(s,eta,c):
    return(np.sin(eta)/(np.cosh(s)-np.cos(eta)))

def zz(s,eta,c):
    return(np.sinh(s)/(np.cosh(s)-np.cos(eta)))
