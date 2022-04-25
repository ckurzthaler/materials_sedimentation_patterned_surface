################################################################################
# bispherical solution for sphere translating perpendicular to a rigid wall
###############################################################################

# import packages and mathematical functions
from math_functions import *

def coefui(n, s, alpha):
    ai = -np.sinh(alpha)**4*n*(n+1.)*(2.*n+1.)/\
        (np.sqrt(2.)*(4.*np.sinh((n+1./2.)*alpha)**2-(2.*n+1.)**2*np.sinh(alpha)**2))
    ci = -ai
    bi = np.sinh(alpha)**2*n*(n+1.)/(np.sqrt(2.)*(2.*n-1.))*\
        ((2.*np.sinh((2.*n+1.)*alpha)+(2.*n+1.)*np.sinh(2.*alpha))\
        /(4.*np.sinh((n+1./2.)*alpha)**2-(2.*n+1.)**2*np.sinh(alpha)**2)-1.)
    di = -(2.*n-1.)/(2.*n+3.)*bi
    coef = np.array([ai,bi,ci,di])
    coef[np.abs(coef)<10**(-14)]=0
    return(coef[0]*np.cosh((n-1./2.)*s)+coef[1]*np.sinh((n-1./2.)*s)+\
            coef[2]*np.cosh((n+3./2.)*s)+coef[3]*np.sinh((n+3./2.)*s))

def stream_fun(U,a,alpha,s,eta,N):
    tmp = 0.
    for i in range(0,N):
        tmp+= coefui(i,s,alpha)*gegenbauer_p(np.cos(eta),i+1)
    return(U*a**2/(np.cosh(s)-np.cos(eta))**(3./2.)*tmp)

def coefui_s(n,s,alpha):
    ai = -np.sinh(alpha)**4*n*(n+1.)*(2.*n+1.)/\
        (np.sqrt(2.)*(4.*np.sinh((n+1./2.)*alpha)**2-(2.*n+1.)**2*np.sinh(alpha)**2))
    ci = -ai
    bi = np.sinh(alpha)**2*n*(n+1.)/(np.sqrt(2.)*(2.*n-1.))*\
        ((2.*np.sinh((2.*n+1.)*alpha)+(2.*n+1.)*np.sinh(2.*alpha))\
        /(4.*np.sinh((n+1./2.)*alpha)**2-(2.*n+1.)**2*np.sinh(alpha)**2)-1.)
    di = -(2.*n-1.)/(2.*n+3.)*bi
    coef = np.array([ai,bi,ci,di])
    coef[np.abs(coef)<10**(-14)]=0.
    return(coef[0]*(n-1./2.)*np.sinh((n-1./2.)*s)+coef[1]*(n-1./2.)*np.cosh((n-1./2.)*s)+\
            coef[2]*(n+3./2.)*np.sinh((n+3./2.)*s)+coef[3]*(n+3./2.)*np.cosh((n+3./2.)*s))

def norm(N,alpha):
    tmp=0.0
    for n in range(1,N):
        ai = -np.sinh(alpha)**4*n*(n+1.)*(2.*n+1.)/\
            (np.sqrt(2.)*(4.*np.sinh((n+1./2.)*alpha)**2-(2.*n+1.)**2*np.sinh(alpha)**2))
        ci = -ai
        bi = np.sinh(alpha)**2*n*(n+1.)/(np.sqrt(2.)*(2.*n-1.))*\
            ((2.*np.sinh((2.*n+1.)*alpha)+(2.*n+1.)*np.sinh(2.*alpha))\
            /(4.*np.sinh((n+1./2.)*alpha)**2-(2.*n+1.)**2*np.sinh(alpha)**2)-1.)
        di = -(2.*n-1.)/(2.*n+3.)*bi
        coef = np.array([ai,bi,ci,di])
        coef[np.abs(coef)<10**(-14)]=0.
        tmp+= coef[0]+coef[1]+coef[2]+coef[3]
    return(tmp)

def force_z(alpha,N):
    tmp =0
    for i in range(1,N):
        tmp += i*(i+1.)/((2.*i-1.)*(2.*i+3.))*\
        ((2.*np.sinh((2.*i+1)*alpha)+(2.*i+1)*np.sinh(2.*alpha))/\
        (4.*np.sinh((i+1./2.)*alpha)**2-(2.*i+1.)**2*np.sinh(alpha)**2)-1)
    return(tmp*4./3.*np.sinh(alpha))


def stream_fun_s(U,a,alpha,s,eta,N):
    tmp = 0.
    tmp1=0.
    for i in range(0,N):
        tmp+= coefui(i,s,alpha)*gegenbauer_p(np.cos(eta),i+1)
        tmp1+= coefui_s(i,s,alpha)*gegenbauer_p(np.cos(eta),i+1)
    t= -3./2.*np.sinh(s)/np.sqrt(np.cosh(s)-np.cos(eta))**5
    t1 = 1./np.sqrt(np.cosh(s)-np.cos(eta))**3
    return(U*a**2*(t*tmp+t1*tmp1))

def stream_fun_eta(U,a,alpha,s,eta,N):
    tmp = 0.
    tmp1=0.
    for i in range(0,N):
        tmp+= coefui(i,s,alpha)*gegenbauer_p(np.cos(eta),i+1)
        tmp1+= coefui(i,s,alpha)*gegenbauer_p_diff(i+1,1,np.cos(eta))*(-np.sin(eta))
        #print tmp
    t= -3./2.*np.sin(eta)/np.sqrt(np.cosh(s)-np.cos(eta))**5
    t1 = 1./np.sqrt(np.cosh(s)-np.cos(eta))**3
    return(U*a**2*(t*tmp+t1*tmp1))

def vel_r(U,a,alpha,r,s,eta,N,c):
    dsz = (1.-np.cos(eta)*np.cosh(s))/c
    detaz=-np.sin(eta)*np.sinh(s)/c
    return(1./r*(dsz*stream_fun_s(U,a,alpha,s,eta,N)+detaz*stream_fun_eta(U,a,alpha,s,eta,N)))



def vel_z(U,a,alpha,r,s,eta,N,c):
    detar = (-1.+np.cos(eta)*np.cosh(s))/c
    dsr=-np.sin(eta)*np.sinh(s)/c
    return(-1./r*(dsr*stream_fun_s(U,a,alpha,s,eta,N)+detar*stream_fun_eta(U,a,alpha,s,eta,N)))

## calculation of the stress tensor requires second order derivatives of the stream besselfunction

# second order metric terms
def srr(s,eta,c):
    tmp = (np.cos(eta)-np.cos(2.*eta)*np.cosh(s))*np.sinh(s)
    return(tmp/c**2)

def szz(s,eta,c):
    tmp = (-np.cos(eta)+np.cos(2.*eta)*np.cosh(s))*np.sinh(s)
    return(tmp/c**2)

def etarr(s,eta,c):
    tmp = (2.*np.cosh(s)*np.sin(eta)-np.cosh(2.*s)*np.sin(2.*eta))/2.
    return(tmp/c**2)

def etazz(s,eta,c):
    tmp = (-2.*np.cosh(s)*np.sin(eta)+np.cosh(2.*s)*np.sin(2.*eta))/2.
    return(tmp/c**2)

def srz(s,eta,c):
    tmp = (-2.*np.cosh(s)*np.sin(eta)+np.cosh(2.*s)*np.sin(2.*eta))/2.
    return(tmp/c**2)

def etarz(s,eta,c):
    tmp = (np.cos(eta)-np.cos(2.*eta)*np.cosh(s))*np.sinh(s)
    return(tmp/c**2)

# second order stream functions
def coefui_ss(n,s,alpha):
    ai = -np.sinh(alpha)**4*n*(n+1.)*(2.*n+1.)/\
        (np.sqrt(2.)*(4.*np.sinh((n+1./2.)*alpha)**2-(2.*n+1.)**2*np.sinh(alpha)**2))
    ci = -ai
    bi = np.sinh(alpha)**2*n*(n+1.)/(np.sqrt(2.)*(2.*n-1.))*\
        ((2.*np.sinh((2.*n+1.)*alpha)+(2.*n+1.)*np.sinh(2.*alpha))\
        /(4.*np.sinh((n+1./2.)*alpha)**2-(2.*n+1.)**2*np.sinh(alpha)**2)-1.)
    di = -(2.*n-1.)/(2.*n+3.)*bi
    coef = np.array([ai,bi,ci,di])
    coef[np.abs(coef)<10**(-14)]=0.

    return(coef[0]*(n-1./2.)**2*np.cosh((n-1./2.)*s)+coef[1]*(n-1./2.)**2*np.sinh((n-1./2.)*s)+\
            coef[2]*(n+3./2.)**2*np.cosh((n+3./2.)*s)+coef[3]*(n+3./2.)**2*np.sinh((n+3./2.)*s))

def stream_fun_ss(U,a,alpha,s,eta,N):
    tmp = 0.
    tmp1 = 0.
    tmp2 = 0.
    tmp3 = 0.
    for i in range(0,N):
        tmp+= coefui_s(i,s,alpha)*gegenbauer_p(np.cos(eta),i+1) #deriv of tmp2
        tmp1+= coefui_ss(i,s,alpha)*gegenbauer_p(np.cos(eta),i+1) #deriv of tmp3
        tmp2+=coefui(i,s,alpha)*gegenbauer_p(np.cos(eta),i+1) #same as 1st order
        tmp3+=coefui_s(i,s,alpha)*gegenbauer_p(np.cos(eta),i+1) #same as 1st order

    t= -3./2.*np.sinh(s)/np.sqrt(np.cosh(s)-np.cos(eta))**5 #same as 1st order
    t1 = 1./np.sqrt(np.cosh(s)-np.cos(eta))**3 #same as 1st order
    t2 = (6.*(np.cos(eta)-np.cosh(s))*np.cosh(s)+15.*np.sinh(s)**2)/(4.*np.sqrt(np.cosh(s)-np.cos(eta))**7)#deriv of t
    t3 = -3./2.*np.sinh(s)/np.sqrt(np.cosh(s)-np.cos(eta))**5#deriv of t1
    return(U*a**2*(t*tmp+t1*tmp1+t2*tmp2+t3*tmp3))

def stream_fun_ee(U,a,alpha,s,eta,N):
    tmp = 0.
    tmp1 = 0.
    tmp2 = 0.
    tmp3 = 0.
    for i in range(0,N):
        tmp+= coefui(i,s,alpha)*gegenbauer_p_diff(i+1,1,np.cos(eta))*(-np.sin(eta))  #deriv of tmp2
        tmp1+= coefui(i,s,alpha)*(gegenbauer_p_diff(i+1,1,np.cos(eta))*(-np.cos(eta))\
                                    +gegenbauer_p_diff(i+1,2,np.cos(eta))*np.sin(eta)**2) #deriv of tmp3
        tmp2+= coefui(i,s,alpha)*gegenbauer_p(np.cos(eta),i+1) # same as 1st order
        tmp3+= coefui(i,s,alpha)*gegenbauer_p_diff(i+1,1,np.cos(eta))*(-np.sin(eta)) #same as 1st order
        #print tmp
    t = -3./2.*np.sin(eta)/np.sqrt(np.cosh(s)-np.cos(eta))**5 # same as 1st order
    t1 = 1./np.sqrt(np.cosh(s)-np.cos(eta))**3 # same as 1st order
    t2 = 15.*np.sin(eta)*np.sinh(s)/(4.*np.sqrt(np.cosh(s)-np.cos(eta))**7)# deriv of t
    t3 = -3./2.*np.sin(eta)/np.sqrt(np.cosh(s)-np.cos(eta))**5# deriv of t1
    return(U*a**2*(t*tmp+t1*tmp1+t2*tmp2+t3*tmp3))

def stream_fun_es(U,a,alpha,s,eta,N):
    tmp = 0.
    tmp1 = 0.
    tmp2 = 0.
    tmp3 = 0.
    for i in range(0,N):
        tmp+= coefui_s(i,s,alpha)*gegenbauer_p(np.cos(eta),i+1) # deriv of tmp2 wrt s
        tmp1+= coefui_s(i,s,alpha)*gegenbauer_p_diff(i+1,1,np.cos(eta))*(-np.sin(eta)) # deriv of tmp3 wrt s
        tmp2+= coefui(i,s,alpha)*gegenbauer_p(np.cos(eta),i+1) # same as 1st order
        tmp3+= coefui(i,s,alpha)*gegenbauer_p_diff(i+1,1,np.cos(eta))*(-np.sin(eta)) #same as 1st order
    t = -3./2.*np.sin(eta)/np.sqrt(np.cosh(s)-np.cos(eta))**5 # same as 1st order
    t1 = 1./np.sqrt(np.cosh(s)-np.cos(eta))**3 # same as 1st order
    t2 = (6.*(np.cos(eta)-np.cosh(s))*np.cosh(s)+15.*np.sinh(s)**2)/(4.*np.sqrt(np.cosh(s)-np.cos(eta))**7)#deriv of t
    t3 = -3./2.*np.sinh(s)/np.sqrt(np.cosh(s)-np.cos(eta))**5#deriv of t1 wrt s
    return(U*a**2*(t*tmp+t1*tmp1+t2*tmp2+t3*tmp3))

# consistency check: change order of derivatives
def stream_fun_se(U,a,alpha,s,eta,N):
    tmp = 0.
    tmp1 = 0.
    tmp2 = 0.
    tmp3 = 0.
    for i in range(0,N):
        tmp+= coefui(i,s,alpha)*gegenbauer_p_diff(i+1,1,np.cos(eta))*(-np.sin(eta)) # deriv of tmp2 wrt eta
        tmp1+= coefui_s(i,s,alpha)*gegenbauer_p_diff(i+1,1,np.cos(eta))*(-np.sin(eta)) # deriv of tmp3 wrt eta
        tmp2+= coefui(i,s,alpha)*gegenbauer_p(np.cos(eta),i+1) # same as 1st order
        tmp3+= coefui_s(i,s,alpha)*gegenbauer_p(np.cos(eta),i+1) #same as 1st order
        #print tmp
    t = -3./2.*np.sin(eta)/np.sqrt(np.cosh(s)-np.cos(eta))**5 # same as 1st order
    t1 = 1./np.sqrt(np.cosh(s)-np.cos(eta))**3 # same as 1st order
    t2 = 15.*np.sin(eta)*np.sinh(s)/(4.*np.sqrt(np.cosh(s)-np.cos(eta))**7)#deriv of t wrt eta
    t3 = -3./2.*np.sin(eta)/np.sqrt(np.cosh(s)-np.cos(eta))**5#deriv of t1 wrt eta
    return(U*a**2*(t*tmp+t1*tmp1+t2*tmp2+t3*tmp3))


# stress tensor:
# dimensionless R = r/c, Z = z/c,...
def sigmazr(alpha,R,s,eta,N):
    U=1.
    a=1.
    c=1.
    dsz = (1.-np.cos(eta)*np.cosh(s))
    detaz=-np.sin(eta)*np.sinh(s)
    detar = (-1.+np.cos(eta)*np.cosh(s))
    dsr=-np.sin(eta)*np.sinh(s)
    psi_r = (dsr*stream_fun_s(U,a,alpha,s,eta,N)+detar*stream_fun_eta(U,a,alpha,s,eta,N))
    psi_rr=  (stream_fun_ss(U,a,alpha,s,eta,N)*dsr + stream_fun_es(U,a,alpha,s,eta,N)*detar)*dsr+\
             stream_fun_s(U,a,alpha,s,eta,N)*srr(s,eta,c)+\
             stream_fun_eta(U,a,alpha,s,eta,N)*etarr(s,eta,c)+\
             (stream_fun_ee(U,a,alpha,s,eta,N)*detar + stream_fun_es(U,a,alpha,s,eta,N)*dsr)*detar

    psi_zz = (stream_fun_ss(U,a,alpha,s,eta,N)*dsz + stream_fun_es(U,a,alpha,s,eta,N)*detaz)*dsz+\
            stream_fun_s(U,a,alpha,s,eta,N)*szz(s,eta,c)+\
            stream_fun_eta(U,a,alpha,s,eta,N)*etazz(s,eta,c)+\
            (stream_fun_ee(U,a,alpha,s,eta,N)*detaz + stream_fun_es(U,a,alpha,s,eta,N)*dsz)*detaz
    tmp = -1./R*psi_rr+(1./R**2)*psi_r+1./R*psi_zz # for s=0: only last term (vr) contributes
    return(tmp)
