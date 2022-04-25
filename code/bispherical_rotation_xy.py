###############################################################################
# bispherical solution for sphere rotating parallel to a rigid, no-slip wall
###############################################################################

# import packages and mathematical functions
from math_functions import *

def kn(n,alpha):
    return((n+1./2.)*1./sp.tanh((n+1./2.)*alpha)-1./sp.tanh(alpha))

def coefR(alpha,N):
    ki=np.zeros(N+1)
    bi=np.zeros(N)
    ai=np.zeros((N,N))
    ki[N]=kn(N,alpha)
    for n in range(1,N+1):
        ki[n-1]=kn(n-1,alpha)
        bi[n-1]=-np.sqrt(2.)*np.exp(-(n+0.5)*alpha)/((2.*n+1.)*np.sinh(alpha))*\
                ((2.*n+1.)**2*(np.exp(alpha)/(2.*n-1.)+np.exp(-alpha)/(2.*n+3.))/(np.sinh((n+0.5)*alpha))-\
                (2.*n-1.)/(np.sinh((n-0.5)*alpha))-(2.*n+3.)/np.sinh((n+1.5)*alpha))
    i=N
    ai[N-1,N-1]=-i/(2.*i+1.)*((2.*i-1.)*ki[(i-1)]-(2.*i-3.)*ki[i])-\
                (i+1.)/(2.*i+1.)*((2.*i+5.)*ki[i])
    ai[N-1,N-2]=(i-1)/(2.*i-1.)*((2.*i-1.)*ki[(i-1)]-(2.*i-3.)*ki[i])
    for i in range(1,N):
        ai[i-1,i-1]=-i/(2.*i+1.)*((2.*i-1.)*ki[(i-1)]-(2.*i-3.)*ki[i])-\
                (i+1.)/(2*i+1.)*((2*i+5.)*ki[i]-(2*i+3.)*ki[(i+1)])
        ai[i-1,i-2]=(i-1.)/(2.*i-1.)*((2.*i-1.)*ki[(i-1)]-(2.*i-3.)*ki[i])
        ai[i-1,i]=(i+2.)/(2.*i+3.)*((2.*i+5.)*ki[i]-(2.*i+3.)*ki[(i+1)])
    return(ai,bi)

def coefbiR(ai):
    N=len(ai)
    bi=np.zeros(N)
    bi[0]=-3.*ai[0]+3.*ai[1]
    for i in range(2,N-1):
        bi[i-1]=(i-1.)*ai[i-2]-(2.*i+1.)*ai[i-1]+(i+2.)*ai[i]
    return(bi)

def coefdiR(ai):
    N=len(ai)
    di=np.zeros(N)
    di[0]=ai[0]
    di[1]=3.*ai[1]
    for i in range(2,N-1):
        di[i]=-1./2.*(i-1.)*i*ai[i-2]+1./2.*(i+1.)*(i+2.)*ai[i]
    return(di)

def coefciR(ai,alpha):
    N=len(ai)
    ci=np.zeros(N)
    ki=np.zeros(N)
    lai = np.zeros(N)
    for n in range(0,N):
        ki[n]=kn(n,alpha)
        lai[n]=-1./np.sqrt(2.)*(np.exp(-(n-0.5)*alpha)/(2.*n-1.)-np.exp(-(n+1.5)*alpha)/(2.*n+3.))
    ci[0]=-2.*ki[1]*(-ai[0]+3.*ai[1]/5.)+4.*lai[1]/(np.sinh(alpha)*np.sinh((0.5+1.)*alpha))
    for i in range(2,N-1):
        ci[i-1]=-2.*ki[i]*((i-1.)*ai[i-2]/(2.*i-1.)-ai[i-1]+(i+2.)*ai[i]/(2.*i+3.))+\
                4.*lai[i]/(np.sinh(alpha)*np.sinh((i+0.5)*alpha))
    return(ci)

def coefeiR(ai,alpha):
    N=len(ai)
    ei=np.zeros(N)
    ki=np.zeros(N)
    lai = np.zeros(N)
    for n in range(0,N):
        ki[n]=kn(n,alpha)
        lai[n]=-1./np.sqrt(2.)*(np.exp(-(n-0.5)*alpha)/(2.*n-1.)-np.exp(-(n+1.5)*alpha)/(2.*n+3.))
    ei[0]=(np.sqrt(2.)*np.exp(-alpha/2.)-lai[0]/np.sinh(alpha))/(sp.sinh(alpha/2.))-\
            ki[0]*2.*ai[0]/3.
    ei[1]=(3.*np.sqrt(2.)*np.exp(-alpha*(3./2.))-lai[1]/np.sinh(alpha))/(sp.sinh(alpha*(3./2.)))-\
            ki[1]*ai[1]*2.*3./5.
    for i in range(2,N-1):
        ei[i]=((2.*i+1.)*np.sqrt(2.)*np.exp(-alpha*(i+1./2.))-lai[i]/np.sinh(alpha))/(np.sinh(alpha*(i+1./2.)))+\
                ki[i]*((i-1.)*i/(2.*i-1.)*ai[i-2]-(i+1.)*(i+2.)*ai[i]/(2.*i+3))
    return(ei)

def coeffiR(ai):
    N=len(ai)
    fi=np.zeros(N)
    for i in range(2,N-1):
        fi[i-2]=1./2.*(ai[i-2]-ai[i])
    return(fi)

def coefgiR(ai,alpha):
    N=len(ai)
    gi=np.zeros(N)
    ki=np.zeros(N)
    lai = np.zeros(N)
    for n in range(0,N):
        ki[n]=kn(n,alpha)
        lai[n]=-1./np.sqrt(2.)*(np.exp(-(n-0.5)*alpha)/(2.*n-1)-np.exp(-(n+1.5)*alpha)/(2.*n+3.))
    for i in range(2,N-1):
        gi[i-2]=-ki[i]*(ai[i-2]/(2.*i-1.)-ai[i]/(2.*i+3.))-4.*lai[i]/np.sinh(alpha)/np.sinh(alpha*(i+0.5))
    return(gi)

## hydrodynamic force and torque on sphere
def force_xR(alpha,en,cn):
    tmp = en[0]
    for i in range(1, len(en)):
        tmp += en[i]+i*(i+1.)*cn[i-1]
    return(-tmp*np.sinh(alpha)**2.*np.sqrt(2.)/6.)

def torque_yR(alpha,an,bn,cn,dn,en):
    tmp = -(en[0])*(1.-1./np.tanh(alpha))*(2.+np.exp(-alpha))-(dn[0])*(1.-1./np.tanh(alpha))*(2.-np.exp(-alpha))
    for i in range(1, len(en)):
        tmp += (2.+np.exp(-(2.*i+1)*alpha))*(i*(i+1.)*(2.*an[i-1]+cn[i-1]/np.tanh(alpha))-(2.*i+1.-1./np.tanh(alpha))*en[i])+\
                (2.-np.exp(-(2.*i+1)*alpha))*(i*(i+1.)*bn[i-1]/np.tanh(alpha)-(2.*i+1.-1./np.tanh(alpha))*dn[i])
    return(-(tmp*np.sinh(alpha)**3/(4.*np.sqrt(2.))-1.)/3.)


## calculation of velocity field
def w1(s,eta,an):
    tmp=0.0
    N=len(an)
    for i in range(1,N):
        tmp= tmp+ an[i-1]*sp.sinh((i+1./2.)*s)*Pndiff(i,1,np.cos(eta)) #*besselfunction
    return(np.sqrt(sp.cosh(s)-np.cos(eta))*np.sin(eta)*tmp)

def Q1(s,eta,bn,cn):
    tmp=0.0
    N=len(bn)
    for i in range(1,N+1):
        tmp= tmp+ (bn[i-1]*sp.cosh((i+1./2.)*s)+cn[i-1]*sp.sinh((i+1./2.)*s))*Pndiff(i,1,np.cos(eta))
    return(np.sqrt(sp.cosh(s)-np.cos(eta))*np.sin(eta)*tmp)

def U0(s,eta,dn,en):
    tmp=0.0
    N=len(dn)
    for i in range(0,N):
        tmp =tmp+ (dn[i]*sp.cosh((i+1./2.)*s)+en[i]*sp.sinh((i+1./2.)*s))*Pn(i,np.cos(eta))
    return(np.sqrt(sp.cosh(s)-np.cos(eta))*tmp)

def U2(s,eta,fn,gn):
    tmp=0.0
    N=len(fn)
    for i in range(2,N+2):
        tmp= tmp + (fn[i-2]*sp.cosh((i+1./2.)*s)+gn[i-2]*sp.sinh((i+1./2.)*s))*Pndiff(i,2,np.cos(eta))
    return(np.sqrt(sp.cosh(s)-np.cos(eta))*(np.sin(eta))**2*tmp)

## calculation of the gradient in z-direction

def Q1_sz(s,eta,bn,cn,c):
    tmp=0.0
    tmp1=0.0
    N=len(bn)
    for i in range(1,N+1):
        tmp= tmp+ (bn[i-1]*np.cosh((i+1./2.)*s)+cn[i-1]*np.sinh((i+1./2.)*s))*Pndiff(i,1,np.cos(eta))
        tmp1 = tmp1+(bn[i-1]*(i+1./2.)*sp.sinh((i+1./2.)*s)+cn[i-1]*(i+1./2.)*sp.cosh((i+1./2.)*s))*Pndiff(i,1,np.cos(eta))
    t1 = np.sin(eta)*(1.-np.cos(eta)*np.cosh(s))*np.sqrt(np.cosh(s)-np.cos(eta))/float(c)
    t = np.sin(eta)*np.sinh(s)*(1.-np.cos(eta)*np.cosh(s))/(2.*np.sqrt(sp.cosh(s)-np.cos(eta))*c)
    return(t1*tmp1+t*tmp)

def Q1_etaz(s,eta,bn,cn,c):
    tmp=0.0
    tmp1=0.0
    N=len(bn)
    for i in range(1,N+1):
        tmp= tmp+(bn[i-1]*sp.cosh((i+1./2.)*s)+cn[i-1]*sp.sinh((i+1./2.)*s))*Pndiff(i,1,np.cos(eta))
        tmp1=tmp1+(bn[i-1]*sp.cosh((i+1./2.)*s)+cn[i-1]*sp.sinh((i+1./2.)*s))*Pndiff(i,2,np.cos(eta))*(-np.sin(eta))

    t =(2.*np.cos(eta)**2-2.*np.cos(eta)*np.cosh(s)-np.sin(eta)**2)*np.sin(eta)*np.sinh(s)/(2.*c*np.sqrt(np.cosh(s)-np.cos(eta)))
    t1 = -np.sin(eta)**2*np.sinh(s)*np.sqrt(sp.cosh(s)-np.cos(eta))/c
    return(t1*tmp1 +t*tmp)

def U0_sz(s,eta,dn,en,c):
    tmp=0.0
    tmp1=0.0
    N=len(dn)
    for i in range(0,N):
        tmp =tmp+ (dn[i]*sp.cosh((i+1./2.)*s)+en[i]*sp.sinh((i+1./2.)*s))*Pn(i,np.cos(eta))
        tmp1 = tmp1 +(dn[i]*(i+1./2.)*sp.sinh((i+1./2.)*s)+en[i]*(i+1./2.)*sp.cosh((i+1./2.)*s))*Pn(i,np.cos(eta))
    t = np.sinh(s)*(1.-np.cos(eta)*np.cosh(s))/(2.*c*np.sqrt(sp.cosh(s)-np.cos(eta)))
    t1 = np.sqrt(sp.cosh(s)-np.cos(eta))*(1.-np.cos(eta)*np.cosh(s))/c
    return(t*tmp+t1*tmp1)

def U0_etaz(s,eta,dn,en,c):
    tmp=0.0
    tmp1=0.0
    N=len(dn)
    for i in range(0,N):
        tmp =tmp+ (dn[i]*sp.cosh((i+1./2.)*s)+en[i]*sp.sinh((i+1./2.)*s))*Pn(i,np.cos(eta))
        tmp1 = tmp1 +(dn[i]*sp.cosh((i+1./2.)*s)+en[i]*sp.sinh((i+1./2.)*s))*Pndiff(i,1,np.cos(eta))*(-np.sin(eta))
    t = -np.sinh(s)*np.sin(eta)**2/(2.*c*np.sqrt(sp.cosh(s)-np.cos(eta)))
    t1 = -np.sqrt(sp.cosh(s)-np.cos(eta))*np.sinh(s)*np.sin(eta)/c
    return(t*tmp+t1*tmp1)

def U2_sz(s,eta,fn,gn,c):
    tmp=0.0
    tmp1=0.0
    N=len(fn)
    for i in range(2,N+2):
        tmp1= tmp1 + (fn[i-2]*(i+1./2.)*np.sinh((i+1./2.)*s)+gn[i-2]*(i+1./2.)*np.cosh((i+1./2.)*s))*Pndiff(i,2,np.cos(eta))
        tmp= tmp + (fn[i-2]*sp.cosh((i+1./2.)*s)+gn[i-2]*sp.sinh((i+1./2.)*s))*Pndiff(i,2,np.cos(eta))
    t=np.sin(eta)**2*np.sinh(s)*(1-np.cos(eta)*np.cosh(s))/(2.*c*np.sqrt(sp.cosh(s)-np.cos(eta)))
    t1=np.sqrt(np.cosh(s)-np.cos(eta))*(np.sin(eta))**2*(1.-np.cos(eta)*np.cosh(s))/c
    return(t1*tmp1+t*tmp)

def U2_etaz(s,eta,fn,gn,c):
    tmp=0.0
    tmp1=0.0
    N=len(fn)
    for i in range(2,N+2):
        tmp= tmp + (fn[i-2]*sp.cosh((i+1./2.)*s)+gn[i-2]*sp.sinh((i+1./2.)*s))*Pndiff(i,2,np.cos(eta))
        tmp1=tmp1+(fn[i-2]*sp.cosh((i+1./2.)*s)+gn[i-2]*sp.sinh((i+1./2.)*s))*Pndiff(i,3,np.cos(eta))*(-np.sin(eta))
    t = -np.sin(eta)*np.sinh(s)*(4.*np.cos(eta)*(np.cosh(s)-np.cos(eta))*np.sin(eta)+np.sin(eta)**3)/(2.*c*np.sqrt(sp.cosh(s)-np.cos(eta)))
    t1 = -np.sqrt(sp.cosh(s)-np.cos(eta))*(np.sin(eta))**2*np.sin(eta)*np.sinh(s)/c
    return(t*tmp+t1*tmp1)

def U0_eta(s,eta,dn,en):
    tmp=0.0
    tmp1=0.0
    N=len(dn)
    for i in range(0,N):
        tmp =tmp+ (dn[i]*sp.cosh((i+1./2.)*s)+en[i]*sp.sinh((i+1./2.)*s))*Pn(i,np.cos(eta))
        tmp1 = tmp1 +(dn[i]*sp.cosh((i+1./2.)*s)+en[i]*sp.sinh((i+1./2.)*s))*Pndiff(i,1,np.cos(eta))*(-np.sin(eta))
    t = np.sin(eta)/(2.*np.sqrt(sp.cosh(s)-np.cos(eta)))
    t1 = np.sqrt(sp.cosh(s)-np.cos(eta))
    return(t*tmp+t1*tmp1)

def U0_s(s,eta,dn,en):
    tmp=0.0
    tmp1=0.0
    N=len(dn)
    for i in range(0,N):
        tmp =tmp+ (dn[i]*sp.cosh((i+1./2.)*s)+en[i]*sp.sinh((i+1./2.)*s))*Pn(i,np.cos(eta))
        tmp1 = tmp1 +(dn[i]*(i+1./2.)*sp.sinh((i+1./2.)*s)+en[i]*(i+1./2.)*sp.cosh((i+1./2.)*s))*Pn(i,np.cos(eta))
    t = np.sinh(s)/(2.*np.sqrt(sp.cosh(s)-np.cos(eta)))
    t1 = np.sqrt(sp.cosh(s)-np.cos(eta))
    return(t*tmp+t1*tmp1)

def U2_s(s,eta,fn,gn):
    tmp=0.0
    tmp1=0.0
    N=len(fn)
    for i in range(2,N+2):
        tmp1= tmp1 + (fn[i-2]*(i+1./2.)*np.sinh((i+1./2.)*s)+gn[i-2]*(i+1./2.)*np.cosh((i+1./2.)*s))*Pndiff(i,2,np.cos(eta))
        tmp= tmp + (fn[i-2]*sp.cosh((i+1./2.)*s)+gn[i-2]*sp.sinh((i+1./2.)*s))*Pndiff(i,2,np.cos(eta))
    t=np.sinh(s)*np.sin(eta)**2/(2.*np.sqrt(sp.cosh(s)-np.cos(eta)))
    t1=np.sqrt(np.cosh(s)-np.cos(eta))*(np.sin(eta))**2
    return(t1*tmp1+t*tmp)

def U2_eta(s,eta,fn,gn):
    tmp=0.0
    tmp1=0.0
    N=len(fn)
    for i in range(2,N+2):
        tmp= tmp + (fn[i-2]*sp.cosh((i+1./2.)*s)+gn[i-2]*sp.sinh((i+1./2.)*s))*Pndiff(i,2,np.cos(eta))
        tmp1=tmp1+(fn[i-2]*sp.cosh((i+1./2.)*s)+gn[i-2]*sp.sinh((i+1./2.)*s))*Pndiff(i,3,np.cos(eta))*(-np.sin(eta))
    t = (4.*np.cos(eta)*(np.cosh(s)-np.cos(eta))*np.sin(eta)+np.sin(eta)**3)/(2.*np.sqrt(sp.cosh(s)-np.cos(eta)))
    t1 = np.sqrt(sp.cosh(s)-np.cos(eta))*(np.sin(eta))**2
    return(t*tmp+t1*tmp1)

def Q1_s(s,eta,bn,cn):
    tmp=0.0
    tmp1=0.0
    N=len(bn)
    for i in range(1,N+1):
        tmp= tmp+ (bn[i-1]*np.cosh((i+1./2.)*s)+cn[i-1]*np.sinh((i+1./2.)*s))*Pndiff(i,1,np.cos(eta))
        tmp1 = tmp1+(bn[i-1]*(i+1./2.)*sp.sinh((i+1./2.)*s)+cn[i-1]*(i+1./2.)*sp.cosh((i+1./2.)*s))*Pndiff(i,1,np.cos(eta))
    t1 = np.sqrt(sp.cosh(s)-np.cos(eta))*np.sin(eta)
    t = np.sinh(s)*np.sin(eta)/(2.*np.sqrt(sp.cosh(s)-np.cos(eta)))
    return(t1*tmp1+t*tmp)

def Q1_eta(s,eta,bn,cn):
    tmp=0.0
    tmp1=0.0
    N=len(bn)
    for i in range(1,N+1):
        tmp= tmp+(bn[i-1]*sp.cosh((i+1./2.)*s)+cn[i-1]*sp.sinh((i+1./2.)*s))*Pndiff(i,1,np.cos(eta))
        tmp1=tmp1+(bn[i-1]*sp.cosh((i+1./2.)*s)+cn[i-1]*sp.sinh((i+1./2.)*s))*Pndiff(i,2,np.cos(eta))*(-np.sin(eta))

    t =(-2.*np.cos(eta)**2+2.*np.cos(eta)*np.cosh(s)+np.sin(eta)**2)/(2.*np.sqrt(np.cosh(s)-np.cos(eta)))
    t1 = np.sin(eta)*np.sqrt(sp.cosh(s)-np.cos(eta))
    return(t1*tmp1 +t*tmp)

def w1_s(s, eta, an):
    tmp=0.0
    tmp1 = 0.0
    N=len(an)
    for i in range(1,N):
        tmp= tmp+ an[i-1]*np.sinh((i+1./2.)*s)*Pndiff(i,1,np.cos(eta))
        tmp1 = tmp + an[i-1]*(i+1./2.)*np.cosh((i+1./2.)*s)*Pndiff(i,1,np.cos(eta))
    t = np.sin(eta)*np.sinh(s)/(2.*np.sqrt(sp.cosh(s)-np.cos(eta)))
    t1 = np.sqrt(sp.cosh(s)-np.cos(eta))*np.sin(eta)
    return(t1*tmp1 + tmp*t)

def w1_eta(s, eta, an):
    tmp=0.0
    tmp1 = 0.0
    N=len(an)
    for i in range(1,N):
        tmp= tmp+ an[i-1]*np.sinh((i+1./2.)*s)*Pndiff(i,1,np.cos(eta))
        tmp1 = tmp + an[i-1]*np.sinh((i+1./2.)*s)*Pndiff(i,2,np.cos(eta))*(-np.sin(eta))

    t = (-2.*np.cos(eta)**2+2.*np.cos(eta)*np.cosh(s)+np.sin(eta)**2)/(2.*np.sqrt(np.cosh(s)-np.cos(eta)))
    t1 = np.sqrt(sp.cosh(s)-np.cos(eta))*np.sin(eta)
    return(t1*tmp1+t*tmp)

## inverse jacobi matrix for bispherical coordinates (s,eta)-(r,z)

def dsz(c,s,eta):
    return((1.-np.cos(eta)*np.cosh(s))/c)

def detaz(c,s,eta):
    return(-np.sin(eta)*np.sinh(s)/c)

def dsr(c,s,eta):
    return(-np.sin(eta)*np.sinh(s)/c)

def detar(c,s,eta):
    return((-1.+np.cos(eta)*np.cosh(s))/c)
