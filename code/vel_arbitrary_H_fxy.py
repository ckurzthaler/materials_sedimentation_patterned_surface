from math_functions import *
import scipy.integrate
from time import time

from bispherical_translation_xy import *
from bispherical_translation_z import *
from bispherical_rotation_xy import *

def vel_H_phi0(h0, a, eps, xs, ys, H, k, phi0):

    Ns=40
    #bispherical coordinates
    alpha=np.arccosh((h0+a)/a)
    c=a*np.sinh(alpha)
    # calculation of coefficients
    # TRANSLATION
    Nc=700
    an,bn=coef(alpha,Nc)
    a1 = np.linalg.solve(an, bn)
    n1 = 80
    a1 = a1[0:n1]
    b1 = coefbi(a1)
    c1 = coefci(a1,alpha)
    d1 = coefdi(a1)
    e1 = coefei(a1,alpha)
    f1 = coeffi(a1)
    g1 = coefgi(a1,alpha)
    # ROTATION
    anR,bnR=coefR(alpha,Nc)
    a1R = np.linalg.solve(anR, bnR)
    a1R = a1R[0:n1]
    b1R = coefbiR(a1R)
    c1R = coefciR(a1R,alpha)
    d1R = coefdiR(a1R)
    e1R = coefeiR(a1R,alpha)
    f1R = coeffiR(a1R)
    g1R = coefgiR(a1R,alpha)

    # integration in R
    R = np.linspace(0.0001,50,500)
    dR = R[1]-R[0]

    # bispherical coordinates
    R1=np.sqrt((c*R)**2)
    Q= np.sqrt((((R1**2)+(c**2))**2))
    s=0.
    eta = np.arccos((((R1**2)-(c**2))/(Q)))

    # rotation induced by translation
    lyR = torque_yR(alpha,a1R,b1R,c1R,d1R,e1R)
    lyT = torque_y(alpha,a1,b1,c1,d1,e1)
    omega = lyT/lyR*c/a

    # hydrodynamic force for parallel motion
    fy_T = force_y(alpha,e1,c1)
    fy_R = force_xR(alpha,e1R,c1R)
    fy = fy_T-lyT/lyR*fy_R

    U = 1./fy

    # hydrodynamic force for perpendicular motion
    rperp_f = force_z(alpha,Ns)

    # integration kernel (R,alpha dependent only)

    # velocity gradient (r,z components)
    ur_z =(uxr0_z(s,eta,R,b1,c1,d1,e1,f1,g1) + omega*uxr0_z(s,eta,R,b1R,c1R,d1R,e1R,f1R,g1R))
    utheta_z = uxtheta0_z(s,eta,R,b1,c1,d1,e1,f1,g1)+omega*uxtheta0_z(s,eta,R,b1R,c1R,d1R,e1R,f1R,g1R)

    # stress-components
    # x-direction
    intR_x1 =R*(ur_z*ur_z+utheta_z*utheta_z)
    intR_x2 =R*(ur_z*ur_z-utheta_z*utheta_z)

    # y-direction
    z=0.
    sigmay_ztheta = sigmay0_ztheta(s,eta,z,R,a1,b1,c1,d1,e1,f1,g1) + omega*sigmay0_ztheta(s,eta,z,R,a1R,b1R,c1R,d1R,e1R,f1R,g1R) # vorher -
    sigmay_zr = sigmay0_zr(s,eta,z,R,a1,b1,c1,d1,e1,f1,g1) + omega*sigmay0_zr(s,eta,z,R,a1R,b1R,c1R,d1R,e1R,f1R,g1R)
    intR_y_1 =R*(ur_z*sigmay_zr+utheta_z*sigmay_ztheta)
    intR_y_2 =R*(ur_z*sigmay_zr-utheta_z*sigmay_ztheta)

    # z-direction
    sigmaz_zr = -sigmazr(alpha,R,s,eta,Ns)
    intR_z =ur_z*sigmaz_zr*R

    # integration of phi
    phi = np.linspace(0., 2.*np.pi, 100)
    tmp_R_x = 0.
    tmp_R_y = 0.
    tmp_R_z = 0.

    # Compute integrals
    for i in range(0,len(R)):
        # x-component
        tmp_phi_x1 = H(k, xs+R[i]*c*np.cos(phi), ys+R[i]*c*np.sin(phi))
        tmp_phi_x2 = np.cos(2.*phi-phi0)*H(k, xs+R[i]*c*np.cos(phi), ys+R[i]*c*np.sin(phi))
        int_phi_i_x1 = (scipy.integrate.trapz(tmp_phi_x1, phi))/2.
        int_phi_i_x2 = (scipy.integrate.trapz(tmp_phi_x2, phi))/2.
        tmp_R_x += int_phi_i_x1*intR_x1[i]*np.cos(phi0)+int_phi_i_x2*intR_x2[i]

        #y-component
        tmp_phi_y1 = np.sin(2.*phi-phi0)*H(k, xs+R[i]*c*np.cos(phi), ys+R[i]*c*np.sin(phi))
        int_phi_i_y1 = (scipy.integrate.trapz(tmp_phi_y1, phi))/2.
        tmp_R_y += int_phi_i_y1*intR_y_1[i]+int_phi_i_x1*intR_y_2[i]*np.sin(phi0)

        # z-component
        tmp_phi_z = np.cos(phi-phi0)*H(k, xs+R[i]*c*np.cos(phi), ys+R[i]*c*np.sin(phi))
        int_phi_i_z= scipy.integrate.trapz(tmp_phi_z, phi)
        tmp_R_z += int_phi_i_z*intR_z[i]

    # VELOCITIES
    u_x = -eps*U/(6.*a*np.pi*fy)*tmp_R_x*dR
    u_y = -eps*U/(6.*a*np.pi*fy)*tmp_R_y*dR
    u_z = -eps*U*a/(6.*np.pi*c**2*rperp_f)*tmp_R_z*dR


    return(np.array([U*np.cos(phi0)+u_x, U*np.sin(phi0)+u_y, u_z, U, u_x, u_y]))
