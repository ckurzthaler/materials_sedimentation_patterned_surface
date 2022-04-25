from math_functions import *
import scipy.integrate

# eps is a length
# integrand lubrication theory (without angular integral)
def lub_int_r_y(R):
    return(R**3*(1.+4.*R**2)/((1.+R**2)**4))
def lub_int_r_perp(R):
    return(R**2*(1.+7.*R**2)/((1.+R**2)**4))
def lub_int_r_para_1(R):
    return(R*(1.+8.*R**2+25.*R**4)/(1.+R**2)**4)
def lub_int_r_para_2(R):
    return(6.*R**3*(1.+4.*R**2)/(1.+R**2)**4)


def vel_H_lub(h0, a, eps, xs, ys, H, k):
    l = np.sqrt(2.*a*h0)

    # hydrodynamic torques due to translation and rotation
    lyT = (1./10.+43./250.*h0/a)*np.log(2.*a/h0)-0.26227
    lyR = 2./5.*np.log(h0/a)-0.3817

    # hydrodynamic forces due to translation and rotation
    fyT = ((8./15.+64./375.*h0/a)*np.log(2.*a/h0)+0.58461)
    fyR = -2./15.*np.log(h0/a)-0.2526
    fy = fyT-lyT/lyR*fyR

    U = 1./fy
    omega = -lyT/lyR

    # integration in R
    R = np.linspace(0.0001,50,500)
    dR = R[1]-R[0]

    # stress times gradient: r-components
    intR_z = lub_int_r_perp(R)
    intR_y = lub_int_r_y(R)
    intR_x1 = lub_int_r_para_1(R)
    intR_x2 = lub_int_r_para_2(R)

    # integration of phi
    phi = np.linspace(0., 2.*np.pi, 100)
    tmp_R_x = 0.
    tmp_R_y = 0.
    tmp_R_z = 0.

    # Compute integrals
    for i in range(0,len(R)):
        # x-component
        tmp_phi_x1 = H(k, xs+R[i]*l*np.cos(phi), ys+R[i]*l*np.sin(phi))
        tmp_phi_x2 = np.cos(2.*phi)*H(k, xs+R[i]*l*np.cos(phi), ys+R[i]*l*np.sin(phi))
        int_phi_i_x1 = (scipy.integrate.trapz(tmp_phi_x1, phi))
        int_phi_i_x2 = (scipy.integrate.trapz(tmp_phi_x2, phi))
        tmp_R_x += int_phi_i_x1*intR_x1[i]+int_phi_i_x2*intR_x2[i]

        #y-component
        tmp_phi_y = np.sin(2.*phi)*H(k, xs+R[i]*l*np.cos(phi), ys+R[i]*l*np.sin(phi))
        int_phi_i_y = (scipy.integrate.trapz(tmp_phi_y, phi))
        tmp_R_y += int_phi_i_y*intR_y[i]

        #eff_con_y[i]=tmp_R_y

        # z-component
        tmp_phi_z = np.cos(phi)*H(k, xs+R[i]*l*np.cos(phi), ys+R[i]*l*np.sin(phi))
        int_phi_i_z= scipy.integrate.trapz(tmp_phi_z, phi)
        tmp_R_z += int_phi_i_z*intR_z[i]

    u_x = -4.*eps*U/(75.*h0*np.pi*fy)*(1.-omega)**2*tmp_R_x*dR
    u_y = -8.*eps*U/(25.*h0*np.pi*fy)*(1.-omega)*(1.+omega)*tmp_R_y*dR
    u_z = 4.*eps*U*(1.-omega)/(5.*np.pi*l)*tmp_R_z*dR

    return(np.array([U+u_x, u_y, u_z, U, u_x, u_y]))
