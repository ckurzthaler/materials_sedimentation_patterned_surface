## surface structures
# 2 stripes of sinusoidal, rectangular, triangular, and sawtooth shape
from math_functions import *

# sinusoidal surface
def H_cos(la, xs, ys):
    N = 2
    if(type(xs)==np.float64):
        if((xs+ys)> N*la):
            return(0.)
        elif((xs+ys)<0):
            return(0.)
        else:
            return(np.sin(2.*np.pi*(xs+ys)/la))
    else:
        tmp = np.array([])
        for j in range(0, len(xs)):
            if((xs[j]+ys[j]) > N*la):
                tmp = np.append(tmp,0.)
            elif((xs[j]+ys[j]) < 0):
                tmp = np.append(tmp,0.)
            else:
                tmp = np.append(tmp,np.sin(2.*np.pi*(xs[j]+ys[j])/la))
        return(tmp)

# rectangular surface
def H_rect(la, xs, ys):
    N=2
    if(type(xs)==np.float64):
        for n in range(0,N):
            if(((xs+ys)>= n*la) and (xs+ys)<(n+1/2.)*la):
                return(1.)
            if(((xs+ys)>= (n+1/2.)*la) and (xs+ys)<(n+1)*la):
                return(-1.)
        else:
            return(0.)
    else:
        tmp = xs*0.
        for j in range(0, len(xs)):
            for n in range(0,N):
                if(((xs[j]+ys[j])> n*la) and (xs[j]+ys[j])<(n+1/2.)*la):
                    tmp[j] = 1.
                elif(((xs[j]+ys[j])> (n+1/2.)*la) and (xs[j]+ys[j])<(n+1.)*la):
                    tmp[j] = -1.
        return(tmp)

# sawtooth surface
def H_sawtooth(la, xs, ys):
    if(type(xs)==np.float64):
        if((xs+ys) > 0 and (xs+ys) < la/2):
            return(2*(xs+ys)/la)
        elif((xs+ys) > la/2 and (xs+ys)<3.*la/2.):
            return((2*(xs+ys)-2*la)/la)
        elif((xs+ys) > 3*la/2 and (xs+ys)<2.*la):
            return((2*(xs+ys)-4*la)/la)
        else:
            return(0)
    else:
        tmp = np.array([])
        for j in range(0, len(xs)):
            if((xs[j]+ys[j]) > 0 and (xs[j]+ys[j]) < la/2):
                tmp = np.append(tmp,2*(xs[j]+ys[j])/la)
            elif((xs[j]+ys[j]) > la/2 and (xs[j]+ys[j])<3.*la/2.):
                    tmp = np.append(tmp,(2*(xs[j]+ys[j])-2*la)/la)
            elif((xs[j]+ys[j]) > 3*la/2 and (xs[j]+ys[j])<2.*la):
                    tmp = np.append(tmp,(2*(xs[j]+ys[j])-4*la)/la)
            else:
                tmp = np.append(tmp,0)
        return(tmp)

# triangular surface
def H_triang(la, xs, ys):
    if(type(xs)==np.float64):
        if((xs+ys) > 0 and (xs+ys)<= la/4.):
            return(4*(xs+ys)/la)
        elif((xs+ys) > (la/4.) and (xs+ys) <= (3.*la/4.)):
            return(4*(-xs-ys)/la+2)
        elif((xs+ys) > (3.*la/4.) and (xs+ys) <= (5.*la/4.)):
            return(4*(xs+ys)/la-4)
        elif((xs+ys) > (5.*la/4.) and (xs+ys) <= (7.*la/4.)):
            return(4*(-xs-ys)/la+6)
        elif((xs+ys) > (7.*la/4.) and (xs+ys) <= (2.*la)):
            return(4*(+xs+ys)/la-8)
        else:
            return(0.)
    else:
        tmp = np.array([])
        for j in range(0, len(xs)):
            if((xs[j]+ys[j]) > 0 and (xs[j]+ys[j])<= la/4.):
                tmp = np.append(tmp,4*(xs[j]+ys[j])/la)
            elif((xs[j]+ys[j]) > (la/4.) and (xs[j]+ys[j]) <= (3.*la/4.)):
                tmp = np.append(tmp,4*(-xs[j]-ys[j])/la+2)
            elif((xs[j]+ys[j]) > (3.*la/4.) and (xs[j]+ys[j]) <= (5.*la/4.)):
                tmp = np.append(tmp,4*(xs[j]+ys[j])/la-4)
            elif((xs[j]+ys[j]) > (5.*la/4.) and (xs[j]+ys[j]) <= (7.*la/4.)):
                tmp = np.append(tmp,4*(-xs[j]-ys[j])/la+6)
            elif((xs[j]+ys[j]) > (7.*la/4.) and (xs[j]+ys[j]) <= (2.*la)):
                tmp = np.append(tmp,4*(+xs[j]+ys[j])/la-8)
            else:
                tmp = np.append(tmp, 0.)
        return(tmp)

## Test: plot surfaces
plot_xz = 1
plot_xy = 1
if(plot_xz==1):
    la = 1
    x_vec = np.linspace(-3*la,5*la,1000)
    y_vec = np.linspace(-3*la,5*la,1000)*0

    z_vec = H_cos(la, x_vec,y_vec)
    plt.plot(x_vec/la, z_vec, label = 'sinusoidal')
    z_vec = H_rect(la, x_vec,y_vec)
    plt.plot(x_vec/la, z_vec, label = 'rectangular')
    z_vec = H_sawtooth(la, x_vec,y_vec)
    plt.plot(x_vec/la, z_vec, label = 'sawtooth')
    z_vec = H_triang(la, x_vec,y_vec)
    plt.plot(x_vec/la, z_vec, label = 'triangular')
    plt.legend(loc = 'upper right')
    plt.xlabel('$x/\lambda$')
    plt.ylabel('H(x,y)')
    plt.title('side view')
    plt.show()

if(plot_xy==1):
    la = 1.
    x = np.arange(-2.*la,2*la,0.01)
    y = x
    Y,X = meshgrid(y, x)

    Z = np.cos(X+Y)*0.
    for i in range(0, len(x)):
        for j in range(0,len(y)):
            Z[i,j] = H_cos(la,x[i],y[j])

    plt.contourf(Y/la, X/la, Z, 20, cmap=cm.binary_r, alpha=0.6)
    plt.colorbar(format='%.2f', label = 'surface height $H(x,y)$',ticks=np.array([-1, -0.5, 0,0.5, 1]))
    plt.xlabel('$y/\lambda$')
    plt.ylabel('$x/\lambda$')
    plt.gca().invert_yaxis()
    plt.title('top view: contour plot')
    plt.show()
