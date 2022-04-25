##################################################################
## Compute trajectory of a sphere moving near a structured surface
##################################################################
# import packages, functions for surface structures, and
# roughness-induced velocities:
from surface_structures import *
from vel_arbitrary_H_fxy import *
from vel_arbitrary_H_lub import *
##################################################################
## Parameters
# particle radius (length scale)
a = 1.5875
# surface wavelength (length scale)
la = 6.
# surface amplitude (length scale)
eps = 0.15
# initial particle-surface distance (length scale)
h_in = 0.2
# direction of the force (here: gravity)
phi0 = 0.
##################################################################
## simulation
## non-dimensional parameters
# time scales rescaled by tau = la/U (U is settling velicity in free space)
# time steps
dt = 0.02
N=1000
# length scales rescaled by la (xs, ys are dimensionless)
xs = -4.
ys  = 0
h0 = h_in/la

print('time, xs, ys, h, ux, uy, uz, u0, ux1, uy1')
for i in range(0,N):
    # exact, bispherical solution for the flow fields
    if (h0>=0.1/la):
        u = vel_H_phi0(h0*la, a, eps, xs*la, ys*la, H_rect, la, 0.)
    # for particles very close to surface we use the lubrication approximation
    if (h0<=0.1/la):
        u = vel_H_lub(h0*la, a, eps, xs*la, ys*la, H_rect, la)
    print((i*dt), xs,ys,h0, u[0], u[1], u[2], u[3], u[4],u[5])
    # propagate
    xs+=dt*u[0]
    ys+=dt*u[1]
    h0+=dt*u[2]
