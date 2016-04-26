# Solve the first order two dimensional pde that derived from the boundary condition at eta = 0
#
# U_t(xi,0,t) = (r - q)S(xi) 1/(dSdxi) U_xi(xi,0,t) - r*U(xi,0,t) + kappa*theta* 2/dvdeta* U_eta(xi,0,t)
#

%matplotlib inline
import numpy
import matplotlib.pyplot as plt

# Some function needed because of the change of variables
S       = lambda xi: K + numpy.sinh( c1*xi + c2*(1-xi) )
dSdxi   = lambda xi: numpy.cosh( c1*xi + c2*(1-xi) ) * (c1 - c2)
d2Sdxi2 = lambda xi: numpy.sinh( c1*xi + c2*(1-xi) ) * (c1 - c2)**2
dvdeta  = lambda eta: numpy.cosh( d*eta ) * d * beta
d2vdv2  = lambda eta: numpy.cosh( d*eta ) * d**2 * beta

# Set grid
# Domain (t x \xi x \eta) : [0, T] x [0, 1] x [0,1]
t_final = 
N = 100
M = 100
L = 100

xi  = numpy.linspace(0, 1, N+2) 
eta = numpy.linspace(0, 1, M+2)
t   = numpy.linspace(0, t_final, L)

delta_xi  =  xi[1] -  xi[0]
delta_eta = eta[1] - eta[0]
delta_t   =   t[1] -   t[0]


S_max = 100.0
S_min = 0.0
v_max = 

alpha = 1.0
beta  = 1.0
K     = 
c1    = numpy.arcsinh( (S_max - K)/alpha )
c2    = numpy.arcsinh( (S_min - K)/alpha )
d     = numpy.arcsinh( v_max/beta )

# r     -- domestic interest rate
# q     -- foreign interest rate
# kappa -- mean reversion speed
# theta -- long run variance 
rho   = 1.0
sigma = 1.0
r     = 1.0
q     = 1.0
kappa = 1.0
theta = 1.0

xi  =  xi[1:-1]
eta = eta[1:-1]
