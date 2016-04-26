%matplotlib inline
import numpy
import matplotlib.pyplot as plt

# Method of lines
# Second order approximation in space; Backward Euler in time

# Heston pde
# U(S, v, t), the europeon call option price
# v     -- variance
# s     -- stock price

# Change of variable
# S --> S(\xi)
# v --> v(\eta)
# alpha -- 
# beta  --
# Range of stock price
S_max = 100.0
S_min = 0.0
v_max = 

alpha = 1.0
beta  = 1.0
K     = 
c1    = numpy.arcsinh( (S_max - K)/alpha )
c2    = numpy.arcsinh( (S_min - K)/alpha )
d     = numpy.arcsinh( v_max/beta )
# Set Coefficient 
# (Citation: https://drive.google.com/drive/u/1/folders/0B97aYusZC-S0NkJ5OURlbnNCZlk)
# rho   -- instantaneous correlation between the Brownian motions
# sigma -- volatility of variance
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

# Set grid
# Domain (t x \xi x \eta) : 
#
# [0, T] x [0, 1] x [0,1]
#
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

# Some function needed because of the change of variables
S       = lambda xi: K + numpy.sinh( c1*xi + c2*(1-xi) )
dSdxi   = lambda xi: numpy.cosh( c1*xi + c2*(1-xi) ) * (c1 - c2)
d2Sdxi2 = lambda xi: numpy.sinh( c1*xi + c2*(1-xi) ) * (c1 - c2)**2
dvdeta  = lambda eta: numpy.cosh( d*eta ) * d * beta
d2vdv2  = lambda eta: numpy.cosh( d*eta ) * d**2 * beta

xi  =  xi[1:-1]
eta = eta[1:-1]
# Boundary condition 
# For left bc i.e. xi[0] 
_l = 2.0 / ( 1.0 + 0.5 * delta_xi * d2Sdxi2(xi[1]) )
l_ = - ( dSdxi(xi[1]) - 0.5 * d2Sdxi2(xi[1]) ) / ( dSdxi(xi[1]) + 0.5 * d2Sdxi2(xi[1]) )

# For right bc i.e. xi[-1]
_r = 2.0 / ( 1.0 - 0.5 * delta_xi * d2Sdxi2(xi[1]) )
r_ = - ( dSdxi(xi[-1]) + 0.5 * d2Sdxi2(xi[-1]) ) / ( dSdxi(xi[-1]) - 0.5 * d2Sdxi2(xi[-1]) )

# For upper bc i.e eta[-1] (Move to rhs)
# rhs[i] = -a[i, M] * S(xi[i-1]) + f[i,M] * S(xi[i]) + a[i,j] * S(xi[i+1])
rhs = numpy.empty(N*M)
rhs[0:M-1] =  

# For lower bc (requires us to first solve a firstorder two-dimensional PDE)
