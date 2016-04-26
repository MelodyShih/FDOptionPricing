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
S_max = 
S_min =

alpha = 
beta  = 

# Set Coefficient 
# (Citation: https://drive.google.com/drive/u/1/folders/0B97aYusZC-S0NkJ5OURlbnNCZlk)
# rho   -- instantaneous correlation between the Brownian motions
# sigma -- volatility of variance
# r     -- domestic interest rate
# q     -- foreign interest rate
# kappa -- mean reversion speed
# theta -- long run variance 

rho   =
sigma = 
r     = 
q     = 
kappa = 
theta =  

# Set grid
# Domain (t x \xi x \eta) : 
#
# [0, T] x [0, 1] x [0,1]
#
t_final = 
n = 100
m = 100
time_steps = 100

delta_xi  = numpy.linspace(0, 1, n) 
delta_ita = numpy.linspace(0, 1, m)
delta_t   = numpy.linspace(0, t_final, time_steps)

# Some function needed because of the change of variables
function DSDxi = 