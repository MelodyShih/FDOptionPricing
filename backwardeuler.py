import numpy
import scipy.sparse as sparse
import scipy.sparse.linalg as linalg
import matplotlib.pyplot as plt

# Method of lines
# Second order approximation in space; Backward Euler in time
# AU^{k+1} = U^{k} + rhs

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
K     = 100

S_max = 200
S_min = 0
v_max = 1.0

alpha = 1.0 
beta  = 1.0 
c1    = numpy.arcsinh( (S_max - K)/alpha )
c2    = numpy.arcsinh( (S_min - K)/alpha )
d1     = numpy.arcsinh( v_max/beta )

# Set Coefficient 
# (Citation: https://drive.google.com/drive/u/1/folders/0B97aYusZC-S0NkJ5OURlbnNCZlk)
# rho   -- instantaneous correlation between the Brownian motions
# sigma -- volatility of variance
# r     -- domestic interest rate
# q     -- foreign interest rate
# kappa -- mean reversion speed
# theta -- long run variance 

rho   = 0.6
sigma = 0.04
r     = 0.01
q     = 0.04
kappa = 3
theta = 1.0 #???

# Set grid
# Domain (t x \xi x \eta) : 
#
# [0, T] x [0, 1] x [0,1]
#  xi[j] = \Delta xi  * j, j = 0, ..., N+1
# eta[j] = \Delta eta * j, j = 0, ..., M+1
#   t[j] = \Delta t   * j, j = 0, ..., L-1
#

t_final = 1.0
N = 100
M = 100
L = 100

xi  = numpy.linspace(0, 1, N+2) 
eta = numpy.linspace(0, 1, M+2)
t   = numpy.linspace(0, t_final, L)

xi  =  xi[1:-1] # Need this
eta = eta[1:-1] # Need this

delta_xi  =  xi[1] -  xi[0]
delta_eta = eta[1] - eta[0]
delta_t   =   t[1] -   t[0]

# Some function needed because of the change of variables
S       = lambda xi: K + numpy.sinh( c1*xi + c2*(1-xi) )
v       = lambda eta: beta*numpy.sinh( d1*eta )
dSdxi   = lambda xi: numpy.cosh( c1*xi + c2*(1-xi) ) * (c1 - c2)
d2Sdxi2 = lambda xi: numpy.sinh( c1*xi + c2*(1-xi) ) * (c1 - c2)**2
dvdeta  = lambda eta: numpy.cosh( d1*eta ) * d1 * beta
d2vdeta2  = lambda eta: numpy.cosh( d1*eta ) * d1**2 * beta



# Allocate solution vector U
U = numpy.empty([N+2, M+2])
U_interior = U[1:-1, 1:-1]

# Note: row-wise reshape, 
# a = [[1,2,3]
#      [4,5,6]
#      [7,8,9]]
# a.reshape(9) = [1,2,3,4,5,6,7,8,9]
#
# U_hat = [U[1,1], U[1,2], ..., U[1, M], U[2,1], ..., U[2,M], ..., U[N,1], ..., U[N,M]]

U_interior_old = numpy.zeros(N*M)
U_interior_new = numpy.zeros(N*M)

# Construct A
# Size: MN * MN
a = lambda i,j: delta_t*rho*sigma/(4.0*delta_xi*delta_eta)*v(eta[j])*S(xi[i])/(dSdxi(xi[i]) * dvdeta(eta[j]))
b = lambda i,j: sigma**2 * delta_t/(2.0*delta_eta**2)*v(eta[j])/dvdeta(eta[i])**2 - \
				delta_t/(2.0*delta_eta)*(kappa*(theta-v(eta[j]))/dvdeta(eta[j]) - 0.5*sigma**2*v(eta[j])*d2vdeta2(eta[j])/dvdeta(eta[j])**3)
c = lambda i,j: delta_t/delta_xi**2*v(eta[j])*S(xi[i])**2/dSdxi(xi[i])**2 - \
				delta_t/(2.0*delta_xi)*((r-q)*S(xi[i])/dSdxi(xi[i]) - 0.5*v(eta[j])*S(xi[i])**2*d2Sdxi2(xi[i])/dSdxi(xi[i])**3)
d = lambda i,j: 1+r*delta_t+delta_t/delta_xi**2*v(eta[j])*S(xi[i])**2/dSdxi(xi[i])**2+sigma**2*delta_t/delta_eta**2*v(eta[j])/dvdeta(eta[i])**2
e = lambda i,j: delta_t/delta_xi**2*v(eta[j])*S(xi[i])**2/dSdxi(xi[i])**2 + \
				delta_t/(2*delta_xi)*((r-q)*S(xi[i])/dSdxi(xi[i]) - 0.5*v(eta[j])*S(xi[i])**2*d2Sdxi2(xi[i])/dSdxi(xi[i])**3)
f = lambda i,j: sigma**2*delta_t/(2*delta_eta**2)*v(eta[j])/dvdeta(eta[i])**2 + \
				delta_t/(2.0*delta_eta)*(kappa*(theta - v(eta[j]))/dvdeta(eta[j]) - 0.5*sigma**2*v(eta[j])*d2vdeta2(eta[j])/dvdeta(eta[j])**3)

# Allocate memory
A = numpy.zeros([M*N, M*N])
rhs = numpy.zeros(N*M)


# Initialize the U vector
for i in range(N):
	U_interior_old[i*M: (i+1)*M] = numpy.maximum(numpy.ones(M) * xi[i] - K , 0.0)

# Boundary condition 
# For bottom bc xi = 0 
_l = 2.0 * dSdxi(xi[0]) / ( dSdxi(xi[0]) + 0.5 * delta_xi * d2Sdxi2(xi[0]) )
l_ = - ( dSdxi(xi[0]) - 0.5* delta_xi * d2Sdxi2(xi[0]) ) / ( dSdxi(xi[0]) + 0.5* delta_xi* d2Sdxi2(xi[0]) )

# For upper bc xi = 1
_r = - ( dSdxi(xi[-1]) + 0.5 * delta_xi * d2Sdxi2(xi[-1]) ) / ( dSdxi(xi[-1]) - 0.5 * delta_xi * d2Sdxi2(xi[-1]) )
r_ = 2.0 * dSdxi(xi[-1]) / ( dSdxi(xi[-1]) - 0.5 * delta_xi * d2Sdxi2(xi[-1]) )

# For right bc eta = 1
# rhs[i] = -a[i, M] * S(xi[i-1]) + f[i,M] * S(xi[i]) + a[i,j] * S(xi[i+1])

for irow in range(M*N):
    # U_hat[irow] == U_interior[i,j]
    i = irow/M
    j = numpy.mod(irow, N)

    if( i == 0 ):
        if (j == 0):
            A[irow, irow + M]     = -e(i,j)
            A[irow, irow + M + 1] = -a(i,j)
            A[irow, irow]         =  d(i,j)
            A[irow, irow + 1]     = -f(i,j)
        elif (j == M-1):
            A[irow, irow + M - 1] =  a(i,j)
            A[irow, irow + M]     = -e(i,j)
            A[irow, irow - 1]     = -b(i,j)
            A[irow, irow]         =  d(i,j)
        else: # Add term for bc condition
            A[irow, irow + M - 1] =  a(i,j) - a(i,j)*l_
            A[irow, irow + M]     = -e(i,j) - c(i,j)*l_
            A[irow, irow + M + 1] = -a(i,j) + a(i,j)*l_
            A[irow, irow - 1]     = -b(i,j) - a(i,j)*_l
            A[irow, irow]         =  d(i,j) - c(i,j)*_l
            A[irow, irow + 1]     = -f(i,j) + a(i,j)*_l
    elif (i == N-1):
        if (j == 0):
            A[irow, irow - M ]    = -c(i,j)
            A[irow, irow - M + 1] =  a(i,j)
            A[irow, irow]         =  d(i,j)
            A[irow, irow + 1]     = -f(i,j)
        elif (j == M-1):
            A[irow, irow - M - 1] = -a(i,j) 
            A[irow, irow - M ]    = -c(i,j)
            A[irow, irow - 1]     = -b(i,j)
            A[irow, irow]         =  d(i,j)
        else:
            A[irow, irow - M - 1] = -a(i,j) + a(i,j)*_r
            A[irow, irow - M ]    = -c(i,j) - e(i,j)*_r
            A[irow, irow - M + 1] =  a(i,j) - a(i,j)*_r
            A[irow, irow - 1]     = -b(i,j) + a(i,j)*r_
            A[irow, irow]         =  d(i,j) - e(i,j)*r_
            A[irow, irow + 1]     = -f(i,j) - a(i,j)*r_
    else:
        if (j == 0):   
            A[irow, irow - M ]    = -e(i,j)
            A[irow, irow - M + 1] = -a(i,j)
            A[irow, irow]         =  d(i,j)
            A[irow, irow + 1]     = -f(i,j)
            A[irow, irow + M]     = -c(i,j)
            A[irow, irow + M + 1] =  a(i,j)
        elif (j == M-1):
            A[irow, irow - M - 1] = -a(i,j) 
            A[irow, irow - M ]    = -c(i,j)
            A[irow, irow - 1]     = -b(i,j)
            A[irow, irow]         =  d(i,j)
            A[irow, irow + M - 1] =  a(i,j)
            A[irow, irow + M]     = -e(i,j)
            rhs[irow] = -a(i,j)*S(xi[i-1]) + f(i,j)*S(xi[i]) + a(i,j)*S(xi[i+1])
        else:
            A[irow, irow - M - 1] = -a(i,j) 
            A[irow, irow - M ]    = -c(i,j)
            A[irow, irow - M + 1] =  a(i,j)
            A[irow, irow - 1]     = -b(i,j)
            A[irow, irow]         =  d(i,j)
            A[irow, irow + 1]     = -f(i,j)
            A[irow, irow + M - 1] =  a(i,j)
            A[irow, irow + M]     = -e(i,j)
            A[irow, irow + M + 1] = -a(i,j)

for timestep in t:
	# For each time step, we need to first solve the pde for lower boundary condition
	# lower bc
	# Construct matrix A_bc for boundary condition
	alpha_bc = numpy.multiply( (r-q) * delta_t/(2.0 * delta_xi), numpy.divide(S(xi), dSdxi(xi)) ) # vector
	beta_bc = kappa*theta*delta_t / ( 2*delta_eta*dvdeta(0) )

	maindiag = numpy.ones(N)*(1+r*delta_t)
	maindiag[0] += alpha_bc[0] * _l
	maindiag[-1] += alpha_bc[-1] * r_

	upperdiag = -alpha_bc
	upperdiag[0] += alpha_bc[0] * l_

	lowerdiag = numpy.zeros(N)
	lowerdiag[:-1] = alpha_bc[1:]
	lowerdiag[-2] -= alpha_bc[-1]*_r

	A_bc = sparse.spdiags([lowerdiag, maindiag, upperdiag], [-1, 0, 1], N, N).tocsr()

	#initialize U_bc vector
	U_bc = numpy.zeros(N)
	U_bc = numpy.maximum(0, S(xi) - K)

	#rhs (Need U vector in previous time step)
	rhs_bc = numpy.zeros(N)
	for i in range(len(rhs_bc)):
	    rhs_bc[i] = (-3.0 * beta_bc + 1)*U_bc[i] + 4.0*U_interior_old[i*M] - beta*U_interior_old[i*M]
	U_bc = linalg.spsolve(A_bc, rhs_bc)

	# After solve the lower boundary condition, update the corresponding rhs
	for i in range(N):
		if i==0 or i==N-1:
			continue
		rhs[i*M] = a(i,0) * U_bc[i-1] + b(i,0) * U_bc[i] - a(i,0) * U_bc[i+1]

	U_interior_new = linalg.spsolve(A, U_interior_old + rhs)
	U_interior_old = numpy.copy(U_interior_new)

U_interior = numpy.reshape(U_interior_new, [M, N])
U[1:-1, 1:-1] = U_interior
# Need to fill in boundary value







