######################################## DISCRETIZATION OF THE AUGMENTED SYSTEM
###############################################################################
# introducing src/DMSfunctions.jl

sNp = 0.25;           # desired terminal position sNp = s^[N+1] 
#
t = 10.0              # time horizon
N = 8                 # discretization of shooting grid
NN = 10               # steps inbetween shotting grids
nx = 4                # state dimension
nu = 1                # input dimension 
ncon = 2              # number of gernal inequality constraints
solver = ODE.RK4()    # numerical integration solver
# ODE.RK4() , ODE.Tsit5() 

dim  = ProblemDimension(t,N,NN,nx,nu,ncon)
prob = DMSproblem(sNp,solver,dim)

################# WARM START
############################
xa0 = [0. ,0., 0., -0.01, 0.];   # augmented inital condition
u0 = [0.005]; 
xs,us = WSxua(xa0,u0,dim,prob)

###### plot warm start
plotDMS(xs,us,dim,prob,"Warm Start DMS")

############## DMS FUNCTIONS
############################
k = dim.N
# solving the inital value problem 
prob.xIVP(xs[k],us[k],prob)
# local sensitivity of the IVP solution
prob.∂xIVP(xs[k],us[k],prob)
# cost of the DMS problem
prob.ϑa(xs[end],us[end],prob)
# general constraints DMS problem
prob.cg(xs[k],us[k],par,dim)[1] #residual of maximum acceleration 
# Jacobian of general constraints 
Array(prob.∂cg(xs[k],us[k],par,dim))