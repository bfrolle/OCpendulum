################ DISCRETIZATION OF STATES AND INPUTS - SHOOTING GRID DEFINITION
###############################################################################
# introducing src/ShootingGrid.jl

t = 10.               # time horizon
N = 2                 # discretization of shooting grid
NN = 10               # steps inbetween shotting grids
nx = 4                # state dimension
nu = 1                # input dimension 
ncon = 2              # number of gernal inequality constraints
solver = ODE.RK4()    # numerical integration solver
# ODE.Exp4() , ODE.RK4() , ODE.Tsit5() , ODE.QNDF() , CVODE_BDF()

dim = ProblemDimension(t,N,NN,nx,nu,ncon)

################# WARM START
############################
x0 = [0., 0., -0.01, 0.];   # inital condition
u0 = [0.005];               # model input... carrier accelearation
xs,us = WSxu(NlModel!,x0,u0,solver,dim)

# ###### plot warm start
plotSG(NlModel!,xs,us,solver,dim,"WS Pendulum");

#### SHOOTING GRID FUNCTIONS
############################
k = dim.N
# solving the inital value problem 
IVP(NlModel!,xs[k],us[k],solver,dim.dt,dim.dτ)
# solving the inital value problem and return intermediate states 
t,xs_fine = IVPt(NlModel!,xs[k],us[k],solver,dim.dt,dim.dτ)
# local sensitivity of the IVP solution
∂IVP(NlModel!,xs[k],us[k],solver,dim.dt,dim.dτ)