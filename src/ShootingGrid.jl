### DISCRETIZATION OF STATES AND INPUTS - SHOOTING GRID DEFINITION
# https://github.tik.uni-stuttgart.de/ac121730/NOOCjulia/wiki/04-Direct-Multiple-Shooting

TimeSeries = Array{AbstractArray,1} #data type for the shooting grid 
######## PROBLEM DIMENSION
##########################
mutable struct DIM
    # shooting grid:
    t::AbstractFloat  # time horizon
    dt::AbstractFloat # time discretization
    dτ::AbstractFloat # intermediate time discretization
    N::Integer     # time grid, number of intervals
    NN::Integer    # steps inbetween shooting intervals
    # system dimensions:
    nx::Integer    # number of states
    nu::Integer    # number of inputs
    nxu::Integer   # total number of states and inputs nxu = nx + nu
    # augmented system:
    nxa::Integer   # number of states of the augmented system nxa = nx+1
    nxua::Integer  # total number of states and inputs of the augmented system nxua = nx+1+nu
    # direct multiple shooting
    ncon::Integer  # number of general constraints  
    Nopt::Integer  # totoal number of optimization variables
    Ncs::Integer   # total number of continuity constraints
    Nc0::Integer   # number of constriants for initial condition and terminal input
    Ncg::Integer   # total number of gneral constraints
    Ncb::Integer   # total number of box constraints
    Ncon::Integer  # total number of constraints
end

"""
    # ProblemDimension(t,N,NN,nx,nu,ncon)::DIM

    Initialize problem dimensions for direct multiple shooting 
    the structure of type DIM contains:
    # shooting grid:
    t...     time horizon
    dt...    time discretization
    dτ...    intermediate time discretization
    N...     time grid, number of intervals dt = t/N
    NN...    steps inbetween shooting intervals dτ = dt/NN
    # system dimensions:
    nx...    number of states
    nu...    number of inputs
    nxu....  total number of states and inputs nxu = nx + nu
    # augmented system:
    nxa...   number of states of the augmented system nxa = nx+1
    nxua...  total number of states and inputs of the augmented system nxua = nx+1+nu
    # direct multiple shooting
    ncon...  number of general constraints  
    Nopt...  totoal number of optimization variables
    Ncs...   total number of continuity constraints
    Nc0...   number of constriants for initial condition and terminal input
    Ncg...   total number of gneral constraints
    Ncb...   total number of box constraints
    Ncon...  total number of constraints
"""
function ProblemDimension(t,N,NN,nx,nu,ncon)::DIM
    # shooting grid:
    dt = t/N
    dτ = dt/NN
    # system dimensions:
    nxu = nx + nu
    # augmented system:
    nxa = nx+1
    nxua = nxa + nu
    # direct Multiple Shooting
    Nopt = (N+1)*(nxa+nu)
    Ncs  = N*nxa
    Nc0  = nxa+nu
    Ncg  = (N+1)*(ncon)
    Ncb  = 2*N*(nxa+nu)
    Ncon = Ncs + Nc0 + Ncg + Ncb
    
    return DIM(t,dt,dτ,N,NN,
    nx,nu,nxu,
    nxa,nxua,
    ncon,Nopt,Ncs,Nc0,Ncg,Ncb,Ncon)
end

#### INITIAL VALUE PROBLEM
##########################
"""
    # IVP(f,x0,u0,solver,dt,stepsize)
    solves the initial value problem for the shooting interval of length dt 
            
            dx = f(x,u)
            x0 = x0
            u  = u0 ∀ t ∈ (0,dt)

    and returns the solution of system f at the end of the shooting interval.
    
    # ODE solver options
    if the numerical integration solver does not have adaptivity, a fixed stepsize must be specified.
    if the solver has adaptivity, stepsize sets the initial stepsize.
    solver examples: adatptive → RK4(), Tsit5() QNDF(), fixed stepsize → Exp4().
"""
function IVP(f,x0,u0,solver,dt,stepsize)
    prob = ODE.ODEProblem(f,x0,(0,dt),u0)     
    # only return final step of integration:
    sol  = ODE.solve(prob, solver, dt = stepsize, save_everystep=false)[end]  
    
    return sol[:]
end

"""
    # IVPt(f,x0,u0,solver,dt,stepsize)
    solves the initial value problem for the shooting interval of length dt 
            
            dx = f(x,u)
            x0 = x0
            u  = u0 ∀ t ∈ (0,dt)
    
    Returns the tuple of the time vector t (containing all intermediate time steps)
    and the corresponding solution x of system f.
    
    # ODE solver options
    if the numerical integration solver does not have adaptivity, a fixed stepsize must be specified.
    if the solver has adaptivity, stepsize sets the initial stepsize.
    solver examples: adatptive → RK4(), Tsit5() QNDF(), fixed stepsize → Exp4().
"""
function IVPt(f,x0,u0,solver,dt,stepsize)
    prob = ODE.ODEProblem(f,x0,(0,dt),u0)              # define initial value problem (IVP)
    sol  = ODE.solve(prob, solver, dt = stepsize)      # save all steps of integration
    
    return sol.t[:], sol[:]
end

"""
    # ∂IVP(f,x0,u0,solver,dt)
    solves the local sensitivity ODE  
            
            d{∂x/∂x0} = ∂f/∂x ∂x/∂x0 
            d{∂x/∂u0} = ∂f/∂x ∂x/∂u0 + ∂f/∂u

    of the initial value problem

            dx = f(x,u)
            x0 = x0
            u  = u0 ∀ t ∈ (0,dt)

    returns the local sensitivity [∂x/∂x0 | ∂x/∂u0] at the end of the shooting interval dt.
    
    # ODE solver options
    if the numerical integration solver does not have adaptivity, a fixed stepsize must be specified.
    if the solver has adaptivity, stepsize sets the initial stepsize.
    solver examples: adatptive → RK4(), Tsit5() QNDF(), fixed stepsize → Exp4().
"""
function ∂IVP(f,x0,u0,solver,dt,stepsize)
    nx = length(x0)
    nu = length(u0) 
    return ForwardDiff.jacobian(xu -> IVP(f,xu[1:nx],xu[nx+1:end],solver,dt,stepsize), vcat(x0,u0))
end
# SOLVE ODES WITH concret_solve VIA SENSITIVITY ALGORITHMS
# sensitvity algortihms see https://docs.juliadiffeq.org/stable/analysis/sensitivity/#Sensitivity-Algorithms-1

############### WARM START
##########################
"""
    # WSu(u0,dim::DIM)
    initialize a TimeSeries for the input trajectory u on 
    the shooting grid specified via the structure dim::DIM.
    It is assumed that the trajectory is constant over the
    entire time horizon. 

    returns the TimeSeries us.
"""
function WSu(u0,dim::DIM)::TimeSeries

    us = TimeSeries(undef,dim.N+1)

    for k in 1:dim.N
         us[k] = u0
    end
    
    us[dim.N+1] = us[dim.N]
    return us
end

"""
    # WSu(u0,dim::DIM)
    initialize a TimeSeries for the input trajectory u and
    the state trajecotry x on the shooting grid specified 
    via the the structure dim::DIM.
    It is assumed that the trajectory u is constant over the
    entire time horizon. 

    returns the TimeSeries xs,us.
"""
function WSxu(f,x0,u0,solver,dim::DIM)
    us = WSu(u0,dim)
    xs = TimeSeries(undef,dim.N+1)
    xs[1] = deepcopy(x0);
    for k in 1:dim.N
        xs[k+1] =  IVP(f,xs[k],us[k],solver,dim.dt,dim.dτ)
    end

    return xs,us
end