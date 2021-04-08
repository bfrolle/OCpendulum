### DIRECT MULTIPLE SHOOTING PORBLEM FORMULATION

#variable notation
# x,u                   ... original state and input
# xa = vcat(xf,x,u)     ... augmented state
# xk = x^[k]            ... state at discrete time k
# xkp = x^[k+1]         ... state at discrete time k+1

mutable struct DMSprob
    # cost and system functions
    ϑa::Function            # DMS terminal constraint 
    xIVP::Function          # IVP solution of a single shooting interval
    ∂xIVP::Function         # local sensitivita of the IVP solution
    cg::Function            # general constraints
    ∂cg::Function           # Jacobian of general constraints
    # cost parmaeter
    wL::AbstractFloat       # weight for xf = ∫L (L...Lagrang term)
    wx::AbstractArray       # weights for Lagrang term x'Q x with Q diagonal
    wu::AbstractArray       # weights for Lagrang term u'R u with R diagonal
    xdes::AbstractArray     # desired terminal state x ∈ R^nx
    # box constraints
    xalb::AbstractArray     # lower bound of augmanted state ∈ R^nx+1
    xaub::AbstractArray     # upper bound of augmented state ∈ R^nx+1
    # ODE solver
    solver::Any             # ODE Solver
    dt::AbstractFloat       # time discretization
    dτ::AbstractFloat       # intermediate time discretization
    # optimizer attributes: 
    max_iter::Integer       # maximum number of iterations 
    print_level::Integer    # define verbosity \ print level 
end

"""
    # DMSproblem(sNp,solver,dim::DIM)::DMSprob
    initialize DMS problem in dependence of the desired final position of the pendulum
    sNp = s^[N+1], the selected ODE solver, and the shooting grid defined within the 
    structure dim::DIM.

    The structure contains 
    # cost and system functions
    ϑ::Function                 DMS terminal constraint → Cost(xaNp,uNp,prob::DMSprob)
    xIVP::Function              IVP solution of a single shooting interval → xIVP(xak,uk,prob::DMSprob)
    ∂xIVP::Function             local sensitivita of the IVP solution → ∂xIVP(xa0,u0,prob::DMSprob)
    cg::Function                general constraints
    ∂cg::Function               Jacobian of general constraints
    # cost parmaeter
    wL...                       weight for xf = ∫L (L...Lagrang term)
    wx...                       weights for Lagrang term x'Q x with Q diagonal
    wu...                       weights for Lagrang term u'R u with R diagonal
    xdes...                     desired terminal state x ∈ R^nx
    # box constraints
    xalb...                     lower bound of augmanted state ∈ R^nx+1
    xaub...                     upper bound of augmented state ∈ R^nx+1
    # ODE solver
    solver...                   ODE Solver
    dt...                       time discretization
    dτ...                       intermediate time discretization
    # optimizer attributes: 
    max_iter...                 maximum number of iterations 
    print_level...              define verbosity of solver output 

"""
function DMSproblem(sNp,solver,dim::DIM)::DMSprob
    # cost parmaeter
    wL =   0.1                              # weight for xf = ∫L (L...Lagrang term)
    wx =   [10.,   1., 1., 1.]              # weights for Lagrang term x'Q x with Q diagonal
    wu =   [1.]                             # weights for Lagrang term u'R u with R diagonal
    xdes = [sNp, 0., 0.,  0.]               # desired terminal state x ∈ R^nx
    # box constraints
    dϕ_max = -40*π/180
    xalb = [0,  par.s_min, par.ds_min,  dϕ_max,-Inf]   # lower bound of augmanted state ∈ R^nx+1
    xaub = [Inf,par.s_max, par.ds_max, -dϕ_max, Inf]   # upper bound of augmented state ∈ R^nx+1 
    # optimizer attributes:
    max_iter    = 500
    print_level = 5

    return DMSprob(Cost,xIVP,∂xIVP,GeneralConstraint,∂GeneralConstraint,
    wL,wx,wu,xdes,
    xalb,xaub,
    solver,dim.dt,dim.dτ,
    max_iter,print_level)
end
######### COST FORMULATION
##########################
"""
    # LagrangeTerm(x,u,prob::DMSprob) 

    the Lagrange term of the original problem formulation
    at discrete time k is chosen to be of the quadratic form

        L^[k] = ||xdes - x^[k]||wx + ||u^[k]||wu

    wich is equal to the sum of the weighted (wx) distance of the desired
    state xades to x^[k] and the weighted (wu) norm of the input u^[k].

    xdes, wx and wu are defined within the structure prob::DMSprob.
"""
function LagrangeTerm(x,u,prob::DMSprob)  
     return sum(prob.wx[i]*(x[i]-prob.xdes[i])^2 for i in 1:dim.nx) + sum(prob.wu[i]*u[i]^2 for i in 1:dim.nu)
end

"""
    # Cost(xaNp,uNp,prob::DMSprob) 

    the terminal cost is chosen as

        wL*xf + ||xdes - x^[N+1]||wx + ||u^[N+1]||wu

    and is equal to the sum of the weighted (wx) distance of the desired
    state xdes to the final state xaNp = x^[N+1], the weighted (wu) norm 
    of the final input uNp = u^[N+1], and the weighted integral of the
    Lagrange term xf = ∫L.

    xdes, wL, wx and wu are defined within the structure prob::DMSprob.
"""
function Cost(xaNp,uNp,prob::DMSprob)
    return prob.wL*xaNp[1] + sum(prob.wx[i]*(xaNp[i+1]-prob.xdes[i])^2 for i in 1:dim.nx) + sum(prob.wu[i]*uNp[i]^2 for i in 1:dim.nu)
end
##### AUGMENTED ODE SYSTEM
##########################
"""
    # AugmentedSystem!(dxa,xa,up,t)

    state space model of the augmented system

        dxa = [L(x,u) | f(x,u)], 
        xa0 = [0 | x0],

    where L is the Lagrangian term and f:R^nx x R^nu → R^nx is the original state space model.
    
    The input up = [u | prob::DMSprob] contains both the model input u and the parameter structure
    with the information of the weights for the Lagrangian term. 
"""
function AugmentedSystem!(dxa,xa,up,t)
    original_state_ind = 2:dim.nxa #origianal state index
    # Lagrange term
    dxa[1] = LagrangeTerm(xa[original_state_ind],up[1],up[2]) #augmented state xf
    # nonlinear system model
    dxa[original_state_ind] = NlModel!(dxa[original_state_ind],xa[original_state_ind],up[1],t)
    return dxa
end

"""
    # xIVP(xak,uk,prob::DMSprob)
    solves the initial value problem for the shooting interval of length dt
    for the augmented system 
            
        dxa = [L(x,u) | f(x,u)]
        xa0 = [0  | x0]
        u  = u0 for t ∈ (0,dt)

    and returns the solution of the system AugmentedSystem! at the end of
    the shooting interval.
"""
function xIVP(xak,uk,prob::DMSprob)
    return IVP(AugmentedSystem!,xak,[uk,prob],prob.solver,prob.dt,prob.dτ)
end

"""
    ∂xIVP(xa0,u0,prob::DMSprob)
    solves the local sensitivity ODE for the augmented system fa = [L | f]

            d{∂xa/∂x0} = ∂fa/∂xa ∂xa/∂xa0 
            d{∂xa/∂u0} = ∂fa/∂xa ∂xa/∂u0 + ∂fa/∂u

    according to the initial value problem

            dxa = [L(x,u) | f(x,u)]
            xa0 = [0  | x0]
            u  = u0 for t ∈ (0,dt)

    wehre L is the Lagrangian term and f is the original state space model
    returns the local sensitivity [∂xa/∂xa0 | ∂xa/∂u0] at the end of the shooting interval.
"""
function ∂xIVP(xa0,u0,prob::DMSprob) 
    return ForwardDiff.jacobian(xu -> xIVP(xu[1:dim.nxa],xu[dim.nxa+1:dim.nxua],prob) , vcat(xa0,u0))
end
############## CONSTRAINTS 
##########################
function GeneralConstraint(xak,uk,par,dim)
    rcg = zeros(Float64,dim.ncon)
    # mges*dssmax = F_max - Ffric(ds)
    ddsmax = (par.F_max - Friction(xak[3]))/par.mges

    rcg[1] = ddsmax - uk[1]
    rcg[2] = uk[1] + ddsmax

    return rcg
end
function ∂GeneralConstraint(xak,uk,par,dim)
    row = [1,1,2,2]
    col = [3,6,3,6]

    ∂xa3_ddsmax = -∂Friction(xak[3])/par.mges
    val = [∂xa3_ddsmax,-1. , ∂xa3_ddsmax, 1.]

    Jcg = sparse(row,col,val,dim.ncon,dim.nxa+dim.nu)

    return Jcg
end

############### WARM START 
##########################
"""
    # WSxua(x0,u0,dim::DIM,prob::DMSprob)
    initialize a TimeSeries for the input trajectory u and
    the augmented state trajecotry x aon the shooting grid specified 
    via the the structure dim::DIM.
    It is assumed that the trajectory u is constant over the
    entire time horizon. 

    returns the TimeSeries xs,us.
"""
function WSxua(x0,u0,dim::DIM,prob::DMSprob)
    us = WSu(u0,dim)
    xs = TimeSeries(undef,dim.N+1)
    xs[1] = deepcopy(x0);
    for k in 1:dim.N
        xs[k+1] =  xIVP(xs[k],us[k],prob)
    end

    return xs,us
end
