### DYNAMIC MODELLING OF THE PENDULUM
# https://github.tik.uni-stuttgart.de/ac121730/NOOCjulia/wiki/03-Pendulum

######## PROBLEM PARMAETER
##########################
mutable struct param
    g::AbstractFloat   #gravitational constnat [m/s^2]
    l::AbstractFloat   #pendulum length [m]
    mp::AbstractFloat  #pendulum mass [kg]
    mw::AbstractFloat  #carriage mass [kg]
    mges::AbstractFloat # overall mass [kg]
    Jp::AbstractFloat  #inertia penedlum [kg.m^2]
    J::AbstractFloat   #overall inertia [kg.m^2]
    s_max::AbstractFloat  # maximum platform vertical travel distance [m]
    s_min::AbstractFloat  # miminimum platform vertical travel distance [m]
    ds_max::AbstractFloat  # maximum carriage speed [m/s]
    ds_min::AbstractFloat  # miminimum carriage speed [m/s]
    F_max::AbstractFloat  # maximum linear axis force [N]
    F_min::AbstractFloat  # minimum linear axis force [N]
end

"""
    # ParamInit()
    Initialize pendelum parameter structure:

    g ...       gravitational constnat [m/s^2]
    l ...       pendulum length [m]
    mp ...      pendulum mass [kg]
    mw ...      carriage mass [kg]
    mges...     mp + mw
    Jp...       inertia penedlum [kg.m^2]
    J...        Jp + mp*l^2 overall inertia [kg.m^2]
    s_max...    maximum platform vertical travel distance [m]
    s_min...    miminimum platform vertical travel distance [m]
    ds_max...   maximum carriage speed [m/s]
    ds_min...   miminimum carriage speed [m/s]
    F_max...    maximum linear axis force [N]
    F_min...    minimum linear axis force [N]

"""
function ParamInit()::param
    g = 9.81   #gravitational constnat [m/s^2]
    l = 0.25   #pendulum length [m]
    mp = 0.31  #pendulum mass [kg]
    mw = 2.99  #carriage mass [kg]
    mges = mp + mw
    Jp = mp*l^2/3 #inertia penedlum [kg.m^2]
    J  = Jp + mp*l^2 #overall inertia [kg.m^2]
    s_max = 0.35 # maximum platform vertical travel distance [m]
    s_min = -0.35 # miminimum platform vertical travel distance [m]
    ds_max = 4. # maximum carriage speed [m/s]
    ds_min = -4 # miminimum carriage speed [m/s]
    F_max = 350. # maximum linear axis force [N]
    F_min = -350. # minimum linear axis force [N]

    return param(g,l,mp,mw,mges,Jp,J,s_max,s_min,ds_max,ds_min,F_max,F_min)
end

######## STATE SPACE MODEL
##########################
""" 
    # NlModel!(dx,x,u,t)
    Nonlinear input affine state space model of the pendelum

    dx = f(x) + g(x)*u 

    with state

    x = [s, ds, ϕ, dϕ], 
    s  ...  carrier position,
    ds ...  carrier speed ds/dt,
    ϕ  ...  pendulum angle,
    dϕ ...  pendulum angular speed dϕ/dt,  

    the input u is equal to the carrier acceleration dds.
"""
function NlModel!(dx,x,u,t)
    # state x = [s ds ϕ dϕ]
    # input u = dds

    # input affine nonlinear state space model
    dx[1] = x[2]
    dx[2] = u[1]
    dx[3] = x[4]
    dx[4] = 3/4/par.l*(par.g*sin(x[3]) + cos(x[3])*u[1])

    return dx
end

function Friction(ds)
    α = [0.25, 0.5, 0.01]*5
    β = [100,  1, 100]

    return α[1]*(tanh(β[1]*ds) - tanh(β[2]*ds)) +  α[2]*tanh(β[3]*ds) + α[2]*ds
end
∂Friction(ds) = ForwardDiff.derivative(Friction,ds)
