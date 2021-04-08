################# SOLVING THE INITAIAL VALUE PROBLEM FOR THE NONLINEAR PENDULUM 
###############################################################################
# introducing src/Pendulum.jl

x0 = [0., 0., -0.01, 0.]                         # inital condition
tspan = (0,10.)                                  # time span   
u  = [0.001]                                     # constant acceleration input 
############# FIXED STEPSIZE
############################
probNL = ODE.ODEProblem(NlModel!,x0,tspan,u)     # explicit Runge-Kutta 4th order with fixed stepsize
solNL  = ODE.solve(probNL, ODE.Exp4(), dt=0.01)  # solve the IVP with the canonical Runge-Kutta 4th order mehtod  
# PLOT IVP
plot(solNL, layout = (4,1), window_title = "Nonlinear Pendulum")
########## ADAPTIVE STEPSIZE
############################
solNLad  = ODE.solve(probNL, ODE.RK4())   # solve the IVP with the canonical Runge-Kutta 4th order method 
# PLOT IVP
plot!(solNLad, layout = (4,1))
gui()
############ COMPARE SOLVERS
############################
#list of all available solvers see  https://docs.juliadiffeq.org/latest/solvers/ode_solve/

tspan = (0,1e2)
prob = ODE.ODEProblem(NlModel!,x0,tspan,u) 

tRK4ex    = @elapsed ODE.solve(prob, ODE.Exp4(), dt = 0.01)   #explicit Runge-Kutta 4th order with fixed stepsize
tRK4ad    = @elapsed ODE.solve(prob, ODE.RK4())               #canonical Runge-Kutta 4th order mehtod with adaptive stepsize
tode45    = @elapsed ODE.solve(prob, ODE.Tsit5())             #ode45
tode15s   = @elapsed ODE.solve(prob, ODE.QNDF())              #ode15s
tCVODE    = @elapsed ODE.solve(prob, CVODE_BDF())             #sundial solver

println("Cumputational time for a time horizon of 100sec:")
println("RK4 explict:     $(tRK4ex*1e3)msec")
println("RK4 adaptive:    $(tRK4ad*1e3)msec")
println("Tsit5 (ode45):   $(tode45*1e3)msec")
println("QNDF (ode15s):   $(tode15s*1e3)msec")
println("CVODE (sundial): $(tCVODE*1e3)msec")