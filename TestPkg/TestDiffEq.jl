using DifferentialEquations, Plots

## Example lorenz attractor
function lorenz!(dx,x,p,t)        #system function with state x and paramtervector p 
 dx[1] = p[1]*(x[2]-x[1])
 dx[2] = x[1]*(p[2]-x[3]) - x[2]
 dx[3] = x[1]*x[2] - p[3]*x[3]
end

x0 = [1.0,0.0,0.0]                  # initial state
tspan = (0.0,100.0)                 # time span
p = [10.0,28.0,8/3]                 # parameter vector 


prob = ODEProblem(lorenz!,x0,tspan,p) # define initial value problem (IVP)
sol = solve(prob)                     # solve the IVP  

## Plot sotuion using the DifferentialEquations plot recipe vars = (t,x1,x2,x3,...,xn) (starting with index 0)
p0 = plot(sol,vars=(1,2,3); window_title = "Phase Diagram", xlabel="x1", ylabel="x2", zlabel="x3", label ="")
gui(p0)

p1 = plot(sol,vars=(0,1), tspan=(0.0,40.0), label = "x1(t)")
p2 = plot(sol,vars=(0,2), tspan=(0.0,40.0),  label = "x2(t)")
p3 = plot(sol,vars=(0,3), tspan=(0.0,40.0), label = "x3(t)")
pt = plot(p1,p2,p3; layout=(3,1), window_title = "States")
gui(pt)