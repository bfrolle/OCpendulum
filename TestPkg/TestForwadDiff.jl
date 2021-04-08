using Plots
plotly()
using ForwardDiff

function f0(x)                      #define a real valued function f0:R^n -> R
 val = sin(x[1]) + cos(x[2])
 return val
end

x0 = [pi/2, 0]
println("Evaluation point is x0 = $(x0)")

println("The value of f0 at x0 is $(f0(x0) )")

df0(x) = ForwardDiff.gradient(f0,x) #define the gradient function df0:R^n -> R^n 
println("The gradient of f0 at x0 is equal to $(df0(x0) )")

Hf0(x) = ForwardDiff.hessian(f0,x) # define the Hessian function of f0 df0:R^n -> R^(n x n) 
println("The hessian of f0 at x0 is $(Hf0(x0))")