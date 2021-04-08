using JuMP, Ipopt

# initialize optimization model with ipopt
m = Model()
set_optimizer(m, Ipopt.Optimizer)
set_optimizer_attributes(m, "print_level" => 0, "max_iter" => 100)

# define optimization variables
x0 = [0, 0] # initial guess

@variable(m, x[k = 1:2], start = x0[k])

# define constraints
@NLconstraint(m, con, (x[1]-1)^2 + (x[2]+1)^3 + exp(-x[1]) <= 1)

# set objective
@NLobjective(m, Min, (x[1]-3)^3 + (x[2]-4)^2)

# solve problem
optimize!(m)

println("The solution of the problem")
println()
println(m)
println()
println("is - minimum cost: $(objective_value(m)) at point x[1]: $(value(x[1])) and x[2]: $(value(x[2]))")