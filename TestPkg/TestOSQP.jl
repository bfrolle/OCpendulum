using JuMP, OSQP

# initialize optimization model with OSQP
m = Model()
set_optimizer(m, OSQP.Optimizer)


# define optimization variables
x0 = [0, 0] # initial guess
@variable(m, x[k = 1:2], start = x0[k])

# define constriant
A = [1  1; 1 0; 0 1]
lb = [1, 0, 0]
ub = [1, 0.7, 0.7]
@constraint(m, clb[i = 1:3], sum(A[i,j]*x[j] for j in 1:2) <= ub[i])
@constraint(m, cub[i = 1:3], sum(A[i,j]*x[j] for j in 1:2) >= lb[i])

# define objective
P = [4 1; 1 2]
q = [1, 1]
@expression(m, exp[i = 1:2], sum(P[i,j]*x[j] for j in 1:2))
@objective(m, Min, sum(exp[i]*x[i] for i in 1:2) + sum(q[i]*x[i] for i in 1:2))

# solve problem
optimize!(m)

println("The solution of the problem")
println()
println(m)
println()
println("is - minimum cost: $(objective_value(m)) at point x[1]: $(value(x[1])) and x[2]: $(value(x[2]))")