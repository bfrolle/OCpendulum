######################################## DISCRETIZATION OF THE AUGMENTED SYSTEM
###############################################################################
# introducing src/DMSfunctions.jl

sNp = 0.3;            # desired terminal position sNp = s^[N+1] 
#
t = 0.5               # time horizon
N = 6                 # discretization of shooting grid
NN = 10               # steps inbetween shotting grids
nx = 4                # state dimension
nu = 1                # input dimension 
ncon = 2              # number of gernal inequality constraints
solver = ODE.Tsit5()  # numerical integration solver
# ODE.RK4() , ODE.Tsit5() 

dim  = ProblemDimension(t,N,NN,nx,nu,ncon)
prob = DMSproblem(sNp,solver,dim)

################# WARM START
############################
xa0 = [0. ,0., 0., 0., 0.];   # augmented inital condition
u0 = [0.005]; 
xs,us = WSxua(xa0,u0,dim,prob)
###### plot warm start
plotDMS(xs,us,dim,prob,"Warm Start DMS")


####### adjust weights
prob.wL = 0.5;
prob.wx = [100.,1.,1.,1.]
prob.wu = [1.]
prob.max_iter = 800;

######## INIT JuMP PROBLEM
##########################
model = Model() 

######### DEFINE VARIABLES
##########################
# http://www.juliaopt.org/JuMP.jl/v0.21/variables/#Variable-containers-1
# http://www.juliaopt.org/JuMP.jl/v0.21/variables/#Start-values-1
# http://www.juliaopt.org/JuMP.jl/v0.21/variables/#Variable-bounds-1  
@variable(model, prob.xalb[i] <= xa[k=1:dim.N+1,i=1:dim.nx+1] <= prob.xaub[i], start = deepcopy(xs[k][i]))
@variable(model, u[k=1:dim.N+1],   start = deepcopy(us[k][1]))

@variable(model, 0 <= sn[k=1:dim.N+1]) #positive slack variable
@variable(model, 0 <= sp[k=1:dim.N+1]) #positive slack variable

for k in 1:dim.N+1
       if us[k][1] >= 0
              set_start_value(sp[k],us[k][1]) 
              set_start_value(sn[k],0.) 
       else
              set_start_value(sp[k],0.) 
              set_start_value(sn[k],us[k][1])
       end
end
              
####### REGISTER FUNCTIONS
##########################
# http://www.juliaopt.org/JuMP.jl/v0.21/nlp/#User-defined-Functions-1
ϑa(xa1,xa2,xa3,xa4,xa5,u)   = prob.ϑa(vcat(xa1,xa2,xa3,xa4,xa5),[u],prob)

x_k1(xa1,xa2,xa3,xa4,xa5,u) = prob.xIVP(vcat(xa1,xa2,xa3,xa4,xa5),[u],prob)[1]
x_k2(xa1,xa2,xa3,xa4,xa5,u) = prob.xIVP(vcat(xa1,xa2,xa3,xa4,xa5),[u],prob)[2]
x_k3(xa1,xa2,xa3,xa4,xa5,u) = prob.xIVP(vcat(xa1,xa2,xa3,xa4,xa5),[u],prob)[3]
x_k4(xa1,xa2,xa3,xa4,xa5,u) = prob.xIVP(vcat(xa1,xa2,xa3,xa4,xa5),[u],prob)[4]
x_k5(xa1,xa2,xa3,xa4,xa5,u) = prob.xIVP(vcat(xa1,xa2,xa3,xa4,xa5),[u],prob)[5]
JuMP.register(model, :ϑa,   dim.nxua,   ϑa, autodiff=true)
JuMP.register(model, :x_k1, dim.nxua, x_k1, autodiff=true)
JuMP.register(model, :x_k2, dim.nxua, x_k2, autodiff=true)
JuMP.register(model, :x_k3, dim.nxua, x_k3, autodiff=true)
JuMP.register(model, :x_k4, dim.nxua, x_k4, autodiff=true)
JuMP.register(model, :x_k5, dim.nxua, x_k5, autodiff=true)


c1g(xa3,u) = prob.cg(vcat(0.,0.,xa3,0.,0.),[u],par,dim)[1]
c2g(xa3,u) = prob.cg(vcat(0.,0.,xa3,0.,0.),[u],par,dim)[2]
function ∇c1g!(g,xa3,u)
       g[1] = prob.∂cg(vcat(0.,0.,xa3,0.,0.),[u],par,dim)[1,3]
       g[2] = prob.∂cg(vcat(0.,0.,xa3,0.,0.),[u],par,dim)[1,6]
       return g
end
function ∇c2g!(g,xa3,u)
       g[1] = prob.∂cg(vcat(0.,0.,xa3,0.,0.),[u],par,dim)[2,3]
       g[2] = prob.∂cg(vcat(0.,0.,xa3,0.,0.),[u],par,dim)[2,6]
       return g
end
JuMP.register(model, :c1g,  2, c1g, ∇c1g!)
JuMP.register(model, :c2g,  2, c2g, ∇c2g!)

######## SETUP CONSTRAINTS
##########################
#http://www.juliaopt.org/JuMP.jl/v0.21/constraints/#Constraint-containers-1
##### CONTINUITY
@NLconstraint(model, cont_x1[k = 1:dim.N], x_k1(xa[k,1],xa[k,2],xa[k,3],xa[k,4],xa[k,5],u[k]) == xa[k+1,1])
@NLconstraint(model, cont_x2[k = 1:dim.N], x_k2(xa[k,1],xa[k,2],xa[k,3],xa[k,4],xa[k,5],u[k]) == xa[k+1,2])
@NLconstraint(model, cont_x3[k = 1:dim.N], x_k3(xa[k,1],xa[k,2],xa[k,3],xa[k,4],xa[k,5],u[k]) == xa[k+1,3])
@NLconstraint(model, cont_x4[k = 1:dim.N], x_k4(xa[k,1],xa[k,2],xa[k,3],xa[k,4],xa[k,5],u[k]) == xa[k+1,4])
@NLconstraint(model, cont_x5[k = 1:dim.N], x_k5(xa[k,1],xa[k,2],xa[k,3],xa[k,4],xa[k,5],u[k]) == xa[k+1,5])


##### INITAL CONDITION AND TERMINAL INPUT
@constraint(model, x0[i = 1:dim.nxa], xa[1,i] == xs[1][i])
@constraint(model, uN,  u[dim.N+1] == u[dim.N])

##### GENERAL CONSTRAINT
@NLconstraint(model, gen1[k = 1:dim.N+1], c1g(xa[k,3],u[k]) >= 0)
@NLconstraint(model, gen2[k = 1:dim.N+1], c2g(xa[k,3],u[k]) >= 0)

##### L1 REGULARIZATION
@constraint(model, L1reg[k = 1:dim.N+1], u[k] == sp[k] - sn[k])

######### DEFINE OBJECTIVE
##########################
# http://www.juliaopt.org/JuMP.jl/v0.21/objective/
Np = dim.N+1
@NLobjective(model, Min, ϑa(xa[Np,1],xa[Np,2],xa[Np,3],xa[Np,4],xa[Np,5],u[Np]) + 1e-3*sum(sp[k]+sn[k] for k in 1:Np))
 
############ SOLVE PROBLEM
##########################
# http://www.juliaopt.org/JuMP.jl/v0.21/solvers/
set_optimizer(model, Ipopt.Optimizer)     # select solver
set_optimizer_attributes(model, "print_level" => prob.print_level, "max_iter" => prob.max_iter) # set solver paramter

optimize!(model) # solve problem

xopt = TimeSeries(undef,dim.N+1)
uopt = TimeSeries(undef,dim.N+1)
##### GET SOLUTION
 for k in 1:dim.N+1
        xopt[k] = value.(xa[k,:])
        uopt[k] = [value(u[k])]
 end


###### plot solution
plotDMS(xopt,uopt,dim,prob,"Solution")

