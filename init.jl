import DifferentialEquations
 const ODE = DifferentialEquations
using Plots
using SparseArrays
using ForwardDiff
using Sundials
using Ipopt
using JuMP
plotly()

##### Include Models
include("src/Pendulum.jl")
par = ParamInit();  # initialize parameter structure

##### Include Shooting Grid Functions
include("src/ShootingGrid.jl")

##### Funtions for Direct Multiple Shooting
include("src/DMSfunctions.jl")

##### Include Plot Functions
include("src/PlotRecipe.jl")

