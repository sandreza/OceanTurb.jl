# tests to see if the modification to KPP does reasonable things or not
# the first test just compares the new model output to the old model output
include("hook.jl")
using OceanTurb, Plots
st = pwd()
les = OceananigansData(st * "/nocurv.jld2")

# set parameters
# parameters = KPP.Parameters( CSL = 𝑪[1], CNL = 𝑪[2], Cb_T = 𝑪[3], CKE = 𝑪[4])
parameters = KPP.Parameters()
# Build the model with a Backward Euler timestepper
N = 16
L = les.L
model = KPP.Model(N = N, L = L, stepper=:BackwardEuler)
# Get grid if necessary

# get average of initial condition of LES
T⁰ = avg(les.T[:,1], N)
# set equal to initial condition of parameterization
model.solution.T[1:N] = copy(T⁰)
# Set boundary conditions
model.bcs.T.top = FluxBoundaryCondition(les.top_T)
model.bcs.T.bottom = GradientBoundaryCondition(les.bottom_T)
# set aside memory
Nt = length(les.t)
𝒢 = zeros(N, Nt)

# loop the model
Δt = 60 * 10
for i in 1:Nt
    t = les.t[i]
    run_until!(model, Δt, t)
    @. 𝒢[:,i] = model.solution.T[1:N]
end

for i in 1:10:Nt
    p1 = plot(𝒢[:,i], model.grid.zc, xlims = (19,20), legend = false)
    display(p1)
end
