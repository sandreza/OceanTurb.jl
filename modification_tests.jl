# tests to see if the modification to KPP does reasonable things or not
# the first test just compares the new model output to the old model output

using OceanTurb, Plots
include("hook.jl")
st = pwd()
les = OceananigansData(st * "/nocurv.jld2")
# set parameters
𝑪 = [0.11803164331592443, 3.7246545857676954, 0.35191154207167974, 6.225750233165317]
parameters = KPP.Parameters( CSL = 𝑪[1], CNL = 𝑪[2], Cb_T = 𝑪[3], CKE = 𝑪[4])
parameters = KPP.Parameters()
# Build the model with a Backward Euler timestepper
constants = Constants(Float64; α = les.α , β = les.β, ρ₀= les.ρ, cP=les.cᵖ, f=les.f⁰, g=les.g)
model = KPP.Model(N=N, L=les.L, stepper=:BackwardEuler, constants = constants, parameters = parameters)
# Get grid if necessary

zp = collect(model.grid.zc)


# get average of initial condition of LES
T⁰ = avg(les.T⁰, N)
# set equal to initial condition of parameterization
model.solution.T[1:N] = copy(T⁰)
# Set boundary conditions
model.bcs.T.top = FluxBoundaryCondition(les.top_T)
model.bcs.T.bottom = GradientBoundaryCondition(les.bottom_T)
# set aside memory
subsample = 1:10:length(les.t)
time_index = subsample
Nt = length(les.t[time_index])
𝒢 = zeros(N, Nt)

# loop the model
ti = collect(time_index)
for i in 1:Nt
    t = les.t[ti[i]]
    run_until!(model, Δt, t)
    @. 𝒢[:,i] = model.solution.T[1:N]
end

for i in 1:10:Nt
    p1 = plot(𝒢[:, i], zp, legend = false)
    display(p1)
end
