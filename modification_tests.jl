# tests to see if the modification to KPP does reasonable things or not
# the first test just compares the new model output to the old model output

using OceanTurb, Plots

# set parameters
parameters = KPP.Parameters( CSL = 𝑪[1], CNL = 𝑪[2], Cb_T = 𝑪[3], CKE = 𝑪[4])
# Build the model with a Backward Euler timestepper
constants = Constants(Float64; α = les.α , β = les.β, ρ₀= les.ρ, cP=les.cᵖ, f=les.f⁰, g=les.g)
model = KPP.Model(N=N, L=les.L, stepper=:BackwardEuler, constants = constants, parameters = parameters)
# Get grid if necessary
if grid != 1
    zp = collect(model.grid.zc)
    @. grid  = zp
end
# get average of initial condition of LES
T⁰ = avg(les.T⁰, N)
# set equal to initial condition of parameterization
model.solution.T[1:N] = copy(T⁰)
# Set boundary conditions
model.bcs.T.top = FluxBoundaryCondition(les.top_T)
model.bcs.T.bottom = GradientBoundaryCondition(les.bottom_T)
# set aside memory
if subsample != 1
    time_index = subsample
else
    time_index = 1:length(les.t)
end
Nt = length(les.t[time_index])
𝒢 = zeros(N, Nt)

# loop the model
ti = collect(time_index)
for i in 1:Nt
    t = les.t[ti[i]]
    run_until!(model, Δt, t)
    @. 𝒢[:,i] = model.solution.T[1:N]
end
