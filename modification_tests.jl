# tests to see if the modification to KPP does reasonable things or not
# the first test just compares the new model output to the old model output

using OceanTurb, Plots
include("hook.jl")
st = pwd()
les = OceananigansData(st * "/nocurv2.jld2")
curved = true
flexible = true
if flexible
    NN = sqrt(les.α * les.g * les.bottom_T)
    default_𝑪 = [0.005741998337334633, 3.629207116893695, 1.1392751590144323, 0.0, 0.4, NN]
    power = 1.0
    parameters = KPP.Parameters( CSL = default_𝑪[1], CNL = default_𝑪[2], Cb_T = default_𝑪[3], CKE = default_𝑪[4], CKE2 = default_𝑪[5], CKE3 = default_𝑪[6], CKE4 = power)
    println("using flexible parameters")
else
    println("using other parameters")
    # set parameters
    𝑪 = [0.11803164331592443, 3.7246545857676954, 0.35191154207167974, 7.225750233165317]
     NN = sqrt(les.α * les. g * les.bottom_T)
    parameters = KPP.Parameters( CSL = 𝑪[1], CNL = 𝑪[2], Cb_T = 𝑪[3], CKE = 𝑪[4])
end
N = 16
# parameters = KPP.Parameters()
# Build the model with a Backward Euler timestepper
constants = Constants(Float64; α = les.α , β = les.β, ρ₀= les.ρ, cP=les.cᵖ, f=les.f⁰, g=les.g)
model = KPP.Model(N=N, L=les.L, stepper=:BackwardEuler, constants = constants, parameters = parameters)
# Get grid if necessary

zp = collect(model.grid.zc)


# get average of initial condition of LES
T⁰ = avg(les.T⁰, N)
if curved
    T⁰ = avg(T[:,1], N)
end
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
Δt = 60 * 10
ti = collect(time_index)
for i in 1:Nt
    t = les.t[ti[i]]
    run_until!(model, Δt, t)
    @. 𝒢[:,i] = model.solution.T[1:N]
end

for i in 1:10:Nt
    p1 = scatter(𝒢[:, i], zp, legend = false)
    if curved
        plot!(T[:,ti[i]], z)
    else
        plot!(les.T[:,ti[i]], les.z)
    end
    display(p1)
end
