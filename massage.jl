st = pwd()
filename = st * "/" * "ocean_convection_punch_profiles.jld2"
les_data = jldopen(filename, "r")
les_keys = keys(les_data)
timeseries_keys = keys(les_data["timeseries"]["t"])

# hold the entries for easy constructor creation
container = []

# size of arrays
Nz = length(collect(les_data["grid"]["zC"]))
Nt = length(timeseries_keys)


# create arrays
#Initial Conditions
T⁰ = zeros(Nz)
S⁰ = zeros(Nz)
U⁰ = zeros(Nz)
V⁰ = zeros(Nz)
#Timeseries
T = zeros(Nz, Nt)
S = zeros(Nz, Nt)
U = zeros(Nz, Nt)
V = zeros(Nz, Nt)
t = zeros(Nt)

#Second Order Statistics
wT = zeros(Nz, Nt)
wS = zeros(Nz, Nt)
uu = zeros(Nz, Nt)
vv = zeros(Nz, Nt)
ww = zeros(Nz, Nt)

# grab arrays
for j in 1:Nt
    # Fields
    key = timeseries_keys[j]
    @. T[:,j] = les_data["timeseries"]["T"][key][2:(end-1)]
    @. S[:,j] = les_data["timeseries"]["S"][key][2:(end-1)]
    @. U[:,j] = les_data["timeseries"]["u"][key][2:(end-1)]
    @. V[:,j] = les_data["timeseries"]["v"][key][2:(end-1)]
    # Second Order Statistics
    @. wT[:,j] = les_data["timeseries"]["wT"][key][2:(end-1)]
    @. wS[:,j] = les_data["timeseries"]["wS"][key][2:(end-1)]
    @. uu[:,j] = les_data["timeseries"]["uu"][key][2:(end-1)]
    @. vv[:,j] = les_data["timeseries"]["vv"][key][2:(end-1)]
    @. ww[:,j] = les_data["timeseries"]["ww"][key][2:(end-1)]

    t[j] = les_data["timeseries"]["t"][key]
end
# Set initial Conditions
@. T⁰ = T[:,1]
@. S⁰ = S[:,1]
@. U⁰ = U[:,1]
@. V⁰ = V[:,1]

# Push initial conditions current stuff into container
push!(container, T⁰, S⁰, V⁰, U⁰)
# Push fields into container
push!(container, T, S, U, V)
# Push second order statistics into container
push!(container, wT, wS, uu, vv, ww)

# Now grab parameter
ρ = les_data["parameters"]["density"]
α = les_data["buoyancy"]["equation_of_state"]["α"]
β = les_data["buoyancy"]["equation_of_state"]["β"]
cᵖ = les_data["parameters"]["specific_heat_capacity"]
f⁰ = les_data["coriolis"]["f"]
g = les_data["buoyancy"]["gravitational_acceleration"]
L = les_data["grid"]["Lz"]

# Push parameters to container
push!(container, ρ, α, β, cᵖ, f⁰, g, L)

# grab domain data
z = collect(les_data["grid"]["zC"])

# push
push!(container, t, z)

# now grab boundary condition data
top_T = les_data["boundary_conditions"]["top"]["FT"]
top_S = 0.0
top_U = les_data["boundary_conditions"]["top"]["Fu"]
top_V = 0.0
#bottom boundary condition data, see string for type of boundary condtion
bottom_T = les_data["boundary_conditions"]["bottom"]["dTdz"]
bottom_S = 0.0
bottom_U = 0.0
bottom_V = 0.0

# push to container
push!(container, top_T, top_S, top_U, top_V, bottom_T, bottom_S, bottom_U, bottom_V)
