module TKE

using OceanTurb

import ..OceanTurb: oncell, onface
import .KPP: ∂B∂z, isunstable, ωτ, ωb
import .ModularKPP: AbstractModularKPPModel

const nsol = 5
@solution U V T S e

minuszero(args...) = -0

Base.@kwdef struct TKEParameters{T} <: AbstractParameters
        CLz :: T = 0.4    # Dissipation parameter
        CLb :: T = Inf    # Dissipation parameter
         Cτ :: T = 400.0  # Dissipation parameter

        CDe :: T = 2.0    # Dissipation parameter

       CK_U :: T = 0.1    # Diffusivity parameter
       CK_T :: T = 0.1    # Diffusivity parameter
       CK_e :: T = 0.1    # Diffusivity parameter

    Ca_unst :: T = -100.0
    Cb_unst :: T = -0.2

    Ca_stab :: T = 2.7
    Cb_stab :: T = -1.0

        KU₀ :: T = 1e-6   # Interior viscosity for velocity
        KT₀ :: T = 1e-7   # Interior diffusivity for temperature
        KS₀ :: T = 1e-9   # Interior diffusivity for salinity
        Ke₀ :: T = 1e-6   # Interior diffusivity for salinity
end

mutable struct State{T} <: FieldVector{6, T}
    Fu :: T
    Fv :: T
    Fθ :: T
    Fs :: T
    Fb :: T
     h :: T
    function State(T=Float64)
        new{T}(0, 0, 0, 0, 0, 0)
    end
end

"""
    update_state!(model)

Update the top flux conditions and mixing depth for `model`
and store in `model.state`.
"""
function update_state!(m)
    m.state.Fu = getbc(m, m.bcs.U.top)
    m.state.Fv = getbc(m, m.bcs.V.top)
    m.state.Fθ = getbc(m, m.bcs.T.top)
    m.state.Fs = getbc(m, m.bcs.S.top)
    m.state.Fb = m.constants.g * (m.constants.α * m.state.Fθ - m.constants.β * m.state.Fs)
    m.state.h  = ModularKPP.mixing_depth(m)
    return nothing
end

struct Model{TKE, H, TS, G, T, S, BC} <: AbstractModel{TS, G, T}
          clock :: Clock{T}
           grid :: G
    timestepper :: TS
       solution :: S
            bcs :: BC
            tke :: TKE
    mixingdepth :: H
      constants :: Constants{T}
          state :: State{T}
end

function Model(; N=10, L=1.0,
            grid = UniformGrid(N, L),
       constants = Constants(),
             tke = TKEParameters(),
     mixingdepth = ModularKPP.LMDMixingDepth(),
         stepper = :ForwardEuler,
    )

    solution = Solution((CellField(grid) for i=1:nsol)...)

    bcs = (
        U = DefaultBoundaryConditions(eltype(grid)),
        V = DefaultBoundaryConditions(eltype(grid)),
        T = DefaultBoundaryConditions(eltype(grid)),
        S = DefaultBoundaryConditions(eltype(grid)),
        e = DefaultBoundaryConditions(eltype(grid))
    )

    Kϕ = (U=KU, V=KV, T=KT, S=KS, e=Ke)
    Rϕ = (U=RU, V=RV, T=RT, S=RS, e=Re)
    Lϕ = (U=minuszero, V=minuszero, T=minuszero, S=minuszero, e=Le)
    eq = Equation(K=Kϕ, R=Rϕ, L=Lϕ, update=update_state!)
    lhs = OceanTurb.build_lhs(solution)

    timestepper = Timestepper(stepper, eq, solution, lhs)

    return Model(Clock(), grid, timestepper, solution, bcs, tke,
                 mixingdepth, constants, State())
end

# Note: to increase readability, we use 'm' to refer to 'model' in function
# definitions below.
#

@inline oncell(f::Function, m, i) = (f(m, i) + f(m, i+1)) / 2
@inline onface(f::Function, m, i) = (f(m, i) + f(m, i-1)) / 2

@inline function mixing_time(m, i)
    #N = sqrt( max(0, ∂B∂z(m, i)) )
    #τb = m.tke.CLb / N

    @inbounds ũ = sqrt( max(0.0, onface(m.solution.e, i)) )
    τu = -m.tke.CLz * m.grid.zf[i] / ũ
    τu = isnan(τu) ? Inf : τu

    τΔ = Δf(m.grid, i) / ũ

    #τu < m.tke.Cτ && @show τu
    τ = 1 / (1/m.tke.Cτ + 1/τu)

    #@show τ

    return max(τ, τΔ)
    #return m.tke.Cτ
    #return min(m.tke.Cτ, τu)
    #return @inbounds min(m.tke.Cτ, τu, τb)
end

"Return the turbuent velocity scale associated with convection."
@inline ωb(Fb, h) = abs(h * Fb)^(1/3)
@inline ωb(m::AbstractModel) = ωb(m.state.Fb, m.state.h)

mixing_length(z, Cκ, Ca, Cb, Fb, Fu::T) where T =
    -Cκ * z * max(zero(T), one(T) - Ca * Cκ * Fb * z / abs(Fu)^T(1.5))^Cb

@inline function mixing_length(m, i)
    if m.state.Fb > 0 # unstable
        Ls = mixing_length(m.grid.zf[i], m.tke.CLz, m.tke.Ca_unst, m.tke.Cb_unst, m.state.Fb, m.state.Fu)
    else
        Ls = mixing_length(m.grid.zf[i], m.tke.CLz, m.tke.Ca_stab, m.tke.Cb_stab, m.state.Fb, m.state.Fu)
    end

    τ = m.state.h / ωb(m)
    Le = τ * maxsqrt(m.solution.e, i)

    return 1 / (1/Ls + 1/Le)
end

@inline ∂B∂z(m, i) = ∂B∂z(m.solution.T, m.solution.S, m.constants.g, m.constants.α,
                            m.constants.β, i)

@inline production(m, i) = KU(m, i) * (∂z(m.solution.U, i)^2 + ∂z(m.solution.V, i)^2)

@inline buoyancy_flux(m, i) = - m.constants.g * (
      m.constants.α * KT(m, i) * ∂z(m.solution.T, i)
    - m.constants.β * KS(m, i) * ∂z(m.solution.S, i)
    )

#@inline dissipation(m, i) = @inbounds m.tke.CDe * m.solution.e[i] / mixing_time(m, i)

@inline dissipation(m, i) =
    @inbounds m.tke.CDe * maxsqrt(m.solution.e, i)^3 / mixing_length(m, i)

#
# Equation entry
#

maxsqrt(ϕ, i) = sqrt(max(0.0, ϕ[i]))
@inline K(m, i) = @inbounds mixing_length(m, i) * onface(maxsqrt, m.solution.e, i)

@inline KU(m, i) = m.tke.KU₀ + m.tke.CK_U * K(m, i)
@inline KT(m, i) = m.tke.KT₀ + m.tke.CK_T * K(m, i)
@inline KS(m, i) = m.tke.KS₀ + m.tke.CK_T * K(m, i)
@inline Ke(m, i) = m.tke.Ke₀ + m.tke.CK_e * K(m, i)

const KV = KU

@inline RU(m, i) = @inbounds   m.constants.f * m.solution.V[i]
@inline RV(m, i) = @inbounds - m.constants.f * m.solution.U[i]
@inline RT(m, i) = 0
@inline RS(m, i) = 0

@inline Re(m, i) = oncell(production, m, i) + oncell(buoyancy_flux, m, i) - dissipation(m, i)
@inline Le(m, i) = -0 #@inbounds m.tke.CDe / mixing_time(m, i)

end # module
