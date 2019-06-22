using OceanTurb, Printf

using OceanTurb.ModularKPP: HoltslagDiffusivity, ROMSMixingDepth

@use_pyplot_utils

usecmbright()

modelsetup = (N=128, L=128, stepper=:BackwardEuler)

#name = "Free convection
#    \\small{with \$ \\overline{w b} |_{z=0} = 10^{-8} \\, \\mathrm{m^2 \\, s^{-3}}\$}"

name = "Windy convection"
#name = "Stable wind, Holtslag vs Large et al" #-driven mixing" #Free convection"
#name = "Neutral wind"

Fb = 1e-8
Fu = -1e-4
N² = 5e-6
Δt = 10minute
Δi = 8hour

        cvmix = ModularKPP.Model(; modelsetup...)

     holtslag = ModularKPP.Model(; modelsetup...,
                                   diffusivity = HoltslagDiffusivity())

         roms = ModularKPP.Model(; modelsetup...,
                                   mixingdepth = ROMSMixingDepth())

holtslag_roms = ModularKPP.Model(; modelsetup...,
                                   diffusivity = HoltslagDiffusivity(),
                                   mixingdepth = ROMSMixingDepth())

# Initial condition and fluxes
dTdz = N² / (cvmix.constants.α * cvmix.constants.g)
T₀(z) = 20 + dTdz*z

#models = (cvmix, holtslag, roms, holtslag_roms)
models = (cvmix, roms)

for model in models
    model.solution.T = T₀

    Fθ = Fb / (model.constants.α * model.constants.g)
    model.bcs.U.top = FluxBoundaryCondition(Fu)

    model.bcs.T.top = FluxBoundaryCondition(Fθ)
    model.bcs.T.bottom = GradientBoundaryCondition(dTdz)

end

fig, axs = subplots()

removespines("top", "right")
xlabel("Temperature \$ \\, {}^\\circ \\mathrm{C} \$")
ylabel(L"z \, \mathrm{(m)}")

for i = 1:5
    if i > 1
        for model in models
            run_until!(model, Δt, (i-1) * Δi)
        end
    end

    @printf("""
        t : %.1f hours

          mixing depths
          =============

               cvmix h : %.2f
            holtslag h : %.2f
                roms h : %.2f
    holtslag + roms  h : %.2f
    \n""", time(cvmix)/hour,
    cvmix.state.h,
    holtslag.state.h,
    roms.state.h,
    holtslag_roms.state.h,
    )

    if i == 1
        vlabel = "CVMix KPP"
        hlabel = "Holtslag \$K\$-profile and CVMix mixing depth model"
        rlabel = "CVMix \$K\$-profile and ROMS mixing depth model"
        mlabel = "Holtslag \$K\$-profile and ROMS mixing depth model"
    else
        vlabel = ""
        hlabel = ""
        rlabel = ""
        mlabel = ""
    end

    if i == 1
        tlabel = text(cvmix.solution.T[end], 0.5,
            @sprintf("\$ t = %.0f \$ hours", time(cvmix)/hour),
            verticalalignment="bottom", horizontalalignment="center", color=defaultcolors[i])
    else
        tlabel = text(maximum(cvmix.solution.T.data)-0.001, -cvmix.state.h,
            @sprintf("\$ t = %.0f \$ hours", time(cvmix)/hour),
            verticalalignment="bottom", horizontalalignment="left", color=defaultcolors[i])
    end

    plot(cvmix.solution.T,          "-",  color=defaultcolors[i], label=vlabel, alpha=0.8, markersize=1.5)
    #plot(holtslag.solution.T,       "-.", color=defaultcolors[i], label=hlabel, alpha=0.8, markersize=1.5)
    plot(roms.solution.T,           "--",  color=defaultcolors[i], label=rlabel, alpha=0.8, markersize=1.5)
    #plot(holtslag_roms.solution.T,  ":", color=defaultcolors[i], label=mlabel, alpha=0.8, markersize=1.5)
end

title(name)
legend(fontsize=10)
gcf()

#savefig("free_convection_intermodel.png", dpi=480)
#savefig("$name.png", dpi=480)
