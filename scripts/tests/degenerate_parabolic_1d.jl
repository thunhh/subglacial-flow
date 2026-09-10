# Exact solution for a degenerate parabolic equation in Cartesian and radial coordinates
# This family of solutions describes flow of lava or ice, water flow in connected cavities or porous fluid flow in deformable medium
# In this file, we assume the initial volume of fluid in the center of the domain
# There is no external potential gradient, i.e. the bed is flat
# Unless the simulation time is large enough, the fluid domain support is compact
#
# Equations:
# 
# ∂H/∂t + ∂q/∂x = 0
# q = -Hᵅ|∂H/∂x|ᵝ⁻² ∂H/∂x
#
# | Case                 |  α  |  β  |
# |----------------------|-----|-----|
# | 1. Darcy-Weisbach    | 5/4 | 3/2 |
# | 2. Newtonian viscous |  3  |  2  |
# | 3. Glen's power law  |  5  |  4  |
# 

using CairoMakie
# arithmetic average
@views avx(H) = 0.5 .* (H[1:(end-1)] .+ H[2:end])
# L2 norm
@views nrm(V) = V .* V
# front position
function front(t, α, β, H₀, R₀, d)
    # notation
    γ = α + β - 2
    k = d * γ + β
    # initial time
    t₀ = R₀ ^ β / (k * (β/γ)^(β-1) * H₀^γ)
    # self-similar variables
    T = 1 + t/t₀
    return R₀ * T^inv(k)
end
# flow profile
function profile(x, t, α, β, H₀, R₀, d)
    # notation
    γ = α + β - 2
    l = β / (β - 1)
    a = (β - 1) / γ
    k = d * γ + β
    # initial time
    t₀ = R₀ ^ β / (k * (β/γ) ^ (β-1) * H₀^γ)
    # self-similar variables
    T = 1 + t / t₀
    H = H₀ * T^(-d/k)
    R = R₀ * T^inv(k)
    return H * max(1 - (abs(x) / R)^l, 0)^a
end
# main function
@views function degenerate_parabolic_1d(; nx=200)
    # physics
    lx = 10.0 # domain length (diameter for axisymmetric)
    H₀ = 1.0  # initial thickness
    R₀ = 1.0  # initial radius
    α  = 3 #5/4  # 5/4 for Darcy-Weisbach, 3 for Newtonian viscous flow, 5 for ice flow
    β  = 2 #3/2  # 3/2 for Darcy-Weisbach, 2 Newtonian viscous flow, 4 for ice flow
    d  = 2 #1   # 1 for cross-section, 2 for axisymmetric
    tₑ = 1.0 #10.0 # total time of the simulation
    # numerics
    nvis = 500
    # preprocessing
    dx = lx / nx
    xv = LinRange(-lx/2, lx/2, nx+1)
    xc = 0.5 .* (xv[1:(end-1)] .+ xv[2:end])
    dt = dx^2 / H₀^α / 4.1 # time step is empirical for simplicity
    # arrays
    H  = zeros(nx)   # flow thickness
    ∇H = zeros(nx-1) # surface gradient
    q  = zeros(nx+1) # flux
    Hᵉ = zeros(nx)   # exact (asymptotic) profile
    # initial conditions
    @. Hᵉ = profile(xc, 0.0, α, β, H₀, R₀, d)
    @. H = Hᵉ
    # front position
    Rsᵉ = Point2f[(0.0, front(0.0, α, β, H₀, R₀, d))]
    Rs = copy(Rsᵉ)
    # visualisation
    fig = Figure()
    axs = (Axis(fig[1, 1]; title="Flow profile", xlabel="x", ylabel="H"),
           Axis(fig[2, 1]; title="Front position", xlabel="t", ylabel="xᶠ"))
    plt = (lines!(axs[1], xc, H; color=:blue, label="initial"),
           lines!(axs[1], xc, H; color=:red, label="numerical"),
           lines!(axs[1], xc, Hᵉ; color=:black, linestyle=:dash, label="exact"),
           lines!(axs[2], Rs; color=:red, label="numerical"),
           lines!(axs[2], Rsᵉ; color=:black, linestyle=:dash, label="exact"))
    axislegend(axs[1])
    axislegend(axs[2]; position=:rb)
    display(fig)
    # time loop
    tₙ = 0.0 # current time
    it = 0   # time step
    while tₙ < tₑ
        @. ∇H = (H[2:end] - H[1:(end-1)]) / dx
        if d == 1 # cross-section flow
            @. q[2:(end-1)] = -$avx(H^α) * (nrm(∇H) + eps()) ^ ((β-2)/2) * ∇H
            @. H -= dt * (q[2:end] - q[1:(end-1)]) / dx
        elseif d == 2 # axisymmetric flow
            @. q[2:(end-1)] = -xv[2:(end-1)] * $avx(H^α) * (nrm(∇H) + eps()) ^ ((β-2)/2) * ∇H
            @. H -= dt / xc * (q[2:end] - q[1:(end-1)]) / dx
        else
            error("d must be 1 or 2, got $d")
        end
        tₙ += dt
        it += 1
        if it % nvis == 0
            # exact profile
            @. Hᵉ = profile(xc, tₙ, α, β, H₀, R₀, d)
            Rᵉ = front(tₙ, α, β, H₀, R₀, d)
            push!(Rsᵉ, Point2f(tₙ, Rᵉ))
            # numerical profile
            for ifr in (nx-1):-1:1
                ϵ = 1e-3H₀
                if H[ifr] > ϵ && H[ifr+1] <= ϵ
                    R = xc[ifr]
                    push!(Rs, Point2f(tₙ, R))
                end
            end
            plt[2][2] = H
            plt[3][2] = Hᵉ
            plt[4][1] = Rs
            plt[5][1] = Rsᵉ
            display(fig)
        end
    end
    return
end

degenerate_parabolic_1d()
