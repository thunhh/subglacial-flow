using CairoMakie
using Printf
using Statistics

function grounding_lines_1d()
    # physics
    lx  = 100e3
    k   = 0.001
    ρⁱg = 910.0 * 9.81
    ρʷg = 1000.0 * 9.81
    alpha = 5/4
    betha = 3/2
    # numerics
    nx   = 100
    nt   = 50000 #100000
    nvis = 5000 # 1000
    # preprocessing
    dx = lx / (nx - 1)
    xn = LinRange(0, lx, nx)
    # arrays
    H   = zeros(nx)
    B   = zeros(nx)
    h   = zeros(nx)
    φ   = zeros(nx)
    ∇φ  = zeros(nx - 1)
    epsi = zeros(nx - 1)
    A   = zeros(nx - 1)
    a   = zeros(nx - 1)
    q   = zeros(nx + 1)
    σnn = zeros(nx)
    # initialisation
    # H - ice thickness
    @. H = 4000.0 - xn / 1e2
    @. H[xn>9lx/10] = 0
    # B - bed elevation
    @. B = 1.4e3 + 0.2e3 * sin(6π * xn / lx) - xn / 1e2
    # σnn - overburden pressure (of ice sheet)
    @. σnn = ρⁱg * H
    # figure
    fig = Figure(; size=(600, 600))
    axs = (Axis(fig[1, 1]; ylabel="z [m]"),
           Axis(fig[2, 1]; ylabel="φ [bar]"),
           Axis(fig[3, 1]; xlabel="x [km]", ylabel="∇φ [bar/km]"))
    plt = (band!(axs[1], xn ./ 1e3, zeros(nx), B; color=:brown),
           band!(axs[1], xn ./ 1e3, B, B .+ h; color=:blue),
           band!(axs[1], xn ./ 1e3, B .+ h, B .+ h .+ H; color=:lightblue),
           lines!(axs[2], xn ./ 1e3, φ ./ 1e5),
           lines!(axs[3], xn[1:end-1] ./ 1e3, ∇φ ./ 1e2))
    display(fig)
    # time loop
    tcur = 0.0
    for it in 1:nt
        # (overburden) hydraulic potential = overburden pressure + elevation potential + water pressure
        @. φ  = σnn + ρʷg * (B + h)
        @. ∇φ = (φ[2:end] - φ[1:end-1]) / dx
        # regularizor more complex
        @. epsi = 1e-5 * mean(abs.((φ[2:end] - φ[1:end-1]) ./ (dx)))
        @. epsi = 1e-5
        @. A = (abs((φ[2:end] - φ[1:end-1]) / dx)^2 + epsi^2 )^(betha/2 - 1) * (φ[2: end] - φ[1:end-1]) / dx
        # @. A = abs(((φ[2:end] - φ[1:end-1]) / dx))^(betha - 2) * (φ[2:end] - φ[1:end-1]) / dx

        @. a = alpha * max(h[1:end-1], h[2:end])^(alpha - 1)
         
        # Darcy(-Weisbach) water flux
        @. q[2:end-1] = -k * 0.5 * (h[1:end-1]^alpha + h[2:end]^alpha) * A - k * a * 0.5 * abs(A) * (h[2:end] - h[1:end-1])
        # advective time step
        dta = dx / k / maximum(abs, ∇φ) / 2.1

        # diffusive time step
        dtd = dx^2 / (k * ρʷg * maximum(h)^alpha * (max(maximum(abs.(∇φ)), 1e-5)^(betha - 2))) / 10 #/ 2.1
        # dtd = dx^2 / (k * ρʷg * maximum(h)^alpha * maximum(abs.(∇φ))^(betha - 2)) / 10 # / 2.1

        dt = min(dta, dtd) 
        # update water sheet thickness using explicit euler scheme
        # y' approx. by (q[2:end] - q[1:end-1]) / dx forward differnces
        @. h -= dt * (q[2:end] - q[1:end-1]) / dx
        h[end] = 4.2e3
        # update plot
        if it % nvis == 0
            @printf(" t = %.1f d, dt [adv] = %1.3e d, dt [dif] = %1.3e d\n", tcur / 3600 / 24, dta / 3600 / 24, dtd / 3600 / 24)

            plt[2][3] = B .+ h
            plt[3][2] = B .+ h
            plt[3][3] = B .+ h .+ H
            plt[4][2] = φ ./ 1e5
            plt[5][2] = ∇φ ./ 1e2
            display(fig)
        end
        tcur += dt
    end
    return
end

time1 = time()
grounding_lines_1d()
time2 = time()
elapsed_time = time2 - time1
print(elapsed_time)

println("Press Enter to quit")
readline()

