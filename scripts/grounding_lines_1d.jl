using CairoMakie
using Printf
using Statistics
using LinearAlgebra


function grounding_lines_1d()
    # physics
    lx  = 100e3
    k   = 0.001
    ρⁱg = 910.0 * 9.81
    ρʷg = 1000.0 * 9.81
    # numerics
    nx   = 100
    # nt   = 100000 #100000
    nvis = 1000 # 1000
    min_DT = Inf

    t_end = 1e5 #1e6 #1.0     # total simulation time
    dt = 10 #20 #11.5 for the 1D script
    nt = Int(ceil(t_end/dt))

    # preprocessing
    dx = lx / (nx - 1)
    xn = LinRange(0, lx, nx)
    # arrays
    H   = zeros(nx)
    B   = zeros(nx)
    h   = zeros(nx)
    φ   = zeros(nx)
    ∇φ  = zeros(nx - 1)
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
        # Darcy(-Weisbach) water flux
        # @. q[2:end-1] = -k * 0.5 * (h[1:end-1] + h[2:end]) * ∇φ - k * 0.5 * abs(∇φ) * (h[2:end] - h[1:end-1])
        @. q[2:end-1] = -k * 0.5 * (h[1:end-1] + h[2:end]) * ∇φ - k * 0.5 * abs(∇φ) * (h[2:end] - h[1:end-1])

        # advective time step
        dta = dx / k / maximum(abs, ∇φ) / 2.1
        # diffusive time step
        dtd = dx^2 / (k * ρʷg * maximum(h)) / 2.1
        # dt = min(dta, dtd)
        min_DT = min(min_DT, dt)
        # update water sheet thickness using explicit euler scheme
        # y' approx. by (q[2:end] - q[1:end-1]) / dx forward differnces
        @. h -= dt * (q[2:end] - q[1:end-1]) / dx
        h[end] = 4.2e3
        # update plot
        if it % nvis == 0
            @printf(" t = %.1f d, dt [adv] = %1.3e d, dt [dif] = %1.3e d\n", tcur / 3600 / 24, dta / 3600 / 24, dtd / 3600 / 24)
            println(min_DT)

            plt[2][3] = B .+ h
            plt[3][2] = B .+ h
            plt[3][3] = B .+ h .+ H
            plt[4][2] = φ ./ 1e5
            plt[5][2] = ∇φ ./ 1e2
            display(fig)
        end
        tcur += dt
    end
    

    ## statistics for scaling the error
    h0_rms = norm(h[2:end-1]) / sqrt(length(h[2:end-1]))
    h0_max = maximum(h[2:end-1])
    println("h0_rms = ", h0_rms)
    println("h0_max = ", h0_max)
    phi_x = @views diff(σnn .+ ρʷg .* (B .+ h)) ./ dx

    Gphi_rms = norm(phi_x) / sqrt(length(phi_x))
    h0 = norm(h[2:end-1]) / sqrt(length(h[2:end-1]))

    Lchar = lx   # or another relevant horizontal scale

    Rchar = k * h0 * Gphi_rms / Lchar

    println("Gphi_rms = ", Gphi_rms)
    println("h0 = ", h0)
    println("Rchar = ", Rchar)

    return
end

time1 = time()
grounding_lines_1d()
time2 = time()
elapsed_time = time2 - time1
print(elapsed_time)

println("Press Enter to quit")
readline()

