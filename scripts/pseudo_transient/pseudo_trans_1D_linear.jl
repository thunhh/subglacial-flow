using CairoMakie
using Printf
using LinearAlgebra

function compute_flux!(h, q, φ, ∇φ, k, dx)
    @. ∇φ = (φ[2:end] - φ[1:end-1]) / dx
    @. q[2:end-1] = -k * 0.5 * (h[1:end-1] + h[2:end]) * ∇φ - k * 0.5 * abs(∇φ) * (h[2:end] - h[1:end-1])
    return
end
    
function compute_update!(h, q, h_old, dt, dτ, dτ_ρ, dx)
    # h = (h + dτ * (h_old / dt - (q[2:end] - q[1:end-1]) / dx)) / (1 + dτ/dt)
    @. h = (h + dτ_ρ * (h_old / dt - (q[2:end] - q[1:end-1]) / dx)) / (1 + dτ_ρ/dt)

    return
end

function check_res!(Resh, h, h_old, q, dt, dx)
    @. Resh = - (h - h_old) / dt - (q[2:end] - q[1:end-1]) / dx  
    return
end  

function pseudo_1D_lin()
    # physics
    lx  = 100e3
    k   = 0.001
    ρⁱg = 910.0 * 9.81
    ρʷg = 1000.0 * 9.81
    # numerics
    nx   = 100
    nvis = 20000 # 1000
    nvistot = 100
    tol  = 1e-8
    maxiter = 1e5
    t_end = 1e5 #1e6 #1.0     # total simulation time
    dt = 20 #12 for the 1D script
    nt = Int(ceil(t_end/dt))
    D = 1

    # preprocessing
    dx = lx / (nx - 1)
    xn = LinRange(0, lx, nx)

    # provisorisches dτ
    CFL    = 0.99       # CFL number
    # Derived numerics
    dx     = lx / nx      # grid size
    Vpdτ   = CFL * dx
    Re     = π + sqrt(π^2 + (lx^2 / D / dt)) # Numerical Reynolds number
    dτ_ρ  = lx / Vpdτ / Re
    println("dτ_ρ = ", dτ_ρ)

    dτ = 1

    dτ_ρ = 0.01
    

    # arrays
    H   = zeros(nx)
    B   = zeros(nx)
    h   = zeros(nx)
    h_old = zeros(nx)
    φ   = zeros(nx)
    ∇φ  = zeros(nx - 1)
    q   = zeros(nx + 1)
    σnn = zeros(nx)
    Resh = zeros(nx)
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

    t = 0.0
    it = 0
    ittot = 0

    # physical time loop
    while it < nt
        iter = 0
        err = 2 * tol

        # pseudo-transient time loop
        while err > tol && iter < maxiter
            @. φ  = σnn + ρʷg * (B + h)
            compute_flux!(h, q, φ, ∇φ, k, dx)
            compute_update!(h, q, h_old, dt, dτ, dτ_ρ, dx)
            h[end] = 4.2e3
            iter += 1

            if iter % nvis == 0
                check_res!(Resh, h, h_old, q, dt, dx)
                err = norm(Resh) / sqrt(length(Resh)) # still need to understand this criteria
            end
        end

        if it % nvistot == 0
            println("t = ", t, ", physical step = ", it, ", pseudo iterations = ", iter)
            plt[2][3] = B .+ h
            plt[3][2] = B .+ h
            plt[3][3] = B .+ h .+ H
            plt[4][2] = φ ./ 1e5
            plt[5][2] = ∇φ ./ 1e2
            display(fig)
        end
        ittot += iter
        it += 1
        t += dt
        # update h
        h_old .= h
        if isnan(err) error("NaN") end
    end

    @printf("Total time = %1.2f, time steps = %d, nx = %d, iterations tot = %d \n", round(t_end, sigdigits=2), it, nx, ittot)   

    # figure
    plt[2][3] = B .+ h
    plt[3][2] = B .+ h
    plt[3][3] = B .+ h .+ H
    plt[4][2] = φ ./ 1e5
    plt[5][2] = ∇φ ./ 1e2
    display(fig)

    return xn, h
end

x, h = pseudo_1D_lin()

println("Press Enter to quit")
readline()
