using CairoMakie
using Printf
using LinearAlgebra

function compute_flux!(h, q_h, q_v, φ, ∇φ_h, ∇φ_v, k, dx, dy)
    @. ∇φ_h = (φ[2:end, :] - φ[1:end-1, :]) / dx
    @. ∇φ_v = (φ[:, 2:end] - φ[:, 1:end-1]) / dy
    # Darcy(-Weisbach) water flux
    @. q_h[2:end-1, :] = -k * 0.5 * (h[1:end-1, :] + h[2:end, :]) * ∇φ_h - k * 0.5 * abs(∇φ_h) * (h[2:end, :] - h[1:end-1, :])
    @. q_v[:, 2:end-1] = -k * 0.5 * (h[:, 1:end-1] + h[:, 2:end]) * ∇φ_v - k * 0.5 * abs(∇φ_v) * (h[:, 2:end] - h[:, 1:end-1])
    return
end
    
function compute_update!(h, q_h, q_v, h_old, dt, dτ, dτ_ρ, dx, dy)
    # h = (h + dτ * (h_old / dt - (q[2:end] - q[1:end-1]) / dx)) / (1 + dτ/dt)
    @. h = (h + dτ_ρ * (h_old / dt - ((q_h[2:end, :] - q_h[1:end-1, :]) / dx + (q_v[:, 2:end] - q_v[:, 1:end-1]) / dy))) / (1 + dτ_ρ/dt)

    return
end

function check_res!(Resh, h, h_old, q_h, q_v, dt, dx, dy)
    @. Resh = - (h - h_old) / dt - ((q_h[2:end, :] - q_h[1:end-1, :]) / dx + (q_v[:, 2:end] - q_v[:, 1:end-1]) / dy)
    return
end  

function pseudo_2D_lin()
    # physics
    lx  = 100e3
    ly  = 100e3
    k   = 0.001
    ρⁱg = 910.0 * 9.81
    ρʷg = 1000.0 * 9.81
    # numerics
    nx   = 100
    ny   = 100
    nvis = 20000 # 1000
    nvistot = 100
    tol  = 1e-8
    maxiter = 1e5
    t_end = 1e5 #1e6 #1.0     # total simulation time
    dt = 20 #12 for the 1D script
    nt = Int(ceil(t_end/dt))
    D = 1   # linear case

    # preprocessing
    dx = lx / (nx - 1)
    dy = ly / (ny - 1)
    xn = LinRange(0, lx, nx)
    yn = LinRange(0, ly, ny)

    # provisorisches dτ
    CFL    = 0.99       # CFL number
    # Derived numerics
    Vpdτ   = CFL * min(dx, dy)
    Re     = π + sqrt(π^2 + (max(lx, ly)^2 / D / dt)) # Numerical Reynolds number
    dτ_ρ   = Vpdτ * max(lx, ly) / D / Re
    println("dτ_ρ = ", dτ_ρ)

    dτ = 1
    dτ_ρ = 0.01
    

    # arrays
    H   = zeros(nx, ny)
    B   = zeros(nx, ny)
    h   = zeros(nx, ny)
    h_old = zeros(nx, ny)
    φ   = zeros(nx, ny)
    ∇φ_h  = zeros(nx - 1, ny)
    ∇φ_v  = zeros(nx, ny - 1)
    q_h   = zeros(nx + 1, ny)
    q_v   = zeros(nx, ny + 1)
    σnn = zeros(nx, ny)
    Resh = zeros(nx, ny)
    # initialisation
    # H - ice thickness
    H .= 4000.0 .- reshape(xn, :, 1) ./ 1e2
    H[xn .> 9lx/10, :] .= 0
    # B - bed elevation
    B .= reshape(1.4e3 .+ 0.2e3 .* sin.(6π .* xn ./ lx) .- xn ./ 1e2, nx, 1)
    # σnn - overburden pressure (of ice sheet)
    @. σnn = ρⁱg * H

    # figure
    fig = Figure(size = (800, 600))

    ax = Axis3(fig[1, 1],
        xlabel = "x [km]",
        ylabel = "y [km]",
        zlabel = "z [m]"
    )

    X = xn ./ 1e3
    Y = yn ./ 1e3

    p1 = surface!(ax, X, Y, B,
        color=:brown, transparency=true, alpha=0.8)

    p2 = surface!(ax, X, Y, B .+ h,
        color=:blue, transparency=true, alpha=0.8)

    p3 = surface!(ax, X, Y, B .+ h .+ H,
        color=:lightblue, transparency=true, alpha=0.8)

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
            compute_flux!(h, q_h, q_v, φ, ∇φ_h, ∇φ_v, k, dx, dy)
            compute_update!(h, q_h, q_v, h_old, dt, dτ, dτ_ρ, dx, dy)
            h[end, :] .= 4.2e3

            iter += 1

            if iter % nvis == 0
                check_res!(Resh, h, h_old, q_h, q_v, dt, dx, dy)
                err = norm(Resh) / sqrt(length(Resh)) # still need to understand this criteria
            end
        end

        if it % nvistot == 0
            println("t = ", t, ", physical step = ", it, ", pseudo iterations = ", iter)
            p1[3] = B
            p2[3] = B .+ h
            p3[3] = B .+ h .+ H
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

x, h = pseudo_2D_lin()

println("Press Enter to quit")
readline()
