using CairoMakie
using Printf
using LinearAlgebra
using Enzyme

@views function compute_phi!(φ, h, ρⁱg, ρʷg, H, B)
    @. φ  = ρⁱg * H + ρʷg * (B + h)
    return
end

@views function compute_flux!(h, q_h, q_v, φ, ∇φ_h, ∇φ_v, k, dx, dy)
    @. ∇φ_h = (φ[2:end, :] - φ[1:end-1, :]) / dx
    @. ∇φ_v = (φ[:, 2:end] - φ[:, 1:end-1]) / dy
    # Darcy(-Weisbach) water flux
    @. q_h[2:end-1, :] = -k * 0.5 * (h[1:end-1, :] + h[2:end, :]) * ∇φ_h - k * 0.5 * abs(∇φ_h) * (h[2:end, :] - h[1:end-1, :])
    @. q_v[:, 2:end-1] = -k * 0.5 * (h[:, 1:end-1] + h[:, 2:end]) * ∇φ_v - k * 0.5 * abs(∇φ_v) * (h[:, 2:end] - h[:, 1:end-1])
    return
end
    
@views function compute_update!(h, q_h, q_v, h_old, dt, dτ, dx, dy)
    @. h = (h + dτ * (h_old / dt - ((q_h[2:end, :] - q_h[1:end-1, :]) / dx + (q_v[:, 2:end] - q_v[:, 1:end-1]) / dy))) / (1 + dτ/dt)
    return
end

function update_h!(h, h_old, φ, ∇φ_h, ∇φ_v, q_h, q_v, k, ρⁱg, ρʷg, H, B, dt, dτ, dx, dy)
    compute_phi!(φ, h, ρⁱg, ρʷg, H, B)
    compute_flux!(h, q_h, q_v, φ, ∇φ_h, ∇φ_v, k, dx, dy)
    compute_update!(h, q_h, q_v, h_old, dt, dτ, dx, dy)
end

function residual!(r, h, h_old, φ, ∇φ_h, ∇φ_v, q_h, q_v, k, ρⁱg, ρʷg, H, B, dt, dx, dy)
    compute_phi!(φ, h, ρⁱg, ρʷg, H, B)
    compute_flux!(h, q_h, q_v, φ, ∇φ_h, ∇φ_v, k, dx, dy)
    
    @. r = -(h[2:end-1, 2:end-1] - h_old[2:end-1, 2:end-1]) / dt - ((q_h[3:end-1, 2:end-1] - q_h[2:end-2, 2:end-1]) / dx + (q_v[2:end-1, 3:end-1] - q_v[2:end-1, 2:end-2]) / dy)
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
    nvistot = 100
    tol  = 1e-12
    tol_change = 1e-14
    maxiter = 1e3 #10 #1e5
    t_end = 60480 #1e5 #1e6 #1.0     # total simulation time
    dt = 12
    dτ = 1
    nt = Int(ceil(t_end/dt))
    epsi = 1e-2

    # preprocessing
    dx = lx / (nx - 1)
    dy = ly / (ny - 1)
    xn = LinRange(0, lx, nx)
    yn = LinRange(0, ly, ny)

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
    r = zeros(nx-2, ny-2)

    # arrays for autodiff
    h_k = zeros(nx,ny)
    h̄ = zeros(nx, ny)
    r̄ = zeros(nx - 2, ny - 2)
    φ_dev = zeros(nx, ny)
    ∇φ_h_dev = zeros(nx - 1, ny)
    ∇φ_v_dev = zeros(nx, ny - 1)
    q_h_dev = zeros(nx + 1, ny)
    q_v_dev = zeros(nx, ny + 1)
    b = zeros(nx - 2, ny - 2)

    # initialisation
    # H - ice thickness
    H .= 4000.0 .- reshape(xn, :, 1) ./ 1e2
    H[xn .> 9lx/10, :] .= 0
    # B - bed elevation
    B .= reshape(1.4e3 .+ 0.2e3 .* sin.(6π .* xn ./ lx) .- xn ./ 1e2, nx, 1)

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
        rel_change = 2 * tol_change

        # pseudo-transient time loop
        # while (err > tol || rel_change > tol_change) && iter < maxiter   
        while rel_change > tol_change && iter < maxiter   
            dτ = min(dx^2, dy^2) / (max(maximum(k * h[2:end-1, 2:end-1]), epsi)) / ρʷg / 4 / 2
            h_k .= h
            
            update_h!(h, h_old, φ, ∇φ_h, ∇φ_v, q_h, q_v, k, ρⁱg, ρʷg, H, B, dt, dτ, dx, dy)
            h[end, :] .= 4.2e3

            iter += 1

            rel_change = norm(h_k[2:end-1,2:end-1] - h[2:end-1,2:end-1]) / norm(h[2:end-1,2:end-1])

            if iter % 5 == 0
                
                h̄ .= h
                r̄ .= 0

                φ_dev .= 0
                ∇φ_h_dev .= 0
                ∇φ_v_dev .= 0
                q_h_dev .= 0
                q_v_dev .= 0

                Enzyme.autodiff(set_runtime_activity(Enzyme.Forward), residual!, Const, Duplicated(r, r̄), Duplicated(h, h̄), Const(h_old), Duplicated(φ, φ_dev), Duplicated(∇φ_h, ∇φ_h_dev), Duplicated(∇φ_v, ∇φ_v_dev), Duplicated(q_h, q_h_dev), Duplicated(q_v, q_v_dev), Const(k), Const(ρⁱg), Const(ρʷg), Const(H), Const(B), Const(dt), Const(dx), Const(dy))
                @. b = r - r̄
                # err = norm(r) / norm(b)
                err = norm(r, Inf) / norm(b, Inf)
            end
        end

        if it % nvistot == 0
            println("t = ", t, ", physical step = ", it, ", pseudo iterations = ", iter)
            println("dτ = ",dτ)
            println("error at final pseudo transient step = ", err)
            println("rel change in h_k = ", rel_change)

            p1[3] = B
            p2[3] = B .+ h
            p3[3] = B .+ h .+ H
            display(fig)

            h_less = h[h .< -eps()]
            println("H below 0: ", h_less)
            # @assert all(h.>= - eps());
        end
        ittot += iter
        it += 1
        t += dt
        # update h
        h_old .= h
        # if isnan(err) error("NaN") end
    end

    @printf("Total time = %1.2f, time steps = %d, nx = %d, iterations tot = %d \n", round(t_end, sigdigits=2), it, nx, ittot)   

    # figure
    p1[3] = B
    p2[3] = B .+ h
    p3[3] = B .+ h .+ H
    display(fig)

    return xn, h
end

x, h = pseudo_2D_lin()

println("Press Enter to quit")
readline()
