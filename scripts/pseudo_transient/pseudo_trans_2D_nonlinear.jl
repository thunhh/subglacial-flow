using CairoMakie
using Printf
using LinearAlgebra
using Enzyme

function smoothmax(a, b)
    return 0.5 * (a + b + sqrt((a - b)^2 + eps()))
end

@views function compute_phi!(φ, h, ρⁱg, ρʷg, H, B)
    @. φ  = ρⁱg * H + ρʷg * (B + h)
    return
end

@views function compute_flux!(h, q_h, q_v, φ, ∇φ_h, ∇φ_v, A_h, A_v, cdiff_x, cdiff_y, a_h, a_v, k, alpha, betha, dx, dy)
    @. ∇φ_h = (φ[2:end, :] - φ[1:end-1, :]) / dx
    @. ∇φ_v = (φ[:, 2:end] - φ[:, 1:end-1]) / dy

    # with gradients at both interfaces, more accurate for non-uniform grids
    @. cdiff_y[:, 2:end-1] = ((φ[2:end, 2:end-1] - φ[2:end, 1:end-2])/dy + (φ[2:end, 3:end] - φ[2:end, 2:end-1])/dy + (φ[1:end-1, 3:end] - φ[1:end-1, 2:end-1])/dy + (φ[1:end-1, 2:end-1] - φ[1:end-1, 1:end-2])/dy)/4
    @. cdiff_x[2:end-1, :] = ((φ[3:end, 2:end] - φ[2:end-1, 2:end])/dx + (φ[2:end-1, 2:end] - φ[1:end-2, 2:end])/dx + (φ[3:end, 1:end-1] - φ[2:end-1, 1:end-1])/dx + (φ[2:end-1, 1:end-1] - φ[1:end-2, 1:end-1])/dx)/4

    # defines edges as central differnces require ghost cells (phi is set to zero)
    @. cdiff_y[:, 1] = ((φ[2:end, 2])/dy + (φ[1:end-1, 2])/dy)/4
    @. cdiff_y[:, end] = - ((φ[2:end, end-1])/dy + ( φ[1:end-1, end-1])/dy)/4
    @. cdiff_x[1, :] = ((φ[2, 2:end])/dx + (φ[2, 1:end-1])/dx)/4
    @. cdiff_x[end, :] = - ((φ[end-1, 2:end])/dx + (φ[end-1, 1:end-1])/dx)/4

    @. A_v = (((φ[2:end, :] - φ[1:end-1, :])/dx)^2 + cdiff_y^2 + eps())^(betha/2 - 1) * (φ[2:end, :] - φ[1:end-1, :]) / dx
    @. A_h = (cdiff_x^2 + ((φ[:, 2:end] - φ[:, 1:end-1])/dy)^2 + eps())^(betha/2 - 1) * (φ[:, 2:end] - φ[:, 1:end-1]) / dy

    # @. a_v = alpha * max(h[1:end-1, :], h[2:end, :])^(alpha - 1)
    # @. a_h = alpha * max(h[:, 1:end-1], h[:, 2:end])^(alpha - 1)

    # use differentiabe function to compute max
    @. a_v = alpha * smoothmax(h[1:end-1, :], h[2:end, :])^(alpha - 1)
    @. a_h = alpha * smoothmax(h[:, 1:end-1], h[:, 2:end])^(alpha - 1)
        
    # Darcy(-Weisbach) water flux
    @. q_v[2:end-1, :] = -k * 0.5 * (h[1:end-1, :]^alpha + h[2:end, :]^alpha) * A_v - k * a_v * 0.5 * abs(A_v) * (h[2:end, :] - h[1:end-1, :])
    @. q_h[:, 2:end-1] = -k * 0.5 * (h[:, 1:end-1]^alpha + h[:, 2:end]^alpha) * A_h - k * a_h * 0.5 * abs(A_h) * (h[:, 2:end] - h[:, 1:end-1])

    return
end
    
@views function compute_update!(h, q_h, q_v, h_old, dt, dτ, dx, dy)
    @. h = (h + dτ * (h_old / dt - ((q_v[2:end, :] - q_v[1:end-1, :]) / dx + (q_h[:, 2:end] - q_h[:, 1:end-1]) / dy))) / (1 + dτ/dt)
    return
end

function update_h!(h, h_old, φ, ∇φ_h, ∇φ_v, q_h, q_v, A_h, A_v, cdiff_x, cdiff_y, a_h, a_v, k, ρⁱg, ρʷg, H, B, alpha, betha, dt, dτ, dx, dy)
    compute_phi!(φ, h, ρⁱg, ρʷg, H, B)
    compute_flux!(h, q_h, q_v, φ, ∇φ_h, ∇φ_v, A_h, A_v, cdiff_x, cdiff_y, a_h, a_v, k, alpha, betha, dx, dy)
    compute_update!(h, q_h, q_v, h_old, dt, dτ, dx, dy)
end

function residual!(r, h, h_old, φ, ∇φ_h, ∇φ_v, q_h, q_v, A_h, A_v, cdiff_x, cdiff_y, a_h, a_v, k, ρⁱg, ρʷg, H, B, alpha, betha, dt, dx, dy)
    compute_phi!(φ, h, ρⁱg, ρʷg, H, B)
    compute_flux!(h, q_h, q_v, φ, ∇φ_h, ∇φ_v, A_h, A_v, cdiff_x, cdiff_y, a_h, a_v, k, alpha, betha, dx, dy)
    ((q_v[2:end, :] - q_v[1:end-1, :]) / dx + (q_h[:, 2:end] - q_h[:, 1:end-1]) / dy)
    
    @. r = -(h[2:end-1, 2:end-1] - h_old[2:end-1, 2:end-1]) / dt - ((q_v[3:end-1, 2:end-1] - q_v[2:end-2, 2:end-1]) / dx + (q_h[2:end-1, 3:end-1] - q_h[2:end-1, 2:end-2]) / dy)
    return
end


function pseudo_2D_lin()
    # physics
    lx  = 100e3
    ly  = 100e3
    k   = 0.001
    ρⁱg = 910.0 * 9.81
    ρʷg = 1000.0 * 9.81
    alpha = 5/4
    betha = 3/2
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
    C_CFL = 1/2

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
    q_h   = zeros(nx, ny + 1)
    q_v   = zeros(nx + 1, ny)
    A_h   = zeros(nx, ny - 1)
    A_v   = zeros(nx - 1, ny)
    cdiff_y = zeros(nx - 1, ny)
    cdiff_x = zeros(nx, ny - 1)
    a_h   = zeros(nx, ny - 1)
    a_v   = zeros(nx - 1, ny)
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
    A_h_dev   = zeros(nx, ny - 1)
    A_v_dev   = zeros(nx - 1, ny)
    cdiff_y_dev = zeros(nx - 1, ny)
    cdiff_x_dev = zeros(nx, ny - 1)
    a_h_dev   = zeros(nx, ny - 1)
    a_v_dev   = zeros(nx - 1, ny)
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
            # D_eff = k .* h[1:end-1, 1:end-1].^alpha .* ρʷg .* (max.(∇φ_v[1:end-1, :], ∇φ_h[:, 1:end-1]).^2 .+ eps()).^((betha - 2)/2)
            # D_eff = k .* h[1:end-1, 1:end-1].^alpha .* (max.(∇φ_v[1:end-1, :], ∇φ_h[:, 1:end-1]).^2 .+ eps()).^((betha - 2)/2)

            # # dτ = C_CFL* min(dx^2, dy^2) / (max(maximum(k * h[2:end-1, 2:end-1]), epsi)) / ρʷg / 4 
            # dτ = C_CFL * min(dx^2, dy^2) / max(maximum(D_eff), epsi)/ ρʷg / 4
            # dτ = min(dx^2, dy^2) / (k * ρʷg * maximum(h)) / 2.1
            dτ = 3
            h_k .= h
            
            update_h!(h, h_old, φ, ∇φ_h, ∇φ_v, q_h, q_v, A_h, A_v, cdiff_x, cdiff_y, a_h, a_v, k, ρⁱg, ρʷg, H, B, alpha, betha, dt, dτ, dx, dy)
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
                A_h_dev   .= 0
                A_v_dev   .= 0
                cdiff_y_dev .= 0
                cdiff_x_dev .= 0
                a_h_dev   .= 0
                a_v_dev   .= 0

                # Enzyme.autodiff(
                #     set_runtime_activity(Enzyme.Forward),
                #     residual!,
                #     Const,
                    
                #     Duplicated(r, r̄),
                #     Duplicated(h, h̄),
                    
                #     Const(h_old),
                    
                #     Duplicated(φ, φ_dev),
                #     Duplicated(∇φ_h, ∇φ_h_dev),
                #     Duplicated(∇φ_v, ∇φ_v_dev),
                #     Duplicated(q_h, q_h_dev),
                #     Duplicated(q_v, q_v_dev),
                #     Duplicated(A_h, A_h_dev),
                #     Duplicated(A_v, A_v_dev),
                #     Duplicated(cdiff_x, cdiff_x_dev),
                #     Duplicated(cdiff_y, cdiff_y_dev),
                #     Duplicated(a_h, a_h_dev),
                #     Duplicated(a_v, a_v_dev),
                 
                #     Const(k),
                #     Const(ρⁱg),
                #     Const(ρʷg),
                #     Const(H),
                #     Const(B),
                #     Const(alpha),
                #     Const(betha),
                #     Const(dt),
                #     Const(dx),
                #     Const(dy)
                # )

                # @. b = r - r̄
                # # err = norm(r) / norm(b)
                # err = norm(r, Inf) / norm(b, Inf)
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
