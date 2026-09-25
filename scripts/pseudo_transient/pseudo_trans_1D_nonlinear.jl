using CairoMakie
using Printf
using LinearAlgebra
using Enzyme

@views function smoothmax(h)
    a = h[1:end-1]
    b = h[2:end]
    m = a .+ b .+ sqrt.((a .- b).^2 .+ eps())
    return m
end

@views function compute_phi!(φ, h, ρⁱg, ρʷg, H, B)
    @. φ  = ρⁱg * H + ρʷg * (B + h)
    return
end

@views function compute_flux!(h, q, φ, ∇φ, A, a, k, alpha, betha, dx)
    @. ∇φ = (φ[2:end] - φ[1:end-1]) / dx
    @. A = (∇φ^2 + eps())^(betha/2 - 1) * ∇φ
    # with smoothmax to have differentiable function
    @. a = alpha * $smoothmax(h)^(alpha - 1)

    @. q[2:end-1] = -k * 0.5 * (h[1:end-1]^alpha + h[2:end]^alpha) * A - k * a * 0.5 * abs(A) * (h[2:end] - h[1:end-1])
    return
end
    
@views function compute_update!(h, q, h_old, dt, dτ, dx)
    @. h[2:end-1] = (h[2:end-1] + dτ * (h_old[2:end-1] / dt - (q[3:end-1] - q[2:end-2]) / dx)) / (1 + dτ/dt)
    return
end

function update_h!(h, h_old, φ, ∇φ, q, A, a, k, ρⁱg, ρʷg, H, B, alpha, betha, dt, dτ, dx)

    compute_phi!(φ, h, ρⁱg, ρʷg, H, B)
    compute_flux!(h, q, φ, ∇φ, A, a, k, alpha, betha, dx)
    compute_update!(h, q, h_old, dt, dτ, dx)

    @. r = -(h[2:end-1] - h_old[2:end-1]) / dt - $diff(q[2:end-1]) / dx
    return
end

function residual!(h, h_old, φ, ∇φ, q, A, a, r, k, ρⁱg, ρʷg, H, B, alpha, betha, dt, dx)
    compute_phi!(φ, h, ρⁱg, ρʷg, H, B)
    compute_flux!(h, q, φ, ∇φ, A, a, k, alpha, betha, dx)
    
    @. r = -(h[2:end-1] - h_old[2:end-1]) / dt - $diff(q[2:end-1]) / dx    
    return
end

@views function compute_entropy_dissipation(h, ∇φ, k, alpha, betha, dx)
    # h evaluated at the midpoint of each interval
    h_inter = 0.5 .* (h[1:end-1] .+ h[2:end])
    # approch integral with midpoint rule
    I = k * dx * sum(h_inter.^alpha .* abs.(∇φ).^betha)

    return I
end

function pseudo_1D_lin()
    # physics
    lx  = 100e3
    k   = 0.001
    ρⁱg = 910.0 * 9.81
    ρʷg = 1000.0 * 9.81
    alpha = 5/4
    betha = 3/2
    # numerics
    nx   = 100
    nvis = 100 # 1000
    nvistot = 100
    tol  = 1e-12
    tol_change = 1e-14
    maxiter = 1e5
    t_end = 1e5    # total simulation time
    dt = 40 #20 #11.5 for the 1D script
    nt = Int(ceil(t_end/dt))
    epsi = 1e-2
    D_eff = 1

    # preprocessing
    dx = lx / (nx - 1)
    xn = LinRange(0, lx, nx)
    

    # arrays
    H   = zeros(nx)
    B   = zeros(nx)
    h   = zeros(nx)
    h_old = zeros(nx)
    φ   = zeros(nx)
    ∇φ  = zeros(nx - 1)
    q   = zeros(nx + 1)
    r = zeros(nx - 2)
    A   = zeros(nx - 1)
    a   = zeros(nx - 1)

    I_history = zeros(nt)
    t_history = zeros(nt)

    # arrays for autodiff
    h_k = zeros(nx)
    h̄ = zeros(nx)
    r̄ = zeros(nx - 2)
    φ_dev = zeros(nx)
    ∇φ_dev = zeros(nx-1)
    q_dev = zeros(nx+1)
    A_dev = zeros(nx - 1)
    a_dev = zeros(nx - 1)
    b = zeros(nx - 2)

    # initialisation
    # H - ice thickness
    @. H = 4000.0 - xn / 1e2
    @. H[xn>9lx/10] = 0
    # B - bed elevation
    @. B = 1.4e3 + 0.2e3 * sin(6π * xn / lx) - xn / 1e2

    dτ = 1

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
        rel_change = 2 * tol

        # pseudo-transient time loop
        while (err > tol || rel_change > tol_change) && iter < maxiter            
            D_eff = k .* h[1:end-1].^alpha .* ρʷg .* (∇φ.^2 .+ eps()).^((betha - 2)/2)
            dτ = dx^2 / 2.1 / max(maximum(D_eff), epsi) / 2    
        
            h_k .= h
            h̄ .= h
            r̄ .= 0
            φ_dev .= 0 
            ∇φ_dev .= 0
            q_dev .= 0
            A_dev .= 0
            a_dev .= 0
            update_h!(h, h_old, φ, ∇φ, q, A, a, r, k, ρⁱg, ρʷg, H, B, alpha, betha, dt, dτ, dx)
            h[end] = 4.2e3
            iter += 1
            
            if iter % 10 == 0
                Enzyme.autodiff(set_runtime_activity(Enzyme.Forward), update_h!, Const, Duplicated(h, h̄), Const(h_old), Duplicated(φ, φ_dev), Duplicated(∇φ, ∇φ_dev), Duplicated(q, q_dev), Duplicated(A, A_dev), Duplicated(a, a_dev), Duplicated(r, r̄), Const(k), Const(ρⁱg), Const(ρʷg), Const(H), Const(B), Const(alpha), Const(betha), Const(dt), Const(dτ), Const(dx))
                @. b = r - r̄
                # err = norm(r) / norm(b)
                err = norm(r, Inf) / norm(b, Inf)
                rel_change = norm(h_k[2:end-1] - h[2:end-1]) / norm(h[2:end-1])
            end
        end
        # update phi with cnverged h
        compute_phi!(φ, h, ρⁱg, ρʷg, H, B)
        @. ∇φ = (φ[2:end] - φ[1:end-1]) / dx

        # cmpute entrpy diddipatin
        I = compute_entropy_dissipation(h, ∇φ, k, alpha, betha, dx)


        if it % nvistot == 0
            println("t = ", t, ", physical step = ", it, ", pseudo iterations = ", iter)
            println("dτ = ",dτ)
            println("error at final pseudo transient step = ", err)
            println("rel change in h_k = ", rel_change)
            println("Entropy dissipation I = ", I)

            I_history[it + 1] = I
            t_history[it + 1] = t


            plt[2][3] = B .+ h
            plt[3][2] = B .+ h
            plt[3][3] = B .+ h .+ H
            plt[4][2] = φ ./ 1e5
            plt[5][2] = ∇φ ./ 1e2
            display(fig)

            h_less = h[h .< -eps()]
            println("H below 0: ", h_less)
            
        end
        ittot += iter
        it += 1
        t += dt
        # update h
        h_old .= h
        # if isnan(err) error("NaN") end
    end

    @printf("Total time = %1.2f, time steps = %d, nx = %d, iterations tot = %d \n", round(t_end, sigdigits=2), it, nx, ittot)   
    println(D_eff)
    # figure
    plt[2][3] = B .+ h
    plt[3][2] = B .+ h
    plt[3][3] = B .+ h .+ H
    plt[4][2] = φ ./ 1e5
    plt[5][2] = ∇φ ./ 1e2
    display(fig)

    fig_I = Figure(; size=(700, 400))

    ax_I = Axis(
        fig_I[1, 1];
        xlabel = "time",
        ylabel = "entropy dissipation I(t)"
    )

    lines!(
        ax_I,
        t_history,
        I_history;
        linewidth = 2
    )

    display(fig_I)

    return xn, h

    
end

x, h = pseudo_1D_lin()

println("Press Enter to quit")
readline()
