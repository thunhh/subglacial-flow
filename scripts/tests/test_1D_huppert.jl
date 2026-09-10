using CairoMakie
using SpecialFunctions


# flow profile
function profile(x, t, Q)
    # rho_a   = 1.204     # density air
    # rho_w   = 1000     # density water
    # g       = 9.81
    # g_      = (rho_w - rho_a) / g
    # v       = 1.004e-6  # viscosity for 20° C water
    # # Q       = 0.6 #1.2 * pi * 5^2        # water volume
  
    v   = 1
    g_  = 1
    # Q   = 1

    Q = 1.5222758135166778
    eta_n   = (1/5 * (3/10)^(1/3) * pi^(1/2) * gamma(1/3) * gamma(5/6))^(-3/5)
    t_0 = 0.314269
    T   = t_0 + t

    # L = eta_n * (g_ * Q^3 * T / 3 / v)^(1/5)
    # println("front = ", L)
    # L = 1.0
    # qq = ((L/eta_n)^5 * 3 * v / g_ / T)^(1/3)
    # println("initial volume = ",qq)

    # arrays
    nx = length(x)
    eta = zeros(nx)
    y   = zeros(nx)
    phi = zeros(nx)
    h   = zeros(nx)

    # similarity solution
    @. eta  = (1/3 * g_ * Q^3 / v)^(-1/5) * x * T^(-1/5)
    # println(eta_n)
    # println(maximum(eta))
    # println(minimum(eta))
    # println(minimum(abs.(eta)))
    @. y    = min(abs(eta/eta_n), 1.0)
    # println("y =", y)
    # @. y    = max(y, 1.0)

    # println(minimum(y))
    # println(maximum(y))
    @. phi  = (3/10)^(1/3) * (1 - y^2)^(1/3)
    @. h    = eta_n^(2/3) * (3*Q^2*v/g_)^(1/5) * T^(-1/5) * phi

    return h
end

# main function
@views function main()
    # physics
    lx = 10.0 # domain length (diameter for axisymmetric)
    Q_0 = 2.0  # initial volume
    t_e = 10.0 # total time of the simulation
    ρʷg = 1.0 #1000.0 * 9.81
    alpha = 3
    c = 1
    v = 1
    g_ = 1
    g = 1

    # numerics
    nvis = 500
    # preprocessing
    nx = 200
    ny = nx
    ly = lx
    dx = lx / nx
    dy = dx
    xv = LinRange(-lx/2, lx/2, nx+1)
    yv = LinRange(-ly/2, ly/2, ny+1)
    xc = 0.5 .* (xv[1:(end-1)] .+ xv[2:end])
    yc = 0.5 .* (yv[1:(end-1)] .+ yv[2:end])
    epsi = eps()    # regularizer

    φ   = zeros(nx, ny)
    ∇φ_h  = zeros(nx - 1, ny)
    ∇φ_v  = zeros(nx, ny - 1)
    q_h   = zeros(nx, ny + 1)
    q_v   = zeros(nx + 1, ny)
    σnn = zeros(nx, ny)

    epsi = 1e-5
    
    A_h   = zeros(nx, ny - 1)
    A_v   = zeros(nx - 1, ny)
    cdiff_y = zeros(nx - 1, ny)
    cdiff_x = zeros(nx, ny - 1)
    a_h   = zeros(nx, ny - 1)
    a_v   = zeros(nx - 1, ny)
    σnn = zeros(nx, ny)

    h  = zeros(nx, ny)   # flow thickness
    hᵉ = zeros(nx)   # exact (asymptotic) profile

    # initial conditions
    hᵉ = profile(xc, 0.0, Q_0)
    h = repeat(hᵉ, 1, ny)

    # for now bullshit, put front later
    tt = LinRange(0, t_e, nvis)

    Rs  = zeros(nvis)
    Rsᵉ = zeros(nvis)

    # visualisation
    fig = Figure()
    axs = (Axis(fig[1, 1]; title="Flow profile", xlabel="x", ylabel="H"),
           Axis(fig[2, 1]; title="Front position", xlabel="t", ylabel="xᶠ"))
    plt = (lines!(axs[1], xc, h[:, 100]; color=:blue, label="initial"),
           lines!(axs[1], xc, h[:, 100]; color=:red, label="numerical"),
           lines!(axs[1], xc, hᵉ; color=:black, linestyle=:dash, label="exact"),
            lines!(axs[2], tt, Rs; color=:red, label="numerical"),
            lines!(axs[2], tt, Rsᵉ; color=:black, linestyle=:dash, label="exact"))
    axislegend(axs[1])
    axislegend(axs[2]; position=:rb)
    display(fig)
    save("initial.png", fig)

    t_n = 0.0 # current time
    it = 0 # time iteration
    while t_n < t_e
        @. φ  = ρʷg * h + epsi
        @. ∇φ_h = (φ[2:end, :] - φ[1:end-1, :]) / dx
        @. ∇φ_v = (φ[:, 2:end] - φ[:, 1:end-1]) / dy

        # central differences
        @. cdiff_y[:, 2:end-1] = ((φ[2:end, 3:end] - φ[2:end, 1:end-2])/dy + (φ[1:end-1, 3:end] - φ[1:end-1, 1:end-2])/dy)/4
        @. cdiff_x[2:end-1, :] = ((φ[3:end, 2:end] - φ[1:end-2, 2:end])/dx + (φ[3:end, 1:end-1] - φ[1:end-2, 1:end-1])/dx)/4

        # defines edges as central differnces require ghost cells (phi is set to zero)
        @. cdiff_y[:, 1] = ((φ[2:end, 2])/dy + (φ[1:end-1, 2])/dy)/4
        @. cdiff_y[:, end] = - ((φ[2:end, end-1])/dy + ( φ[1:end-1, end-1])/dy)/4
        @. cdiff_x[1, :] = ((φ[2, 2:end])/dx + (φ[2, 1:end-1])/dx)/4
        @. cdiff_x[end, :] = - ((φ[end-1, 2:end])/dx + (φ[end-1, 1:end-1])/dx)/4

        @. A_v = (φ[2:end, :] - φ[1:end-1, :]) / dx
        @. A_h = (φ[:, 2:end] - φ[:, 1:end-1]) / dy

        # f is monotone function -> max of derivative of f is derivative of f at max (h) on intervall
        @. a_v = alpha * max(h[1:end-1, :], h[2:end, :])^(alpha - 1)
        @. a_h = alpha * max(h[:, 1:end-1], h[:, 2:end])^(alpha - 1)

        # Darcy(-Weisbach) water flux
        @. q_v[2:end-1, :] = -c * 0.5 * (h[1:end-1, :]^alpha + h[2:end, :]^alpha) * A_v - c * a_v * 0.5 * abs(A_v) * (h[2:end, :] - h[1:end-1, :])
        @. q_h[:, 2:end-1] = -c * 0.5 * (h[:, 1:end-1]^alpha + h[:, 2:end]^alpha) * A_h - c * a_h * 0.5 * abs(A_h) * (h[:, 2:end] - h[:, 1:end-1])

        # diffusive time step
        hmax = maximum(h)
        dt = 3 * v * dx^2 * dy^2 / (2 * g * hmax^alpha * (dx^2 + dy^2)) / 10

        # update water sheet thickness using explicit euler scheme
        @. h -= dt * ((q_v[2:end, :] - q_v[1:end-1, :]) / dx + (q_h[:, 2:end] - q_h[:, 1:end-1]) / dy)
                
        t_n += dt
        it  += 1

        if it % nvis == 0
            # exact profile (similarity solution)
            hᵉ = profile(xc, t_n, Q_0)

            # update plot
            plt[2][2] = h[:, 100]
            plt[3][2] = hᵉ
            plt[4][1] = Rs
            plt[5][1] = Rsᵉ
            display(fig)
        end
    end
    return
end

main()

