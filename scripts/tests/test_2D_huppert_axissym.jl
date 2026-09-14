using CairoMakie
using SpecialFunctions
using Printf

ETA_N   = (1/5 * (3/10)^(1/3) * pi^(1/2) * gamma(1/3) * gamma(5/6))^(-3/5)
V       = 1.004e-6
RHO_A   = 1.204     # density air
RHO_W   = 1000     # density water
G       = 9.81
G_      = (RHO_W - RHO_A) / G
Q_0       = 1.0301548313168503

function front(t, q)
    t_0 = 1e-8
    T = t + t_0
    return ETA_N * (G_ * Q_0^3 * T / 3 / V)^(1/5) 
end

# flow profile
function profile(x, t)
    rho_a   = 1.204     # density air
    rho_w   = 1000     # density water
    g       = 9.81
    g_      = (rho_w - rho_a) / g
    v       = 1.004e-6  # viscosity for 20° C water
    # Q       = 0.6 #1.2 * pi * 5^2        # water volume
  
    # v   = 1
    # g_  = 1
    # Q   = 1

    Q = 1.0301548313168503 #0.2219401304692966 #0.047815551619325485
    eta_n   = (1/5 * (3/10)^(1/3) * pi^(1/2) * gamma(1/3) * gamma(5/6))^(-3/5)
    t_0 = 1e-8 #0.0001
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
    t_e = 10.0 # total time of the simulation

    # water/air variables
    ρⁱg = 910.0 * 9.81
    ρʷg = 1000.0 * 9.81
    alpha = 3
    betha = 2
    g = 9.81
    v = 1e-6
    rho_a   = 1.204     # density air
    rho_w   = 1000     # density water
    g_      = (rho_w - rho_a) / g
    c = g_/v/3

    # simplified problem variables
    # ρʷg = 1.0 #1000.0 * 9.81
    # alpha = 3
    # c = 1
    # v = 1
    # g_ = 1
    # g = 1

    # numerics
    nvis = 2000
    # preprocessing
    nx = 500 #200
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
    hᵉ = profile(xc, 0.0)
    h = repeat(hᵉ, 1, ny)
    H_0 = maximum(hᵉ)


    # initial front position
    Rsᵉ = Point2f[(0.0, front(0.0, Q_0))] 
    Rs = copy(Rsᵉ)  
    Rsee = copy(Rsᵉ)

    # visualisation
    fig = Figure()
    
    axs = (Axis(fig[1, 1]; title="Flow profile", xlabel="x", ylabel="H"),
           Axis(fig[2, 1]; title="Front position", xlabel="t", ylabel="xᶠ"))
    axs[1].title = "Flow profile - initial conditions"

    plt = (lines!(axs[1], xc, h[:, 100]; color=:blue, label="initial"),
           lines!(axs[1], xc, h[:, 100]; color=:red, label="numerical"),
           lines!(axs[1], xc, hᵉ; color=:black, linestyle=:dash, label="exact"),
            lines!(axs[2], Rs; color=:red, label="numerical"),
            lines!(axs[2], Rsᵉ; color=:black, linestyle=:dash, label="exact from formula"),
            lines!(axs[2], Rsee; color=:black, label="exact from vector"))
    axislegend(axs[1], labelsize=10)
    # axislegend(axs[2]; position=:lt, labelsize=10)
    axislegend(
    axs[2];
    position=:lt,
    labelsize=10,
    patchsize=(12, 8),
    padding=(3, 3, 3, 3),
    rowgap=0
)
    display(fig)
    save("initial.png", fig)

    t_n = 0.0 # current time
    it = 0 # time iteration
    while t_n < t_e
        # @. φ  = ρʷg * h + epsi
        # @. ∇φ_h = (h[2:end, :] - h[1:end-1, :]) / dx
        # @. ∇φ_v = (h[:, 2:end] - h[:, 1:end-1]) / dy

        # # central differences
        # @. cdiff_y[:, 2:end-1] = ((h[2:end, 3:end] - h[2:end, 1:end-2])/dy + (h[1:end-1, 3:end] - h[1:end-1, 1:end-2])/dy)/4
        # @. cdiff_x[2:end-1, :] = ((h[3:end, 2:end] - h[1:end-2, 2:end])/dx + (h[3:end, 1:end-1] - h[1:end-2, 1:end-1])/dx)/4

        # # defines edges as central differnces require ghost cells (phi is set to zero)
        # @. cdiff_y[:, 1] = ((h[2:end, 2])/dy + (h[1:end-1, 2])/dy)/4
        # @. cdiff_y[:, end] = - ((h[2:end, end-1])/dy + (h[1:end-1, end-1])/dy)/4
        # @. cdiff_x[1, :] = ((h[2, 2:end])/dx + (h[2, 1:end-1])/dx)/4
        # @. cdiff_x[end, :] = - ((h[end-1, 2:end])/dx + (h[end-1, 1:end-1])/dx)/4

        # @. A_v = (((φ[2:end, :] - φ[1:end-1, :])/dx)^2 + cdiff_y^2 + epsi^2)^(betha/2 - 1) * (φ[2:end, :] - φ[1:end-1, :]) / dx
        # @. A_h = (cdiff_x^2 + ((φ[:, 2:end] - φ[:, 1:end-1])/dy)^2 + epsi^2)^(betha/2 - 1) * (φ[:, 2:end] - φ[:, 1:end-1]) / dy


        @. A_v = (h[2:end, :] - h[1:end-1, :]) / dx
        @. A_h = (h[:, 2:end] - h[:, 1:end-1]) / dy

        # f is monotone function -> max of derivative of f is derivative of f at max (h) on intervall
        @. a_v = alpha * max(h[1:end-1, :], h[2:end, :])^(alpha - 1)
        @. a_h = alpha * max(h[:, 1:end-1], h[:, 2:end])^(alpha - 1)

        # Darcy(-Weisbach) water flux
        @. q_v[2:end-1, :] = -c * 0.5 * (h[1:end-1, :]^alpha + h[2:end, :]^alpha) * A_v - c * a_v * 0.5 * abs(A_v) * (h[2:end, :] - h[1:end-1, :])
        @. q_h[:, 2:end-1] = -c * 0.5 * (h[:, 1:end-1]^alpha + h[:, 2:end]^alpha) * A_h - c * a_h * 0.5 * abs(A_h) * (h[:, 2:end] - h[:, 1:end-1])

        # diffusive time step
        hmax = maximum(h)
        dt = 1e-12 #1e-10 #1e-7

        # update water sheet thickness using explicit euler scheme
        @. h -= dt * ((q_v[2:end, :] - q_v[1:end-1, :]) / dx + (q_h[:, 2:end] - q_h[:, 1:end-1]) / dy)
        
        t_n += dt
        it  += 1

        if it % nvis == 0
            # exact profile (similarity solution)
            hᵉ = profile(xc, t_n)

            # estimate front
            Rᵉ = front(t_n, Q_0)
            push!(Rsᵉ, Point2f(t_n, Rᵉ))
            
            for ifr in (nx-1):-1:1
                ϵ = 1e-3H_0
                if h[ifr] > ϵ && h[ifr+1] <= ϵ
                    R = xc[ifr]
                    push!(Rs, Point2f(t_n, R))                
                end
                if hᵉ[ifr] > ϵ && hᵉ[ifr+1] <= ϵ
                    Ree = xc[ifr]
                    push!(Rsee, Point2f(t_n, Ree))                
                end
            end

            # update plot
            axs[1].title = "Flow profile —  t = $(@sprintf("%.3e", t_n))"
            plt[2][2] = h[:, 100]
            plt[3][2] = hᵉ
            plt[4][1] = Rs
            plt[5][1] = Rsᵉ
            plt[6][1] = Rsee
            display(fig)
        end
    end
    return
end

main()

