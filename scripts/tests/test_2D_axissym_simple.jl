using CairoMakie
using SpecialFunctions
using Printf

ETA_N   = (2^10 / (3^4 * pi^3))^(1/8)
V       = 1
G_      = 1
Q_0     = pi 

function front(t, q)
    t_0 = 1e-8
    T = t + t_0
    return ETA_N * (G_ * Q_0^3 * T / 3 / V)^(1/8) 
end

# flow profile
function profile(x, t)
    rho_a   = 1.204     # density air
    rho_w   = 1000     # density water
    g       = 9.81
    g_      = (rho_w - rho_a) / g
    v       = 1.004e-6  # viscosity for 20° C water
    Q = Q_0 
    t_0 = 1e-8 #0.0001
    T   = t_0 + t

    eta = (1/3 * g_ * Q^3 / v)^(-1/8) .* x .* T^(-1/8)
    y = min.(abs.(eta ./ ETA_N), 1.0)
    phi = (3/16)^(1/3) * (1 .- y.^2).^(1/3)
    h = ETA_N^(2/3) * (3 * Q^2 * v / g_)^(1/4) * T^(-1/4) .* phi
    return h
end

function profile_deg(x, t, α, β, H₀, R₀, d)
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
@views function main()
    # physics
    lx = 10.0 # domain length (diameter for axisymmetric)
    t_e = 10.0 # total time of the simulation

    # for deg profile
    H_0 = 1.0  # initial thickness
    R_0 = 1.0  # initial radius
    alpha  = 3 #5/4  # 5/4 for Darcy-Weisbach, 3 for Newtonian viscous flow, 5 for ice flow
    betha  = 2 #3/2  # 3/2 for Darcy-Weisbach, 2 Newtonian viscous flow, 4 for ice flow
    d  = 2   # 1 for cross-section, 2 for axisymmetric

    # water/air variables
    c = 1 # for Huppert set to 1/3

    # numerics
    nvis = 2000
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
    h_deg = zeros(nx)

    # initial conditions
    r = sqrt.(xc.^2 .+ yc'.^2)
    # re = abs.(xc)
    j = argmin(abs.(yc))
    #  j = 80 # test also outside the middle axis
    re = r[:, j]
    @. h_deg = profile_deg(re, 0.0, alpha, betha, H_0, R_0, d)
    @. h = profile_deg(r, 0.0, alpha, betha, H_0, R_0, d)

    fig1 = Figure(size=(900, 700))
    ax = Axis3(
        fig1[1, 1],
        xlabel="x",
        ylabel="y",
        zlabel="h(x,y,0)",
        title="Similarity solution at t = 0"
    )

    surface!(
        ax, xc, yc, h;
        colormap=:viridis
    )

    display(fig1)

    
    hᵉ = profile(re, 0.0)
    # initial front position
    Rsᵉ = Point2f[(0.0, front(0.0, Q_0))] 
    Rs = copy(Rsᵉ)  
    Rsee = copy(Rsᵉ)

    # visualisation
    fig = Figure()

    axs = (Axis(fig[1, 1]; title="Flow profile", xlabel="x", ylabel="H"),)
    axs[1].title = "Flow profile - initial conditions"

    plt = (lines!(axs[1], xc, h[:, j]; color=:blue, label="initial"),
           lines!(axs[1], xc, h[:, j]; color=:red, label="numerical"),
           lines!(axs[1], xc, h_deg; color=:green, linestyle=:dash, label="deg sol"))
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
        dt = 1e-4 

        # update water sheet thickness using explicit euler scheme
        @. h -= dt * ((q_v[2:end, :] - q_v[1:end-1, :]) / dx + (q_h[:, 2:end] - q_h[:, 1:end-1]) / dy)

        
        t_n += dt
        it  += 1

        if it % nvis == 0
            # exact profile (similarity solution)
            hᵉ = profile(xc, t_n)
            @. h_deg = profile_deg(re, t_n, alpha, betha, H_0, R_0, d)

            # estimate front
            Rᵉ = front(t_n, Q_0)
            push!(Rsᵉ, Point2f(t_n, Rᵉ))
            
            for ifr in (nx-1):-1:1
                ϵ = 1e-3H_0
                if h[ifr, j] > ϵ && h[ifr+1, j] <= ϵ
                    R = xc[ifr]
                    push!(Rs, Point2f(t_n, R))                
                end
                if hᵉ[ifr] > ϵ && hᵉ[ifr+1] <= ϵ
                    Ree = xc[ifr]
                    push!(Rsee, Point2f(t_n, Ree))                
                end
            end

            axs[1].title = "Flow profile —  t = $(@sprintf("%.3e", t_n))"
            plt[2][2] = h[:, j]
            plt[3][2] = h_deg
            display(fig)
        end
    end
    return
end

main()

