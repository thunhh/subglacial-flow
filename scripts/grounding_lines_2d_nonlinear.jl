using CairoMakie
using Printf
using Statistics

function grounding_lines_1d()
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
    nt   = 50000 #100000
    nvis = 5000 # 1000
    # preprocessing
    dx = lx / (nx - 1)
    dy = ly / (ny - 1)
    xn = LinRange(0, lx, nx)
    yn = LinRange(0, ly, ny)
    # arrays
    H   = zeros(nx, ny)
    B   = zeros(nx, ny)
    h   = zeros(nx, ny)
    # without ghost cells
    φ   = zeros(nx, ny)

    # with ghost cells
    # φ   = zeros(nx + 2, ny + 2)

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
    # initialisation
    # H - ice thickness
    # H .= 4000.0 .- reshape(xn, :, 1) ./ 1e2
    yfactor = 0.2 * sin.(range(0, π, length=ny))
    H .= (4000.0 .- reshape(xn, :, 1) ./ 1e2) .* reshape(yfactor, 1, :)
    H[xn .> 9lx/10, :] .= 0
    # B - bed elevation
    @. B = 1.4e3 + 0.2e3 * sin(6π * xn / lx) - xn / 1e2
    # σnn - overburden pressure (of ice sheet)
    @. σnn = ρⁱg * H
    # figure
    # figure 3D
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
    # time loop
    tcur = 0.0
    for it in 1:nt
        # (overburden) hydraulic potential = overburden pressure + elevation potential + water pressure
        # with ghost cells
        # @. φ[2:end-1, 2:end-1]  = σnn + ρʷg * (B + h)
        @. φ  = σnn + ρʷg * (B + h)
        @. ∇φ_h = (φ[2:end, :] - φ[1:end-1, :]) / dx
        @. ∇φ_v = (φ[:, 2:end] - φ[:, 1:end-1]) / dy
        # regularizor more complex
        #@. epsi = 1e-5 * mean(abs.((φ[2:end] - φ[1:end-1]) ./ (dx)))


        # central differences
        # @. cdiff_y[:, 2:end-1] = ((φ[2:end, 3:end] - φ[2:end, 1:end-2])/dy + (φ[1:end-1, 3:end] - φ[1:end-1, 1:end-2])/dy)/4
        # @. cdiff_x[2:end-1, :] = ((φ[3:end, 2:end] - φ[1:end-2, 2:end])/dx + (φ[3:end, 1:end-1] - φ[1:end-2, 1:end-1])/dx)/4

        # with gradients at both interfaces, more accurate for non-uniform grids
        @. cdiff_y[:, 2:end-1] = ((φ[2:end, 2:end-1] - φ[2:end, 1:end-2])/dy + (φ[2:end, 3:end] - φ[2:end, 2:end-1])/dy + (φ[1:end-1, 3:end] - φ[1:end-1, 2:end-1])/dy + (φ[1:end-1, 2:end-1] - φ[1:end-1, 1:end-2])/dy)/4
        @. cdiff_x[2:end-1, :] = ((φ[3:end, 2:end] - φ[2:end-1, 2:end])/dx + (φ[2:end-1, 2:end] - φ[1:end-2, 2:end])/dx + (φ[3:end, 1:end-1] - φ[2:end-1, 1:end-1])/dx + (φ[2:end-1, 1:end-1] - φ[1:end-2, 1:end-1])/dx)/4

        # defines edges as central differnces require ghost cells (phi is set to zero)
        @. cdiff_y[:, 1] = ((φ[2:end, 2])/dy + (φ[1:end-1, 2])/dy)/4
        @. cdiff_y[:, end] = - ((φ[2:end, end-1])/dy + ( φ[1:end-1, end-1])/dy)/4
        @. cdiff_x[1, :] = ((φ[2, 2:end])/dx + (φ[2, 1:end-1])/dx)/4
        @. cdiff_x[end, :] = - ((φ[end-1, 2:end])/dx + (φ[end-1, 1:end-1])/dx)/4

        @. A_v = (((φ[2:end, :] - φ[1:end-1, :])/dx)^2 + cdiff_y^2 + epsi^2)^(betha/2 - 1) * (φ[2:end, :] - φ[1:end-1, :]) / dx
        @. A_h = (cdiff_x^2 + ((φ[:, 2:end] - φ[:, 1:end-1])/dy)^2 + epsi^2)^(betha/2 - 1) * (φ[:, 2:end] - φ[:, 1:end-1]) / dy

        
        @. a_v = alpha * max(h[1:end-1, :], h[2:end, :])^(alpha - 1)
        @. a_h = alpha * max(h[:, 1:end-1], h[:, 2:end])^(alpha - 1)
         
        # Darcy(-Weisbach) water flux
        @. q_v[2:end-1, :] = -k * 0.5 * (h[1:end-1, :]^alpha + h[2:end, :]^alpha) * A_v - k * a_v * 0.5 * abs(A_v) * (h[2:end, :] - h[1:end-1, :])
        @. q_h[:, 2:end-1] = -k * 0.5 * (h[:, 1:end-1]^alpha + h[:, 2:end]^alpha) * A_h - k * a_h * 0.5 * abs(A_h) * (h[:, 2:end] - h[:, 1:end-1])

        # advective time step
        #dta = dx / k / maximum(abs, ∇φ) / 2.1

        # diffusive time step
        # dtd = dx^2 / (k * ρʷg * maximum(h)^alpha * (max(maximum(abs.(∇φ)), 1e-5)^(betha - 2))) / 10 #/ 2.1
        # dtd = dx^2 / (k * ρʷg * maximum(h)^alpha * maximum(abs.(∇φ))^(betha - 2)) / 10 # / 2.1
        dta = min(dx, dy) / k / max(maximum(abs, ∇φ_h), maximum(abs, ∇φ_v)) / 2.1
        # diffusive time step
        dtd = min(dx^2, dy^2) / (k * ρʷg * maximum(h)) / 2.1

        dt = min(dta, dtd) 
        # update water sheet thickness using explicit euler scheme
        # y' approx. by (q[2:end] - q[1:end-1]) / dx forward differnces
        # @. h -= dt * (q[2:end] - q[1:end-1]) / dx
        @. h -= dt * ((q_v[2:end, :] - q_v[1:end-1, :]) / dx + (q_h[:, 2:end] - q_h[:, 1:end-1]) / dy)
        # h[end] = 4.2e3
        h[end, :] .= 4.2e3

        # update plot
        if it % nvis == 0
            @printf(" t = %.1f d, dt [adv] = %1.3e d, dt [dif] = %1.3e d\n", tcur / 3600 / 24, dta / 3600 / 24, dtd / 3600 / 24)

            p1[3] = B
            p2[3] = B .+ h
            p3[3] = B .+ h .+ H

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

