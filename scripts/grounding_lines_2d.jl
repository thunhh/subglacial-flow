using CairoMakie
using Printf

function grounding_lines_2d()
    # physics
    lx  = 100e3
    ly  = 100e3
    k   = 0.001
    ρⁱg = 910.0 * 9.81
    ρʷg = 1000.0 * 9.81
    # numerics
    nx   = 100
    ny   = 100
    nt   = 100000 #100000
    nvis = 20000 # 1000 # 20000
    # preprocessing
    dx = lx / (nx - 1)
    dy = ly / (ny - 1)
    xn = LinRange(0, lx, nx)
    yn = LinRange(0, ly, ny)
    # arrays
    H   = zeros(nx, ny)
    B   = zeros(nx, ny)
    h   = zeros(nx, ny)
    φ   = zeros(nx, ny)
    ∇φ_h  = zeros(nx - 1, ny)
    ∇φ_v  = zeros(nx, ny - 1)
    q_h   = zeros(nx + 1, ny)
    q_v   = zeros(nx, ny + 1)
    σnn = zeros(nx, ny)
    # initialisation
    # H - ice thickness
    #H .= reshape(@. 4000.0 - xn / 1e2, :, 1)
    H .= 4000.0 .- reshape(xn, :, 1) ./ 1e2
    H[xn .> 9lx/10, :] .= 0
    # B - bed elevation
    B .= reshape(1.4e3 .+ 0.2e3 .* sin.(6π .* xn ./ lx) .- xn ./ 1e2, nx, 1)
    # σnn - overburden pressure (of ice sheet)
    @. σnn = ρⁱg * H
    # figure 2D
    fig = Figure(; size=(600, 600))
    axs = (Axis(fig[1, 1]; ylabel="z [m]"),
           Axis(fig[2, 1]; ylabel="φ [bar]"),
           Axis(fig[3, 1]; xlabel="x [km]", ylabel="∇φ [bar/km]"))
    plt = (band!(axs[1], xn ./ 1e3, zeros(nx), B[:,1]; color=:brown),
           band!(axs[1], xn ./ 1e3, B[:,1], B[:,1] .+ h[:,1]; color=:blue),
           band!(axs[1], xn ./ 1e3, B[:,1] .+ h[:,1], B[:,1] .+ h[:,1] .+ H[:,1]; color=:lightblue),
           lines!(axs[2], xn ./ 1e3, φ[:,1] ./ 1e5),
           lines!(axs[3], xn[1:end-1] ./ 1e3, ∇φ_h[:,1] ./ 1e2))
    display(fig)

    # fig = Figure(size = (800, 600))
    # ax = Axis3(fig[1, 1],
    #     xlabel = "x [km]",
    #     ylabel = "y [km]",
    #     zlabel = "z [m]"
    # )

    # X = xn ./ 1e3          # length nx
    # Y = yn ./ 1e3          # length ny

    # surface!(ax, xn, yn, B')
    # surface!(ax, xn, yn, (B .+ h)')
    # surface!(ax, xn, yn, (B .+ h .+ H)')
    # display(fig)

    # fig = Figure(size = (800, 600))

    # ax = Axis3(fig[1, 1],
    #     xlabel = "x [km]",
    #     ylabel = "y [km]",
    #     zlabel = "z [m]"
    # )

    # X = xn ./ 1e3
    # Y = yn ./ 1e3

    # # Convert interfaces to 2D arrays
    # Z0 = zeros(nx, ny)
    # Z1 = B
    # Z2 = B .+ h
    # Z3 = B .+ h .+ H

    # surface!(ax, X, Y, Z0, color=:brown, transparency=true, alpha=0.8)
    # surface!(ax, X, Y, Z1, color=:blue, transparency=true, alpha=0.8)
    # surface!(ax, X, Y, Z2, color=:lightblue, transparency=true, alpha=0.8)
    # surface!(ax, X, Y, Z3, color=:cyan, transparency=true, alpha=0.8)

    # figure 3D
    # fig = Figure(size = (800, 600))

    # ax = Axis3(fig[1, 1],
    #     xlabel = "x [km]",
    #     ylabel = "y [km]",
    #     zlabel = "z [m]"
    # )

    # X = xn ./ 1e3
    # Y = yn ./ 1e3

    # p1 = surface!(ax, X, Y, B,
    #     color=:brown, transparency=true, alpha=0.8)

    # p2 = surface!(ax, X, Y, B .+ h,
    #     color=:blue, transparency=true, alpha=0.8)

    # p3 = surface!(ax, X, Y, B .+ h .+ H,
    #     color=:lightblue, transparency=true, alpha=0.8)


    # display(fig)

    # time loop
    tcur = 0.0
    #     # (overburden) hydraulic potential = overburden pressure + elevation potential + water pressure
    #     @. φ  = σnn + ρʷg * (B + h)
    #     @. ∇φ = (φ[2:end] - φ[1:end-1]) / dx
    #     # Darcy(-Weisbach) water flux
    #     @. q[2:end-1] = -k * 0.5 * (h[1:end-1] + h[2:end]) * ∇φ - k * 0.5 * abs(∇φ) * (h[2:end] - h[1:end-1])
    #     # advective time step
    #     dta = dx / k / maximum(abs, ∇φ) / 2.1
    #     # diffusive time step
    #     dtd = dx^2 / (k * ρʷg * maximum(h)) / 2.1
    #     dt = min(dta, dtd) 
    #     # update water sheet thickness using explicit euler scheme
    #     # y' approx. by (q[2:end] - q[1:end-1]) / dx forward differnces
    #     @. h -= dt * (q[2:end] - q[1:end-1]) / dx
    #     h[end] = 4.2e3
    #     # update plot
    #     if it % nvis == 0
    #         @printf(" t = %.1f d, dt [adv] = %1.3e d, dt [dif] = %1.3e d\n", tcur / 3600 / 24, dta / 3600 / 24, dtd / 3600 / 24)

    #         plt[2][3] = B .+ h
    #         plt[3][2] = B .+ h
    #         plt[3][3] = B .+ h .+ H
    #         plt[4][2] = φ ./ 1e5
    #         plt[5][2] = ∇φ ./ 1e2
    #         display(fig)
    #     end
    #     tcur += dt
    # end

    # time loop
    # tcur = 0.0
    for it in 1:nt
        # (overburden) hydraulic potential = overburden pressure + elevation potential + water pressure
        @. φ  = σnn + ρʷg * (B + h)
        @. ∇φ_h = (φ[2:end, :] - φ[1:end-1, :]) / dx
        @. ∇φ_v = (φ[:, 2:end] - φ[:, 1:end-1]) / dy
        # Darcy(-Weisbach) water flux
        @. q_h[2:end-1, :] = -k * 0.5 * (h[1:end-1, :] + h[2:end, :]) * ∇φ_h - k * 0.5 * abs(∇φ_h) * (h[2:end, :] - h[1:end-1, :])
        @. q_v[:, 2:end-1] = -k * 0.5 * (h[:, 1:end-1] + h[:, 2:end]) * ∇φ_v - k * 0.5 * abs(∇φ_v) * (h[:, 2:end] - h[:, 1:end-1])
        # advective time step
        dta = min(dx, dy) / k / max(maximum(abs, ∇φ_h), maximum(abs, ∇φ_v)) / 2.1
        # diffusive time step
        dtd = min(dx^2, dy^2) / (k * ρʷg * maximum(h)) / 2.1
        dt = min(dta, dtd) 
        # update water sheet thickness using explicit euler scheme
        # y' approx. by (q[2:end] - q[1:end-1]) / dx forward differnces
        # @. h -= dt * (q_h[2:end, :] - q_h[1:end-1, :] + q_v[:, 2:end] - q_v[:, 1:end-1]) / dx / dy
        @. h -= dt * ((q_h[2:end, :] - q_h[1:end-1, :]) / dx + (q_v[:, 2:end] - q_v[:, 1:end-1]) / dy)
        h[end, :] .= 4.2e3
        # update plot
        # if it % nvis == 0
        #     @printf(" t = %.1f d, dt [adv] = %1.3e d, dt [dif] = %1.3e d\n", tcur / 3600 / 24, dta / 3600 / 24, dtd / 3600 / 24)

        #     p1[3] = B
        #     p2[3] = B .+ h
        #     p3[3] = B .+ h .+ H

        #     display(fig)
        # end
        # tcur += dt

    end
    return
end

time1 = time()
grounding_lines_2d()
time2 = time()
elapsed_time = time2 - time1
print(elapsed_time)

println("Press Enter to quit")
readline()

