using CairoMakie
using Printf
using Statistics
using Serialization

function grounding_lines_1d()
    # physics
    lx  = 30
    ly  = 30
    k   = 0.001
    ρⁱg = 910.0 * 9.81
    ρʷg = 1000.0 * 9.81
    alpha = 3
    betha = 2
    g = 9.81
    v = 1e-6
    c = g/v/3
    # numerics
    nx   = 301 #100
    ny   = 301 #100
    nt   = 400000 #300000 #10000 #100000
    nvis = 30000 #30000 # 1000
    T =  0.0005055718055700804

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
    # set ice thickness and surface to zero
    # initialize h with instantanius droplet
    mid_x = nx ÷ 2
    mid_y = ny ÷ 2
    idx_x = 1:nx
    idx_y = 1:ny

    mask = zeros(nx, ny)
    dist = sqrt.((idx_x .- mid_x).^2 .+ (idx_y' .- mid_y).^2)
    mask = dist .< 50

    h = 0.12 .* mask #2000.0 .* mask


    # other potential initial condition -> droplet form
    # R0 = 10.0
    # H0 = 2000.0

    # h = H0 .* (1 .- (dist ./ R0).^2) .* mask
    # σnn - overburden pressure (of ice sheet)
    @. σnn = ρⁱg * H

    # figure
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

    # figure 2D
    fig = Figure(; size=(600, 600))
    axs = (Axis(fig[1, 1]; ylabel="z [m]"),
           Axis(fig[2, 1]; ylabel="φ [bar]"),
           Axis(fig[3, 1]; xlabel="x [km]", ylabel="∇φ [bar/km]"))
    plt = (band!(axs[1], xn ./ 1e3, zeros(nx), B[:,mid_y]; color=:brown),
           band!(axs[1], xn ./ 1e3, B[:,mid_y], B[:,mid_y] .+ h[:,mid_y]; color=:blue),
           band!(axs[1], xn ./ 1e3, B[:,mid_y] .+ h[:,mid_y], B[:,mid_y] .+ h[:,mid_y] .+ H[:,mid_y]; color=:lightblue),
           lines!(axs[2], xn ./ 1e3, φ[:,mid_y] ./ 1e5),
           lines!(axs[3], xn[1:end-1] ./ 1e3, ∇φ_h[:,mid_y] ./ 1e2))
    display(fig)

    # time loop
    tcur = 0.0
    for it in 1:nt
        # (overburden) hydraulic potential = overburden pressure + elevation potential + water pressure
        # with ghost cells
        # @. φ[2:end-1, 2:end-1]  = σnn + ρʷg * (B + h)
        @. φ  = σnn + ρʷg * (B + h) + epsi
        @. ∇φ_h = (φ[2:end, :] - φ[1:end-1, :]) / dx
        @. ∇φ_v = (φ[:, 2:end] - φ[:, 1:end-1]) / dy
        # regularizor more complex
        #@. epsi = 1e-5 * mean(abs.((φ[2:end] - φ[1:end-1]) ./ (dx)))


        # central differences
        @. cdiff_y[:, 2:end-1] = ((φ[2:end, 3:end] - φ[2:end, 1:end-2])/dy + (φ[1:end-1, 3:end] - φ[1:end-1, 1:end-2])/dy)/4
        @. cdiff_x[2:end-1, :] = ((φ[3:end, 2:end] - φ[1:end-2, 2:end])/dx + (φ[3:end, 1:end-1] - φ[1:end-2, 1:end-1])/dx)/4

        # with gradients at both interfaces, more accurate for non-uniform grids
        # @. cdiff_y[:, 2:end-1] = ((φ[2:end, 2:end-1] - φ[2:end, 1:end-2])/dy + (φ[2:end, 3:end] - φ[2:end, 2:end-1])/dy + (φ[1:end-1, 3:end] - φ[1:end-1, 2:end-1])/dy + (φ[1:end-1, 2:end-1] - φ[1:end-1, 1:end-2])/dy)/4
        # @. cdiff_x[2:end-1, :] = ((φ[3:end, 2:end] - φ[2:end-1, 2:end])/dx + (φ[2:end-1, 2:end] - φ[1:end-2, 2:end])/dx + (φ[3:end, 1:end-1] - φ[2:end-1, 1:end-1])/dx + (φ[2:end-1, 1:end-1] - φ[1:end-2, 1:end-1])/dx)/4

        # defines edges as central differnces require ghost cells (phi is set to zero)
        @. cdiff_y[:, 1] = ((φ[2:end, 2])/dy + (φ[1:end-1, 2])/dy)/4
        @. cdiff_y[:, end] = - ((φ[2:end, end-1])/dy + ( φ[1:end-1, end-1])/dy)/4
        @. cdiff_x[1, :] = ((φ[2, 2:end])/dx + (φ[2, 1:end-1])/dx)/4
        @. cdiff_x[end, :] = - ((φ[end-1, 2:end])/dx + (φ[end-1, 1:end-1])/dx)/4

        
        # @. A_v = (((φ[2:end, :] - φ[1:end-1, :])/dx)^2 + cdiff_y^2 + epsi^2)^(betha/2 - 1) * (φ[2:end, :] - φ[1:end-1, :]) / dx
        # @. A_h = (cdiff_x^2 + ((φ[:, 2:end] - φ[:, 1:end-1])/dy)^2 + epsi^2)^(betha/2 - 1) * (φ[:, 2:end] - φ[:, 1:end-1]) / dy

        @. A_v = (φ[2:end, :] - φ[1:end-1, :]) / dx
        @. A_h = (φ[:, 2:end] - φ[:, 1:end-1]) / dy

        # with ghost cells
        # @. A_v = (((φ[3:end-1, 2:end] - φ[2:end-2, 2:end])/dx)^2 + (I)^2 + epsi^2)^(betha/2 - 1) * (φ[3:end-1, 2:end] - φ[2:end-2, 2:end]) / dx
        # @. A_h = ((II)^2 + ((φ[2:end, 3:end-1] - φ[2:end, 2:end-2])/dy)^2 + epsi^2)^(betha/2 - 1) * (φ[2:end, 3:end-1] - φ[:, 2:end-2]) / dy

        # f is monotone function -> max of derivative of f is derivative of f at max (h) on intervall
        @. a_v = alpha * max(h[1:end-1, :], h[2:end, :])^(alpha - 1)
        @. a_h = alpha * max(h[:, 1:end-1], h[:, 2:end])^(alpha - 1)
         
        # Darcy(-Weisbach) water flux
        @. q_v[2:end-1, :] = -c * 0.5 * (h[1:end-1, :]^alpha + h[2:end, :]^alpha) * A_v - c * a_v * 0.5 * abs(A_v) * (h[2:end, :] - h[1:end-1, :])
        @. q_h[:, 2:end-1] = -c * 0.5 * (h[:, 1:end-1]^alpha + h[:, 2:end]^alpha) * A_h - c * a_h * 0.5 * abs(A_h) * (h[:, 2:end] - h[:, 1:end-1])

        # advective time step
        #dta = dx / k / maximum(abs, ∇φ) / 2.1

        # diffusive time step
        # dtd = dx^2 / (k * ρʷg * maximum(h)^alpha * (max(maximum(abs.(∇φ)), 1e-5)^(betha - 2))) / 10 #/ 2.1
        # dtd = dx^2 / (k * ρʷg * maximum(h)^alpha * maximum(abs.(∇φ))^(betha - 2)) / 10 # / 2.1
        # dta = min(dx, dy) / c / max(maximum(abs, ∇φ_h), maximum(abs, ∇φ_v)) / 2.1 
        # # diffusive time step
        # dtd = min(dx^2, dy^2) / (c * ρʷg * maximum(h)) / 2.1
        # dt = min(dta, dtd) 

        # time stepping    
        hmax = maximum(h)

        gradmax = max(
            maximum(abs, ∇φ_h),
            maximum(abs, ∇φ_v)
        )

        dta = min(dx, dy) / (c * hmax^alpha * max(gradmax, 1e-12)) / 10
        dtd = min(dx^2, dy^2) / (c * ρʷg * hmax^alpha) / 10
        dt = min(dta, dtd)


        # update water sheet thickness using explicit euler scheme
        h_old = copy(h)
        @. h -= dt * ((q_v[2:end, :] - q_v[1:end-1, :]) / dx + (q_h[:, 2:end] - q_h[:, 1:end-1]) / dy)
        
        # figure 3D
        # update plot
        # if it % nvis == 0
        #     @printf(" t = %.1f s, dt [adv] = %1.3e s, dt [dif] = %1.3e s\n", tcur, dta, dtd)
        #     println(maximum(h))
        #     p1[3] = B
        #     p2[3] = B .+ h
        #     p3[3] = B .+ h .+ H

        #     display(fig)
        #     change = maximum(abs.(h .- h_old))
        #     relative_change = change / max(maximum(abs.(h_old)), 1e-12)

        #     @printf("it = %d, max(h) = %.6e, change = %.6e, rel = %.6e\n",
        #             it, maximum(h), change, relative_change)

        #     end
        #     tcur += dt

            # figure 2D
        if it % nvis == 0
            @printf(" t = %.5f , dt [adv] = %1.3e , dt [dif] = %1.3e \n", tcur, dta , dtd)
            println(maximum(h))
            println(minimum(h))

            plt[2][3] = B[:,mid_y] .+ h[:,mid_y]
            plt[3][2] = B[:,mid_y] .+ h[:,mid_y]
            plt[3][3] = B[:,mid_y] .+ h[:,mid_y] .+ H[:,mid_y]
            plt[4][2] = φ[:,mid_y] ./ 1e5
            plt[5][2] = ∇φ_h[:,mid_y] ./ 1e2
            display(fig)
        end
        tcur += dt

        if tcur >= T
            serialize("h_huppert_2D.jls", h)
            # read fie via v = deserialize("vector.jls")
            plt[2][3] = B[:,mid_y] .+ h[:,mid_y]
            plt[3][2] = B[:,mid_y] .+ h[:,mid_y]
            plt[3][3] = B[:,mid_y] .+ h[:,mid_y] .+ H[:,mid_y]
            plt[4][2] = φ[:,mid_y] ./ 1e5
            plt[5][2] = ∇φ_h[:,mid_y] ./ 1e2
            display(fig)
            break
        end

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

