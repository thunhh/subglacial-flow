using CairoMakie
using Printf

function grounding_lines_1d()
    # physics
    lx  = 100e3
    k   = 0.001
    ρⁱg = 910.0 * 9.81
    ρʷg = 1000.0 * 9.81
    # numerics
    nx   = 100
    nt   = 100000 #100000
    nvis = 10000 # 1000
    # preprocessing
    dx = lx / (nx - 1)
    xn = LinRange(0, lx, nx)
    # arrays
    H   = zeros(nx)
    B   = zeros(nx)
    h   = zeros(nx)
    hh  = zeros(nx + 2) # add 2 ghost cells
    φ   = zeros(nx)
    ∇φ  = zeros(nx - 1)
    φ_1st   = zeros(nx)
    ∇φ_1st  = zeros(nx - 1)
    q   = zeros(nx + 1) # also boundary fluxes (interior would be nx - 1)
    q_1st = zeros(nx + 1)
    σnn = zeros(nx)
    sl  = ones(nx)
    h_neg = zeros(nx+2)
    h_pos  = zeros(nx+2)
    limiter = zeros(nx) 
    # initialisation
    # H - ice thickness
    @. H = 4000.0 - xn / 1e2
    @. H[xn>9lx/10] = 0
    # B - bed elevation
    @. B = 1.4e3 + 0.2e3 * sin(6π * xn / lx) - xn / 1e2
    # σnn - overburden pressure (of ice sheet)
    @. σnn = ρⁱg * H
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
    # time loop
    tcur = 0.0
    for it in 1:nt
        # (overburden) hydraulic potential = overburden pressure + elevation potential + water pressure
        @. φ  = σnn + ρʷg * (B + hh[2:end-1])
        @. ∇φ = (φ[2:end] - φ[1:end-1]) / dx

        @. φ_1st  = σnn + ρʷg * (B + h)
        @. ∇φ_1st = (φ_1st[2:end] - φ_1st[1:end-1]) / dx

        
        @assert hh[2:end-1] ≈ h # are equal
        @assert φ_1st ≈ φ # are not equal
        @assert ∇φ_1st ≈ ∇φ # are equal

        # slope limiter
        # dif1 = h[2:end-1] - h[1:end-2]
        # dif2 = h[3:end] - h[2:end-1]
        # @. sl = ifelse(dif2 == 0, 1.0, dif1/dif2)

        # @. sl = (h[2:end-1] - h[1:end-2]) / (h[3:end] - h[2:end-1]) -> leads to NaNs
        # @. limiter = (sl + abs(sl)) / (1 + abs(sl))
        # für h_neg(i = 1) bräuchte ich i = 0 für limiter, existiert nicht, deshalb bei 2 starten und end -1 enden
        h_neg = copy(hh)
        h_pos = copy(hh)
        @. h_neg[2:end-1] += 0.5 * limiter[1:end] * (hh[3:end] - hh[2:end-1])
        @. h_pos[2:end-1] -= 0.5 * limiter[1:end] * (hh[3:end] - hh[2:end-1])

        @assert h_neg[2:end-1] ≈ h
        @assert h_pos[2:end-1] ≈ h

        # 2nd order scheme for Darcy(-Weisbach) water flux
        # update interior fluxes
        @. q[2:end-1] = -k * 0.5 * (h_neg[2:end-2] + h_pos[3:end-1]) * ∇φ - k * 0.5 * abs(∇φ) * (h_pos[3:end-1] - h_neg[2:end-2])

        # Darcy(-Weisbach) water flux
        @. q_1st[2:end-1] = -k * 0.5 * (h[1:end-1] + h[2:end]) * ∇φ - k * 0.5 * abs(∇φ) * (h[2:end] - h[1:end-1])
        @assert q_1st ≈ q  # assert fails

        # advective time step
        dta = dx / k / maximum(abs, ∇φ) / 2.1
        # diffusive time step
        dtd = dx^2 / (k * ρʷg * maximum(hh)) / 2.1
        dt = min(dta, dtd) 

        dta_1st = dx / k / maximum(abs, ∇φ_1st) / 2.1
        # diffusive time step
        dtd_1st = dx^2 / (k * ρʷg * maximum(h)) / 2.1
        dt_1st = min(dta_1st, dtd_1st) 

        dt = 1.356e-4 * 3600 * 24
        dt_1st = 1.356e-4 * 3600 *24

        # update water sheet thickness using explicit euler scheme
        # y' approx. by (q[2:end] - q[1:end-1]) / dx forward differnces
        @. hh[2:end-1] -= dt * (q[2:end] - q[1:end-1]) / dx
        @. h -= dt_1st * (q_1st[2:end] - q_1st[1:end-1]) / dx
        # before φ
        @assert isapprox(
            dt * (q[2:end] - q[1:end-1]) / dx,
            dt_1st * (q_1st[2:end] - q_1st[1:end-1]) / dx,
            rtol=1e-12
        )

        @assert hh[2:end-1] ≈ h


        # boundary conditions -> how to def. ghost cells?
        hh[end-1] = 4.2e3
        hh[end] = hh[end-1]
        #hh[1] = hh[2]
        h[end] = 4.2e3

        
        
   

        # update plot
        if it % nvis == 0
            @printf(" t = %1.7f d, dt [adv] = %1.3e d, dt [dif] = %1.3e d\n", tcur / 3600 / 24, dta / 3600 / 24, dtd / 3600 / 24)
            # compare the qs
            #@assert isapprox(q_1st, q, atol=1e-14)
            println("iteration", it)
            plt[2][3] = B .+ hh[2:end-1]
            plt[3][2] = B .+ hh[2:end-1]
            plt[3][3] = B .+ hh[2:end-1] .+ H
            plt[4][2] = φ ./ 1e5
            plt[5][2] = ∇φ ./ 1e2
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

