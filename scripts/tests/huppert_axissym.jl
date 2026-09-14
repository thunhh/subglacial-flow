using Printf
using SpecialFunctions
using Plots
using Serialization




function huppert_axissym()
    # physics
    lx      = 10 #100e3
    rho_a   = 1.204     # density air
    rho_w   = 1000     # density water
    g       = 9.81
    g_      = (rho_w - rho_a) / g
    v       = 1.004e-6  # viscosity for 20° C water
    # eta_n   = (1/56 * (3/10)^(1/3) * pi^(1/2) * gamma(1/3) * gamma(5/6))^(-3/5)
    eta_n   = (2^10 / (3^4 * pi^3))^(1/8)
    q       = 2 * pi * 5^2 #0.6 * pi * 5^2       # water volume -> wieso hmax > 1.2

    # numerics
    nx  = 1000
    t   = lx^8 * 3 * v / (eta_n^8) / g_ / (q^3)
    println(t)

    # preprocessing
    dx  = lx / (nx - 1)
    xn  = LinRange(0, lx, nx)

    # arrays
    eta = zeros(nx)
    y   = zeros(nx)
    phi = zeros(nx)
    h   = zeros(nx)

    # similarity solution 
    @. eta  = (1/3 * g_ * q^3 / v)^(-1/8) * xn * t^(-1/8)
    println(maximum(eta))
    @. y    = eta/eta_n
    y[end]  = min(1,y[end])
    print(maximum(y))
    @. phi  = (3/16)^(1/3) * (1 - y^2)^(1/3)
    @. h    = eta_n^(2/3) * (3*q*v/g_)^(1/4) * t^(-1/4) * phi

    serialize("h_huppert_selfsim.jls", h)


    println("h_max = ", maximum(h))
    
    p = plot(y, phi, xlabel="y", ylabel="ϕ")
    display(p)
    savefig("huppert_y_phi_axsym.png")

    p = plot(xn, phi, xlabel="x", ylabel="ϕ")
    display(p)
    savefig("huppert_x_phi_axsym.png")


    p = plot(xn, h, xlabel="x", ylabel="h")
    display(p)
    savefig("huppert_x_h_axsym.png")

    # compute volume
    # I = sum(diff(xn) .* (h[2:end] .+ h[1:end-1]) ./ 2)
    I = sum(diff(xn) .* ((xn[2:end].*h[2:end]) + (xn[1:end-1].*h[1:end-1])) ./ 2)
    println("Initial Volume Q = ", q)
    println("Volume of solution = ", 2* pi*I)


end

huppert_axissym()
    
# T = 0.0005055718055700804 for lx = 10
# T = 5.055718055700805e-9 for lx = 1