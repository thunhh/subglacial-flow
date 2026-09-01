using Printf
using SpecialFunctions
using Plots



function huppert_analytic()
    # physics
    lx      = 10 #100e3
    rho_a   = 1.204     # density air
    rho_w   = 1000     # density water
    g       = 9.81
    g_      = (rho_w - rho_a) / g
    v       = 1.004e-6  # viscosity for 20° C water
    eta_n   = (1/5 * (3/10)^(1/3) * pi^(1/2) * gamma(1/3) * gamma(5/6))^(-3/5)
    q       = 0.6 #1.2 * pi * 5^2        # water volume

    # numerics
    nx  = 100
    t   = 0.0005055718055700804

    # preprocessing
    dx  = lx / (nx - 1)
    xn  = LinRange(0, lx, nx)

    # arrays
    eta = zeros(nx)
    y   = zeros(nx)
    phi = zeros(nx)
    h   = zeros(nx)

    # similarity solution
    @. eta  = (1/3 * g_ * q^3 / v)^(-1/5) * xn * t^(-1/5)
    @. y    = eta/eta_n
    @. phi  = (3/10)^(1/3) * (1 - y^2)^(1/3)
    @. h    = eta_n^(2/3) * (3*q^2*v/g_)^(1/5) * t^(-1/5) * phi

    T = (eta_n / lx)^(-5) * 3 * (v/g_/(q^3))
    TT= (lx / eta_n)^5 * 3 * v /g_ / q / q / q

    println(y)
    println(phi)
    println(h)
    println(eta_n)
    println(T)
    println(TT)
    println("h_max = ", maximum(h))
    
    p = plot(y, phi)
    display(p)
    savefig("huppert.png")

    p = plot(xn, h)
    display(p)
    savefig("huppert_h.png")

    # compute volume
    I = sum(diff(xn) .* (h[2:end] .+ h[1:end-1]) ./ 2)
    println("Initial Volume Q = ", q)
    println("Volume of solution = ", I)



end

huppert_analytic()
    
# T = 0.0005055718055700804 for lx = 10
# T = 5.055718055700805e-9 for lx = 1