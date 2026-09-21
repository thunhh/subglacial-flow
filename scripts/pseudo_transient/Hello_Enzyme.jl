using Enzyme
using LinearAlgebra

@views avx(A) = 0.5 .* (A[1:end-1] .+ A[2:end])

@views function residual!(r, h, hᵒ, k, ρⁱg, ρʷg, H, B, dt, dx)
    @. r = (h[2:end-1] - hᵒ[2:end-1]) / dt - $diff(k * $avx(h) * $diff(ρⁱg * H + ρʷg * (B + h)) / dx) / dx
    # @. r = (h[2:end-1] - hᵒ[2:end-1]) / dt - $diff(k * $diff(ρⁱg * H + ρʷg * (B + h)) / dx) / dx
    return
end

function main(nx)
    x = LinRange(0, 1, nx)
    h   = @. sin(4π * x)
    hᵒ  = @. sin(4π * x + 0.01)
    H   = @. 0.5sin(3π * x)
    B   = @. 0.5cos(3π * x)
    k   = @. 0.5exp(-(10*(x[1:end-1]-0.5))^2)
    ρⁱg = rand()
    ρʷg = rand()
    dx  = 1 / (nx-1)
    dt  = 0.01
    r   = zeros(nx-2)
    b   = zeros(nx-2)

    r̄ = make_zero(r)
    h̄ = make_zero(h)

    # r = b - A * h
    # ∂r/∂h = -A
    # b = r + A * h = r - ∂r/∂h * h

    # Enzyme computes ∂r/∂h * h̄ + ∂r/∂k * k̄ + ...

    h̄ .= h
    r̄ .= 0
    residual!(r, h, hᵒ, k, ρⁱg, ρʷg, H, B, dt, dx)
    Enzyme.autodiff(set_runtime_activity(Enzyme.Forward), residual!, Const, Duplicated(r, r̄), Duplicated(h, h̄), Const(hᵒ), Const(k), Const(ρⁱg), Const(ρʷg), Const(H), Const(B), Const(dt), Const(dx))
    @. b = r - r̄

    b_frozen = @. -hᵒ[2:end-1] / dt - $diff(k * $avx(h) * $diff(ρⁱg * H + ρʷg * B) / dx) / dx

    display(b)
    display(b_frozen)
    display(b_frozen - b)

    display(norm(b, Inf))
    display(norm(b_frozen, Inf))
    display(norm(b_frozen - b, Inf)/norm(b, Inf))
    return
end

main(10)

