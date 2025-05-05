"""
    Y_zero(x, lambda, ν, κ, ϵ, ξ₁, λ::CGLParams; tol::Float64 = 1e-11)
    Y_zero(x, lambda, κ, ϵ, ξ₁, Q_hat, λ::CGLParams; tol::Float64 = 1e-11)

Compute the solution to the ODE on the interval ``[0, ξ₁]``. Returns a
vector with two complex values, where the first is the value at `ξ₁`
and the second is the derivative.
"""
function Y_zero(
    x::Acb,
    lambda::Acb,
    ν::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    tol::Float64 = 1e-11,
) where {Arb}
    return Y_zero_capd(x, lambda, ν, κ, ϵ, ξ₁, λ; tol)
end

function Y_zero(
    x::ComplexF64,
    lambda::ComplexF64,
    ν::ComplexF64,
    κ::Float64,
    ϵ::Float64,
    ξ₁::Float64,
    λ::CGLParams{Float64};
    tol::Float64 = 1e-11,
)
    Q_hat = CGL2.Q_hat_zero_float_curve(real(ν), imag(ν), κ, ϵ, ξ₁, λ)

    return Y_zero(x, lambda, κ, ϵ, ξ₁, Q_hat, λ)
end

function Y_zero(
    x::ComplexF64,
    lambda::ComplexF64,
    κ::Float64,
    ϵ::Float64,
    ξ₁::Float64,
    Q_hat,
    λ::CGLParams{Float64};
    tol::Float64 = 1e-11,
)
    return Y_zero_float(x, lambda, κ, ϵ, ξ₁, Q_hat, λ; tol)
end

"""
    Y_zero_jacobian(x, lambda, ν, κ, ϵ, ξ₁, λ::CGLParams; tol::Float64 = 1e-11)
    Y_zero_jacobian(x, lambda, κ, ϵ, ξ₁, Q_hat, λ::CGLParams; tol::Float64 = 1e-11)

This function computes the Jacobian of [`Y_zero`](@ref) w.r.t. the
parameters `x` and `lambda`.
"""
function Y_zero_jacobian(
    x::Acb,
    lambda::Acb,
    ν::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    tol::Float64 = 1e-11,
) where {Arb}
    return Y_zero_jacobian_capd(x, lambda, ν, κ, ϵ, ξ₁, λ; tol)
end

function Y_zero_jacobian(
    x::ComplexF64,
    lambda::ComplexF64,
    ν::ComplexF64,
    κ::Float64,
    ϵ::Float64,
    ξ₁::Float64,
    λ::CGLParams{Float64};
    tol::Float64 = 1e-11,
)
    Q_hat = CGL2.Q_hat_zero_float_curve(real(ν), imag(ν), κ, ϵ, ξ₁, λ)

    return Y_zero_jacobian(x, lambda, κ, ϵ, ξ₁, Q_hat, λ)
end

function Y_zero_jacobian(
    x::ComplexF64,
    lambda::ComplexF64,
    κ::Float64,
    ϵ::Float64,
    ξ₁::Float64,
    Q_hat,
    λ::CGLParams{Float64};
    tol::Float64 = 1e-11,
)
    return Y_zero_jacobian_float(x, lambda, κ, ϵ, ξ₁, Q_hat, λ; tol)
end
