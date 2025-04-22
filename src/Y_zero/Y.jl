"""
    Y_zero(x, lambda, κ, ϵ, ξ₁, Q_hat, λ::CGLParams; tol::Float64 = 1e-11)

Compute the solution to the ODE on the interval ``[0, ξ₁]``. Returns a
vector with two complex values, where the first is the value at `ξ₁`
and the second is the derivative.
"""
function Y_zero(
    x::Complex{T},
    lambda::Complex{T},
    κ::T,
    ϵ::T,
    ξ₁::T,
    Q_hat,
    λ::CGLParams{T};
    tol::Float64 = 1e-11,
) where {T}
    if T == Arb
        return Y_zero_capd(x, lambda, κ, ϵ, ξ₁, Q_hat, λ; tol)
    else
        return Y_zero_float(x, lambda, κ, ϵ, ξ₁, Q_hat, λ; tol)
    end
end
