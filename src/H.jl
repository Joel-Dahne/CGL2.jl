"""
H(x, c, κ, ϵ, ξ₁, Q_hat_ξ₁, λ::CGLParams)

Compute
```
H(ξ) = Y_0(ξ) - Y_inf(ξ)
```
where `Y_0` is given by [`Y_zero`](@ref) and `Y_inf` by
[`Y_infinity`](@ref). This function returns a vector with two complex
vectors, the first is the value of `H` at `ξ₁` and the second is the
derivative.
"""
function H(
    x::Complex{T},
    c::SVector{2,Complex{T}},
    lambda::Complex{T},
    κ::T,
    ϵ::T,
    ξ₁::T,
    Q_hat,
    Q_hat_ξ₁::Complex{T},
    λ::CGLParams{T},
) where {T}
    Y_0 = Y_zero(x, lambda, κ, ϵ, ξ₁, Q_hat, λ)
    Y_inf = Y_infinity(c, lambda, κ, ϵ, ξ₁, Q_hat_ξ₁, λ)

    return Y_0 - Y_inf
end
