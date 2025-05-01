"""
H(x, c, ν, κ, ϵ, ξ₁, λ::CGLParams)

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
    x::Acb,
    c::SVector{2,Acb},
    lambda::Acb,
    ν::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
)
    Y_0 = Y_zero(x, lambda, ν, κ, ϵ, ξ₁, λ)
    Y_inf = Y_infinity(c, lambda, ν, κ, ϵ, ξ₁, λ)

    return Y_0 - Y_inf
end

function H(
    x::Complex{T},
    c::SVector{2,Complex{T}},
    lambda::Complex{T},
    ν::Complex{T},
    κ::T,
    ϵ::T,
    ξ₁::T,
    λ::CGLParams{T},
) where {T}
    Q_hat = CGL2.Q_hat_zero_float_curve(real(ν), imag(ν), κ, ϵ, ξ₁, λ)
    Q_hat_ξ₁ = complex(Q_hat(ξ₁)[1:2]...)

    return H(x, c, lambda, κ, ϵ, ξ₁, Q_hat, Q_hat_ξ₁, λ)
end

"""
H(x, c, κ, ϵ, ξ₁, Q_hat, Q_hat_ξ₁, λ::CGLParams)

Version of `H` which uses precomputed values for `Q_hat` and `Q_hat_ξ₁`.
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
