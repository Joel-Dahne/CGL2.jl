"""
    H(ν_real, ν_imag, μ, γ, κ, ϵ, ξ₁, λ::CGLParams)

Compute
```
H(ξ) = Q(ξ) - Q_hat(ξ)
```
where `Q` is the backwards self-similar solution and `Q_hat` is the
forward self-similar solution.
"""
function H(
    ν_real::T,
    ν_imag::T,
    μ::T,
    γ::Complex{T},
    κ::T,
    ϵ::T,
    ξ₁::T,
    λ::CGLParams{T},
) where {T}
    Q, _ = Q_zero(μ, κ, ϵ, ξ₁, λ)
    Q_hat, _ = Q_hat_zero(_complex(ν_real, ν_imag), κ, ϵ, ξ₁, λ)

    res = Q - Q_hat

    return SVector(real(res), imag(res))
end
