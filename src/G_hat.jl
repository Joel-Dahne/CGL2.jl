"""
    G_hat(ν, γ₁, γ₂, κ, ϵ, ξ₁, λ::CGLParams)

Compute
```
G_hat(ξ) = Q_hat_0(ξ) - Q_hat_inf(ξ)
```
where `Q_hat_0` is given by [`Q_hat_zero`](@ref) and `Q_hat_inf` by
[`Q_hat_infinity`](@ref). This function returns a vector with two
complex values, the first is the value of `G_hat` at `ξ₁` and the
second is its derivative.
"""
function G_hat(
    ν::Complex{T}, # TODO: Should support Acb later on
    γ₁::Complex{T},
    γ₂::Complex{T},
    κ::T,
    ϵ::T,
    ξ₁::T,
    λ::CGLParams{T},
) where {T}
    Q_hat_0, Q_hat_0_dξ = Q_hat_zero(ν, κ, ϵ, ξ₁, λ)
    Q_hat_inf, Q_hat_inf_dξ = Q_hat_infinity(γ₁, γ₂, κ, ϵ, ξ₁, λ)

    return SVector(Q_hat_0 - Q_hat_inf, Q_hat_0_dξ - Q_hat_inf_dξ)
end
