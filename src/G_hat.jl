"""
    G_hat(ν, γ₁, γ₂, κ, ϵ, ξ₁, Λ::CGLParams)

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
    ν::Union{Complex{T},Acb},
    γ₁::Union{Complex{T},Acb},
    γ₂::Union{Complex{T},Acb},
    κ::T,
    ϵ::T,
    ξ₁::T,
    Λ::CGLParams{T},
) where {T}
    Q_hat_0, Q_hat_0_dξ = Q_hat_zero(ν, κ, ϵ, ξ₁, Λ)
    Q_hat_inf, Q_hat_inf_dξ = Q_hat_infinity(γ₁, γ₂, κ, ϵ, ξ₁, Λ)

    return SVector(Q_hat_0 - Q_hat_inf, Q_hat_0_dξ - Q_hat_inf_dξ)
end

"""
    G_hat_jacobian(ν, γ₁, γ₂, κ, ϵ, ξ₁, Λ::CGLParams)

This function computes the Jacobian of [`G_hat`](@ref) w.r.t. the
parameters `ν` and `γ₂`.
"""
function G_hat_jacobian(
    ν::Union{Complex{T},Acb},
    γ₁::Union{Complex{T},Acb},
    γ₂::Union{Complex{T},Acb},
    κ::T,
    ϵ::T,
    ξ₁::T,
    Λ::CGLParams{T},
) where {T}
    Q_hat_0_J = Q_hat_zero_jacobian(ν, κ, ϵ, ξ₁, Λ)
    Q_hat_inf_J = Q_hat_infinity_jacobian(γ₁, γ₂, κ, ϵ, ξ₁, Λ)

    return SMatrix{2,2}(
        # Derivative w.r.t. ν
        Q_hat_0_J[1, 1],
        Q_hat_0_J[2, 1],
        # Derivative w.r.t. γ₁
        Q_hat_inf_J[1, 1],
        Q_hat_inf_J[2, 1],
    )
end
