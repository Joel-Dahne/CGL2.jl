"""
    G_hat(ν, γ₁, γ₂, κ, ϵ, ξ₁, Λ::CGLParams)

Compute
```
G_hat(ξ) = Q_hat_0(ξ) - Q_hat_inf(ξ)
```
where `Q_hat_0` is given by [`Q_hat_zero`](@ref) and `Q_hat_inf` by
[`Q_hat_infinity`](@ref). This function returns a vector with four
real values, the first two are the real and imaginary values of
`G_hat` at `ξ₁` and the last two the real and imaginary values of its
derivative.
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

    return SVector(
        real(Q_hat_0) - real(Q_hat_inf),
        imag(Q_hat_0) - imag(Q_hat_inf),
        real(Q_hat_0_dξ) - real(Q_hat_inf_dξ),
        imag(Q_hat_0_dξ) - imag(Q_hat_inf_dξ),
    )
end

"""
    G_hat_jacobian(ν, γ₁, γ₂, κ, ϵ, ξ₁, Λ::CGLParams)

This function computes the Jacobian of [`G_hat`](@ref) w.r.t. the real
and imaginary parts of `ν` and `γ₂`.
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

    return SMatrix{4,4}(
        # Derivative w.r.t. real(ν)
        real(Q_hat_0_J[1, 1]),
        imag(Q_hat_0_J[1, 1]),
        real(Q_hat_0_J[2, 1]),
        imag(Q_hat_0_J[2, 1]),
        # Derivative w.r.t. imag(ν)
        real(Q_hat_0_J[1, 2]),
        imag(Q_hat_0_J[1, 2]),
        real(Q_hat_0_J[2, 2]),
        imag(Q_hat_0_J[2, 2]),
        # Derivative w.r.t. real(γ₁)
        real(Q_hat_inf_J[1, 1]),
        imag(Q_hat_inf_J[1, 1]),
        real(Q_hat_inf_J[2, 1]),
        imag(Q_hat_inf_J[2, 1]),
        # Derivative w.r.t. imag(γ₁)
        real(Q_hat_inf_J[1, 2]),
        imag(Q_hat_inf_J[1, 2]),
        real(Q_hat_inf_J[2, 2]),
        imag(Q_hat_inf_J[2, 2]),
    )
end
