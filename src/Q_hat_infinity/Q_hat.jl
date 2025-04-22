"""
    Q_hat_infinity(γ₁, γ₂, κ, ϵ, ξ₁, λ::CGLParams)

Compute the solution to the ODE on the interval ``[ξ₁, ∞)``. Returns a
vector with two complex values, where the first is the value at `ξ₁`
and the second is the derivative.
"""
function Q_hat_infinity(γ₁::Acb, γ₂::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    # FIXME: Implement rigorous enclosure

    # Compute first order approximation of Q_hat
    Q_hat = γ₁ * P_hat(ξ₁, κ, ϵ, λ) + γ₂ * E_hat(ξ₁, κ, ϵ, λ)
    dQ_hat = γ₁ * P_hat_dξ(ξ₁, κ, ϵ, λ) + γ₂ * E_hat_dξ(ξ₁, κ, ϵ, λ)

    Q_hat = add_error(Q_hat, Mag(1e-8))
    dQ_hat = add_error(dQ_hat, Mag(1e-8))

    return SVector(Q_hat, dQ_hat)
end

function Q_hat_infinity(
    γ₁::ComplexF64,
    γ₂::ComplexF64,
    κ::Float64,
    ϵ::Float64,
    ξ₁::Float64,
    λ::CGLParams{Float64},
)
    # Compute first order approximation of Q_hat
    Q_hat = γ₁ * P_hat(ξ₁, κ, ϵ, λ) + γ₂ * E_hat(ξ₁, κ, ϵ, λ)
    dQ_hat = γ₁ * P_hat_dξ(ξ₁, κ, ϵ, λ) + γ₂ * E_hat_dξ(ξ₁, κ, ϵ, λ)

    # TODO: Compute improved approximation

    return SVector(Q_hat, dQ_hat)
end

"""
    Q_hat_infinity_jacobian(γ₁, γ₂, κ, ϵ, ξ₁, λ::CGLParams)

This function computes the Jacobian of [`Q_hat_infinity`](@ref) w.r.t. the
parameter `γ₂`.
"""
function Q_hat_infinity_jacobian(
    γ₁::Acb,
    γ₂::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
)
    # FIXME: Implement rigorous enclosure

    # Compute first order approximation
    Q_hat_dγ₂ = E_hat(ξ₁, κ, ϵ, λ)
    dQ_hat_dγ₂ = E_hat_dξ(ξ₁, κ, ϵ, λ)

    Q_hat_dγ₂ = add_error(Q_hat_dγ₂, Mag(1e-8))
    dQ_hat_dγ₂ = add_error(dQ_hat_dγ₂, Mag(1e-8))

    return SMatrix{2,1}(Q_hat_dγ₂, dQ_hat_dγ₂)
end

function Q_hat_infinity_jacobian(
    γ₁::ComplexF64,
    γ₂::ComplexF64,
    κ::Float64,
    ϵ::Float64,
    ξ₁::Float64,
    λ::CGLParams{Float64},
)
    # Compute first order approximation
    Q_hat_dγ₂ = E_hat(ξ₁, κ, ϵ, λ)
    dQ_hat_dγ₂ = E_hat_dξ(ξ₁, κ, ϵ, λ)

    # TODO: Compute improved approximation

    return SMatrix{2,1}(Q_hat_dγ₂, dQ_hat_dγ₂)
end
