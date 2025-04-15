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
