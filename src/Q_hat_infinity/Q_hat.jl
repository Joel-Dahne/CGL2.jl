function Q_hat_infinity(
    c10_hat::ComplexF64,
    c20_hat::ComplexF64,
    κ::Float64,
    ϵ::Float64,
    ξ₁::Float64,
    λ::CGLParams{Float64},
)
    # Compute first order approximation of Q_hat
    Q_hat = c10_hat * P_hat(ξ₁, κ, ϵ, λ) + c20_hat * E_hat(ξ₁, κ, ϵ, λ)
    dQ_hat = c10_hat * P_hat_dξ(ξ₁, κ, ϵ, λ) + c20_hat * E_hat_dξ(ξ₁, κ, ϵ, λ)

    # TODO: Compute improved approximation

    return SVector(Q_hat, dQ_hat)
end
