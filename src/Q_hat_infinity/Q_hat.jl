"""
    Q_hat_infinity(γ₁, γ₂, κ, ϵ, ξ₁, λ::CGLParams)

Compute the solution to the forward self-similar ODE on the interval
``[ξ₁, ∞)``. Returns a vector with two complex values, where the first
is the value at `ξ₁` and the second is the derivative at `ξ₁`.
"""
function Q_hat_infinity(γ₁::Acb, γ₂::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    (; d, σ) = λ
    c = _c(κ, ϵ, λ)

    # Precompute functions as well as function and norm bounds
    F = FunctionEnclosures_hat(κ, ϵ, ξ₁, λ)
    C = FunctionBounds_hat(κ, ϵ, ξ₁, λ)
    norms = NormBounds_hat(γ₁, γ₂, κ, ϵ, ξ₁, λ, C)

    # Enclosure of Q_hat
    I_E_hat = I_E_hat_enclosure(γ₁, γ₂, κ, ϵ, ξ₁, λ, F, C, norms)
    I_P_hat = zero(Acb)

    Q_hat = γ₁ * F.P_hat + γ₂ * F.E_hat + F.P_hat * I_E_hat + F.E_hat * I_P_hat

    # Enclosure of Q_hat_dξ
    I_E_hat_dξ = -F.J_E_hat * abs(Q_hat)^2σ * Q_hat
    I_P_hat_dξ = F.J_P_hat * abs(Q_hat)^2σ * Q_hat

    Q_hat_dξ =
        γ₁ * F.P_hat_dξ +
        γ₂ * F.E_hat_dξ +
        F.P_hat_dξ * I_E_hat +
        F.P_hat * I_E_hat_dξ +
        F.E_hat_dξ * I_P_hat +
        F.E_hat * I_P_hat_dξ

    return SVector(Q_hat, Q_hat_dξ)
end

Q_hat_infinity(
    γ₁::ComplexF64,
    γ₂::ComplexF64,
    κ::Float64,
    ϵ::Float64,
    ξ₁::Float64,
    λ::CGLParams{Float64},
) = ComplexF64.(
    Q_hat_infinity(Acb(γ₁), Acb(γ₂), Arb(κ), Arb(ϵ), Arb(ξ₁), CGLParams{Arb}(λ)),
)


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
    (; d, σ) = λ
    c = _c(κ, ϵ, λ)

    # Precompute functions as well as function and norm bounds
    F = FunctionEnclosures_hat(κ, ϵ, ξ₁, λ)
    C = FunctionBounds_hat(κ, ϵ, ξ₁, λ)
    norms = NormBounds_hat(γ₁, γ₂, κ, ϵ, ξ₁, λ, C)

    # Enclosure of Q_hat and Q_hat_dγ₂
    I_E_hat = I_E_hat_enclosure(γ₁, γ₂, κ, ϵ, ξ₁, λ, F, C, norms)
    I_E_hat_dγ₂ = I_E_hat_dγ₂_enclosure(γ₁, γ₂, κ, ϵ, ξ₁, λ, F, C, norms)

    I_P_hat = zero(Acb)
    I_P_hat_dγ₂ = zero(Acb)

    Q_hat = γ₁ * F.P_hat + γ₂ * F.E_hat + F.P_hat * I_E_hat + F.E_hat * I_P_hat
    Q_hat_dγ₂ = F.E_hat + F.P_hat * I_E_hat_dγ₂ + F.E_hat * I_P_hat_dγ₂

    # Enclosure of Q_hat_dξ_dγ₂
    I_E_hat_dξ = -F.J_E_hat * abs(Q_hat)^2σ * Q_hat
    I_P_hat_dξ = F.J_P_hat * abs(Q_hat)^2σ * Q_hat

    I_E_hat_dξ_dγ₂ =
        -F.J_E_hat *
        abs(Q_hat)^(2σ - 2) *
        (2σ * real(conj(Q_hat) * Q_hat_dγ₂) * Q_hat + abs(Q_hat)^2 * Q_hat_dγ₂)
    I_P_hat_dξ_dγ₂ =
        F.J_P_hat *
        abs(Q_hat)^(2σ - 2) *
        (2σ * real(conj(Q_hat) * Q_hat_dγ₂) * Q_hat + abs(Q_hat)^2 * Q_hat_dγ₂)

    Q_hat_dξ_dγ₂ =
        F.E_hat_dξ +
        F.P_hat_dξ * I_E_hat_dγ₂ +
        F.P_hat * I_E_hat_dξ_dγ₂ +
        F.E_hat_dξ * I_P_hat_dγ₂ +
        F.E_hat * I_P_hat_dξ_dγ₂

    return SMatrix{2,1}(Q_hat_dγ₂, Q_hat_dξ_dγ₂)
end
