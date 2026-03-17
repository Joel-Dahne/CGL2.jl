"""
    Q_hat_infinity(γ₁, γ₂, κ, ϵ, ξ₁, λ::CGLParams)

Compute the solution to the ODE on the interval ``[ξ₁, ∞)``. Returns a
vector with two complex values, where the first is the value at `ξ₁`
and the second is the derivative.
"""
function Q_hat_infinity(γ₁::Acb, γ₂::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    (; d, σ) = λ
    _, _, c = _abc(κ, ϵ, λ)
    v = Arb("0") # TODO: Prove that we can take v = 0

    # Precompute functions as well as function and norm bounds
    F = FunctionEnclosures_hat(κ, ϵ, ξ₁, λ)
    C = FunctionBounds_hat(κ, ϵ, ξ₁, λ)
    norms = NormBounds_hat(γ₁, γ₂, κ, ϵ, ξ₁, v, λ, C)

    # Enclosure of Q_hat
    I_E_hat = I_E_hat_enclosure(γ₁, γ₂, κ, ϵ, ξ₁, v, λ, F, C, norms)
    I_P_hat = zero(Acb)

    Q_hat = γ₁ * F.P_hat + γ₂ * F.E_hat + F.P_hat * I_E_hat + F.E_hat * I_P_hat

    # Enclosure of dQ_hat
    I_E_hat_dξ = -F.J_E_hat * abs(Q_hat)^2σ * Q_hat
    I_P_hat_dξ = F.J_P_hat * abs(Q_hat)^2σ * Q_hat

    dQ_hat =
        γ₁ * F.P_hat_dξ +
        γ₂ * F.E_hat_dξ +
        F.P_hat_dξ * I_E_hat +
        F.P_hat * I_E_hat_dξ +
        F.E_hat_dξ * I_P_hat +
        F.E_hat * I_P_hat_dξ

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
    a, b, c = _abc(κ, ϵ, λ)
    (; d, σ) = λ

    P_hat, E_hat = CGL2.P_hat(ξ₁, κ, ϵ, λ), CGL2.E_hat(ξ₁, κ, ϵ, λ)
    P_hat_dξ, E_hat_dξ = CGL2.P_hat_dξ(ξ₁, κ, ϵ, λ), CGL2.E_hat_dξ(ξ₁, κ, ϵ, λ)

    # Approximation of Q_hat
    p_Q_hat = CGL2.p_Q_hat(γ₁, κ, ϵ, λ)
    p_J_E_hat = B_W_hat(κ, ϵ, λ) * c^(a - b)
    I_E_hat = abs(p_Q_hat)^2 * p_Q_hat * p_J_E_hat / abs(-4real(a)) * ξ₁^(-4real(a))
    I_P_hat = zero(I_E_hat)

    Q_hat = γ₁ * P_hat + γ₂ * E_hat + P_hat * I_E_hat + E_hat * I_P_hat

    # Approximation of dQ_hat
    I_E_hat_dξ = -J_E_hat(ξ₁, κ, ϵ, λ) * abs(Q_hat)^2σ * Q_hat
    I_P_hat_dξ = J_P_hat(ξ₁, κ, ϵ, λ) * abs(Q_hat)^2σ * Q_hat

    dQ_hat =
        γ₁ * P_hat_dξ +
        γ₂ * E_hat_dξ +
        P_hat_dξ * I_E_hat +
        P_hat * I_E_hat_dξ +
        E_hat_dξ * I_P_hat +
        E_hat * I_P_hat_dξ

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
    (; d, σ) = λ
    _, _, c = _abc(κ, ϵ, λ)
    v = Arb("0") # TODO: Prove that we can take v = 0

    # Precompute functions as well as function and norm bounds
    F = FunctionEnclosures_hat(κ, ϵ, ξ₁, λ)
    C = FunctionBounds_hat(κ, ϵ, ξ₁, λ)
    norms = NormBounds_hat(γ₁, γ₂, κ, ϵ, ξ₁, v, λ, C)

    # Enclosure of Q_hat and Q_hat_dγ₂
    I_E_hat = I_E_hat_enclosure(γ₁, γ₂, κ, ϵ, ξ₁, v, λ, F, C, norms)
    I_E_hat_dγ₂ = I_E_hat_dγ₂_enclosure(γ₁, γ₂, κ, ϵ, ξ₁, v, λ, F, C, norms)

    I_P_hat = zero(Acb)
    I_P_hat_dγ₂ = zero(Acb)

    Q_hat = γ₁ * F.P_hat + γ₂ * F.E_hat + F.P_hat * I_E_hat + F.E_hat * I_P_hat
    Q_hat_dγ₂ = F.E_hat + F.P_hat * I_E_hat_dγ₂ + F.E_hat * I_P_hat_dγ₂

    I_E_hat_dξ = -F.J_E_hat * abs(Q_hat)^2σ * Q_hat
    I_P_hat_dξ = F.J_P_hat * abs(Q_hat)^2σ * Q_hat

    # Enclosure of dQ_hat_dγ₂
    I_E_hat_dξ_dγ₂ =
        -F.J_E_hat *
        abs(Q_hat)^(2σ - 2) *
        (2σ * real(conj(Q_hat) * Q_hat_dγ₂) * Q_hat + abs(Q_hat)^2 * Q_hat_dγ₂)
    I_P_hat_dξ_dγ₂ =
        F.J_P_hat *
        abs(Q_hat)^(2σ - 2) *
        (2σ * real(conj(Q_hat) * Q_hat_dγ₂) * Q_hat + abs(Q_hat)^2 * Q_hat_dγ₂)

    dQ_hat_dγ₂ =
        F.E_hat_dξ +
        F.P_hat_dξ * I_E_hat_dγ₂ +
        F.P_hat * I_E_hat_dξ_dγ₂ +
        F.E_hat_dξ * I_P_hat_dγ₂ +
        F.E_hat * I_P_hat_dξ_dγ₂

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

    return SMatrix{2,1}(Q_hat_dγ₂, dQ_hat_dγ₂)
end
