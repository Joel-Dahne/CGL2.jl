function p_Q_hat(γ₁, κ, ϵ, λ::CGLParams{T}) where {T}
    a, b, c = CGL2._abc(κ, ϵ, λ)
    return γ₁ * (-c)^-a
end

function C_R_Q_hat_1(
    γ₁::Acb,
    γ₂::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    F::FunctionEnclosures_hat,
    C::FunctionBounds_hat,
    norms::NormBounds_hat,
)
    a, b, c = _abc(κ, ϵ, λ)
    (; d, σ) = λ

    # Requirement of Lemma REF(lemma:Q-hat-leading-term)
    # The requirements from Lemma REF(lemma:I_E_hat-I_P_hat-bounds)
    # are checked in the computation of `C`, where the associated
    # constants are computed.
    @assert (2 / σ - (2σ + 1) * v - d + 2) * ξ₁^-2 < 2real(c)

    return abs(γ₂) * C.E_hat * exp(-real(c) * ξ₁^2) * ξ₁^(2 / σ - (2σ + 1) * v - d + 2) +
           C.T_hat * norms.Q_hat^(2σ + 1)
end

function C_R_Q_hat_2(
    γ₁::Acb,
    γ₂::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    F::FunctionEnclosures_hat,
    C::FunctionBounds_hat,
    norms::NormBounds_hat,
)
    a, b, c = _abc(κ, ϵ, λ)
    (; d, σ) = λ

    return abs(γ₁) * C.R_P_hat * ξ₁^(-(2σ + 1) * v) +
           C_R_Q_hat_1(γ₁, γ₂, κ, ϵ, ξ₁, v, λ, F, C, norms)
end

# This is used by Y_infinity, so doesn't take any precomputed values.
function C_Q_hat(γ₁::Acb, γ₂::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    # TODO: Add checks for parameters
    a, b, c = _abc(κ, ϵ, λ)
    (; d, σ, δ) = λ

    #TODO: Can we take this to be zero?
    v = Arb("0.0")

    C = FunctionBounds_hat(v, κ, ϵ, ξ₁, λ)
    norms = NormBounds_hat(γ₁, γ₂, κ, ϵ, ξ₁, v, λ, C)

    return norms.Q_hat

    # TODO: This is an old version which doesn't need to use v = 0 to
    # work. Depending on if we can use v = 0 or not we might want it.
    p_Q_hat = CGL2.p_Q_hat(γ₁, κ, ϵ, λ)

    C_R_Q_hat =
        abs(γ₁) * C.R_P_hat * ξ₁^(-2 / σ + d - 2) +
        abs(γ₂) * C.E_hat * exp(-real(c) * ξ₁^2) +
        (C.P_hat * C.I_E_hat + C.E_hat * C.I_P_hat * ξ₁^-2) *
        ξ₁^(-2 / σ + (2σ + 1) * v + d - 2) *
        norms.Q_hat^(2σ + 1)

    return abs(p_Q_hat) + C_R_Q_hat * ξ₁^(2 / σ - d)
end

# This is used by Y_infinity, so doesn't take any precomputed values.
function C_Q_hat_dξ(γ₁::Acb, γ₂::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    # TODO: Add checks for parameters
    a, b, c = _abc(κ, ϵ, λ)
    (; d, σ) = λ

    #TODO: Can we take this to be zero?
    v = Arb("0.0")

    C = FunctionBounds_hat(v, κ, ϵ, ξ₁, λ)
    norms = NormBounds_hat(γ₁, γ₂, κ, ϵ, ξ₁, v, λ, C)

    return abs(γ₁) * C.P_hat_dξ +
           abs(γ₂) * C.E_hat_dξ * exp(-real(c) * ξ₁^2) * ξ₁^(2 / σ - d + 2) +
           (
               C.P_hat_dξ * C.I_E_hat +
               C.P_hat * C.J_E_hat +
               C.E_hat_dξ * C.I_P_hat +
               C.E_hat * C.J_P_hat
           ) *
           ξ₁^((2σ + 1) * v - 2) *
           norms.Q_hat^(2σ + 1)
end
