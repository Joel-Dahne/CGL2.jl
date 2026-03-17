"""
    p_Q_hat(γ₁::Acb, κ::Arb, ϵ::Arb, λ::CGLParams{T})

Compute `p_Q_hat` from Lemma REF(lemma:Q-hat-leading-term), giving the
leading asymptotic behavior of `Q_hat`.
"""
function p_Q_hat(γ₁, κ, ϵ, λ::CGLParams{T}) where {T}
    a, b, c = CGL2._abc(κ, ϵ, λ)
    return γ₁ * (-c)^-a
end

"""
    C_R_Q_hat(γ₁, γ₂, κ, ϵ, ξ₁, λ, F, C, norms)

Compute `C_R_Q_hat` from Lemma REF(lemma:Q-hat-leading-term).
"""
function C_R_Q_hat(
    γ₁::Acb,
    γ₂::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    F::FunctionEnclosures_hat,
    C::FunctionBounds_hat,
    norms::NormBounds_hat,
)
    c = _c(κ, ϵ, λ)
    (; d, σ) = λ

    # Requirement of Lemma REF(lemma:Q-hat-leading-term)
    # The requirements from Lemma REF(lemma:I_E_hat-I_P_hat-bounds)
    # are checked in the computation of `C`, where the associated
    # constants are computed.
    @assert -2real(c) + (2 / σ - d + 2) * ξ₁^-2 < 0

    return abs(γ₂) * C.E_hat * exp(-real(c) * ξ₁^2) * ξ₁^(2 / σ - d + 2) +
           C.T_hat * norms.Q_hat^(2σ + 1)
end

"""
    C_Q_hat(γ₁, γ₂, κ, ϵ, ξ₁, λ)

Compute the constant `C_Q_hat` from Lemma REF(lemma:Q-hat-bounds).

This function is used in the computation of `Y_infinity` and for that
reason it doesn't take any precomputed values, but computes them by
itself.
"""
function C_Q_hat(γ₁::Acb, γ₂::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    C = FunctionBounds_hat(κ, ϵ, ξ₁, λ)
    norms = NormBounds_hat(γ₁, γ₂, κ, ϵ, ξ₁, λ, C)
    return norms.Q_hat
end

"""
    C_Q_hat_dξ(γ₁, γ₂, κ, ϵ, ξ₁, λ)

Compute the constant `C_Q_hat_dξ` from Lemma REF(lemma:Q-hat-bounds).

This function is used in the computation of `Y_infinity` and for that
reason it doesn't take any precomputed values, but computes them by
itself.
"""
function C_Q_hat_dξ(γ₁::Acb, γ₂::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    c = _c(κ, ϵ, λ)
    (; d, σ) = λ

    @assert -2real(c) + (2 / σ - d + 2) * ξ₁^-2 < 0

    C = FunctionBounds_hat(κ, ϵ, ξ₁, λ)
    norms = NormBounds_hat(γ₁, γ₂, κ, ϵ, ξ₁, λ, C)

    return abs(γ₁) * C.P_hat_dξ +
           abs(γ₂) * C.E_hat_dξ * exp(-real(c) * ξ₁^2) * ξ₁^(2 / σ - d + 2) +
           (
               C.P_hat_dξ * C.I_E_hat +
               C.P_hat * C.J_E_hat +
               C.E_hat_dξ * C.I_P_hat +
               C.E_hat * C.J_P_hat
           ) *
           ξ₁^-2 *
           norms.Q_hat^(2σ + 1)
end
