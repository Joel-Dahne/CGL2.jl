"""
    p_Q_hat(γ₁::Acb, κ::Arb, ϵ::Arb, Λ::CGLParams{T})

Compute `p_Q_hat` from Lemma REF(lemma:Q-hat-leading-term), giving the
leading asymptotic behavior of `Q_hat`.
"""
function p_Q_hat(γ₁, κ, ϵ, Λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, Λ)
    return γ₁ * (-c)^-a
end

"""
    C_R_Q_hat(γ₁, γ₂, κ, ϵ, ξ₁, Λ, F, C, norms)

Compute `C_R_Q_hat` from Lemma REF(lemma:Q-hat-leading-term).
"""
function C_R_Q_hat(
    γ₁::Acb,
    γ₂::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    Λ::CGLParams{Arb},
    F::FunctionEnclosures_hat,
    C::FunctionBounds_hat,
    norms::NormBounds_hat,
)
    c = _c(κ, ϵ, Λ)
    (; d, σ) = Λ

    # Requirement of Lemma REF(lemma:Q-hat-leading-term)
    # The requirements from Lemma REF(lemma:I_E_hat-I_P_hat-bounds)
    # are checked in the computation of `C`, where the associated
    # constants are computed.

    return abs(γ₂) * C.E_hat * exp(-real(c) * ξ₁^2) * ξ₁^(2 / σ - d + 2) +
           C.T_hat * norms.Q_hat^(2σ + 1)
end
