function I_E_hat_enclosure(
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

    # Requirements of Lemma REF(lemma:I_E_hat-enclosure)
    # The requirements from Lemma REF(lemma:Q-hat-leading-term) are
    # checked internally by the function C_R_Q_hat.
    @assert isone(σ)

    p_Q_hat = CGL2.p_Q_hat(γ₁, κ, ϵ, λ)
    p_J_E_hat = B_W_hat(κ, ϵ, λ) * c^(a - b)
    I_E_hat_main = abs(p_Q_hat)^2 * p_Q_hat * p_J_E_hat * σ / 2 * ξ₁^(-2 / σ)

    R_I_E_hat_1_bound = C.R_J_E_hat / abs(-2 / σ - 2) * ξ₁^(-2 / σ - 2)
    R_I_E_hat_1 = add_error(zero(Acb), R_I_E_hat_1_bound)

    C_R_Q_hat_2 = CGL2.C_R_Q_hat_2(γ₁, γ₂, κ, ϵ, ξ₁, v, λ, F, C, norms)
    R_I_E_hat_2_bound =
        C.J_E_hat *
        (
            3abs(p_Q_hat)^2 * C_R_Q_hat_2 / abs(-2 / σ + (2σ + 1) * v - 2) +
            3abs(p_Q_hat) * C_R_Q_hat_2^2 / abs(-2 / σ + 2(2σ + 1) * v - 4) *
            ξ₁^((2σ + 1) * v - 2) +
            C_R_Q_hat_2^3 / abs(-2 / σ + 3(2σ + 1) * v - 6) * ξ₁^(2(2σ + 1) * v - 4)
        ) *
        ξ₁^(-2 / σ + (2σ + 1) * v - 2)
    R_I_E_hat_2 = add_error(zero(Acb), R_I_E_hat_2_bound)

    return I_E_hat_main + abs(p_Q_hat)^2 * p_Q_hat * R_I_E_hat_1 + R_I_E_hat_2
end

function I_E_hat_dγ₂_enclosure(
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
    (; σ) = λ

    # Requirements of Lemma REF(lemma:I_E_hat-I_P_hat-dgamma-bounds)
    # are checked in the computation of `C`, where the associated
    # constants are computed.

    I_E_hat_dγ₂_bound =
        (2σ + 1) * C.I_E_hat * ξ₁^((2σ + 1) * v - 2) * norms.Q_hat^2σ * norms.Q_hat_dγ₂

    return add_error(zero(Acb), I_E_hat_dγ₂_bound)
end
