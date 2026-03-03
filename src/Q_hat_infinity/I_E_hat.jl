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
    c = _c(κ, ϵ, λ)
    (; d, σ) = λ

    @assert v > 0

    C_I_E_hat = C.J_E_hat / abs((2σ + 1) * v - 2)

    exponent = (2σ + 1) * v - 2

    I_E_hat_bound = C_I_E_hat * ξ₁^exponent * norms.Q_hat^(2σ + 1)

    return add_error(zero(Acb), I_E_hat_bound)
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
    c = _c(κ, ϵ, λ)
    (; d, σ) = λ

    @assert v > 0

    C_I_E_hat = C.J_E_hat / abs((2σ + 1) * v - 2)

    exponent = (2σ + 1) * v - 2

    I_E_hat_dγ₂_bound =
        (2σ + 1) * C_I_E_hat * ξ₁^exponent * norms.Q_hat^2σ * norms.Q_hat_dγ₂

    return add_error(zero(Acb), I_E_hat_dγ₂_bound)
end
