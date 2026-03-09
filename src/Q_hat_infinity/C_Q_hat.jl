function C_Q_hat(γ₁::Acb, γ₂::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    # TODO: Add checks for parameters
    a, b, c = _abc(κ, ϵ, λ)
    (; d, σ, δ) = λ

    #TODO: Does this need to be the same as in Q_hat_infinity? Probably not?
    v = Arb("0.0")

    C = FunctionBounds_hat(κ, ϵ, ξ₁, λ)
    norms = NormBounds_hat(γ₁, γ₂, κ, ϵ, ξ₁, v, λ, C)

    p_Q_hat = CGL2.p_Q_hat(γ₁, κ, ϵ, λ)

    C_I_E_hat = C.J_P_hat / abs((2σ + 1) * v - 2)
    C_I_P_hat = C.J_P_hat / 2real(c)
    C_R_Q_hat =
        abs(γ₁) * C.R_P_hat * ξ₁^(-2 / σ + d - 2) +
        abs(γ₂) * C.E_hat * exp(-real(c) * ξ₁^2) +
        (C.P_hat * C_I_E_hat + C.E_hat * C_I_P_hat * ξ₁^-2) *
        ξ₁^(-2 / σ + (2σ + 1) * v + d - 2) *
        norms.Q_hat^(2σ + 1)

    return abs(p_Q_hat) + C_R_Q_hat * ξ₁^(2 / σ - d)
end

function C_Q_hat_dξ(γ₁::Acb, γ₂::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    #TODO: Does this need to be the same as in Q_hat_infinity? Probably not?
    v = Arb("0.0")

    C = FunctionBounds_hat(κ, ϵ, ξ₁, λ)
    norms = NormBounds_hat(γ₁, γ₂, κ, ϵ, ξ₁, v, λ, C)

    # FIXME: This is only an approximation

    (; d, σ, δ) = λ
    Q_hat, Q_hat_dξ = Q_hat_infinity(γ₁, γ₂, κ, ϵ, ξ₁, λ)

    return 1.1abs(Q_hat_dξ) / ξ₁^(-1 / σ - 1)
end
