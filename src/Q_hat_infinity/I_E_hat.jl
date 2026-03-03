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

    @assert v > 0

    C_I_E_hat = C.J_E_hat / abs((2σ + 1) * v - 2)

    exponent = (2σ + 1) * v - 2

    I_E_hat_bound = C_I_E_hat * ξ₁^exponent * norms.Q_hat^(2σ + 1)

    # FIXME: Don't cheat
    #I_E_hat_bound = 1e-3 * I_E_hat_bound
    #return add_error(zero(Acb), I_E_hat_bound)

    # Improved version

    @assert isone(σ)

    p_Q_hat = CGL2.p_Q_hat(γ₁, κ, ϵ, λ)
    p_J_E_hat = B_W_hat(κ, ϵ, λ) * c^(a - b) # PROVE: That this is the right one

    q1 = J_E_hat(ξ₁, κ, ϵ, λ) / ξ₁^(2a - 1)
    q2 = B_W_hat(κ, ϵ, λ) * c^(a - b)
    w1 = J_E_hat(ξ₁, κ, ϵ, λ)
    w2 = J_E_hat(ξ₁, κ, ϵ, λ) - B_W_hat(κ, ϵ, λ) * c^(a - b) * ξ₁^(2a - 1)
    #@show abs(w1) abs(w2)

    I_E_hat_main = abs(p_Q_hat)^2 * p_Q_hat * p_J_E_hat * σ / 2 * ξ₁^(-2 / σ)

    I_E_hat_remainder_1_bound = C.R_J_E_hat / abs(-2 / σ - 2) * ξ₁^(-2 / σ - 2)
    I_E_hat_remainder_1 = add_error(zero(Acb), I_E_hat_remainder_1_bound)

    C_R_Q_hat = Arb(1) # FIXME: Implement this bound
    I_E_hat_remainder_2_bound =
        C.J_E_hat *
        (
            3abs(p_Q_hat)^2 * C_R_Q_hat / abs(-2 / σ + (2σ + 1) * v - 3) +
            3abs(p_Q_hat) * C_R_Q_hat^2 / abs(-2 / σ + 2(2σ + 1) * v - 6) *
            ξ₁^((2σ + 1) * v - 3) +
            C_R_Q_hat^3 / abs(-2 / σ + 3(2σ + 1) * v - 9) * ξ₁^(2(2σ + 1) * v - 6)
        ) *
        ξ₁^(-2 / σ + (2σ + 1) * v - 3)
    I_E_hat_remainder_2 = add_error(zero(Acb), I_E_hat_remainder_2_bound)

    return I_E_hat_main +
           abs(p_Q_hat)^2 * p_Q_hat * I_E_hat_remainder_1 +
           I_E_hat_remainder_2
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
