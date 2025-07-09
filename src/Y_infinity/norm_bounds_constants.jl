function C_T_Y(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C_Y::FunctionBounds_Y,
)
    (; d, σ) = λ
    real_a = real_a12(κ, ϵ)

    @assert v > 0
    @assert real_a < 0

    return 2C_Y.Y_12 * C_I_K_1(v, C_Y) + 2C_Y.Y_34 * C_I_K_2(lambda, κ, v, λ, C_Y)
end

function C_Y_dλ_1(v::Arb, C_Y::FunctionBounds_Y)
    return 2C_Y.Y_12 * exp(Arb(-1)) / v
end

function C_Y_dλ_2(
    lambda::Acb,
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C_Y::FunctionBounds_Y,
)
    return 2exp(Arb(-1)) / v *
           (C_Y.Y_12_dλ * C_I_K_1(v, C_Y) + C_Y.Y_34_dλ * C_I_K_2(lambda, κ, v, λ, C_Y)) *
           ξ₁^(v - 2)
end

function C_Y_dλ_3(
    lambda::Acb,
    κ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C_Y::FunctionBounds_Y,
)
    return 2(
        C_Y.Y_12 *
        (C_I_K_1_dλ_1_1(v, C_Y) * log(ξ₁) + C_I_K_1_dλ_1_2(v, C_Y) + C_I_K_1_dλ_2(v, C_Y)) +
        C_Y.Y_34 * (
            C_I_K_2_dλ_1_1(lambda, κ, v, λ, C_Y) * exp(Arb(-1)) / v +
            (C_I_K_2_dλ_1_2(lambda, κ, v, λ, C_Y) + C_I_K_2_dλ_2(lambda, κ, v, λ, C_Y)) *
            ξ₁^-v
        )
    ) * ξ₁^(v - 2)
end
