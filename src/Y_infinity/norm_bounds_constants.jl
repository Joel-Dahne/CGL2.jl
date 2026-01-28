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

function C_T_12(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C_Y::FunctionBounds_Y_new,
)
    (; d, σ) = λ
    _, _, c = _abc(κ, ϵ, λ)

    @assert v > 0
    @assert real(c) > 0
    return 2max(C_Y.E_1 * C_I_K_1_1(v, C_Y), C_Y.E_2 * C_I_K_1_2(v, C_Y)) +
           2max(
        C_Y.P_1 * C_I_K_2_1(lambda, κ, ϵ, v, λ, C_Y),
        C_Y.P_2 * C_I_K_2_2(lambda, κ, ϵ, v, λ, C_Y),
    ) * ξ₁^-2
end

function C_Y_dλ_1(v::Arb, C_Y::FunctionBounds_Y)
    return 2C_Y.Y_12 * exp(Arb(-1)) / v
end

function C_Y_dλ_1(c::SVector{2,Acb}, v::Arb, C_Y::FunctionBounds_Y_new)
    return 2max(C_Y.E_1_dλ * abs(c[1]), C_Y.E_2_dλ * abs(c[2])) * exp(Arb(-1)) / v
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

function C_Y_dλ_2(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C_Y::FunctionBounds_Y_new,
)
    return 2max(C_Y.E_1_dλ * C_I_K_1_1(v, C_Y), C_Y.E_2_dλ * C_I_K_1_2(v, C_Y)) *
           exp(Arb(-1)) / v * ξ₁^(v - 2) +
           2max(
               C_Y.P_1_dλ * C_I_K_2_1(lambda, κ, ϵ, v, λ, C_Y),
               C_Y.P_2_dλ * C_I_K_2_2(lambda, κ, ϵ, v, λ, C_Y),
           ) *
           log(ξ₁) *
           ξ₁^-4 +
           2max(
               C_Y.E_1 * C_I_K_1_dλ_1_1(ξ₁, v, C_Y),
               C_Y.E_2 * C_I_K_1_dλ_1_2(ξ₁, v, C_Y),
           ) *
           log(ξ₁) *
           ξ₁^-2 +
           2max(
               C_Y.P_1 * C_I_K_2_dλ_1_1(lambda, κ, ϵ, ξ₁, v, λ, C_Y),
               C_Y.P_2 * C_I_K_2_dλ_1_2(lambda, κ, ϵ, ξ₁, v, λ, C_Y),
           ) *
           log(ξ₁) *
           ξ₁^-4
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

function C_Y_dλ_3(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C_Y::FunctionBounds_Y_new,
)
    return 2max(C_Y.E_1 * C_I_K_1_dλ_2_1(v, C_Y), C_Y.E_2 * C_I_K_1_dλ_2_2(v, C_Y)) *
           ξ₁^-2 +
           2max(
        C_Y.P_1 * C_I_K_2_dλ_2_1(lambda, κ, ϵ, v, λ, C_Y),
        C_Y.P_2 * C_I_K_2_dλ_2_2(lambda, κ, ϵ, v, λ, C_Y),
    ) * ξ₁^-4
end
