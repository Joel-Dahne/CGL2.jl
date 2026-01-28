struct FunctionBounds_Y
    J_N::Arb
    E_1::Arb
    E_2::Arb
    P_1::Arb
    P_2::Arb
    E_1_dξ::Arb
    E_2_dξ::Arb
    P_1_dξ::Arb
    P_2_dξ::Arb
    K_1_1::Arb
    K_1_2::Arb
    K_2_1::Arb
    K_2_2::Arb
    E_1_dλ::Arb
    E_2_dλ::Arb
    P_1_dλ::Arb
    P_2_dλ::Arb
    E_1_dλ_dξ::Arb
    E_2_dλ_dξ::Arb
    P_1_dλ_dξ::Arb
    P_2_dλ_dξ::Arb
    J_E_1::Arb
    J_E_2::Arb
    J_P_1::Arb
    J_P_2::Arb
    J_E_1_dλ::Arb
    J_E_2_dλ::Arb
    J_P_1_dλ::Arb
    J_P_2_dλ::Arb
    K_1_dλ_1::Arb
    K_1_dλ_2::Arb
    K_2_dλ_1::Arb
    K_2_dλ_2::Arb

    function FunctionBounds_Y(
        lambda::Acb,
        γ₁::Acb,
        γ₂::Acb,
        κ::Arb,
        ϵ::Arb,
        ξ₁::Arb,
        λ::CGLParams{Arb};
        include_dλ::Bool = false,
    )
        a, b, c = _abc(κ, ϵ, λ)
        CU = UBounds(a - lambda / 2κ, b, -c, ξ₁, include_da = include_dλ)
        CU_conj = UBounds(conj(a) - lambda / 2κ, b, -conj(c), ξ₁, include_da = include_dλ)

        C = new(
            C_J_N(lambda, γ₁, γ₂, κ, ϵ, ξ₁, λ),
            C_E_1(lambda, κ, ϵ, ξ₁, λ, CU),
            C_E_2(lambda, κ, ϵ, ξ₁, λ, CU_conj),
            C_P_1(lambda, κ, ϵ, ξ₁, λ, CU),
            C_P_2(lambda, κ, ϵ, ξ₁, λ, CU_conj),
            C_E_1_dξ(lambda, κ, ϵ, ξ₁, λ, CU),
            C_E_2_dξ(lambda, κ, ϵ, ξ₁, λ, CU_conj),
            C_P_1_dξ(lambda, κ, ϵ, ξ₁, λ, CU),
            C_P_2_dξ(lambda, κ, ϵ, ξ₁, λ, CU_conj),
            indeterminate(κ), # C_K_1_1
            indeterminate(κ), # C_K_2_2
            indeterminate(κ), # C_K_1_1
            indeterminate(κ), # C_K_2_2
            indeterminate(κ), # C_E_1_dλ
            indeterminate(κ), # C_E_2_dλ
            indeterminate(κ), # C_P_1_dλ
            indeterminate(κ), # C_P_2_dλ
            indeterminate(κ), # C_E_1_dλ_dξ
            indeterminate(κ), # C_E_2_dλ_dξ
            indeterminate(κ), # C_P_1_dλ_dξ
            indeterminate(κ), # C_P_2_dλ_dξ
            indeterminate(κ), # J_E_1
            indeterminate(κ), # J_E_2
            indeterminate(κ), # J_P_1
            indeterminate(κ), # J_P_2
            indeterminate(κ), # J_E_1_dλ
            indeterminate(κ), # J_E_2_dλ
            indeterminate(κ), # J_P_1_dλ
            indeterminate(κ), # J_P_2_dλ
            indeterminate(κ), # K_1_dλ_1
            indeterminate(κ), # K_2_dλ_2
            indeterminate(κ), # K_1_dλ_1
            indeterminate(κ), # K_2_dλ_2
        )

        C.J_E_1[] = C_J_E_1(lambda, κ, ϵ, ξ₁, λ, C)
        C.J_E_2[] = C_J_E_2(lambda, κ, ϵ, ξ₁, λ, C)
        C.J_P_1[] = C_J_P_1(lambda, κ, ϵ, ξ₁, λ, C)
        C.J_P_2[] = C_J_P_2(lambda, κ, ϵ, ξ₁, λ, C)

        C.K_1_1[] = C_K_1_1(lambda, κ, ϵ, ξ₁, λ, C)
        C.K_1_2[] = C_K_1_2(lambda, κ, ϵ, ξ₁, λ, C)
        C.K_2_1[] = C_K_2_1(lambda, κ, ϵ, ξ₁, λ, C)
        C.K_2_2[] = C_K_2_2(lambda, κ, ϵ, ξ₁, λ, C)

        if include_dλ
            C.E_1_dλ[] = C_E_1_dλ(lambda, κ, ϵ, ξ₁, λ, CU)
            C.E_2_dλ[] = C_E_2_dλ(lambda, κ, ϵ, ξ₁, λ, CU_conj)

            C.P_1_dλ[] = C_P_1_dλ(lambda, κ, ϵ, ξ₁, λ, CU)
            C.P_2_dλ[] = C_P_2_dλ(lambda, κ, ϵ, ξ₁, λ, CU_conj)

            C.E_1_dλ_dξ[] = C_E_1_dλ_dξ(lambda, κ, ϵ, ξ₁, λ, CU)
            C.E_2_dλ_dξ[] = C_E_2_dλ_dξ(lambda, κ, ϵ, ξ₁, λ, CU_conj)

            C.P_1_dλ_dξ[] = C_P_1_dλ_dξ(lambda, κ, ϵ, ξ₁, λ, CU)
            C.P_2_dλ_dξ[] = C_P_2_dλ_dξ(lambda, κ, ϵ, ξ₁, λ, CU_conj)

            C.J_E_1_dλ[] = C_J_E_1_dλ(lambda, κ, ϵ, ξ₁, λ, C)
            C.J_E_2_dλ[] = C_J_E_2_dλ(lambda, κ, ϵ, ξ₁, λ, C)
            C.J_P_1_dλ[] = C_J_P_1_dλ(lambda, κ, ϵ, ξ₁, λ, C)
            C.J_P_2_dλ[] = C_J_P_2_dλ(lambda, κ, ϵ, ξ₁, λ, C)

            C.K_1_dλ_1[] = C_K_1_dλ_1(lambda, κ, ϵ, ξ₁, λ, C)
            C.K_1_dλ_2[] = C_K_1_dλ_2(lambda, κ, ϵ, ξ₁, λ, C)
            C.K_2_dλ_1[] = C_K_2_dλ_1(lambda, κ, ϵ, ξ₁, λ, C)
            C.K_2_dλ_2[] = C_K_2_dλ_2(lambda, κ, ϵ, ξ₁, λ, C)
        end

        return C
    end
end

function C_J_N(lambda::Acb, γ₁::Acb, γ₂::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    (; σ, δ) = λ

    # FIXME
    Q_hat, _ = Q_hat_infinity(γ₁, γ₂, κ, ϵ, ξ₁, λ)
    a, b = reim(Q_hat)
    C_ab = 1.1max(abs(a), abs(b)) / ξ₁^(-1 / σ)

    return 2^(σ - 1) *
           C_ab^2σ *
           (
               1 +
               abs(δ) * (1 + 2σ) +
               2σ * (1 + abs(δ)) +
               max(1 + 2σ * abs(δ), 1 + abs(δ) * (1 + 2σ))
           )
end

function C_E_1(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU::UBounds)
    a, b, c = _abc(κ, ϵ, λ)
    return CU.U_bma_b * abs(c^(-b + a - lambda / 2κ))
end

function C_E_2(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU_conj::UBounds)
    a, b, c = _abc(κ, ϵ, λ)
    return CU_conj.U_bma_b * abs(conj(c)^(-b + conj(a) - lambda / 2κ))
end

function C_P_1(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU::UBounds)
    a, b, c = _abc(κ, ϵ, λ)
    return CU.U_a_b * abs((-c)^(-a + lambda / 2κ))
end

function C_P_2(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU_conj::UBounds)
    a, b, c = _abc(κ, ϵ, λ)
    return CU_conj.U_a_b * abs(conj(-c)^(-conj(a) + lambda / 2κ))
end

function C_E_1_dξ(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU::UBounds)
    a, b, c = _abc(κ, ϵ, λ)
    return 2(CU.U_bma_b * abs(c) + CU.U_dz_bma_b * ξ₁^-2) * abs(c^(-b + a - lambda / 2κ))
end

function C_E_2_dξ(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU_conj::UBounds)
    a, b, c = _abc(κ, ϵ, λ)
    return 2(CU_conj.U_bma_b * abs(conj(c)) + CU_conj.U_dz_bma_b * ξ₁^-2) *
           abs(conj(c)^(-b + conj(a) - lambda / 2κ))
end

function C_P_1_dξ(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU::UBounds)
    a, b, c = _abc(κ, ϵ, λ)
    return 2CU.U_dz_a_b * abs((-c)^(-a + lambda / 2κ))
end

function C_P_2_dξ(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU_conj::UBounds)
    a, b, c = _abc(κ, ϵ, λ)
    return 2CU_conj.U_dz_a_b * abs((-conj(c))^(-conj(a) + lambda / 2κ))
end

function C_E_1_dλ(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU::UBounds)
    a, b, c = _abc(κ, ϵ, λ)
    return inv(2κ) *
           CU.U_da_bma_b *
           (2 + abs(log(c)) / log(ξ₁)) *
           abs(c^(-b + a - lambda / 2κ))
end

function C_E_2_dλ(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU_conj::UBounds)
    a, b, c = _abc(κ, ϵ, λ)
    return inv(2κ) *
           CU_conj.U_da_bma_b *
           (2 + abs(log(conj(c))) / log(ξ₁)) *
           abs(conj(c)^(-b + conj(a) - lambda / 2κ))
end

function C_P_1_dλ(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU::UBounds)
    a, b, c = _abc(κ, ϵ, λ)
    return inv(2κ) *
           CU.U_da_a_b *
           (2 + abs(log(-c)) / log(ξ₁)) *
           abs((-c)^(-a + lambda / 2κ))
end

function C_P_2_dλ(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU_conj::UBounds)
    a, b, c = _abc(κ, ϵ, λ)
    return inv(2κ) *
           CU_conj.U_da_a_b *
           (2 + abs(log(-conj(c))) / log(ξ₁)) *
           abs((-conj(c))^(-conj(a) + lambda / 2κ))
end

function C_E_1_dλ_dξ(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU::UBounds)
    a, b, c = _abc(κ, ϵ, λ)
    return inv(κ) *
           (CU.U_da_a_b * abs(c) + CU.U_da_dz_a_b * ξ₁^-2) *
           (2 + abs(log(c)) / log(ξ₁)) *
           abs(c^(-b + a - lambda / 2κ))
end

function C_E_2_dλ_dξ(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    CU_conj::UBounds,
)
    a, b, c = _abc(κ, ϵ, λ)
    return inv(κ) *
           (CU_conj.U_da_a_b * abs(conj(c)) + CU_conj.U_da_dz_a_b * ξ₁^-2) *
           (2 + abs(log(conj(c))) / log(ξ₁)) *
           abs(conj(c)^(-b + conj(a) - lambda / 2κ))
end

function C_P_1_dλ_dξ(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU::UBounds)
    a, b, c = _abc(κ, ϵ, λ)
    return inv(κ) *
           CU.U_da_dz_a_b *
           (2 + abs(log(-c)) / log(ξ₁)) *
           abs((-c)^(-a + lambda / 2κ))
end

function C_P_2_dλ_dξ(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    CU_conj::UBounds,
)
    a, b, c = _abc(κ, ϵ, λ)
    return inv(κ) *
           CU_conj.U_da_dz_a_b *
           (2 + abs(log(-conj(c))) / log(ξ₁)) *
           abs((-conj(c))^(-conj(a) + lambda / 2κ))
end

function C_J_E_1(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    return abs(B_W_1(lambda, κ, ϵ, λ)) * C.E_1
end

function C_J_E_2(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    return abs(B_W_2(lambda, κ, ϵ, λ)) * C.E_2
end

function C_J_P_1(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    return abs(B_W_1(lambda, κ, ϵ, λ)) * C.P_1
end

function C_J_P_2(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    return abs(B_W_2(lambda, κ, ϵ, λ)) * C.P_2
end

function C_J_E_1_dλ(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    return abs(B_W_1_dλ(lambda, κ, ϵ, λ)) * C.E_1 * inv(log(ξ₁)) +
           abs(B_W_1(lambda, κ, ϵ, λ)) * C.E_1_dλ
end

function C_J_E_2_dλ(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    return abs(B_W_2_dλ(lambda, κ, ϵ, λ)) * C.E_2 * inv(log(ξ₁)) +
           abs(B_W_2(lambda, κ, ϵ, λ)) * C.E_2_dλ
end

function C_J_P_1_dλ(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    return abs(B_W_1_dλ(lambda, κ, ϵ, λ)) * C.P_1 * inv(log(ξ₁)) +
           abs(B_W_1(lambda, κ, ϵ, λ)) * C.P_1_dλ
end

function C_J_P_2_dλ(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    return abs(B_W_2_dλ(lambda, κ, ϵ, λ)) * C.P_2 * inv(log(ξ₁)) +
           abs(B_W_2(lambda, κ, ϵ, λ)) * C.P_2_dλ
end

function C_K_1_1(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    return inv(sqrt(1 + ϵ^2)) * C.J_P_1
end

function C_K_1_2(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    return inv(sqrt(1 + ϵ^2)) * C.J_P_2
end

function C_K_2_1(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    return inv(sqrt(1 + ϵ^2)) * C.J_E_1
end

function C_K_2_2(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    return inv(sqrt(1 + ϵ^2)) * C.J_E_2
end

function C_K_1_dλ_1(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    return inv(sqrt(1 + ϵ^2)) * C.J_P_1_dλ
end

function C_K_1_dλ_2(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    return inv(sqrt(1 + ϵ^2)) * C.J_P_2_dλ
end

function C_K_2_dλ_1(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    return inv(sqrt(1 + ϵ^2)) * C.J_E_1_dλ
end

function C_K_2_dλ_2(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    return inv(sqrt(1 + ϵ^2)) * C.J_E_2_dλ
end
