"""
    FunctionBounds_Y(lambda, γ₁, γ₂, κ, ϵ, ξ₁, λ; include_dλ = false)

This contains all the bounds of functions that are needed in the
enclosure asymptotic expansion of `Y` at infinity.

More precisely it contains the bounds from

- Lemma REF(lemma:P_i_E_i-bounds)
- Lemma REF(lemma:I_K_1-I_K_2-bounds)
- Lemma REF(lemma:bound-J_N)
- Lemma REF(lemma:H-bounds)

If `include_dλ = false` it doesn't include the bounds corresponding to
derivatives in `lambda`.

It checks all the conditions on the parameters that these lemmas
assume. If any of these conditions are not satisfied it will throw an
error. When using this struct the bounds can therefore safely be
assume to hold.
"""
struct FunctionBounds_Y
    J_N::Arb
    J_N_dξ::Arb
    I_N::Arb
    I_N_dξ::Arb
    E_1::Arb
    E_2::Arb
    P_1::Arb
    P_2::Arb
    E_1_dξ::Arb
    E_2_dξ::Arb
    P_1_dξ::Arb
    P_2_dξ::Arb
    exp_E_1_dξ::Arb
    exp_E_2_dξ::Arb
    exp_P_1_dξ::Arb
    exp_P_2_dξ::Arb
    K_1_1::Arb
    K_1_2::Arb
    K_2_1::Arb
    K_2_2::Arb
    K_2_dξ_1::Arb
    K_2_dξ_2::Arb
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
    J_E_1_dξ::Arb
    J_E_2_dξ::Arb
    J_E_1_dλ::Arb
    J_E_2_dλ::Arb
    J_P_1_dλ::Arb
    J_P_2_dλ::Arb
    K_1_dλ_1::Arb
    K_1_dλ_2::Arb
    K_2_dλ_1::Arb
    K_2_dλ_2::Arb
    H_11::Arb
    H_12::Arb
    H_21::Arb
    H_22::Arb
    H_11_dξ::Arb
    H_12_dξ::Arb
    H_21_dξ::Arb
    H_22_dξ::Arb

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
        (; σ, δ) = λ
        a, b, c = _abc(κ, ϵ, λ)

        # This is the only direct condition in Lemma
        # REF(lemma:P_i_E_i-bounds), REF(lemma:I_K_1-I_K_2-bounds) and
        # REF(lemma:H-bounds).
        # The conditions related to the bounds for U are checked by
        # Ubounds.
        ξ₁ > 1 || throw(ArgumentError("ξ₁ > 1 not satisfied"))

        # Add checks for the conditions related to C_J_N once it is
        # fully implemented.
        isone(σ) || throw(ArgumentError("σ = 1 not satisfied"))
        iszero(δ) || throw(ArgumentError("δ = 0 not satisfied"))
        Q_hat, Q_hat_dξ = Q_hat_infinity(γ₁, γ₂, κ, ϵ, ξ₁, λ)
        C_Q_hat = 1.1abs(Q_hat) / ξ₁^-1
        C_Q_hat_dξ = 1.1abs(Q_hat_dξ) / ξ₁^-2

        CU = UBounds(a - lambda / 2κ, b, -c, ξ₁, include_da = include_dλ)
        CU_conj = UBounds(conj(a) - lambda / 2κ, b, -conj(c), ξ₁, include_da = include_dλ)

        C = new(
            C_J_N(lambda, γ₁, γ₂, κ, ϵ, ξ₁, λ),
            C_J_N_dξ(lambda, γ₁, γ₂, κ, ϵ, ξ₁, λ),
            C_I_N(lambda, κ, ϵ, ξ₁, λ, C_Q_hat),
            C_I_N_dξ(lambda, κ, ϵ, ξ₁, λ, C_Q_hat, C_Q_hat_dξ),
            C_E_1(lambda, κ, ϵ, ξ₁, λ, CU),
            C_E_2(lambda, κ, ϵ, ξ₁, λ, CU_conj),
            C_P_1(lambda, κ, ϵ, ξ₁, λ, CU),
            C_P_2(lambda, κ, ϵ, ξ₁, λ, CU_conj),
            C_E_1_dξ(lambda, κ, ϵ, ξ₁, λ, CU),
            C_E_2_dξ(lambda, κ, ϵ, ξ₁, λ, CU_conj),
            C_P_1_dξ(lambda, κ, ϵ, ξ₁, λ, CU),
            C_P_2_dξ(lambda, κ, ϵ, ξ₁, λ, CU_conj),
            indeterminate(κ), # C_exp_E_1_dξ
            indeterminate(κ), # C_exp_E_2_dξ
            indeterminate(κ), # C_exp_P_1_dξ
            indeterminate(κ), # C_exp_P_2_dξ
            indeterminate(κ), # C_K_1_1
            indeterminate(κ), # C_K_1_2
            indeterminate(κ), # C_K_2_1
            indeterminate(κ), # C_K_2_2
            indeterminate(κ), # C_K_2_dξ_1
            indeterminate(κ), # C_K_2_dξ_2
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
            indeterminate(κ), # J_E_1_dξ
            indeterminate(κ), # J_E_2_dξ
            indeterminate(κ), # J_E_1_dλ
            indeterminate(κ), # J_E_2_dλ
            indeterminate(κ), # J_P_1_dλ
            indeterminate(κ), # J_P_2_dλ
            indeterminate(κ), # K_1_dλ_1
            indeterminate(κ), # K_2_dλ_2
            indeterminate(κ), # K_1_dλ_1
            indeterminate(κ), # K_2_dλ_2
            indeterminate(κ), # H_11
            indeterminate(κ), # H_12
            indeterminate(κ), # H_21
            indeterminate(κ), # H_22
            indeterminate(κ), # H_11_dξ
            indeterminate(κ), # H_12_dξ
            indeterminate(κ), # H_21_dξ
            indeterminate(κ), # H_22_dξ
        )

        C.exp_E_1_dξ[] = C_exp_E_1_dξ(lambda, κ, ϵ, ξ₁, λ, CU)
        C.exp_E_2_dξ[] = C_exp_E_2_dξ(lambda, κ, ϵ, ξ₁, λ, CU_conj)
        C.exp_P_1_dξ[] = C_exp_P_1_dξ(lambda, κ, ϵ, ξ₁, λ, C)
        C.exp_P_2_dξ[] = C_exp_P_2_dξ(lambda, κ, ϵ, ξ₁, λ, C)

        C.J_E_1[] = C_J_E_1(lambda, κ, ϵ, ξ₁, λ, C)
        C.J_E_2[] = C_J_E_2(lambda, κ, ϵ, ξ₁, λ, C)
        C.J_P_1[] = C_J_P_1(lambda, κ, ϵ, ξ₁, λ, C)
        C.J_P_2[] = C_J_P_2(lambda, κ, ϵ, ξ₁, λ, C)

        C.J_E_1_dξ[] = C_J_E_1_dξ(lambda, κ, ϵ, ξ₁, λ, C)
        C.J_E_2_dξ[] = C_J_E_2_dξ(lambda, κ, ϵ, ξ₁, λ, C)

        C.K_1_1[] = C_K_1_1(lambda, κ, ϵ, ξ₁, λ, C)
        C.K_1_2[] = C_K_1_2(lambda, κ, ϵ, ξ₁, λ, C)
        C.K_2_1[] = C_K_2_1(lambda, κ, ϵ, ξ₁, λ, C)
        C.K_2_2[] = C_K_2_2(lambda, κ, ϵ, ξ₁, λ, C)

        C.K_2_dξ_1[] = C_K_2_dξ_1(lambda, κ, ϵ, ξ₁, λ, C)
        C.K_2_dξ_2[] = C_K_2_dξ_2(lambda, κ, ϵ, ξ₁, λ, C)

        C.H_11[] = C_H_11(lambda, κ, ϵ, ξ₁, λ, C, C_Q_hat)
        C.H_12[] = C_H_12(lambda, κ, ϵ, ξ₁, λ, C, C_Q_hat)
        C.H_21[] = C_H_21(lambda, κ, ϵ, ξ₁, λ, C, C_Q_hat)
        C.H_22[] = C_H_22(lambda, κ, ϵ, ξ₁, λ, C, C_Q_hat)
        C.H_11_dξ[] = C_H_11_dξ(lambda, κ, ϵ, ξ₁, λ, C, C_Q_hat, C_Q_hat_dξ)
        C.H_12_dξ[] = C_H_12_dξ(lambda, κ, ϵ, ξ₁, λ, C, C_Q_hat, C_Q_hat_dξ)
        C.H_21_dξ[] = C_H_21_dξ(lambda, κ, ϵ, ξ₁, λ, C, C_Q_hat, C_Q_hat_dξ)
        C.H_22_dξ[] = C_H_22_dξ(lambda, κ, ϵ, ξ₁, λ, C, C_Q_hat, C_Q_hat_dξ)

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

function C_J_N_dξ(lambda::Acb, γ₁::Acb, γ₂::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    (; σ, δ) = λ
    # TODO: If we assume that σ is one then the factor
    # abs2(Q_hat_ξ)^(σ - 1) doesn't play a role and the derivative is
    # much simpler. Do we need to care about the general case?
    @assert isone(σ)

    # FIXME
    Q_hat, Q_hat_dξ = Q_hat_infinity(γ₁, γ₂, κ, ϵ, ξ₁, λ)
    a, b = reim(Q_hat)
    a_dξ, b_dξ = reim(Q_hat_dξ)
    C_ab = 1.1max(abs(a), abs(b)) / ξ₁^(-1 / σ)
    C_ab_dξ = 1.1max(abs(a_dξ), abs(b_dξ)) / ξ₁^(-1 / σ - 1)

    return 2C_ab *
           C_ab_dξ *
           (
               1 +
               abs(δ) * (1 + 2σ) +
               2σ * (1 + abs(δ)) +
               max(1 + 2σ * abs(δ), 1 + abs(δ) * (1 + 2σ))
           )
end

function C_I_N(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, C_Q_hat::Arb)
    return 3C_Q_hat^2
end

function C_I_N_dξ(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C_Q_hat::Arb,
    C_Q_hat_dξ::Arb,
)
    return 6C_Q_hat * C_Q_hat_dξ
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

function C_exp_E_1_dξ(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU::UBounds)
    a, b, c = _abc(κ, ϵ, λ)
    return 2CU.U_dz_bma_b * abs(c^(-b + a - lambda / 2κ))
end

function C_exp_E_2_dξ(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    CU_conj::UBounds,
)
    a, b, c = _abc(κ, ϵ, λ)
    return 2CU_conj.U_dz_bma_b * abs(conj(c)^(-b + conj(a) - lambda / 2κ))
end

function C_exp_P_1_dξ(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    a, b, c = _abc(κ, ϵ, λ)
    return 2abs(c) * C.P_1 + C.P_1_dξ * ξ₁^-2
end

function C_exp_P_2_dξ(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    a, b, c = _abc(κ, ϵ, λ)
    return 2abs(c) * C.P_2 + C.P_2_dξ * ξ₁^-2
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

function C_J_E_1_dξ(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    (; d) = λ
    # IMPROVE: This bound could likely be improved by about a factor 2
    # by taking into account cancellations between the two terms.
    return abs(B_W_1(lambda, κ, ϵ, λ)) * (C.exp_E_1_dξ + (d - 1) * C.E_1)
end

function C_J_E_2_dξ(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    (; d) = λ
    # IMPROVE: This bound could likely be improved by about a factor 2
    # taking into account cancellations between the two terms.
    return abs(B_W_2(lambda, κ, ϵ, λ)) * (C.exp_E_2_dξ + (d - 1) * C.E_2)
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

function C_K_2_dξ_1(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    return inv(sqrt(1 + ϵ^2)) * C.J_E_1_dξ
end

function C_K_2_dξ_2(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    return inv(sqrt(1 + ϵ^2)) * C.J_E_2_dξ
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

function C_H_11(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
    C_Q_hat::Arb,
)
    return 2inv(sqrt(1 + ϵ)) * C.J_E_1 * C_Q_hat^2
end

function C_H_12(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
    C_Q_hat::Arb,
)
    return inv(sqrt(1 + ϵ)) * C.J_E_1 * C_Q_hat^2
end

function C_H_21(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
    C_Q_hat::Arb,
)
    return inv(sqrt(1 + ϵ)) * C.J_E_2 * C_Q_hat^2
end

function C_H_22(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
    C_Q_hat::Arb,
)
    return 2inv(sqrt(1 + ϵ)) * C.J_E_2 * C_Q_hat^2
end

function C_H_11_dξ(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
    C_Q_hat::Arb,
    C_Q_hat_dξ::Arb,
)
    return 2inv(sqrt(1 + ϵ)) * (C.J_E_1_dξ * C_Q_hat^2 + 2C.J_E_1 * C_Q_hat_dξ * C_Q_hat)
end

function C_H_12_dξ(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
    C_Q_hat::Arb,
    C_Q_hat_dξ::Arb,
)
    return inv(sqrt(1 + ϵ)) * (C.J_E_1_dξ * C_Q_hat^2 + 2C.J_E_1 * C_Q_hat_dξ * C_Q_hat)
end

function C_H_21_dξ(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
    C_Q_hat::Arb,
    C_Q_hat_dξ::Arb,
)
    return inv(sqrt(1 + ϵ)) * (C.J_E_2_dξ * C_Q_hat^2 + 2C.J_E_2 * C_Q_hat_dξ * C_Q_hat)
end

function C_H_22_dξ(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
    C_Q_hat::Arb,
    C_Q_hat_dξ::Arb,
)
    return 2inv(sqrt(1 + ϵ)) * (C.J_E_2_dξ * C_Q_hat^2 + 2C.J_E_2 * C_Q_hat_dξ * C_Q_hat)
end
