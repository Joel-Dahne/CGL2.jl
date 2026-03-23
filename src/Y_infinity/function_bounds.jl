"""
    FunctionBounds_Y(lambda, γ₁, γ₂, κ, ϵ, ξ₁, λ)

Contains the constants involved in asymptotic bounds for functions
that are needed in the enclosure of `Q_hat` at infinity.

More precisely it contains the bounds from

- Lemma REF(lemma:P_i_E_i-bounds)
- Lemma REF(lemma:I_K_1-I_K_2-bounds)
- Lemma REF(lemma:bound-I_N)
- Lemma REF(lemma:H-bounds)

It checks all the conditions on the parameters that these lemmas
assume. If any of these conditions are not satisfied it will throw an
error. When using this struct the bounds can therefore safely be
assume to hold.
"""
struct FunctionBounds_Y
    # Lemma REF(lemma:P_i_E_i-bounds)
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
    # Lemma REF(lemma:bound-K_1-K_2)
    K_1_1::Arb
    K_1_2::Arb
    K_2_1::Arb
    K_2_2::Arb
    K_2_dξ_1::Arb
    K_2_dξ_2::Arb
    K_1_dλ_1::Arb
    K_1_dλ_2::Arb
    K_2_dλ_1::Arb
    K_2_dλ_2::Arb
    # Lemma REF(lemma:bound-I_N)
    I_N::Arb
    I_N_dξ::Arb
    # Lemma REF(lemma:H-bounds)
    H_11::Arb
    H_12::Arb
    H_21::Arb
    H_22::Arb
    H_11_dξ::Arb
    H_12_dξ::Arb
    H_21_dξ::Arb
    H_22_dξ::Arb

    FunctionBounds_Y() = new(
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
    )
end

function FunctionBounds_Y(
    lambda::Acb,
    γ₁::Acb,
    γ₂::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
)
    (; σ, δ) = λ
    a, b, c = _abc(κ, ϵ, λ)

    # This is the only direct condition in Lemma
    # REF(lemma:P_i_E_i-bounds), REF(lemma:I_K_1-I_K_2-bounds) and
    # REF(lemma:H-bounds).
    # The conditions related to the bounds for U are checked by
    # Ubounds.
    ξ₁ > 1 || throw(ArgumentError("ξ₁ > 1 not satisfied"))

    # This arethe only direct condition for Lemma
    # REF(lemma:bound-I_N). Note that the conditions for C_Q_hat and
    # C_Q_hat_dξ and checked by their respective methods.
    isone(σ) || throw(ArgumentError("σ = 1 not satisfied"))
    iszero(δ) || throw(ArgumentError("δ = 0 not satisfied"))

    C = FunctionBounds_Y()

    CU = UBounds(a - lambda / 2κ, b, -c, ξ₁, include_da = true)
    CU_conj = UBounds(conj(a) - lambda / 2κ, b, -conj(c), ξ₁, include_da = true)

    C.E_1[] = C_E_1(lambda, κ, ϵ, ξ₁, λ, CU)
    C.E_2[] = C_E_2(lambda, κ, ϵ, ξ₁, λ, CU_conj)
    C.P_1[] = C_P_1(lambda, κ, ϵ, ξ₁, λ, CU)
    C.P_2[] = C_P_2(lambda, κ, ϵ, ξ₁, λ, CU_conj)

    C.E_1_dξ[] = C_E_1_dξ(lambda, κ, ϵ, ξ₁, λ, CU)
    C.E_2_dξ[] = C_E_2_dξ(lambda, κ, ϵ, ξ₁, λ, CU_conj)
    C.P_1_dξ[] = C_P_1_dξ(lambda, κ, ϵ, ξ₁, λ, CU)
    C.P_2_dξ[] = C_P_2_dξ(lambda, κ, ϵ, ξ₁, λ, CU_conj)

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

    C.K_1_1[] = C_K_1_1(lambda, κ, ϵ, ξ₁, λ, C)
    C.K_1_2[] = C_K_1_2(lambda, κ, ϵ, ξ₁, λ, C)
    C.K_2_1[] = C_K_2_1(lambda, κ, ϵ, ξ₁, λ, C)
    C.K_2_2[] = C_K_2_2(lambda, κ, ϵ, ξ₁, λ, C)

    C.K_1_dλ_1[] = C_K_1_dλ_1(lambda, κ, ϵ, ξ₁, λ, C)
    C.K_1_dλ_2[] = C_K_1_dλ_2(lambda, κ, ϵ, ξ₁, λ, C)
    C.K_2_dλ_1[] = C_K_2_dλ_1(lambda, κ, ϵ, ξ₁, λ, C)
    C.K_2_dλ_2[] = C_K_2_dλ_2(lambda, κ, ϵ, ξ₁, λ, C)

    C.K_2_dξ_1[] = C_K_2_dξ_1(lambda, κ, ϵ, ξ₁, λ, C)
    C.K_2_dξ_2[] = C_K_2_dξ_2(lambda, κ, ϵ, ξ₁, λ, C)

    C_Q_hat = CGL2.C_Q_hat(γ₁, γ₂, κ, ϵ, ξ₁, λ)
    C_Q_hat_dξ = CGL2.C_Q_hat_dξ(γ₁, γ₂, κ, ϵ, ξ₁, λ)
    C.I_N[] = C_I_N(C_Q_hat)
    C.I_N_dξ[] = C_I_N_dξ(C_Q_hat, C_Q_hat_dξ)

    C.H_11[] = C_H_11(lambda, κ, ϵ, ξ₁, λ, C, C_Q_hat)
    C.H_12[] = C_H_12(lambda, κ, ϵ, ξ₁, λ, C, C_Q_hat)
    C.H_21[] = C_H_21(lambda, κ, ϵ, ξ₁, λ, C, C_Q_hat)
    C.H_22[] = C_H_22(lambda, κ, ϵ, ξ₁, λ, C, C_Q_hat)
    C.H_11_dξ[] = C_H_11_dξ(lambda, κ, ϵ, ξ₁, λ, C, C_Q_hat, C_Q_hat_dξ)
    C.H_12_dξ[] = C_H_12_dξ(lambda, κ, ϵ, ξ₁, λ, C, C_Q_hat, C_Q_hat_dξ)
    C.H_21_dξ[] = C_H_21_dξ(lambda, κ, ϵ, ξ₁, λ, C, C_Q_hat, C_Q_hat_dξ)
    C.H_22_dξ[] = C_H_22_dξ(lambda, κ, ϵ, ξ₁, λ, C, C_Q_hat, C_Q_hat_dξ)

    return C
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
    a, b, c = _abc(κ, ϵ, λ)
    a_tilde = b - a + lambda / 2κ
    z₁ = c * ξ₁^2
    n = 5

    S = sum(0:(n-1)) do k
        abs((d - 1) * p_U(k, a_tilde, b, z₁) - 2a_tilde * p_U(k, a_tilde + 1, b + 1, z₁))
    end

    R =
        abs(d - 1) * C_R_U(n, a_tilde, b, z₁) +
        2abs(a_tilde) * C_R_U(n, a_tilde + 1, b + 1, z₁)

    return abs(B_W_1(lambda, κ, ϵ, λ)) * abs(c^-a_tilde) * (S + R * abs(z₁)^-n)
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
    a, b, c = _abc(κ, ϵ, λ)
    c_conj = conj(c)
    a_tilde_conj = b - conj(a) + lambda / 2κ
    z₁_conj = c_conj * ξ₁^2
    n = 5

    S = sum(0:(n-1)) do k
        abs(
            (d - 1) * p_U(k, a_tilde_conj, b, z₁_conj) -
            2a_tilde_conj * p_U(k, a_tilde_conj + 1, b + 1, z₁_conj),
        )
    end

    R =
        abs(d - 1) * C_R_U(n, a_tilde_conj, b, z₁_conj) +
        2abs(a_tilde_conj) * C_R_U(n, a_tilde_conj + 1, b + 1, z₁_conj)

    return abs(B_W_2(lambda, κ, ϵ, λ)) *
           abs(c_conj^-a_tilde_conj) *
           (S + R * abs(z₁_conj)^-n)
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

C_I_N(C_Q_hat::Arb) = 3C_Q_hat^2

C_I_N_dξ(C_Q_hat::Arb, C_Q_hat_dξ::Arb) = 6C_Q_hat * C_Q_hat_dξ

function C_H_11(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
    C_Q_hat::Arb,
)
    return 2inv(sqrt(1 + ϵ^2)) * C.J_E_1 * C_Q_hat^2
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
    return inv(sqrt(1 + ϵ^2)) * C.J_E_1 * C_Q_hat^2
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
    return inv(sqrt(1 + ϵ^2)) * C.J_E_2 * C_Q_hat^2
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
    return 2inv(sqrt(1 + ϵ^2)) * C.J_E_2 * C_Q_hat^2
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
    return 2inv(sqrt(1 + ϵ^2)) * (C.J_E_1_dξ * C_Q_hat^2 + 2C.J_E_1 * C_Q_hat_dξ * C_Q_hat)
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
    return inv(sqrt(1 + ϵ^2)) * (C.J_E_1_dξ * C_Q_hat^2 + 2C.J_E_1 * C_Q_hat_dξ * C_Q_hat)
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
    return inv(sqrt(1 + ϵ^2)) * (C.J_E_2_dξ * C_Q_hat^2 + 2C.J_E_2 * C_Q_hat_dξ * C_Q_hat)
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
    return 2inv(sqrt(1 + ϵ^2)) * (C.J_E_2_dξ * C_Q_hat^2 + 2C.J_E_2 * C_Q_hat_dξ * C_Q_hat)
end
