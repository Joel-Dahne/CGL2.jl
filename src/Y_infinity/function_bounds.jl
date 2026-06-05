"""
    FunctionBounds_Y(lambda, γ₁, γ₂, κ, ϵ, ξ₁, Λ)

Contains the constants involved in asymptotic bounds for functions
that are needed in the enclosure of `Q_hat` at infinity.

More precisely it contains the bounds from

- Lemma REF(lemma:P_i_E_i-bounds)
- Lemma REF(lemma:bound-K_1-K_2)
- Lemma REF(lemma:bound-I_N)
- Lemma REF(lemma:I_K_1-I_K_2-bounds)
- Lemma REF(lemma:Z-fixed-point-bounds)
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
    J_E_1::Arb
    J_E_2::Arb
    J_P_1::Arb
    J_P_2::Arb
    J_E_1_dξ::Arb
    J_E_2_dξ::Arb
    # Lemma REF(lemma:bound-K_1-K_2)
    K_1_1::Arb
    K_1_2::Arb
    K_2_1::Arb
    K_2_2::Arb
    # Lemma REF(lemma:bound-I_N)
    I_N::Arb
    # Lemma REF(lemma:I_K_1-I_K_2-bounds)
    I_K_1_1::Arb
    I_K_1_2::Arb
    I_K_2_1::Arb
    I_K_2_2::Arb
    # Lemma REF(lemma:Z-fixed-point-bounds)
    T_12::Arb
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
    )
end

function FunctionBounds_Y(
    lambda::Acb,
    γ₁::Acb,
    γ₂::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    Λ::CGLParams{Arb};
)
    (; d, σ, δ) = Λ
    a, b, c = _abc(κ, ϵ, Λ)

    # This is the only direct condition in Lemma
    # REF(lemma:P_i_E_i-bounds), REF(lemma:I_K_1-I_K_2-bounds) and
    # REF(lemma:H-bounds).
    # The conditions related to the bounds for U are checked by
    # Ubounds.
    ξ₁ > 1 || throw(ArgumentError("ξ₁ > 1 not satisfied"))

    # These are the only direct condition for Lemma
    # REF(lemma:bound-I_N). Note that the conditions for C_Q_hat is
    # checked by its method.
    isone(σ) || throw(ArgumentError("σ = 1 not satisfied"))
    iszero(δ) || throw(ArgumentError("δ = 0 not satisfied"))

    # These are requirements of Lemmas REF(lemma:I_K_1-I_K_2-bounds)
    exponent = 2 / σ - d - 2real(lambda) / κ - 4
    real(c) > 0 || throw(ArgumentError("real(c) > 0 not satisfied"))
    exponent < 0 || throw(ArgumentError("exponent < 0 not satisfied"))

    # The requirements for Lemma REF(lemma:Z-fixed-point-bounds) are
    # the same as for REF(lemma:I_K_1-I_K_2-bounds).

    C = FunctionBounds_Y()

    CU = UBounds(a - lambda / 2κ, b, -c, ξ₁, include_da = true)
    CU_conj = UBounds(conj(a) - lambda / 2κ, b, -conj(c), ξ₁, include_da = true)

    # Lemma REF(lemma:P_i_E_i-bounds)

    C.E_1[] = C_E_1(lambda, κ, ϵ, ξ₁, Λ, CU)
    C.E_2[] = C_E_2(lambda, κ, ϵ, ξ₁, Λ, CU_conj)
    C.P_1[] = C_P_1(lambda, κ, ϵ, ξ₁, Λ, CU)
    C.P_2[] = C_P_2(lambda, κ, ϵ, ξ₁, Λ, CU_conj)

    C.E_1_dξ[] = C_E_1_dξ(lambda, κ, ϵ, ξ₁, Λ, CU)
    C.E_2_dξ[] = C_E_2_dξ(lambda, κ, ϵ, ξ₁, Λ, CU_conj)
    C.P_1_dξ[] = C_P_1_dξ(lambda, κ, ϵ, ξ₁, Λ, CU)
    C.P_2_dξ[] = C_P_2_dξ(lambda, κ, ϵ, ξ₁, Λ, CU_conj)

    C.exp_E_1_dξ[] = C_exp_E_1_dξ(lambda, κ, ϵ, ξ₁, Λ, CU)
    C.exp_E_2_dξ[] = C_exp_E_2_dξ(lambda, κ, ϵ, ξ₁, Λ, CU_conj)
    C.exp_P_1_dξ[] = C_exp_P_1_dξ(lambda, κ, ϵ, ξ₁, Λ, C)
    C.exp_P_2_dξ[] = C_exp_P_2_dξ(lambda, κ, ϵ, ξ₁, Λ, C)

    C.J_E_1[] = C_J_E_1(lambda, κ, ϵ, ξ₁, Λ, C)
    C.J_E_2[] = C_J_E_2(lambda, κ, ϵ, ξ₁, Λ, C)
    C.J_P_1[] = C_J_P_1(lambda, κ, ϵ, ξ₁, Λ, C)
    C.J_P_2[] = C_J_P_2(lambda, κ, ϵ, ξ₁, Λ, C)

    C.J_E_1_dξ[] = C_J_E_1_dξ(lambda, κ, ϵ, ξ₁, Λ, C)
    C.J_E_2_dξ[] = C_J_E_2_dξ(lambda, κ, ϵ, ξ₁, Λ, C)

    # Lemma REF(lemma:bound-K_1-K_2)

    C.K_1_1[] = C_K_1_1(lambda, κ, ϵ, ξ₁, Λ, C)
    C.K_1_2[] = C_K_1_2(lambda, κ, ϵ, ξ₁, Λ, C)
    C.K_2_1[] = C_K_2_1(lambda, κ, ϵ, ξ₁, Λ, C)
    C.K_2_2[] = C_K_2_2(lambda, κ, ϵ, ξ₁, Λ, C)

    # Lemma REF(lemma:bound-I_N)

    C_Q_hat = CGL2.C_Q_hat(γ₁, γ₂, κ, ϵ, ξ₁, Λ)
    C_Q_hat_dξ = CGL2.C_Q_hat_dξ(γ₁, γ₂, κ, ϵ, ξ₁, Λ)
    C.I_N[] = C_I_N(C_Q_hat)

    # Lemma REF(lemma:I_K_1-I_K_2-bounds)

    C.I_K_1_1[] = C_I_K_1_1(C)
    C.I_K_1_2[] = C_I_K_1_2(C)
    C.I_K_2_1[] = C_I_K_2_1(κ, ϵ, Λ, C)
    C.I_K_2_2[] = C_I_K_2_2(κ, ϵ, Λ, C)

    # Lemma REF(lemma:Z-fixed-point-bounds)

    C.T_12[] = C_T_12(lambda, κ, ϵ, ξ₁, Λ, C)

    # Lemma REF(lemma:H-bounds)

    C.H_11[] = C_H_11(lambda, κ, ϵ, ξ₁, Λ, C, C_Q_hat)
    C.H_12[] = C_H_12(lambda, κ, ϵ, ξ₁, Λ, C, C_Q_hat)
    C.H_21[] = C_H_21(lambda, κ, ϵ, ξ₁, Λ, C, C_Q_hat)
    C.H_22[] = C_H_22(lambda, κ, ϵ, ξ₁, Λ, C, C_Q_hat)
    C.H_11_dξ[] = C_H_11_dξ(lambda, κ, ϵ, ξ₁, Λ, C, C_Q_hat, C_Q_hat_dξ)
    C.H_12_dξ[] = C_H_12_dξ(lambda, κ, ϵ, ξ₁, Λ, C, C_Q_hat, C_Q_hat_dξ)
    C.H_21_dξ[] = C_H_21_dξ(lambda, κ, ϵ, ξ₁, Λ, C, C_Q_hat, C_Q_hat_dξ)
    C.H_22_dξ[] = C_H_22_dξ(lambda, κ, ϵ, ξ₁, Λ, C, C_Q_hat, C_Q_hat_dξ)

    return C
end

# Lemma REF(lemma:P_i_E_i-bounds)

function C_E_1(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, Λ::CGLParams{Arb}, CU::UBounds)
    a, b, c = _abc(κ, ϵ, Λ)
    return CU.U_bma_b * abs(c^(-b + a - lambda / 2κ))
end

function C_E_2(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, Λ::CGLParams{Arb}, CU_conj::UBounds)
    a, b, c = _abc(κ, ϵ, Λ)
    return CU_conj.U_bma_b * abs(conj(c)^(-b + conj(a) - lambda / 2κ))
end

function C_P_1(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, Λ::CGLParams{Arb}, CU::UBounds)
    a, b, c = _abc(κ, ϵ, Λ)
    return CU.U_a_b * abs((-c)^(-a + lambda / 2κ))
end

function C_P_2(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, Λ::CGLParams{Arb}, CU_conj::UBounds)
    a, b, c = _abc(κ, ϵ, Λ)
    return CU_conj.U_a_b * abs(conj(-c)^(-conj(a) + lambda / 2κ))
end

function C_E_1_dξ(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, Λ::CGLParams{Arb}, CU::UBounds)
    a, b, c = _abc(κ, ϵ, Λ)
    return 2(CU.U_bma_b * abs(c) + CU.U_dz_bma_b * ξ₁^-2) * abs(c^(-b + a - lambda / 2κ))
end

function C_E_2_dξ(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, Λ::CGLParams{Arb}, CU_conj::UBounds)
    a, b, c = _abc(κ, ϵ, Λ)
    return 2(CU_conj.U_bma_b * abs(conj(c)) + CU_conj.U_dz_bma_b * ξ₁^-2) *
           abs(conj(c)^(-b + conj(a) - lambda / 2κ))
end

function C_P_1_dξ(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, Λ::CGLParams{Arb}, CU::UBounds)
    a, b, c = _abc(κ, ϵ, Λ)
    return 2CU.U_dz_a_b * abs((-c)^(-a + lambda / 2κ))
end

function C_P_2_dξ(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, Λ::CGLParams{Arb}, CU_conj::UBounds)
    a, b, c = _abc(κ, ϵ, Λ)
    return 2CU_conj.U_dz_a_b * abs((-conj(c))^(-conj(a) + lambda / 2κ))
end

function C_exp_E_1_dξ(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, Λ::CGLParams{Arb}, CU::UBounds)
    a, b, c = _abc(κ, ϵ, Λ)
    return 2CU.U_dz_bma_b * abs(c^(-b + a - lambda / 2κ))
end

function C_exp_E_2_dξ(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    Λ::CGLParams{Arb},
    CU_conj::UBounds,
)
    a, b, c = _abc(κ, ϵ, Λ)
    return 2CU_conj.U_dz_bma_b * abs(conj(c)^(-b + conj(a) - lambda / 2κ))
end

function C_exp_P_1_dξ(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    a, b, c = _abc(κ, ϵ, Λ)
    return 2abs(c) * C.P_1 + C.P_1_dξ * ξ₁^-2
end

function C_exp_P_2_dξ(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    a, b, c = _abc(κ, ϵ, Λ)
    return 2abs(c) * C.P_2 + C.P_2_dξ * ξ₁^-2
end

function C_J_E_1(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    return abs(B_W_1(lambda, κ, ϵ, Λ)) * C.E_1
end

function C_J_E_2(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    return abs(B_W_2(lambda, κ, ϵ, Λ)) * C.E_2
end

function C_J_P_1(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    return abs(B_W_1(lambda, κ, ϵ, Λ)) * C.P_1
end

function C_J_P_2(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    return abs(B_W_2(lambda, κ, ϵ, Λ)) * C.P_2
end

function C_J_E_1_dξ(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    (; d) = Λ
    a, b, c = _abc(κ, ϵ, Λ)
    a_tilde = b - a + lambda / 2κ
    z₁ = c * ξ₁^2
    n = 5

    S = sum(0:(n-1)) do k
        abs((d - 1) * p_U(k, a_tilde, b, z₁) - 2a_tilde * p_U(k, a_tilde + 1, b + 1, z₁))
    end

    R =
        abs(d - 1) * C_R_U(n, a_tilde, b, z₁) +
        2abs(a_tilde) * C_R_U(n, a_tilde + 1, b + 1, z₁)

    return abs(B_W_1(lambda, κ, ϵ, Λ)) * abs(c^-a_tilde) * (S + R * abs(z₁)^-n)
end

function C_J_E_2_dξ(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    (; d) = Λ
    a, b, c = _abc(κ, ϵ, Λ)
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

    return abs(B_W_2(lambda, κ, ϵ, Λ)) *
           abs(c_conj^-a_tilde_conj) *
           (S + R * abs(z₁_conj)^-n)
end

# Lemma REF(lemma:bound-K_1-K_2)

function C_K_1_1(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    return inv(sqrt(1 + ϵ^2)) * C.J_P_1
end

function C_K_1_2(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    return inv(sqrt(1 + ϵ^2)) * C.J_P_2
end

function C_K_2_1(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    return inv(sqrt(1 + ϵ^2)) * C.J_E_1
end

function C_K_2_2(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    return inv(sqrt(1 + ϵ^2)) * C.J_E_2
end

# Lemma REF(lemma:bound-I_N)

C_I_N(C_Q_hat::Arb) = 3C_Q_hat^2

# Lemma REF(lemma:I_K_1-I_K_2-bounds)

function C_I_K_1_1(C_Y::FunctionBounds_Y)
    return C_Y.K_1_1 * C_Y.I_N / 2
end

function C_I_K_1_2(C_Y::FunctionBounds_Y)
    return C_Y.K_1_2 * C_Y.I_N / 2
end

function C_I_K_2_1(κ::Arb, ϵ::Arb, Λ::CGLParams{Arb}, C_Y::FunctionBounds_Y)
    c = _c(κ, ϵ, Λ)
    return C_Y.K_2_1 * C_Y.I_N / (2real(c))
end

function C_I_K_2_2(κ::Arb, ϵ::Arb, Λ::CGLParams{Arb}, C_Y::FunctionBounds_Y)
    c = _c(κ, ϵ, Λ)
    return C_Y.K_2_2 * C_Y.I_N / (2real(c))
end

# Lemma REF(lemma:Z-fixed-point-bounds)

function C_T_12(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    Λ::CGLParams{Arb},
    C_Z::FunctionBounds_Y,
)
    return max(C_Z.E_1 * C_Z.I_K_1_1, C_Z.E_2 * C_Z.I_K_1_2) +
           max(C_Z.P_1 * C_Z.I_K_2_1, C_Z.P_2 * C_Z.I_K_2_2) * ξ₁^-2
end

# Lemma REF(lemma:H-bounds)

function C_H_11(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    Λ::CGLParams{Arb},
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
    Λ::CGLParams{Arb},
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
    Λ::CGLParams{Arb},
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
    Λ::CGLParams{Arb},
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
    Λ::CGLParams{Arb},
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
    Λ::CGLParams{Arb},
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
    Λ::CGLParams{Arb},
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
    Λ::CGLParams{Arb},
    C::FunctionBounds_Y,
    C_Q_hat::Arb,
    C_Q_hat_dξ::Arb,
)
    return 2inv(sqrt(1 + ϵ^2)) * (C.J_E_2_dξ * C_Q_hat^2 + 2C.J_E_2 * C_Q_hat_dξ * C_Q_hat)
end
