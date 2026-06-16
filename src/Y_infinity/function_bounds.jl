"""
    FunctionBounds_Y(λ, γ₁, γ₂, κ, ϵ, ξ₁, Λ)

Contains the constants involved in asymptotic bounds for functions
that are needed in the enclosure of `Y` at infinity.

More precisely it contains the bounds from

- Lemma REF(lemma:P_i_E_i-bounds)
- Lemma REF(lemma:bound-K_1-K_2)
- Lemma REF(lemma:bound-I_N)
- Lemma REF(lemma:I_K_1-I_K_2-bounds)
- Lemma REF(lemma:Z-fixed-point-bounds)

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
    J_E_1::Arb
    J_E_2::Arb
    J_P_1::Arb
    J_P_2::Arb
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
    )
end

function FunctionBounds_Y(
    λ::Acb,
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
    # REF(lemma:P_i_E_i-bounds) and REF(lemma:I_K_1-I_K_2-bounds).
    # The conditions related to the bounds for U are checked by
    # Ubounds.
    @assert ξ₁ > 1

    # These are the only direct conditions for Lemma
    # REF(lemma:bound-I_N). Note that the conditions for C_Q_hat is
    # checked by its method.
    @assert isone(σ)
    @assert iszero(δ)

    # These are requirements of Lemmas REF(lemma:I_K_1-I_K_2-bounds)
    exponent = 2 / σ - d - 2real(λ) / κ - 4
    @assert real(c) > 0
    @assert exponent < 0

    # The requirements for Lemma REF(lemma:Z-fixed-point-bounds) are
    # the same as for REF(lemma:I_K_1-I_K_2-bounds).

    C = FunctionBounds_Y()

    CU = UBounds(a - λ / 2κ, b, -c, ξ₁, include_da = true)
    CU_conj = UBounds(conj(a) - λ / 2κ, b, -conj(c), ξ₁, include_da = true)

    # Lemma REF(lemma:P_i_E_i-bounds)

    C.E_1[] = C_E_1(λ, κ, ϵ, ξ₁, Λ, CU)
    C.E_2[] = C_E_2(λ, κ, ϵ, ξ₁, Λ, CU_conj)
    C.P_1[] = C_P_1(λ, κ, ϵ, ξ₁, Λ, CU)
    C.P_2[] = C_P_2(λ, κ, ϵ, ξ₁, Λ, CU_conj)

    C.J_E_1[] = C_J_E_1(λ, κ, ϵ, ξ₁, Λ, C)
    C.J_E_2[] = C_J_E_2(λ, κ, ϵ, ξ₁, Λ, C)
    C.J_P_1[] = C_J_P_1(λ, κ, ϵ, ξ₁, Λ, C)
    C.J_P_2[] = C_J_P_2(λ, κ, ϵ, ξ₁, Λ, C)

    # Lemma REF(lemma:bound-K_1-K_2)

    C.K_1_1[] = C_K_1_1(λ, κ, ϵ, ξ₁, Λ, C)
    C.K_1_2[] = C_K_1_2(λ, κ, ϵ, ξ₁, Λ, C)
    C.K_2_1[] = C_K_2_1(λ, κ, ϵ, ξ₁, Λ, C)
    C.K_2_2[] = C_K_2_2(λ, κ, ϵ, ξ₁, Λ, C)

    # Lemma REF(lemma:bound-I_N)

    C.I_N[] = C_I_N(γ₁, γ₂, κ, ϵ, ξ₁, Λ)

    # Lemma REF(lemma:I_K_1-I_K_2-bounds)

    C.I_K_1_1[] = C_I_K_1_1(C)
    C.I_K_1_2[] = C_I_K_1_2(C)
    C.I_K_2_1[] = C_I_K_2_1(κ, ϵ, Λ, C)
    C.I_K_2_2[] = C_I_K_2_2(κ, ϵ, Λ, C)

    # Lemma REF(lemma:Z-fixed-point-bounds)

    C.T_12[] = C_T_12(λ, κ, ϵ, ξ₁, Λ, C)

    return C
end

# Lemma REF(lemma:P_i_E_i-bounds)

function C_E_1(λ::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, Λ::CGLParams{Arb}, CU::UBounds)
    a, b, c = _abc(κ, ϵ, Λ)
    return CU.U_bma_b * abs(c^(-b + a - λ / 2κ))
end

function C_E_2(λ::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, Λ::CGLParams{Arb}, CU_conj::UBounds)
    a, b, c = _abc(κ, ϵ, Λ)
    return CU_conj.U_bma_b * abs(conj(c)^(-b + conj(a) - λ / 2κ))
end

function C_P_1(λ::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, Λ::CGLParams{Arb}, CU::UBounds)
    a, b, c = _abc(κ, ϵ, Λ)
    return CU.U_a_b * abs((-c)^(-a + λ / 2κ))
end

function C_P_2(λ::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, Λ::CGLParams{Arb}, CU_conj::UBounds)
    a, b, c = _abc(κ, ϵ, Λ)
    return CU_conj.U_a_b * abs(conj(-c)^(-conj(a) + λ / 2κ))
end

C_J_E_1(λ::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, Λ::CGLParams{Arb}, C::FunctionBounds_Y) =
    abs(B_W_1(λ, κ, ϵ, Λ)) * C.E_1

C_J_E_2(λ::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, Λ::CGLParams{Arb}, C::FunctionBounds_Y) =
    abs(B_W_2(λ, κ, ϵ, Λ)) * C.E_2


C_J_P_1(λ::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, Λ::CGLParams{Arb}, C::FunctionBounds_Y) =
    abs(B_W_1(λ, κ, ϵ, Λ)) * C.P_1

C_J_P_2(λ::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, Λ::CGLParams{Arb}, C::FunctionBounds_Y) =
    abs(B_W_2(λ, κ, ϵ, Λ)) * C.P_2

# Lemma REF(lemma:bound-K_1-K_2)

C_K_1_1(λ::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, Λ::CGLParams{Arb}, C::FunctionBounds_Y) =
    inv(sqrt(1 + ϵ^2)) * C.J_P_1

C_K_1_2(λ::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, Λ::CGLParams{Arb}, C::FunctionBounds_Y) =
    inv(sqrt(1 + ϵ^2)) * C.J_P_2

C_K_2_1(λ::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, Λ::CGLParams{Arb}, C::FunctionBounds_Y) =
    inv(sqrt(1 + ϵ^2)) * C.J_E_1

C_K_2_2(λ::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, Λ::CGLParams{Arb}, C::FunctionBounds_Y) =
    inv(sqrt(1 + ϵ^2)) * C.J_E_2

# Lemma REF(lemma:bound-I_N)

function C_I_N(γ₁::Acb, γ₂::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, Λ::CGLParams{Arb})
    # Compute a bound for the norm of Q_hat
    C = FunctionBounds_hat(κ, ϵ, ξ₁, Λ)
    norms = NormBounds_hat(γ₁, γ₂, κ, ϵ, ξ₁, Λ, C)

    return 3norms.Q_hat^2
end

# Lemma REF(lemma:I_K_1-I_K_2-bounds)

C_I_K_1_1(C_Y::FunctionBounds_Y) = C_Y.K_1_1 * C_Y.I_N / 2
C_I_K_1_2(C_Y::FunctionBounds_Y) = C_Y.K_1_2 * C_Y.I_N / 2

function C_I_K_2_1(κ::Arb, ϵ::Arb, Λ::CGLParams{Arb}, C_Y::FunctionBounds_Y)
    c = _c(κ, ϵ, Λ)
    return C_Y.K_2_1 * C_Y.I_N / (2real(c))
end

function C_I_K_2_2(κ::Arb, ϵ::Arb, Λ::CGLParams{Arb}, C_Y::FunctionBounds_Y)
    c = _c(κ, ϵ, Λ)
    return C_Y.K_2_2 * C_Y.I_N / (2real(c))
end

# Lemma REF(lemma:Z-fixed-point-bounds)

function C_T_12(λ::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, Λ::CGLParams{Arb}, C_Z::FunctionBounds_Y)
    return max(C_Z.E_1 * C_Z.I_K_1_1, C_Z.E_2 * C_Z.I_K_1_2) +
           max(C_Z.P_1 * C_Z.I_K_2_1, C_Z.P_2 * C_Z.I_K_2_2) * ξ₁^-2
end
