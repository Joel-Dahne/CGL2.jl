"""
    FunctionBounds_hat(κ, ϵ, ξ₁, Λ)

Contains the constants involved in asymptotic bounds for functions
that are needed in the enclosure of `Q_hat` at infinity.

More precisely it contains the bounds from

- Lemma REF(lemma:P_hat-E_hat-bounds)
- Lemma REF(lemma:I_E_hat-I_P_hat-bounds)
- Lemma REF(lemma:fixed-point-bounds)

It checks all the conditions on the parameters that these lemmas
assume. If any of these conditions are not satisfied it will throw an
error. When using this struct the bounds can therefore safely be
assume to hold.
"""
struct FunctionBounds_hat
    # Lemma REF(lemma:P_hat-E_hat-bounds)
    P_hat::Arb
    E_hat::Arb
    J_P_hat::Arb
    J_E_hat::Arb
    # Lemma REF(lemma:I_E_hat-I_P_hat-bounds)
    I_E_hat::Arb
    I_P_hat::Arb
    # Lemma REF(lemma:fixed-point-bounds)
    T_hat::Arb

    FunctionBounds_hat() = new(
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
        indeterminate(Arb),
    )
end

function FunctionBounds_hat(κ::Arb, ϵ::Arb, ξ₁::Arb, Λ::CGLParams{Arb})
    (; d, σ) = Λ
    a, b, c = _abc(κ, ϵ, Λ)

    # Requirement of Lemma REF(lemma:P_hat-E_hat-bounds)
    @assert ξ₁ > 1

    # Requirements of Lemma REF(lemma:I_E_hat-I_P_hat-bounds)
    @assert ξ₁ > 1
    @assert real(c) > 0
    @assert -2 / σ + d - 4 < 0
    @assert 2real(c) + (-2 / σ + d - 4) * ξ₁^-2 > 0

    # Requirements of Lemma REF(lemma:Q-hat-fixed-point-bounds) are
    # the same as for Lemma REF(lemma:I_E_hat-I_P_hat-bounds) plus the
    # following one
    @assert 2 / σ - d < 0

    CU_hat = UBounds(a, b, -c, ξ₁)

    C_hat = FunctionBounds_hat()

    C_hat.P_hat[] = C_P_hat(κ, ϵ, ξ₁, Λ, CU_hat)

    C_hat.E_hat[] = C_E_hat(κ, ϵ, ξ₁, Λ, CU_hat)

    C_hat.J_E_hat[] = C_J_E_hat(κ, ϵ, ξ₁, Λ, C_hat)
    C_hat.J_P_hat[] = C_J_P_hat(κ, ϵ, ξ₁, Λ, C_hat)

    C_hat.I_E_hat[] = C_I_E_hat(κ, ϵ, ξ₁, Λ, C_hat)
    C_hat.I_P_hat[] = C_I_P_hat(κ, ϵ, ξ₁, Λ, C_hat)

    C_hat.T_hat[] = C_hat.P_hat * C_hat.I_E_hat + C_hat.E_hat * C_hat.I_P_hat * ξ₁^-2

    return C_hat
end

function C_P_hat(κ::Arb, ϵ::Arb, ξ₁::Arb, Λ::CGLParams{Arb}, CU_hat::UBounds)
    a, b, c = _abc(κ, ϵ, Λ)
    return CU_hat.U_a_b * abs((-c)^-a)
end

function C_E_hat(κ::Arb, ϵ::Arb, ξ₁::Arb, Λ::CGLParams{Arb}, CU_hat::UBounds)
    a, b, c = _abc(κ, ϵ, Λ)
    return CU_hat.U_bma_b * abs(c^(a - b))
end

C_J_P_hat(κ::Arb, ϵ::Arb, ξ₁::Arb, Λ::CGLParams{Arb}, C_hat::FunctionBounds_hat) =
    abs(B_W_hat(κ, ϵ, Λ)) * C_hat.P_hat

C_J_E_hat(κ::Arb, ϵ::Arb, ξ₁::Arb, Λ::CGLParams{Arb}, C_hat::FunctionBounds_hat) =
    abs(B_W_hat(κ, ϵ, Λ)) * C_hat.E_hat

C_I_E_hat(κ::Arb, ϵ::Arb, ξ₁::Arb, Λ::CGLParams{Arb}, C_hat::FunctionBounds_hat) =
    C_hat.J_E_hat / 2

function C_I_P_hat(κ::Arb, ϵ::Arb, ξ₁::Arb, Λ::CGLParams{Arb}, C_hat::FunctionBounds_hat)
    (; d, σ) = Λ
    c = _c(κ, ϵ, Λ)
    return C_hat.J_P_hat / (2real(c) + (-2 / σ + d - 4) * ξ₁^-2)
end
