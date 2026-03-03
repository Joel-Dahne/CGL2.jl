"""
    FunctionBounds_hat(κ, ϵ, ξ₁, λ)

This contains all the bounds of functions that are needed in the
enclosure asymptotic expansion of `Q_hat` at infinity.

More precisely it contains the bounds from

- Lemma REF(lemma:P_hat-E_hat-bounds)

It checks all the conditions on the parameters that these lemmas
assume. If any of these conditions are not satisfied it will throw an
error. When using this struct the bounds can therefore safely be
assume to hold.
"""
struct FunctionBounds_hat
    P_hat::Arb
    E_hat::Arb
    J_P_hat::Arb
    J_E_hat::Arb

    FunctionBounds_hat() =
        new(indeterminate(Arb), indeterminate(Arb), indeterminate(Arb), indeterminate(Arb))
end

function FunctionBounds_hat(κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    # TODO: Add checks for conditions

    a, b, c = _abc(κ, ϵ, λ)
    CU_hat = UBounds(a, b, -c, ξ₁)

    C_hat = FunctionBounds_hat()

    C_hat.P_hat[] = C_P_hat(κ, ϵ, ξ₁, λ, CU_hat)

    C_hat.E_hat[] = C_E_hat(κ, ϵ, ξ₁, λ, CU_hat)

    C_hat.J_P_hat[] = C_J_P_hat(κ, ϵ, ξ₁, λ, C_hat)
    C_hat.J_E_hat[] = C_J_E_hat(κ, ϵ, ξ₁, λ, C_hat)

    return C_hat
end

function C_P_hat(κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU_hat::UBounds)
    a, b, c = _abc(κ, ϵ, λ)
    return CU_hat.U_a_b * abs((-c)^-a)
end

function C_E_hat(κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU_hat::UBounds)
    a, b, c = _abc(κ, ϵ, λ)

    return CU_hat.U_bma_b * abs(c^(a - b))
end

C_J_P_hat(κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, C_hat::FunctionBounds_hat) =
    abs(B_W_hat(κ, ϵ, λ)) * C_hat.P_hat

C_J_E_hat(κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, C_hat::FunctionBounds_hat) =
    abs(B_W_hat(κ, ϵ, λ)) * C_hat.E_hat
