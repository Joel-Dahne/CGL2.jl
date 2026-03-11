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
    P_hat_dξ::Arb
    E_hat::Arb
    E_hat_dξ::Arb
    J_P_hat::Arb
    R_P_hat::Arb
    J_E_hat::Arb
    R_J_E_hat::Arb

    FunctionBounds_hat() = new(
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

function FunctionBounds_hat(κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    # Requirement of Lemma REF(lemma:P_hat-E_hat-bounds)
    @assert ξ₁ > 1

    a, b, c = _abc(κ, ϵ, λ)
    CU_hat = UBounds(a, b, -c, ξ₁)

    C_hat = FunctionBounds_hat()

    C_hat.P_hat[] = C_P_hat(κ, ϵ, ξ₁, λ, CU_hat)
    C_hat.P_hat_dξ[] = C_P_hat_dξ(κ, ϵ, ξ₁, λ, CU_hat)

    C_hat.E_hat[] = C_E_hat(κ, ϵ, ξ₁, λ, CU_hat)
    C_hat.E_hat_dξ[] = C_E_hat_dξ(κ, ϵ, ξ₁, λ, CU_hat)

    C_hat.J_E_hat[] = C_J_E_hat(κ, ϵ, ξ₁, λ, C_hat)
    C_hat.J_P_hat[] = C_J_P_hat(κ, ϵ, ξ₁, λ, C_hat)

    C_hat.R_P_hat[] = C_R_P_hat(κ, ϵ, ξ₁, λ)
    C_hat.R_J_E_hat[] = C_R_J_E_hat(κ, ϵ, ξ₁, λ)

    return C_hat
end

function C_P_hat(κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU_hat::UBounds)
    a, b, c = _abc(κ, ϵ, λ)
    return CU_hat.U_a_b * abs((-c)^-a)
end

function C_P_hat_dξ(κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU_hat::UBounds)
    a, b, c = _abc(κ, ϵ, λ)
    return CU_hat.U_dz_a_b * abs(2(-c)^-a)
end

function C_E_hat(κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU_hat::UBounds)
    a, b, c = _abc(κ, ϵ, λ)
    return CU_hat.U_bma_b * abs(c^(a - b))
end

function C_E_hat_dξ(κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU_hat::UBounds)
    a, b, c = _abc(κ, ϵ, λ)
    return 2(CU_hat.U_bma_b * abs(c) + CU_hat.U_dz_bma_b * ξ₁^-2) * abs(c^(a - b))
end

C_J_P_hat(κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, C_hat::FunctionBounds_hat) =
    abs(B_W_hat(κ, ϵ, λ)) * C_hat.P_hat

C_J_E_hat(κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, C_hat::FunctionBounds_hat) =
    abs(B_W_hat(κ, ϵ, λ)) * C_hat.E_hat

function C_R_P_hat(κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    a, b, c = _abc(κ, ϵ, λ)

    z₁ = -c * ξ₁^2
    n = 20
    S = sum(1:(n-1)) do k
        # p_U(k, a, b, z₁) gives us coefficient with (-z₁)^k in
        # denominator, we want (-z₁)^(k - 1) so multiply by -z₁
        abs(p_U(k, a, b, z₁) * (-z₁))
    end

    return (S + C_R_U(n, a, b, z₁) * abs(z₁)^(-n + 1)) * abs((-c)^(-a - 1))
end

function C_R_J_E_hat(κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    a, b, c = _abc(κ, ϵ, λ)

    z₁ = c * ξ₁^2
    n = 20
    S = sum(1:(n-1)) do k
        # p_U(k, b - a, b, z₁) gives us coefficient with (-z₁)^k in
        # denominator, we want (-z₁)^(k - 1) so multiply by -z₁
        abs(p_U(k, b - a, b, z₁) * (-z₁))
    end

    return abs(B_W_hat(κ, ϵ, λ)) *
           abs(c^(a - b - 1)) *
           (S + C_R_U(n, b - a, b, z₁) * abs(z₁)^(-n + 1))
end
