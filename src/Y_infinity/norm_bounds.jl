"""
    NormBounds_hat(...)

Contains bounds for the norms of `Z` and `Z_dλ`.

The bound for `Z` is based on Lemma REF(lemma:Z-fixed-point-bounds)
and Proposition REF(prop:Y-fixed-point). The bound for `Z_dλ` is based
on Lemma REF(lemma:Z-lambda-fixed-point-bounds).

It checks all the conditions on the parameters that these lemmas
assume. If any of these conditions are not satisfied it will throw an
error. When using this struct the bounds can therefore safely be
assume to hold.
"""
struct NormBounds_Y
    Z::Arb
    Z_dλ::Arb

    NormBounds_Y() = new(indeterminate(Arb), indeterminate(Arb))
end

function NormBounds_Y(
    c::SVector{2,Acb},
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    norms = NormBounds_Y()

    # The requirements for Lemma REF(lemma:Z-fixed-point-bounds) are
    # checked in the computation of `C`, where the associated
    # constants are computed.

    # The requirements for the fixed point in Proposition
    # REF(prop:Z-fixed-point) are checked by the norm_bound_Z
    # function.

    # The requirements for Lemma
    # REF(lemma:Z-lambda-fixed-point-bounds) that are related to
    # Lemmas REF(lemma:I_E_hat-I_P_hat-dgamma-bounds) and
    # REF(lemma:I_K_1_I_K_2-lambda-2-bounds) are checked in the
    # computation of `C`, where the associated constants are computed.
    # The other conditions are checked in the norm_bound_Z_dλ
    # function.

    norms.Z[] = norm_bound_Z(c, lambda, κ, ϵ, ξ₁, v, λ, C)
    norms.Z_dλ[] = norm_bound_Z_dλ(c, lambda, κ, ϵ, ξ₁, v, λ, C, norms)

    return norms
end

"""
    norm_bound_Z(c, lambda, κ, ϵ, ξ₁, v, λ, C)

To apply the fixed point theorem in Proposition
REF(prop:Z-fixed-point) we need to find `ρ` satisfying the inequality

```
max(C_E_1 * abs(c[1]), C_E_2 * abs(c[2])) + C_T_12 * ξ₁^-2 * ρ <= ρ
```

and

```
C_T_12 * ξ₁^-2 < 1
```

The second inequality means that the first inequality is satisfied if
we take

```
ρ = (1 - C_T_12 * ξ₁^-2)^-1 * max(C_E_1 * abs(c[1]), C_E_2 * abs(c[2]))
```
"""
function norm_bound_Z(
    c::SVector{2,Acb},
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    if C.T_12 * ξ₁^-2 < 1
        return inv(1 - C.T_12 * ξ₁^-2) * max(C.E_1 * abs(c[1]), C.E_2 * abs(c[2])) * ξ₁^-v
    else
        throw(ErrorException("could not verify T_12 * ξ₁^-2 < 1"))
    end
end

"""
    norm_bound_Z_dλ(c, lambda, κ, ϵ, ξ₁, v, λ, C)

Compute a bound for the norm of `Z_dλ` based on Lemma
REF(lemma:Z-lambda-fixed-point-bounds).
"""
function norm_bound_Z_dλ(
    c::SVector{2,Acb},
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
    norms_Z::NormBounds_Y,
)
    num =
        max(C.E_1_dλ * abs(c[1]), C.E_2_dλ * abs(c[2])) * exp(Arb(-1)) / v +
        C.Z_dλ_1 * norms_Z.Z
    den = 1 - C.Z_dλ_2
    if den > 0
        return num / den
    else
        throw(ErrorException("could not compute bound for norm of Z_dλ"))
    end
end
