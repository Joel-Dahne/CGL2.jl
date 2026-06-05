"""
    NormBounds_hat(...)

Contains bounds for the norm of `Z`.

The bound is based on Lemma REF(lemma:Z-fixed-point-bounds) and
Proposition REF(prop:Y-fixed-point).

It checks all the conditions on the parameters that these lemmas
assume. If any of these conditions are not satisfied it will throw an
error. When using this struct the bounds can therefore safely be
assume to hold.

We use a struct for this, even though there is only one norm involved,
to keep the same structure as for `Q` and `Q_hat`.
"""
struct NormBounds_Y
    Z::Arb

    NormBounds_Y() = new(indeterminate(Arb))
end

function NormBounds_Y(
    c_0::SVector{2,Acb},
    λ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    norms = NormBounds_Y()

    # The requirements for Lemma REF(lemma:Z-fixed-point-bounds) are
    # checked in the computation of `C`, where the associated
    # constants are computed.

    # The requirements for the fixed point in Proposition
    # REF(prop:Z-fixed-point) are checked by the norm_bound_Z
    # function.

    norms.Z[] = norm_bound_Z(c_0, λ, κ, ϵ, ξ₁, Λ, C)

    return norms
end

"""
    norm_bound_Z(c_0, λ, κ, ϵ, ξ₁, Λ, C)

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

From Proposition REF(prop:Z-fixed-point) we then have that the norm of
`Z` is bounded by `ρ`.
"""
function norm_bound_Z(
    c_0::SVector{2,Acb},
    λ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    if C.T_12 * ξ₁^-2 < 1
        return inv(1 - C.T_12 * ξ₁^-2) * max(C.E_1 * abs(c_0[1]), C.E_2 * abs(c_0[2]))
    else
        throw(ErrorException("could not verify T_12 * ξ₁^-2 < 1"))
    end
end
