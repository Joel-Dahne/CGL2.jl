"""
    NormBounds_hat(...)

Contains bounds for the norms of `Q_hat` and `Q_hat_dγ₂`.

The bound for `Q_hat` is based on Lemma
REF(lemma:Q-hat-fixed-point-bounds) and Proposition
REF(prop:Q-hat-fixed-point). The bound for `Q_hat_dγ₂` is based on
Lemma REF(lemma:Q-hat-dgamma-bound).

It checks all the conditions on the parameters that these lemmas
assume. If any of these conditions are not satisfied it will throw an
error. When using this struct the bounds can therefore safely be
assume to hold.
"""
struct NormBounds_hat
    Q_hat::Arb
    Q_hat_dγ₂::Arb

    NormBounds_hat() = new(indeterminate(Arb), indeterminate(Arb))
end

function NormBounds_hat(
    γ₁::Acb,
    γ₂::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds_hat,
)
    norms = NormBounds_hat()

    # The requirements for Lemma REF(lemma:Q-hat-fixed-point-bounds)
    # are checked in the computation of `C`, where the associated
    # constants are computed.

    # The requirements for the fixed point in Proposition
    # REF(prop:Q-hat-fixed-point) are checked by the norm_bound_Q_hat
    # function.

    # The requirements for Lemma REF(lemma:Q-hat-dgamma-bound) that
    # are the same as for Lemma REF(lemma:Q-hat-fixed-point-bounds)
    # are already checked above. The condition on the numerator is
    # checked in the norm_bound_Q_hat_dγ₂ function.

    norms.Q_hat[] = norm_bound_Q_hat(γ₁, γ₂, κ, ϵ, ξ₁, Λ, C)
    norms.Q_hat_dγ₂[] = norm_bound_Q_hat_dγ₂(κ, ϵ, ξ₁, Λ, C, norms)

    return norms
end

"""
    norm_bound_Q_hat(γ₁, γ₂, κ, ϵ, ξ₁, Λ, C)

To apply the fixed point theorem in Proposition
REF(prop:Q-hat-fixed-point) we need to find `ρ` satisfying the
inequality

```
C_P_hat * abs(γ₁) +
    E_hat * abs(γ₂) * exp(-real(c) * ξ₁^2) * ξ₁^(2 / σ - d) +
    C_T_hat * ξ₁^-2 * ρ^(2σ + 1) <= ρ
```

and

```
2M_σ * C_T_hat * ξ₁^-2 * ρ^2σ < 1
```

The second inequality gives us a direct upper bound for `ρ`, we take
`ρ_bound` to be a value slightly less than this upper bound. For the
first inequality we show that

```
f(ρ) = C_P_hat * abs(γ₁) +
    E_hat * abs(γ₂) * exp(-real(c) * ξ₁^2) * ξ₁^(2 / σ - d) +
    C_T_hat * ξ₁^-2 * ρ^(2σ + 1) - ρ
```

has a unique root on the interval ``0 <= ρ < ρ_bound``. The expression
is always non-negative at `ρ = 0`, so to the right of the root it will
be negative and satisfy the inequality. The zero itself is the
smallest possible `ρ` satisfying the inequality.

To prove that there is a unique root on the interval ``0 <= ρ <
ρ_bound`` we note that

```
f'(ρ) = (2σ + 1) * C_T_hat * ξ₁^-2 * ρ^2σ - 1
```

has a unique root for `ρ > 0`. It follows that that `f(ρ)` has a
unique critical point. If `f(ρ_bound)` is negative it then follows
that `f(ρ)` has a unique root on the interval (possibly at `ρ = 0`).

From Proposition REF(prop:Q-hat-fixed-point) we have that the norm of
`Q_hat` is bounded by `ρ`.
"""
function norm_bound_Q_hat(
    γ₁::Acb,
    γ₂::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds_hat,
)
    c = _c(κ, ϵ, Λ)
    (; d, σ) = Λ

    # Upper bound for ρ from second inequality. We take a value that
    # is strictly lower than this (by eps(Arb)), so that we know that
    # the strict inequality is satisfied.
    ρ_bound = lbound((2C.T_hat * M(σ) * ξ₁^-2)^(-1 / 2σ) - eps(Arb))

    isfinite(ρ_bound) || return indeterminate(Arb)

    # Precompute constants
    w_1 = C.P_hat * abs(γ₁) + C.E_hat * abs(γ₂) * exp(-real(c) * ξ₁^2) * ξ₁^(2 / σ - d)
    w_2 = C.T_hat * ξ₁^-2
    f(ρ) = w_1 + w_2 * ρ^(2σ + 1) - ρ

    # Check that f is negative at the right endpoint.
    Arblib.isnegative(f(Arb(ρ_bound))) || return indeterminate(Arb)

    # We now know that there is a unique root on the interval 0 <= ρ <
    # ρ_bound.

    if Arblib.ispositive(f(Arb(0)))
        # When f(0) is positive we first get a rough enclosure using
        # bisection and then refine it using interval Newton.

        ρ_initial = ArbExtras.refine_root_bisection(f, Arf(0), ρ_bound, rtol = Arb(1e-3))

        return ArbExtras.refine_root(f, Arb(ρ_initial), strict = false)
    else
        # If f(0) contains zero we can't use refine_root_bisection
        # directly since it can determine the sign at the left
        # endpoint. Instead we take a uniform grid on the interval [0,
        # ρ_bound] and return the first value on which f is negative.
        # This doesn't give a tight enclosure, but is good enough
        # since this case doesn't occur in most cases.
        ρs = range(Arb(0), Arb(ρ_bound), 20)[2:(end-1)]

        i = findfirst(ρ -> Arblib.isnegative(f(ρ)), ρs)

        if isnothing(i)
            return Arb(ρ_bound)
        else
            return ρs[i]
        end
    end
end

"""
    norm_bound_Q_hat_dγ₂(γ₁, γ₂, κ, ϵ, ξ₁, Λ, C, norms)

Compute a bound for the norm of `Q_hat_dγ₂` based on Lemma
REF(lemma:Q-hat-dgamma-bound).
"""
function norm_bound_Q_hat_dγ₂(
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    Λ::CGLParams{Arb},
    C::FunctionBounds_hat,
    norms::NormBounds_hat,
)
    c = _c(κ, ϵ, Λ)
    (; d, σ) = Λ

    num = C.E_hat * exp(-real(c) * ξ₁^2) * ξ₁^(2 / σ - d)
    den = 1 - (2σ + 1) * C.T_hat * ξ₁^-2 * norms.Q_hat^2σ

    if Arblib.ispositive(den)
        return num / den
    else
        throw(ErrorException("could not compute bound for norm of Q_hat_dγ₂"))
    end
end
