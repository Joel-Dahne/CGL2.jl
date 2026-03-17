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
    λ::CGLParams{Arb},
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
    # are related to Lemma REF(lemma:I_E_hat-I_P_hat-dgamma-bounds)
    # are checked in the computation of `C`, where the associated
    # constants are computed. The other conditions are checked in the
    # norm_bound_Q_hat_dγ₂ function.

    norms.Q_hat[] = norm_bound_Q_hat(γ₁, γ₂, κ, ϵ, ξ₁, λ, C)
    norms.Q_hat_dγ₂[] = norm_bound_Q_hat_dγ₂(κ, ϵ, ξ₁, λ, C, norms)

    return norms
end

"""
    norm_bound_Q_hat(γ₁, γ₂, κ, ϵ, ξ₁, λ, C)

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

The second inequality gives us a direct upper bound for `ρ`, which we
denote `ρ_bound`. For the first inequality we show that

```
C_P_hat * abs(γ₁) +
    E_hat * abs(γ₂) * exp(-real(c) * ξ₁^2) * ξ₁^(2 / σ - d) +
    C_T_hat * ξ₁^-2 * ρ^(2σ + 1) - ρ
```

has a unique root on the interval ``0 < ρ < ρ_bound``. The expression
is always positive at `ρ = 0`, so to the right of the root it will be
negative and satisfy the inequality. The zero itself is the smallest
possible `ρ` satisfying the inequality.

From Proposition REF(prop:Q-hat-fixed-point) we then have that the
norm of `Q_hat` is bounded by `ρ`.
"""
function norm_bound_Q_hat(
    γ₁::Acb,
    γ₂::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_hat,
)
    c = _c(κ, ϵ, λ)
    (; d, σ) = λ

    # Upper bound for ρ from second inequality
    ρ_bound = (2C.T_hat * M(σ) * ξ₁^-2)^(-1 / 2σ)

    isfinite(ρ_bound) || throw(ErrorException("could not compute bound for norm of Q_hat"))

    f(ρ) =
        C.P_hat * abs(γ₁) +
        C.E_hat * abs(γ₂) * exp(-real(c) * ξ₁^2) * ξ₁^(2 / σ - d) +
        C.T_hat * ξ₁^-2 * abspow(ρ, 2σ + 1) - ρ

    # Isolate roots
    roots, flags = ArbExtras.isolate_roots(f, Arf(0), lbound(ρ_bound))

    if length(roots) == 1 && only(flags)
        # Refine a little bit with bisection. This helps a lot with
        # improving the numerical stability.
        ρ_initial = ArbExtras.refine_root_bisection(f, only(roots)..., rtol = Arb(1e-2))

        ρ = ArbExtras.refine_root(f, Arb(ρ_l_initial))

        return ρ
    else
        throw(ErrorException("could not isolate root in ρ when computing norm of Q_hat"))
    end
end

"""
    norm_bound_Q_hat_dγ₂(γ₁, γ₂, κ, ϵ, ξ₁, λ, C, norms)

Compute a bound for the norm of `Q_hat_dγ₂` based on Lemma
REF(lemma:Q-hat-dgamma-bound).
"""
function norm_bound_Q_hat_dγ₂(
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_hat,
    norms::NormBounds_hat,
)
    c = _c(κ, ϵ, λ)
    (; d, σ) = λ

    num = C.E_hat * exp(-real(c) * ξ₁^2) * ξ₁^(2 / σ - d)
    den = 1 - (2σ + 1) * C.T_hat * ξ₁^-2 * norms.Q_hat^2σ

    if Arblib.ispositive(den)
        return num / den
    else
        throw(ErrorException("could not compute bound for norm of Q_hat_dγ₂"))
    end
end
