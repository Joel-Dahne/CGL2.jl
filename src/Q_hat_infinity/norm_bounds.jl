"""
    NormBounds_hat(...)

Contains bounds for the norms of `Q_hat` and `Q_hat_dγ₂`.

The bound for `Q_hat` is based on Lemmas REF(lemma:fixed-point-bounds)
and REF(prop:fixed-point). The bound for `Q_hat_dγ₂` is based on Lemma
REF(lemma:Q-hat-dgamma-bound).

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
    v::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_hat,
)
    norms = NormBounds_hat()

    # The requirements for Lemma REF(lemma:fixed-point-bounds) are
    # checked in the computation of `C`, where the associated
    # constants are computed.

    # The requirements for the fixed point in Lemma
    # REF(prop:fixed-point) are checked by the norm_bound_Q_hat
    # function.

    # The requirements for Lemma REF(lemma:Q-hat-dgamma-bound) are
    # checked by the norm_bound_Q_hat_dγ₂ function.

    norms.Q_hat[] = norm_bound_Q_hat(γ₁, γ₂, κ, ϵ, ξ₁, v, λ, C)
    norms.Q_hat_dγ₂[] = norm_bound_Q_hat_dγ₂(κ, ϵ, ξ₁, v, λ, C, norms)

    return norms
end

# IMPROVE: Add documentation
function norm_bound_Q_hat(
    γ₁::Acb,
    γ₂::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_hat,
)
    c = _c(κ, ϵ, λ)
    (; d, σ) = λ

    # Upper bounds from second inequality
    ρ_bound = (2C.T_hat * M(σ) * ξ₁^(-2 + 2σ * v))^(-1 / 2σ)

    isfinite(ρ_bound) || throw(ErrorException("could not compute bound for norm of Q_hat"))

    f(ρ) =
        C.P_hat * abs(γ₁) * ξ₁^-v +
        C.E_hat * abs(γ₂) * exp(-real(c) * ξ₁^2) * ξ₁^(2 / σ - d - v) +
        C.T_hat * ξ₁^(-2 + 2σ * v) * abspow(ρ, 2σ + 1) - ρ

    # Isolate roots
    roots, flags = ArbExtras.isolate_roots(f, Arf(0), ubound(ρ_bound))

    if length(roots) == 1 && only(flags)
        # Refine a little bit with bisection. This helps a lot with
        # improving the numerical stability.
        ρ_l_initial = ArbExtras.refine_root_bisection(f, only(roots)..., rtol = Arb(1e-2))

        return ArbExtras.refine_root(f, Arb(ρ_l_initial), strict = false)
    elseif length(roots) == 2 && all(flags)
        ρ_l_initial = ArbExtras.refine_root_bisection(f, roots[1]..., rtol = Arb(1e-2))

        return ArbExtras.refine_root(f, Arb(ρ_l_initial), strict = false)
    else
        throw(ErrorException("could not compute bound for norm of Q_hat"))
    end
end

# IMPROVE: Add documentation
function norm_bound_Q_hat_dγ₂(
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_hat,
    norms::NormBounds_hat,
)
    c = _c(κ, ϵ, λ)
    (; d, σ) = λ

    num = C.E_hat * exp(-real(c) * ξ₁^2) * ξ₁^(2 / σ - d - v)
    den = 1 - (2σ + 1) * C.T_hat * ξ₁^(-2 + 2σ * v) * norms.Q_hat^2σ

    if den > 0
        return num / den
    else
        throw(ErrorException("could not compute bound for norm of Q_hat_dγ₂"))
    end
end
