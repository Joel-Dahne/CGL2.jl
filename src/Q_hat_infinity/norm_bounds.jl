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
    C::FunctionBounds_hat;
    include_dγ₂::Bool = false,
)
    norms = NormBounds_hat()

    # TODO: Add check for lemma requirements

    norms.Q_hat[] = norm_bound_Q_hat(γ₁, γ₂, κ, ϵ, ξ₁, v, λ, C)
    norms.Q_hat_dγ₂[] = norm_bound_Q_hat_dγ₂(κ, ϵ, ξ₁, v, λ, C, norms)

    return norms
end

# IMPROVE: Add documentation and update so that constants are computed
# elsewhere.
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

    @assert -((2σ + 1) * v - 2 / σ + d - 4) * ξ₁^-2 < 2real(c)

    C_I_E_hat = C.J_P_hat / abs((2σ + 1) * v - 2)
    C_I_P_hat = C.J_P_hat / (2real(c) + ((2σ + 1) * v - 2 / σ + d - 4) * ξ₁^-2)

    C_T_hat = C.P_hat * C_I_E_hat + C.E_hat * C_I_P_hat * ξ₁^-2

    # Upper bounds from second inequality
    ρ_bound = (2C_T_hat * M(σ) * ξ₁^(-2 + 2σ * v))^(-1 / 2σ)

    isfinite(ρ_bound) || return indeterminate(ρ_bound)

    f(ρ) =
        C.P_hat * abs(γ₁) * ξ₁^-v +
        C.E_hat * abs(γ₂) * exp(-real(c) * ξ₁^2) * ξ₁^(2 / σ - d - v) +
        C_T_hat * ξ₁^(-2 + 2σ * v) * abspow(ρ, 2σ + 1) - ρ

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
        return indeterminate(ρ_bound)
    end
end

# IMPROVE: Add documentation and update so that constants are computed
# elsewhere.
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

    @assert -((2σ + 1) * v - 2 / σ + d - 4) * ξ₁^-2 < 2real(c)

    C_I_E_hat = C.J_P_hat / abs((2σ + 1) * v - 2)
    C_I_P_hat = C.J_P_hat / (2real(c) + ((2σ + 1) * v - 2 / σ + d - 4) * ξ₁^-2)

    C_T_hat = C.P_hat * C_I_E_hat + C.E_hat * C_I_P_hat * ξ₁^-2

    num = C.E_hat * exp(-real(c) * ξ₁^2) * ξ₁^(2 / σ - d - v)
    den = 1 - (2σ + 1) * C_T_hat * ξ₁^(-2 + 2σ * v) * norms.Q_hat^2σ

    return Arblib.ispositive(den) ? num / den : indeterminate(num)
end
