function _C_T_hat_1(κ::Arb, ϵ::Arb, ξ₁::Arb, v::Arb, λ::CGLParams{Arb})
    (; d, σ) = λ
    @assert (2σ + 1) * v < 2 + 2 / σ - d
    @assert 2 / d < σ

    return C_P_hat(κ, ϵ, ξ₁, λ) * C_J_E_hat(κ, ϵ, ξ₁, λ) / abs((2σ + 1) * v - 2) +
           C_E_hat(κ, ϵ, ξ₁, λ) * C_J_P_hat(κ, ϵ, ξ₁, λ) / abs((2σ + 1) * v - 2 / σ + d - 2)
end


"""
    Q_hat_infinity_fixed_point(γ, κ, ϵ, ξ₁, v, λ::CGLParams)

Consider the fixed point problem that is used in the paper to get a
solution on ``[ξ₁, ∞)``. This function computes `ρ_l, ρ_u` such that
there exists a unique fixed point in the ball of radius `ρ_u` and it
is contained in the ball of radius `ρ_l`.

We have a unique fixed points in a ball of radius `ρ` if
```
C_P_hat * abs(c10_hat) * ξ₁^-v +
    C_E_hat * abs(c20_hat) * exp(-real(c) * ξ₁^2) * ξ₁^-(2 / σ - d - v) +
    C_T_hat_1 * ξ₁^(-2 + 2σ * v) * ρ^(2σ + 1) <= ρ
```
and
```
2C_T_hat_2 * ρ^2σ * ξ₁^(-2 + 2σ * v) < 1
```

The second inequality gives us a direct upper bound for `ρ`. For the
first inequality we find the zeros of
```
C_P_hat * abs(c10_hat) * ξ₁^-v +
    C_E_hat * abs(c20_hat) * exp(-real(c)) * ξ₁^-(2 / σ - d - v) +
    C_T_hat_1 * ξ₁^(-2 + 2σ * v) * ρ^(2σ + 1) - ρ
```
This is always positive at `ρ = 0` so it is enough to find the
smallest root to find `ρ_l`.
"""
function Q_hat_infinity_fixed_point(
    c10_hat::Acb,
    c20_hat::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb};
    throw_on_failure::Bool = false,
)
    (; d, σ) = λ
    @assert v > 0 # Required for the below bounds to be valid

    c = _c(κ, ϵ, λ)
    CP_hat = C_P_hat(κ, ϵ, ξ₁, λ)
    CE_hat = C_E_hat(κ, ϵ, ξ₁, λ)
    C_T_hat_1 = _C_T_hat_1(κ, ϵ, ξ₁, v, λ)
    C_T_hat_2 = M(σ) * C_T_hat_1

    # Upper from second inequality
    ρ_bound = (2C_T_hat_2 * ξ₁^(-2 + 2σ * v))^(-1 / 2σ)

    isfinite(ρ_bound) || return indeterminate(ρ_bound), indeterminate(ρ_bound)

    f(ρ) =
        CP_hat * abs(c10_hat) * ξ₁^-v +
        CE_hat * abs(c20_hat) * exp(-real(c) * ξ₁^2) * ξ₁^-(2 / σ - d - v) +
        +C_T_hat_1 * ξ₁^(-2 + 2σ * v) * abspow(ρ, 2σ + 1) - ρ

    # Isolate roots
    roots, flags = ArbExtras.isolate_roots(f, Arf(0), ubound(ρ_bound))

    if length(roots) == 1 && only(flags)
        # Refine a little bit with bisection. This helps a lot with
        # improving the numerical stability.
        ρ_l_initial = ArbExtras.refine_root_bisection(f, only(roots)..., rtol = Arb(1e-2))

        ρ_l = ArbExtras.refine_root(f, Arb(ρ_l_initial), strict = false)
        ρ_u = ρ_bound
    elseif length(roots) == 2 && all(flags)
        ρ_l_initial = ArbExtras.refine_root_bisection(f, roots[1]..., rtol = Arb(1e-2))
        ρ_u_initial = ArbExtras.refine_root_bisection(f, roots[2]..., rtol = Arb(1e-2))

        ρ_l = ArbExtras.refine_root(f, Arb(ρ_l_initial), strict = false)
        ρ_u = min(ρ_bound, ArbExtras.refine_root(f, Arb(ρ_u_initial)), strict = false)
    elseif throw_on_failure
        if isempty(roots)
            error("could not find any roots when bounding fixed point")
        elseif !all(flags)
            error("could not isolate roots when bounding fixed point")
        else
            error("found more than two roots when bounding fixed point")
        end
    else
        ρ_l = indeterminate(ρ_bound)
        ρ_u = indeterminate(ρ_bound)
    end

    return ρ_l, ρ_u
end
