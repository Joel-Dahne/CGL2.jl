"""
    Q_hat_zero_capd(ν, κ, ϵ, ξ₁, λ::CGLParams; tol::Float64 = 1e-11)
    Q_hat_zero_capd(ν, κ, ϵ, ξ₀, ξ₁, λ::CGLParams; tol::Float64 = 1e-11)

Compute the solution to the ODE on the interval ``[0, ξ₁]``. Returns a
vector with four real values, the first two are the real and imaginary
values at `ξ₁` and the second two are their derivatives.

The solution is computed using the rigorous CAPD integrator.

If `ξ₀` is given then it uses a single Taylor expansion on the
interval `[0, ξ₀]` and CAPD on `[ξ₀, ξ₁]`. For `λ.d != 1` this is
automatically used (`ξ₀ = 1e-2` by default) to handle the removable
singularity at zero.

If the given `ξ₀` gives a non-finite enclosure on `[0, ξ₀]`, then it
tries with half that value. If it fails again it tries to halve it
once more, iterating like this for a maximum of a few times.
"""
Q_hat_zero_capd(
    ν_real::Arb,
    ν_imag::Arb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    tol::Float64 = 1e-11,
) = Q_hat_zero_capd(
    ν_real,
    ν_imag,
    κ,
    ϵ,
    ifelse(isone(λ.d), zero(Arb), Arb(1e-2)),
    ξ₁,
    λ;
    tol,
)

function Q_hat_zero_capd(
    ν_real::Arb,
    ν_imag::Arb,
    κ::Arb,
    ϵ::Arb,
    ξ₀::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    tol::Float64 = 1e-11,
)
    S = BareInterval{Float64}

    Q_hat_ξ₀ = if !iszero(ξ₀)
        @assert 0 < ξ₀ < ξ₁
        # Integrate system on [0, ξ₀] using Taylor expansion at zero
        Q_hat_ξ₀ = Q_hat_zero_taylor(ν_real, ν_imag, κ, ϵ, ξ₀, λ)
        if !all(isfinite, Q_hat_ξ₀)
            iterations = 0
            while !all(isfinite, Q_hat_ξ₀) && iterations < 5
                iterations += 1
                ξ₀ /= 2
                Q_hat_ξ₀ = Q_hat_zero_taylor(ν_real, ν_imag, κ, ϵ, ξ₀, λ)
            end
            iterations == 5 && @debug "Non-finite enclosure for smallest ξ₀" ξ₀
        end
        convert(SVector{4,S}, Q_hat_ξ₀)
    else
        SVector{4,S}(ν_real, ν_imag, bareinterval(0.0), bareinterval(0.0))
    end

    # Integrate system on [ξ₀, ξ₁] using capd.
    # We use the fact that the equation is identical to the one for Q,
    # except for the change in sign for κ and ω.
    Q_hat = _Q_zero_capd(
        Q_hat_ξ₀,
        convert(S, -κ),
        convert(S, ϵ),
        convert(S, ξ₀),
        convert(S, ξ₁),
        CGLParams{S}(λ, ω = -λ.ω);
        tol,
    )

    return Arb.(Q_hat)
end

"""
    Q_hat_zero_jacobian_capd(ν, κ, ϵ, ξ₁, λ::CGLParams; tol::Float64 = 1e-11)
    Q_hat_zero_jacobian_capd(ν, κ, ϵ, ξ₀, ξ₁, λ::CGLParams; tol::Float64 = 1e-11)

This function computes the Jacobian of [`Q_hat_zero_capd`](@ref)
w.r.t. the parameter `ν`.
"""
Q_hat_zero_jacobian_capd(
    ν_real::Arb,
    ν_imag::Arb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    tol::Float64 = 1e-11,
) = Q_hat_zero_jacobian_capd(
    ν_real,
    ν_imag,
    κ,
    ϵ,
    ifelse(isone(λ.d), zero(Arb), Arb(1e-2)),
    ξ₁,
    λ;
    tol,
)

function Q_hat_zero_jacobian_capd(
    ν_real::Arb,
    ν_imag::Arb,
    κ::Arb,
    ϵ::Arb,
    ξ₀::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    tol::Float64 = 1e-11,
)
    S = BareInterval{Float64}

    Q_hat_ξ₀, J_ξ₀ = let
        if !iszero(ξ₀)
            @assert 0 < ξ₀ < ξ₁
            # Integrate system on [0, ξ₀] using Taylor expansion at zero
            Q_hat_ξ₀, J_ξ₀ = Q_hat_zero_jacobian_taylor(ν_real, ν_imag, κ, ϵ, ξ₀, λ)
            if !all(isfinite, Q_hat_ξ₀) && !all(isfinite, J_ξ₀)
                iterations = 0
                while !all(isfinite, Q_hat_ξ₀) && iterations < 5
                    iterations += 1
                    ξ₀ /= 2
                    Q_hat_ξ₀, J_ξ₀ =
                        Q_hat_zero_jacobian_taylor(ν_real, ν_imag, κ, ϵ, ξ₀, λ)
                end
                iterations == 5 && @debug "Non-finite enclosure for smallest ξ₀" ξ₀
            end

            Q_hat_ξ₀ = convert(SVector{4,S}, Q_hat_ξ₀)
            J_ξ₀ = convert(SMatrix{4,2,S}, J_ξ₀)
        else
            Q_hat_ξ₀ = SVector{4,S}(ν_real, ν_imag, bareinterval(0.0), bareinterval(0.0))
            # Empty integration so the only non-zero derivatives are
            # the ones of Q_hat_ξ₀[1] and Q_hat_ξ₀[2] w.r.t. ν_real
            # and ν_imag, which are both 1.
            J_ξ₀ = SMatrix{4,2,S}(
                bareinterval(1.0),
                bareinterval(0.0),
                bareinterval(0.0),
                bareinterval(0.0),
                bareinterval(1.0),
                bareinterval(0.0),
                bareinterval(0.0),
                bareinterval(0.0),
            )
        end

        Q_hat_ξ₀, J_ξ₀
    end

    # Integrate system on [ξ₀, ξ₁] using capd.
    # We use the fact that the equation is identical to the one for Q,
    # except for the change in sign for κ and ω.
    J_ξ₀_ξ₁ = _Q_zero_capd(
        Q_hat_ξ₀,
        convert(S, -κ),
        convert(S, ϵ),
        convert(S, ξ₀),
        convert(S, ξ₁),
        CGLParams{S}(λ, ω = -λ.ω),
        output_jacobian_only_init = Val{true}();
        tol,
    )

    # The Jacobian on the interval [0, ξ₁] is the product of the one
    # on [0, ξ₀] and the one on [ξ₀, ξ₁].
    J = J_ξ₀_ξ₁ * J_ξ₀

    return Arb.(J)
end
