"""
    Q_hat_zero_capd(ν_real, ν_imag, κ, ϵ, ξ₁, Λ::CGLParams; ξ₀, tol)

Compute the solution to the ODE on the interval ``[0, ξ₁]``. Returns a
vector with four real values, the first two are the real and imaginary
values at `ξ₁` and the second two are their derivatives.

The solution is computed using the rigorous CAPD integrator.

If `ξ₀` is non-zero it uses a single Taylor expansion on the interval
`[0, ξ₀]` and the CAPD integrator on `[ξ₀, ξ₁]`. This is needed to
avoid the removable singularity at `ξ = 0` which CAPD cannot handle
directly. For `Λ.d = 1` there is no removable singularity and the
default value is `ξ₀ = 0`, otherwise the default value is `ξ₀ = 1e-2`.

If the given `ξ₀` gives a non-finite enclosure on `[0, ξ₀]`, then it
tries with half that value. If it fails again it tries to halve it
once more, iterating like this for a maximum of a few times.
"""
function Q_hat_zero_capd(
    ν_real::Arb,
    ν_imag::Arb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    Λ::CGLParams{Arb};
    ξ₀::Arb = ifelse(isone(Λ.d), zero(Arb), Arb(1e-2)),
    tol::Float64 = 1e-11,
)
    Q_hat_ξ₀ = if !iszero(ξ₀)
        @assert 0 < ξ₀ < ξ₁
        # Integrate system on [0, ξ₀] using Taylor expansion at zero
        Q_hat_ξ₀ = Q_hat_zero_taylor(ν_real, ν_imag, κ, ϵ, ξ₀, Λ)
        if !all(isfinite, Q_hat_ξ₀)
            iterations = 0
            while !all(isfinite, Q_hat_ξ₀) && iterations < 5
                iterations += 1
                ξ₀ /= 2
                Q_hat_ξ₀ = Q_hat_zero_taylor(ν_real, ν_imag, κ, ϵ, ξ₀, Λ)
            end
            iterations == 5 && @debug "Non-finite enclosure for smallest ξ₀" ξ₀
        end
        Q_hat_ξ₀
    else
        SVector{4,Arb}(ν_real, ν_imag, 0, 0)
    end

    # Integrate system on [ξ₀, ξ₁] using CAPD.
    # We use the fact that the equation is identical to the one for Q,
    # except for the change in sign for κ and ω.
    return Q_hat = _Q_zero_capd(Q_hat_ξ₀, -κ, ϵ, ξ₀, ξ₁, CGLParams(Λ, ω = -Λ.ω); tol)
end

"""
    Q_hat_zero_jacobian_capd(ν_real, ν_imag, κ, ϵ, ξ₁, Λ::CGLParams; ξ₀, tol)

This function computes the Jacobian of [`Q_hat_zero_capd`](@ref)
w.r.t. the parameters `ν_real` and `ν_imag`.

In general it works similarly to [`Q_hat_zero_capd`](@ref).
"""
function Q_hat_zero_jacobian_capd(
    ν_real::Arb,
    ν_imag::Arb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    Λ::CGLParams{Arb};
    ξ₀::Arb = ifelse(isone(Λ.d), zero(Arb), Arb(1e-2)),
    tol::Float64 = 1e-11,
)
    Q_hat_ξ₀, J_ξ₀ = let
        if !iszero(ξ₀)
            @assert 0 < ξ₀ < ξ₁
            # Integrate system on [0, ξ₀] using Taylor expansion at zero
            Q_hat_ξ₀, J_ξ₀ = Q_hat_zero_jacobian_taylor(ν_real, ν_imag, κ, ϵ, ξ₀, Λ)
            if !all(isfinite, Q_hat_ξ₀) && !all(isfinite, J_ξ₀)
                iterations = 0
                while !all(isfinite, Q_hat_ξ₀) && iterations < 5
                    iterations += 1
                    ξ₀ /= 2
                    Q_hat_ξ₀, J_ξ₀ =
                        Q_hat_zero_jacobian_taylor(ν_real, ν_imag, κ, ϵ, ξ₀, Λ)
                end
                iterations == 5 && @debug "Non-finite enclosure for smallest ξ₀" ξ₀
            end
        else
            Q_hat_ξ₀ = SVector{4,Arb}(ν_real, ν_imag, 0, 0)
            # Empty integration so the only non-zero derivatives are
            # the ones of Q_hat_ξ₀[1] and Q_hat_ξ₀[2] w.r.t. ν_real
            # and ν_imag, which are both 1.
            J_ξ₀ = SMatrix{4,2,Arb}(1, 0, 0, 0, 1, 0, 0, 0)
        end

        Q_hat_ξ₀, J_ξ₀
    end

    # Integrate system on [ξ₀, ξ₁] using capd.
    # We use the fact that the equation is identical to the one for Q,
    # except for the change in sign for κ and ω.
    J_ξ₀_ξ₁ = _Q_zero_capd(
        Q_hat_ξ₀,
        -κ,
        ϵ,
        ξ₀,
        ξ₁,
        CGLParams(Λ, ω = -Λ.ω);
        output_jacobian = Val(true),
        include_parameter_derivatives = Val(false),
        tol,
    )

    # The Jacobian on the interval [0, ξ₁] is the product of the one
    # on [0, ξ₀] and the one on [ξ₀, ξ₁].
    return J_ξ₀_ξ₁ * J_ξ₀
end
