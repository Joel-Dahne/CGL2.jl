"""
    Q_hat_zero_float(ν_real, ν_imag, κ, ϵ, ξ₁, Λ::CGLParams; tol::Float64 = 1e-11)

Compute the solution to the ODE on the interval ``[0, ξ₁]``. Returns a
vector with four real values, the first two are the real and imaginary
values at `ξ₁` and the second two are their derivatives.

The solution is computed using [`ODEProblem`](@ref). The computations
are always done in `Float64`. However, for `Arb` input with wide
intervals for `ν`, `κ` and/or `ϵ` it computes it at the corners of the
box they form. This means you still get something that resembles an
enclosure.
"""
function Q_hat_zero_float(ν_real, ν_imag, κ, ϵ, ξ₁, Λ::CGLParams; tol::Float64 = 1e-11)
    # We use the fact that the equation is identical to the one for Q,
    # except for the change in sign for κ and ω.
    prob = ODEProblem{false}(
        cgl_equation_real,
        SVector(ν_real, ν_imag, 0, 0),
        (zero(ξ₁), ξ₁),
        (-κ, ϵ, CGLParams(Λ, ω = -Λ.ω)),
    )

    # Used to exit early in extreme cases. Sometimes the
    # refine_approximation methods en up in a bad region and the
    # solutions are highly oscillating, which makes the solver
    # extremely slow. For that reason we exit early in case the norm
    # of the solution grows to large.
    unstable_check = (dt, u, p, t) -> any(isnan, u) || norm(u) > 100

    sol = solve(
        prob,
        Vern7(),
        abstol = tol,
        reltol = tol,
        save_everystep = false,
        verbose = false;
        unstable_check,
    )

    return sol.u[end]
end

"""
    Q_hat_zero_jacobian_float(ν_real, ν_imag, κ, ϵ, ξ₁, Λ::CGLParams; tol::Float64 = 1e-11)

This function computes the Jacobian of [`Q_hat_zero_float`](@ref) w.r.t.
the parameters `ν_real` and `ν_imag`.
"""
function Q_hat_zero_jacobian_float(
    ν_real,
    ν_imag,
    κ,
    ϵ,
    ξ₁,
    Λ::CGLParams;
    tol::Float64 = 1e-11,
)
    return ForwardDiff.jacobian(SVector(ν_real, ν_imag)) do (ν_real, ν_imag)
        Q_hat_zero_float(ν_real, ν_imag, κ, ϵ, ξ₁, Λ; tol)
    end
end

"""
    Q_hat_zero_float_curve(ν_real, ν_imag, κ, ϵ, ξ₁, Λ::CGLParams; tol::Float64 = 1e-11)

Similar to [`Q_hat_zero_float`](@ref) but returns the whole
solution object given by the ODE solver, instead of just the value at
the final point.
"""
function Q_hat_zero_float_curve(
    ν_real,
    ν_imag,
    κ,
    ϵ,
    ξ₁,
    Λ::CGLParams;
    tol::Float64 = 1e-11,
    saveat = [],
)
    # We use the fact that the equation is identical to the one for Q,
    # except for the change in sign for κ and ω.
    prob = ODEProblem{false}(
        cgl_equation_real,
        SVector(ν_real, ν_imag, 0, 0),
        (zero(ξ₁), ξ₁),
        (-κ, ϵ, CGLParams(Λ, ω = -Λ.ω)),
    )

    sol = solve(prob, Vern7(), abstol = tol, reltol = tol, verbose = false; saveat)

    return sol
end
