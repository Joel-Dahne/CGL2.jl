"""
    Y_zero_float(Q_hat, κ, ϵ, ξ₁, λ::CGLParams; tol::Float64 = 1e-11)

Compute the solution to the ODE on the interval ``[0, ξ₁]``. Returns a
vector with four complex values, the first two are the values at `ξ₁`
and the second two are their derivatives.

The solution is computed using [`ODEProblem`](@ref). The computations
are always done in `ComplexF64`.
"""
function Y_zero_float(x, lambda, κ, ϵ, ξ₁, Q_hat, λ::CGLParams; tol::Float64 = 1e-11)
    prob = ODEProblem{false}(
        cgl_linearization_equation_real,
        SVector{4,ComplexF64}(x, 1, 0, 0),
        (zero(ξ₁), ξ₁),
        (lambda, κ, ϵ, Q_hat, λ),
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
    Y_zero_float_curve(Q_hat, κ, ϵ, ξ₁, λ::CGLParams; tol::Float64 = 1e-11)

Similar to [`Y_hat_zero_float`](@ref) but returns the whole
solution object given by the ODE solver, instead of just the value at
the final point.
"""
function Y_zero_float_curve(
    x,
    lambda,
    κ,
    ϵ,
    ξ₁,
    Q_hat,
    λ::CGLParams;
    y = one(x),
    Y₀_deriv = SVector{2,ComplexF64}(0, 0),
    ξ₀ = zero(ξ₁),
    tol::Float64 = 1e-11,
    saveat = [],
)
    prob = ODEProblem{false}(
        cgl_linearization_equation_real,
        SVector{4,ComplexF64}(x, y, Y₀_deriv[1], Y₀_deriv[2]),
        (ξ₀, ξ₁),
        (lambda, κ, ϵ, Q_hat, λ),
    )

    sol = solve(prob, Vern7(), abstol = tol, reltol = tol, verbose = false; saveat)

    return sol
end
