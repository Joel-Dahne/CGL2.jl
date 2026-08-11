"""
    Y_zero_float(Y₀, λ, κ, ϵ, ξ₁, Q_hat, Λ::CGLParams; tol::Float64 = 1e-11)

Compute the solution to the ODE on the interval ``[0, ξ₁]`` with
initial values given by `Y₀`. Returns a vector with four complex
values, where the first two are the values at `ξ₁` and the last two
are the derivatives.

The solution is computed using [`ODEProblem`](@ref). The computations
are always done in `ComplexF64`.

Note that this, contrary to the CAPD version, uses the complex
formulation of the ODE. This allows for testing that the two versions
agree.
"""
function Y_zero_float(Y₀, λ, κ, ϵ, ξ₁, Q_hat, Λ::CGLParams; tol::Float64 = 1e-11)
    prob = ODEProblem{false}(
        cgl_linearization_equation,
        SVector{4,ComplexF64}(Y₀[1], Y₀[2], 0, 0),
        (zero(ξ₁), ξ₁),
        (λ, κ, ϵ, Q_hat, Λ),
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
        verbose = DiffEqBase.DEVerbosity(SciMLLogging.None());
        unstable_check,
    )

    return sol.u[end]
end

"""
    Y_zero_float_curve(Y₀, λ, κ, ϵ, ξ₁, Q_hat, Λ::CGLParams; tol::Float64 = 1e-11)

Similar to [`Y_hat_zero_float`](@ref) but returns the whole
solution object given by the ODE solver, instead of just the value at
the final point.
"""
function Y_zero_float_curve(
    Y₀,
    λ,
    κ,
    ϵ,
    ξ₁,
    Q_hat,
    Λ::CGLParams;
    Y₀_deriv = SVector{2,ComplexF64}(0, 0),
    ξ₀ = zero(ξ₁),
    tol::Float64 = 1e-11,
    saveat = [],
)
    prob = ODEProblem{false}(
        cgl_linearization_equation,
        SVector{4,ComplexF64}(Y₀[1], Y₀[2], Y₀_deriv[1], Y₀_deriv[2]),
        (ξ₀, ξ₁),
        (λ, κ, ϵ, Q_hat, Λ),
    )

    sol = solve(
        prob,
        Vern7(),
        abstol = tol,
        reltol = tol,
        verbose = DiffEqBase.DEVerbosity(SciMLLogging.None());
        saveat,
    )

    return sol
end
