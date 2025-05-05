"""
    Y_zero_float(x, lambda, κ, ϵ, ξ₁, Q_hat, λ::CGLParams; tol::Float64 = 1e-11)

Compute the solution to the ODE on the interval ``[0, ξ₁]``. Returns a
vector with four complex values, the first two are the values at `ξ₁`
and the second two are their derivatives.

The solution is computed using [`ODEProblem`](@ref). The computations
are always done in `ComplexF64`.
"""
function Y_zero_float(x, lambda, κ, ϵ, ξ₁, Q_hat, λ::CGLParams; tol::Float64 = 1e-11)
    prob = ODEProblem{false}(
        cgl_linearization_equation,
        SVector{4,ComplexF64}(x, 1, 0, 0),
        (zero(ξ₁), ξ₁),
        (lambda, κ, ϵ, Q_hat, λ),
    )
    #return cgl_linearization_equation(SVector{4,ComplexF64}(x, 1, 0, 0), (lambda, κ, ϵ, Q_hat, λ), zero(ξ₁))
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

function Y_zero_float_real(
    x_real,
    x_imag,
    lambda_real,
    lambda_imag,
    κ,
    ϵ,
    ξ₁,
    Q_hat,
    λ::CGLParams;
    tol::Float64 = 1e-11,
)
    prob = ODEProblem{false}(
        cgl_linearization_equation_real,
        SVector(x_real, 1, x_imag, 0, 0, 0, 0, 0),
        (zero(ξ₁), ξ₁),
        (lambda_real, lambda_imag, κ, ϵ, Q_hat, λ),
    )
    #return cgl_linearization_equation_real(SVector(x_real, 1, x_imag, 0, 0, 0, 0, 0), (lambda_real, lambda_imag, κ, ϵ, Q_hat, λ), zero(ξ₁))
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

    # We return the real result. This can be converted to the complex
    # version with.
    #complex.(
    #    SVector(sol.u[end][1], sol.u[end][2], sol.u[end][5], sol.u[end][6]),
    #    SVector(sol.u[end][3], sol.u[end][4], sol.u[end][7], sol.u[end][8]),
    #)

    return sol.u[end]
end

"""
    Y_zero_jacobian_float(x, lambda, κ, ϵ, ξ₁, Q_hat, λ::CGLParams; tol::Float64 = 1e-11)

This function computes the Jacobian of [`Y_zero_float`](@ref) w.r.t.
the parameters `x` and `lambda`.
"""
function Y_zero_jacobian_float(
    x,
    lambda,
    κ,
    ϵ,
    ξ₁,
    Q_hat,
    λ::CGLParams;
    tol::Float64 = 1e-11,
)
    # We compute the Jacobian using Y_zero_float_real since
    # ForwardDiff cannot differentiate through Complex directly.
    J = ForwardDiff.jacobian(
        SVector(real(x), imag(x), real(lambda), imag(lambda)),
    ) do (x_real, x_imag, lambda_real, lambda_imag)
        Y_zero_float_real(x_real, x_imag, lambda_real, lambda_imag, κ, ϵ, ξ₁, Q_hat, λ; tol)
    end

    return SMatrix{4,2}(
        complex(J[1, 1], J[3, 1]), # ∂Y₁ / ∂x
        complex(J[2, 1], J[4, 1]), # ∂Y₂ / ∂x
        complex(J[5, 1], J[7, 1]), # ∂Z₁ / ∂x
        complex(J[6, 1], J[8, 1]), # ∂Z₂ / ∂x
        complex(J[1, 3], J[3, 3]), # ∂Y₁ / ∂λ
        complex(J[2, 3], J[4, 3]), # ∂Y₂ / ∂λ
        complex(J[5, 3], J[7, 3]), # ∂Z₁ / ∂λ
        complex(J[6, 3], J[8, 3]), # ∂Z₂ / ∂λ
    )
end

"""
    Y_zero_float_curve(x, lambda, κ, ϵ, ξ₁, Q_hat, λ::CGLParams; tol::Float64 = 1e-11)

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
        cgl_linearization_equation,
        SVector{4,ComplexF64}(x, y, Y₀_deriv[1], Y₀_deriv[2]),
        (ξ₀, ξ₁),
        (lambda, κ, ϵ, Q_hat, λ),
    )

    sol = solve(prob, Vern7(), abstol = tol, reltol = tol, verbose = false; saveat)

    return sol
end
