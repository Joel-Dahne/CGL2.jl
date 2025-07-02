"""
    Y_zero_float(Y₀, lambda, κ, ϵ, ξ₁, Q_hat, λ::CGLParams; tol::Float64 = 1e-11)

Compute the solution to the ODE on the interval ``[0, ξ₁]`` with
initial values given by `Y₀`. Returns a vector with four complex
values, where the first two are the values at `ξ₁` and the last two
are the derivatives.

The solution is computed using [`ODEProblem`](@ref). The computations
are always done in `ComplexF64`.
"""
function Y_zero_float(Y₀, lambda, κ, ϵ, ξ₁, Q_hat, λ::CGLParams; tol::Float64 = 1e-11)
    prob = ODEProblem{false}(
        cgl_linearization_equation,
        SVector{4,ComplexF64}(Y₀[1], Y₀[2], 0, 0),
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

function Y_zero_float_real(
    Y₀_real,
    Y₀_imag,
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
        SVector(Y₀_real[1], Y₀_real[2], Y₀_imag[1], Y₀_imag[2], 0, 0, 0, 0),
        (zero(ξ₁), ξ₁),
        (lambda_real, lambda_imag, κ, ϵ, Q_hat, λ),
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

    # We return the real result. This can be converted to the complex
    # version with.
    #complex.(
    #    SVector(sol.u[end][1], sol.u[end][2], sol.u[end][5], sol.u[end][6]),
    #    SVector(sol.u[end][3], sol.u[end][4], sol.u[end][7], sol.u[end][8]),
    #)

    return sol.u[end]
end

"""
    Y_zero_derivative_float(x, lambda, κ, ϵ, ξ₁, Q_hat, λ::CGLParams; tol::Float64 = 1e-11)

This function computes the derivative of [`Y_zero_float`](@ref) w.r.t.
the parameter `lambda`.
"""
function Y_zero_derivative_float(
    Y₀,
    lambda,
    κ,
    ϵ,
    ξ₁,
    Q_hat,
    λ::CGLParams;
    tol::Float64 = 1e-11,
)
    # IMPROVE: In principle it should be enough to compute only
    # the derivatives with respect to the real parts since
    # everything is analytic.

    # We compute the Jacobian using Y_zero_float_real since
    # ForwardDiff cannot differentiate through Complex directly.
    J = ForwardDiff.jacobian(
        SVector(real(lambda), imag(lambda)),
    ) do (lambda_real, lambda_imag)
        Y_zero_float_real(real(Y₀), imag(Y₀), lambda_real, lambda_imag, κ, ϵ, ξ₁, Q_hat, λ; tol)
    end

    return SVector(
        complex(J[1, 1], J[3, 1]), # ∂Y₁ / ∂λ
        complex(J[2, 1], J[4, 1]), # ∂Y₂ / ∂λ
        complex(J[5, 1], J[7, 1]), # ∂Z₁ / ∂λ
        complex(J[6, 1], J[8, 1]), # ∂Z₂ / ∂λ
    )
end

"""
    Y_zero_float_curve(Y₀, lambda, κ, ϵ, ξ₁, Q_hat, λ::CGLParams; tol::Float64 = 1e-11)

Similar to [`Y_hat_zero_float`](@ref) but returns the whole
solution object given by the ODE solver, instead of just the value at
the final point.
"""
function Y_zero_float_curve(
    Y₀,
    lambda,
    κ,
    ϵ,
    ξ₁,
    Q_hat,
    λ::CGLParams;
    Y₀_deriv = SVector{2,ComplexF64}(0, 0),
    ξ₀ = zero(ξ₁),
    tol::Float64 = 1e-11,
    saveat = [],
)
    prob = ODEProblem{false}(
        cgl_linearization_equation,
        SVector{4,ComplexF64}(Y₀[1], Y₀[2], Y₀_deriv[1], Y₀_deriv[2]),
        (ξ₀, ξ₁),
        (lambda, κ, ϵ, Q_hat, λ),
    )

    sol = solve(prob, Vern7(), abstol = tol, reltol = tol, verbose = false; saveat)

    return sol
end
