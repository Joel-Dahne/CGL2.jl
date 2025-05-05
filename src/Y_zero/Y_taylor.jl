"""
    Y_zero_taylor(x, lambda, ν, κ, ϵ, ξ₀, λ::CGLParams; degree = 20, enclose_curve = Val{false}())

Compute the solution to the ODE on the interval ``[0, ξ₀]``. Returns a
vector with four complex values, the first two are the values at `ξ₀`
and the second two are their derivatives.

The solution is computed using a Taylor expansion at `ξ = 0`. This
only works well for small values of `ξ₀` and is intended to be used
for handling the removable singularity at `ξ = 0`.

If `enclose_curve` is set to `Val{true}()` then don't return the value
at `ξ₀`, but an enclosure valid for `0 <= ξ <= ξ₀`. In this case it
also returns a second value containing enclosures of the real and
imaginary parts of the second order derivatives.
"""
function Y_zero_taylor(
    x::Acb,
    lambda::Acb,
    ν::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₀::Arb,
    λ::CGLParams{Arb};
    degree = 20,
    enclose_curve::Union{Val{false},Val{true}} = Val{false}(),
)
    # Compute expansion
    Y1, Y2 = cgl_linearization_equation_taylor(
        SVector{2,NTuple{2,Acb}}((x, 0), (1, 0)),
        lambda,
        ν,
        κ,
        ϵ,
        zero(ξ₀),
        λ;
        degree,
    )

    # FIXME: Implement remainder
    remainder, remainder_derivative, remainder_derivative2 = 0, 0, 0
    #_Q_zero_taylor_remainder(a, b, κ, ϵ, ξ₀, λ)

    if enclose_curve isa Val{true}
        Y10, Y11 = Arblib.evaluate2(Y1, Arb((0, ξ₀)))
        Y20, Y21 = Arblib.evaluate2(Y2, Arb((0, ξ₀)))
    else
        Y10, Y11 = Arblib.evaluate2(Y1, ξ₀)
        Y20, Y21 = Arblib.evaluate2(Y2, ξ₀)
    end

    Y10 += remainder
    Y11 += remainder_derivative
    Y20 += remainder
    Y21 += remainder_derivative

    if enclose_curve isa Val{true}
        Y12 = Arblib.derivative(Y1, 2)(Arb((0, ξ₀)))
        Y22 = Arblib.derivative(Y2, 2)(Arb((0, ξ₀)))
        Y12 += remainder_derivative2
        Y22 += remainder_derivative2

        return SVector(Y10, Y20, Y11, Y21), SVector(Y12, Y22)
    else
        return SVector(Y10, Y20, Y11, Y21)
    end
end

"""
    Y_zero_jacobian_taylor(x, lambda, ν, κ, ϵ, ξ₀, λ::CGLParams; degree = 20)

This function computes the Jacobian of [`Y_zero_taylor`](@ref) w.r.t.
the parameters `x` and `lambda`. It also returns the result of
[`Y_zero_taylor`](@ref).
"""
function Y_zero_jacobian_taylor(
    x::Acb,
    lambda::Acb,
    ν::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₀::Arb,
    λ::CGLParams{Arb};
    degree = 20,
)
    # Compute expansion
    Y1, Y2 = cgl_linearization_equation_taylor(
        SVector{2,NTuple{2,Acb}}((x, 0), (1, 0)),
        lambda,
        ν,
        κ,
        ϵ,
        zero(ξ₀),
        λ;
        degree,
    )

    # Compute expansion of derivative w.r.t. x
    # Note that this uses exactly the same equation, just different initial conditions.
    Y1_dx, Y2_dx = cgl_linearization_equation_taylor(
        SVector{2,NTuple{2,Acb}}((1, 0), (0, 0)),
        lambda,
        ν,
        κ,
        ϵ,
        zero(ξ₀),
        λ;
        degree,
    )

    # Compute expansion of derivative w.r.t. lambda
    Y1_dλ, Y2_dλ = cgl_linearization_equation_dλ_taylor(
        SVector{2,NTuple{2,Acb}}((0, 0), (0, 0)),
        Y1,
        Y2,
        lambda,
        ν,
        κ,
        ϵ,
        zero(ξ₀),
        λ;
        degree,
    )

    # FIXME: Implement remainder
    remainder, remainder_derivative = 0, 0
    remainder_dx, remainder_derivative_dx = 0, 0
    remainder_dλ, remainder_derivative_dλ = 0, 0

    Y10, Y11 = Arblib.evaluate2(Y1, ξ₀)
    Y20, Y21 = Arblib.evaluate2(Y2, ξ₀)
    Y10_dx, Y11_dx = Arblib.evaluate2(Y1_dx, ξ₀)
    Y20_dx, Y21_dx = Arblib.evaluate2(Y2_dx, ξ₀)
    Y10_dλ, Y11_dλ = Arblib.evaluate2(Y1_dλ, ξ₀)
    Y20_dλ, Y21_dλ = Arblib.evaluate2(Y2_dλ, ξ₀)

    Y10 += remainder
    Y11 += remainder_derivative
    Y20 += remainder
    Y21 += remainder_derivative
    Y10_dx += remainder_dx
    Y11_dx += remainder_derivative_dx
    Y20_dx += remainder_dx
    Y21_dx += remainder_derivative_dx
    Y10_dλ += remainder_dλ
    Y11_dλ += remainder_derivative_dλ
    Y20_dλ += remainder_dλ
    Y21_dλ += remainder_derivative_dλ

    Y = SVector(Y10, Y20, Y11, Y21)

    J = SMatrix{4,2}(Y10_dx, Y20_dx, Y11_dx, Y21_dx, Y10_dλ, Y20_dλ, Y11_dλ, Y21_dλ)

    return Y, J
end
