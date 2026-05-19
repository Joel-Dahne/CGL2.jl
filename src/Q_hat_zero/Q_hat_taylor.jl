"""
    Q_hat_zero_taylor(ν_real, ν_imag, κ, ϵ, ξ₀, Λ::CGLParams; degree = 20, enclose_curve = Val{false}())

Compute the solution to the ODE on the interval ``[0, ξ₀]``. Returns a
vector with four real values, the first two are the real and imaginary
values at `ξ₀` and the second two are their derivatives.

The solution is computed using a Taylor expansion at `ξ = 0`. This
only works well for small values of `ξ₀` and is intended to be used
for handling the removable singularity at `ξ = 0`.

If `enclose_curve` is set to `Val{true}()` then don't return the value
at `ξ₀`, but an enclosure valid for `0 <= ξ <= ξ₀`. In this case it
also returns a second value containing enclosures of the real and
imaginary parts of the second order derivatives.
"""
function Q_hat_zero_taylor(
    ν_real::Arb,
    ν_imag::Arb,
    κ::Arb,
    ϵ::Arb,
    ξ₀::Arb,
    Λ::CGLParams{Arb};
    degree = 20,
    enclose_curve::Union{Val{false},Val{true}} = Val{false}(),
)
    # Compute expansion
    # We use the fact that the equation is identical to the one for Q,
    # except for the change in sign for κ and ω.
    a, b = cgl_equation_real_taylor(
        SVector{2,NTuple{2,Arb}}((ν_real, 0), (ν_imag, 0)),
        -κ,
        ϵ,
        zero(ξ₀),
        CGLParams(Λ, ω = -Λ.ω);
        degree,
    )

    remainder, remainder_derivative, remainder_derivative2 =
        _Q_zero_taylor_remainder(a, b, -κ, ϵ, ξ₀, CGLParams(Λ, ω = -Λ.ω))

    if enclose_curve isa Val{true}
        a0, a1 = Arblib.evaluate2(a, Arb((0, ξ₀)))
        b0, b1 = Arblib.evaluate2(b, Arb((0, ξ₀)))
    else
        a0, a1 = Arblib.evaluate2(a, ξ₀)
        b0, b1 = Arblib.evaluate2(b, ξ₀)
    end

    a0 += remainder
    a1 += remainder_derivative
    b0 += remainder
    b1 += remainder_derivative

    if enclose_curve isa Val{true}
        a2 = Arblib.derivative(a, 2)(Arb((0, ξ₀)))
        b2 = Arblib.derivative(b, 2)(Arb((0, ξ₀)))
        a2 += remainder_derivative2
        b2 += remainder_derivative2

        return SVector(a0, b0, a1, b1), SVector(a2, b2)
    else
        return SVector(a0, b0, a1, b1)
    end
end

"""
    Q_hat_zero_jacobian_taylor(ν_real, ν_imag, κ, ϵ, ξ₀, Λ::CGLParams; degree = 20, enclose_curve = Val{false}())

This function computes the Jacobian of [`Q_hat_zero_taylor`](@ref)
w.r.t. the parameters `ν_real` and `ν_imag`. It also returns the
result of [`Q_hat_zero_taylor`](@ref).
"""
function Q_hat_zero_jacobian_taylor(
    ν_real::Arb,
    ν_imag::Arb,
    κ::Arb,
    ϵ::Arb,
    ξ₀::Arb,
    Λ::CGLParams{Arb};
    degree = 20,
)
    # Compute expansion
    # We use the fact that the equation is identical to the one for Q,
    # except for the change in sign for κ and ω.
    a, b = cgl_equation_real_taylor(
        SVector{2,NTuple{2,Arb}}((ν_real, 0), (ν_imag, 0)),
        -κ,
        ϵ,
        zero(ξ₀),
        CGLParams(Λ, ω = -Λ.ω);
        degree,
    )

    # Note that cgl_equation_real_dμ_taylor is the ODE given by
    # differentiating with respect to the initial conditions. We can
    # thus use the same ODE for the derivative w.r.t. ν_real and
    # ν_imag. The only difference is the initial conditions to use.

    # Compute expansion of derivative w.r.t. ν_real
    a_dν_real, b_dν_real = cgl_equation_real_dμ_taylor(
        SVector{2,NTuple{2,Arb}}((1, 0), (0, 0)),
        a,
        b,
        -κ,
        ϵ,
        zero(ξ₀),
        CGLParams(Λ, ω = -Λ.ω);
        degree,
    )

    # Compute expansion of derivative w.r.t. ν_real
    a_dν_imag, b_dν_imag = cgl_equation_real_dμ_taylor(
        SVector{2,NTuple{2,Arb}}((0, 0), (1, 0)),
        a,
        b,
        -κ,
        ϵ,
        zero(ξ₀),
        CGLParams(Λ, ω = -Λ.ω);
        degree,
    )

    remainder, remainder_derivative, _ =
        _Q_zero_taylor_remainder(a, b, -κ, ϵ, ξ₀, CGLParams(Λ, ω = -Λ.ω))

    remainder_dν_real, remainder_derivative_dν_real =
        _Q_zero_taylor_remainder_dμ(a_dν_real, b_dν_real, -κ, ϵ, ξ₀, CGLParams(Λ, ω = -Λ.ω))

    remainder_dν_imag, remainder_derivative_dν_imag =
        _Q_zero_taylor_remainder_dμ(a_dν_imag, b_dν_imag, -κ, ϵ, ξ₀, CGLParams(Λ, ω = -Λ.ω))

    a0, a1 = Arblib.evaluate2(a, ξ₀)
    b0, b1 = Arblib.evaluate2(b, ξ₀)
    a0_dν_real, a1_dν_real = Arblib.evaluate2(a_dν_real, ξ₀)
    b0_dν_real, b1_dν_real = Arblib.evaluate2(b_dν_real, ξ₀)
    a0_dν_imag, a1_dν_imag = Arblib.evaluate2(a_dν_imag, ξ₀)
    b0_dν_imag, b1_dν_imag = Arblib.evaluate2(b_dν_imag, ξ₀)

    a0 += remainder
    a1 += remainder_derivative
    b0 += remainder
    b1 += remainder_derivative
    a0_dν_real += remainder_dν_real
    a1_dν_real += remainder_derivative_dν_real
    b0_dν_real += remainder_dν_real
    b1_dν_real += remainder_derivative_dν_real
    a0_dν_imag += remainder_dν_imag
    a1_dν_imag += remainder_derivative_dν_imag
    b0_dν_imag += remainder_dν_imag
    b1_dν_imag += remainder_derivative_dν_imag

    Q_hat = SVector(a0, b0, a1, b1)

    J = SMatrix{4,2}(
        a0_dν_real,
        b0_dν_real,
        a1_dν_real,
        b1_dν_real,
        a0_dν_imag,
        b0_dν_imag,
        a1_dν_imag,
        b1_dν_imag,
    )

    return Q_hat, J
end
