"""
    _Y_zero_taylor_remainder(
        Y1::AcbSeries,
        Y2::AcbSeries,
        lambda::Acb,
        ν::Acb,
        κ::Arb,
        ϵ::Arb,
        ξ₀::Arb,
        λ::CGLParams{Arb},
    )

Compute an enclosure of the remainder term in [`Y_zero_taylor`](@ref).

Let `N` be the degree of the arguments `Y1` and `Y2`. The remainder
term is bounded by finding `r` such that `abs(Y1[n])` and `abs(Y2[n])`
are bounded by `r^n` for `n > N`. The remainder term is then bounded
as
```
abs(sum(Y1[n] * ξ₀^n for n = N+1:Inf)) <= (r * ξ₀)^(N + 1) / (1 - r * ξ₀)
```
and similarly for `Y2`. For the derivatives we instead get the bound
```
abs(sum(n * Y1[n] * ξ₀^(n - 1) for n = N+1:Inf)) <=
    (r * ξ₀)^N * (N + 1 - N * r * ξ₀) / (1 - r * ξ₀)^2
```
and the same for `Y2`.

For details on how we find `r` see lemma:tail-bound-Y in the paper.
"""
function _Y_zero_taylor_remainder(
    Y1::AcbSeries,
    Y2::AcbSeries,
    lambda::Acb,
    ν::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₀::Arb,
    λ::CGLParams{Arb},
)
    @assert Arblib.degree(Y1) == Arblib.degree(Y2)

    indeterminate_result = (indeterminate(κ), indeterminate(κ))

    isfinite(Y1) && isfinite(Y2) || return indeterminate_result

    (; d, ω, σ, δ) = λ

    # Compute expansion of forward solution
    a, b = cgl_equation_real_taylor(
        SVector{2,NTuple{2,Arb}}((real(ν), 0), (imag(ν), 0)),
        -κ,
        ϵ,
        zero(ξ₀),
        CGLParams(λ, ω = -λ.ω);
        degree = Arblib.degree(Y1),
    )

    N = Arblib.degree(Y1)

    if σ == 1
        # r such that abs(a[n]), abs(b[n]), abs(Y1[n]), abs(Y2[n]) < r^n for n > N
        r = let
            # Value of r is a tuning parameter. Lower value gives tighter
            # enclosures but makes it harder to verify the requirements.
            r = inv(16ξ₀)

            # Find C such that abs(a[n]), abs(b[n]), abs(Y1[n]), abs(Y2[n]) < C * r^n for 0 <= k <= N
            C = 1.01max(
                maximum(n -> abs(a[n] / r^n), 0:N),
                maximum(n -> abs(b[n] / r^n), 0:N),
                maximum(n -> abs(Y1[n] / r^n), 0:N),
                maximum(n -> abs(Y2[n] / r^n), 0:N),
            )

            # Find M such that abs(a[n]), abs(b[n]), abs(Y1[n]), abs(Y2[n]) <= r^n for M <= n <= N
            M =
                let M = findlast(
                        n -> !(
                            abs(a[n]) <= r^n &&
                            abs(b[n]) <= r^n &&
                            abs(Y1[n]) <= r^n &&
                            abs(Y2[n]) <= r^n
                        ),
                        0:N,
                    )
                    isnothing(M) ? 0 : M
                end

            # Verify that r, C and M satisfy the requirements
            all(n -> abs(a[n]) <= C * r^n, 0:(M-1)) || return indeterminate_result
            all(n -> abs(a[n]) <= r^n, M:N) || return indeterminate_result
            all(n -> abs(b[n]) <= C * r^n, 0:(M-1)) || return indeterminate_result
            all(n -> abs(b[n]) <= r^n, M:N) || return indeterminate_result
            all(n -> abs(Y1[n]) <= C * r^n, 0:(M-1)) || return indeterminate_result
            all(n -> abs(Y1[n]) <= r^n, M:N) || return indeterminate_result
            all(n -> abs(Y2[n]) <= C * r^n, 0:(M-1)) || return indeterminate_result
            all(n -> abs(Y2[n]) <= r^n, M:N) || return indeterminate_result

            # This is needed for the lemma to apply
            3M < N || return indeterminate_result

            D =
                (1 + abs(ϵ)) / (1 + ϵ^2) * (
                    abs(κ) / (N + d) +
                    (abs(κ / σ) + abs(lambda) + abs(ω)) / ((N + 2) * (N + d)) +
                    3(1 + abs(δ)) * (1 + 6M * C^3 / (N + d))
                )

            D <= r^2 || return indeterminate_result

            r
        end

        @assert 0 < r * ξ₀ < 1

        remainder_bound = (r * ξ₀)^(N + 1) / (1 - r * ξ₀)
        remainder_derivative_bound = (r * ξ₀)^N * (N + 1 - N * r * ξ₀) / (1 - r * ξ₀)^2

        remainder = add_error(Acb(0), remainder_bound)
        remainder_derivative = add_error(Acb(0), remainder_derivative_bound)
    else
        error("No implementation of remainder for σ != 1")
    end

    return remainder, remainder_derivative
end

"""
    _Y_zero_taylor_remainder_dλ(
        Y1::AcbSeries,
        Y2::AcbSeries,
        Y1_dλ::AcbSeries,
        Y2_dλ::AcbSeries,
        lambda::Acb,
        ν::Acb,
        κ::Arb,
        ϵ::Arb,
        ξ₀::Arb,
        λ::CGLParams{Arb},
    )

Compute an enclosure of the remainder term for the derivative w.r.t. λ
in [`Y_zero_jacobian_taylor`](@ref).

It works in the same way as [`_Y_zero_taylor_remainder`](@ref), except
it doesn't return a bound for the second derivative.

For details on how we find `r` see lemma:tail-bound-Y-dlambda in the paper.
"""
function _Y_zero_taylor_remainder_dλ(
    Y1::AcbSeries,
    Y2::AcbSeries,
    Y1_dλ::AcbSeries,
    Y2_dλ::AcbSeries,
    lambda::Acb,
    ν::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₀::Arb,
    λ::CGLParams{Arb},
)
    @assert Arblib.degree(Y1) ==
            Arblib.degree(Y2) ==
            Arblib.degree(Y1_dλ) ==
            Arblib.degree(Y2_dλ)

    # TODO: Implement from here!an

    indeterminate_result = (indeterminate(κ), indeterminate(κ))

    isfinite(Y1) && isfinite(Y2) || return indeterminate_result

    (; d, ω, σ, δ) = λ

    # Compute expansion of forward solution
    a, b = cgl_equation_real_taylor(
        SVector{2,NTuple{2,Arb}}((real(ν), 0), (imag(ν), 0)),
        -κ,
        ϵ,
        zero(ξ₀),
        CGLParams(λ, ω = -λ.ω);
        degree = Arblib.degree(Y1),
    )

    N = Arblib.degree(Y1)

    if σ == 1
        # r such that abs(Y1[n]), abs(Y2[n]) < r^n for n > N
        r = let
            # Value of r is a tuning parameter. Lower value gives tighter
            # enclosures but makes it harder to verify the requirements.
            r = inv(16ξ₀)

            # Find C such that abs(a[n]), abs(b[n]), abs(Y1[n]), abs(Y2[n]),
            # abs(Y1_dλ[n]), abs(Y2_dλ[n]) < C * r^n for 0 <= k <= N
            C = 1.01max(
                maximum(n -> abs(a[n] / r^n), 0:N),
                maximum(n -> abs(b[n] / r^n), 0:N),
                maximum(n -> abs(Y1[n] / r^n), 0:N),
                maximum(n -> abs(Y2[n] / r^n), 0:N),
                maximum(n -> abs(Y1_dλ[n] / r^n), 0:N),
                maximum(n -> abs(Y2_dλ[n] / r^n), 0:N),
            )

            # Find M such that abs(a[n]), abs(b[n]), abs(Y1[n]), abs(Y2[n]),
            # abs(Y1_dλ[n]), abs(Y2_dλ[n]) <= r^n for M <= n <= N
            M =
                let M = findlast(
                        n -> !(
                            abs(a[n]) <= r^n &&
                            abs(b[n]) <= r^n &&
                            abs(Y1[n]) <= r^n &&
                            abs(Y2[n]) <= r^n &&
                            abs(Y1_dλ[n]) <= r^n &&
                            abs(Y2_dλ[n]) <= r^n
                        ),
                        0:N,
                    )
                    isnothing(M) ? 0 : M
                end

            # Verify that r, C and M satisfy the requirements
            all(n -> abs(a[n]) <= C * r^n, 0:(M-1)) || return indeterminate_result
            all(n -> abs(a[n]) <= r^n, M:N) || return indeterminate_result
            all(n -> abs(b[n]) <= C * r^n, 0:(M-1)) || return indeterminate_result
            all(n -> abs(b[n]) <= r^n, M:N) || return indeterminate_result
            all(n -> abs(Y1[n]) <= C * r^n, 0:(M-1)) || return indeterminate_result
            all(n -> abs(Y1[n]) <= r^n, M:N) || return indeterminate_result
            all(n -> abs(Y2[n]) <= C * r^n, 0:(M-1)) || return indeterminate_result
            all(n -> abs(Y2[n]) <= r^n, M:N) || return indeterminate_result
            all(n -> abs(Y1_dλ[n]) <= C * r^n, 0:(M-1)) || return indeterminate_result
            all(n -> abs(Y1_dλ[n]) <= r^n, M:N) || return indeterminate_result
            all(n -> abs(Y2_dλ[n]) <= C * r^n, 0:(M-1)) || return indeterminate_result
            all(n -> abs(Y2_dλ[n]) <= r^n, M:N) || return indeterminate_result


            # This is needed for the lemma to apply
            3M < N || return indeterminate_result

            D =
                (1 + abs(ϵ)) / (1 + ϵ^2) * (
                    abs(κ) / (N + d) +
                    (abs(κ / σ) + abs(lambda) + abs(ω) + 1) / ((N + 2) * (N + d)) +
                    3(1 + abs(δ)) * (1 + 6M * C^3 / (N + d))
                )

            D <= r^2 || return indeterminate_result

            r
        end

        @assert 0 < r * ξ₀ < 1

        remainder_bound = (r * ξ₀)^(N + 1) / (1 - r * ξ₀)
        remainder_derivative_bound = (r * ξ₀)^N * (N + 1 - N * r * ξ₀) / (1 - r * ξ₀)^2

        remainder = add_error(Acb(0), remainder_bound)
        remainder_derivative = add_error(Acb(0), remainder_derivative_bound)
    else
        error("No implementation of remainder for σ != 1")
    end

    return remainder, remainder_derivative
end

"""
    Y_zero_taylor(Y₀, lambda, ν, κ, ϵ, ξ₀, λ::CGLParams; degree = 20)

Compute the solution to the ODE on the interval ``[0, ξ₀]``. Returns a
vector with four complex values, the first two are the values at `ξ₀`
and the second two are their derivatives.

The solution is computed using a Taylor expansion at `ξ = 0`. This
only works well for small values of `ξ₀` and is intended to be used
for handling the removable singularity at `ξ = 0`.
"""
function Y_zero_taylor(
    Y₀::SVector{2,Acb},
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
        SVector{2,NTuple{2,Acb}}((Y₀[1], 0), (Y₀[2], 0)),
        lambda,
        ν,
        κ,
        ϵ,
        zero(ξ₀),
        λ;
        degree,
    )

    remainder, remainder_derivative =
        _Y_zero_taylor_remainder(Y1, Y2, lambda, ν, κ, ϵ, ξ₀, λ)

    Y10, Y11 = Arblib.evaluate2(Y1, ξ₀)
    Y20, Y21 = Arblib.evaluate2(Y2, ξ₀)

    Y10 += remainder
    Y11 += remainder_derivative
    Y20 += remainder
    Y21 += remainder_derivative

    return SVector(Y10, Y20, Y11, Y21)
end

"""
    Y_zero_derivative_taylor(Y₀, lambda, ν, κ, ϵ, ξ₀, λ::CGLParams; degree = 20)

This function computes the derivative of [`Y_zero_taylor`](@ref)
w.r.t. the parameter `lambda`. It also returns the result of
[`Y_zero_taylor`](@ref).
"""
function Y_zero_derivative_taylor(
    Y₀::SVector{2,Acb},
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
        SVector{2,NTuple{2,Acb}}((Y₀[1], 0), (Y₀[2], 0)),
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

    remainder, remainder_derivative =
        _Y_zero_taylor_remainder(Y1, Y2, lambda, ν, κ, ϵ, ξ₀, λ)
    remainder_dλ, remainder_derivative_dλ =
        _Y_zero_taylor_remainder_dλ(Y1, Y2, Y1_dλ, Y2_dλ, lambda, ν, κ, ϵ, ξ₀, λ)

    Y10, Y11 = Arblib.evaluate2(Y1, ξ₀)
    Y20, Y21 = Arblib.evaluate2(Y2, ξ₀)
    Y10_dλ, Y11_dλ = Arblib.evaluate2(Y1_dλ, ξ₀)
    Y20_dλ, Y21_dλ = Arblib.evaluate2(Y2_dλ, ξ₀)

    Y10 += remainder
    Y11 += remainder_derivative
    Y20 += remainder
    Y21 += remainder_derivative
    Y10_dλ += remainder_dλ
    Y11_dλ += remainder_derivative_dλ
    Y20_dλ += remainder_dλ
    Y21_dλ += remainder_derivative_dλ

    Y = SVector(Y10, Y20, Y11, Y21)

    J = SVector(Y10_dλ, Y20_dλ, Y11_dλ, Y21_dλ)

    return Y, J
end
