"""
    _Y_zero_taylor_remainder_check_conditions(M, N, C, r, a_hat, b_hat, Y1, Y2, κ, ϵ, Λ)

Check that `M`,`N`, `C` and `r` satisfy the conditions of Lemma
REF(lemma:tail-bound-Y).
"""
function _Y_zero_taylor_remainder_check_conditions(
    M::Int,
    N::Int,
    C::Arb,
    r::Arb,
    a_hat::ArbSeries,
    b_hat::ArbSeries,
    Y1::AcbSeries,
    Y2::AcbSeries,
    λ::Acb,
    κ::Arb,
    ϵ::Arb,
    Λ::CGLParams{Arb},
)
    # Ensure the series have coefficient up to N computed
    @assert Arblib.degree(a_hat) ==
            Arblib.degree(b_hat) ==
            Arblib.degree(Y1) ==
            Arblib.degree(Y2) ==
            N

    (; d, ω, σ, δ) = Λ
    @assert isone(σ) # The lemma is only valid for σ = 1

    iseven(N) || return false
    3M < N || return false

    all(n -> abs(a_hat[n]) <= C * r^n, 0:(M-1)) || return false
    all(n -> abs(b_hat[n]) <= C * r^n, 0:(M-1)) || return false
    all(n -> abs(a_hat[n]) <= r^n, M:N) || return false
    all(n -> abs(b_hat[n]) <= r^n, M:N) || return false
    all(n -> abs(Y1[n]) <= C * r^n, 0:(M-1)) || return false
    all(n -> abs(Y2[n]) <= C * r^n, 0:(M-1)) || return false
    all(n -> abs(Y1[n]) <= r^n, M:N) || return false
    all(n -> abs(Y2[n]) <= r^n, M:N) || return false

    m = (M + 1) ÷ 2 # This is the same as ceil(M / 2)

    D =
        (1 + abs(ϵ)) / (1 + ϵ^2) * (
            abs(κ) / (N + d) +
            (abs(ω) + 2abs(κ) * abs(λ)) / ((N + 2) * (N + d)) +
            6(1 + abs(δ)) * (
                1 // 8 +
                1 // 2N +
                3m * C * ((1 + 3 // N) // 2(N + d)) +
                3m^2 * C^2 / ((N + 2) * (N + d))
            )
        )

    return D <= r^2
end

"""
    _Y_zero_taylor_remainder(
        Y1::AcbSeries,
        Y2::AcbSeries,
        λ::Acb,
        ν::Acb,
        κ::Arb,
        ϵ::Arb,
        ξ₀::Arb,
        Λ::CGLParams{Arb},
    )

Compute an enclosure of the remainder term in [`Y_zero_taylor`](@ref).

The approach is similar to [`Q_zero_taylor`](@ref) and is based on
Lemma REF(lemma:tail-bound-Y).
"""
function _Y_zero_taylor_remainder(
    Y1::AcbSeries,
    Y2::AcbSeries,
    λ::Acb,
    ν::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₀::Arb,
    Λ::CGLParams{Arb},
)
    @assert Arblib.degree(Y1) == Arblib.degree(Y2)

    (; d, ω, σ, δ) = Λ
    N = Arblib.degree(Y1)

    isone(σ) || error("No implementation of remainder for σ != 1")

    isfinite(Y1) && isfinite(Y2) || return indeterminate(Arb), indeterminate(Arb)

    # Compute expansion for Q_hat.
    # We use the fact that the equation is identical to the one for Q,
    # except for the change in sign for κ and ω.
    a_hat, b_hat = cgl_equation_real_taylor(
        SVector{2,NTuple{2,Arb}}((real(ν), 0), (imag(ν), 0)),
        -κ,
        ϵ,
        zero(ξ₀),
        CGLParams(Λ, ω = -Λ.ω);
        degree = N,
    )

    # Value of r is a tuning parameter. Lower value gives tighter
    # enclosures but makes it harder to verify the requirements.
    r = inv(16ξ₀)
    @assert 0 < r * ξ₀ < 1 # Sanity check

    # Find C such that abs(a[n]), abs(b[n]) < C * r^n for 0 <= k <= N
    C = ubound(
        Arb,
        1.01max(
            maximum(n -> abs(a_hat[n] / r^n), 0:N),
            maximum(n -> abs(b_hat[n] / r^n), 0:N),
            maximum(n -> abs(Y1[n] / r^n), 0:N),
            maximum(n -> abs(Y2[n] / r^n), 0:N),
        ),
    )

    M =
        let M = findlast(
                n -> !(
                    abs(a_hat[n]) <= r^n &&
                    abs(b_hat[n]) <= r^n &&
                    abs(Y1[n]) <= r^n &&
                    abs(Y2[n]) <= r^n
                ),
                0:N,
            )
            isnothing(M) ? 0 : M
        end

    # Check that the conditions of the lemma are satisfied
    if !_Y_zero_taylor_remainder_check_conditions(
        M,
        N,
        C,
        r,
        a_hat,
        b_hat,
        Y1,
        Y2,
        λ,
        κ,
        ϵ,
        Λ,
    )
        return indeterminate(Arb), indeterminate(Arb)
    end

    remainder_bound = (r * ξ₀)^(N + 1) / (1 - r * ξ₀)
    remainder_derivative_bound = r * (r * ξ₀)^N * (N + 1 - N * r * ξ₀) / (1 - r * ξ₀)^2

    remainder = add_error(Arb(0), remainder_bound)
    remainder_derivative = add_error(Arb(0), remainder_derivative_bound)

    return remainder, remainder_derivative
end

"""
    Y_zero_taylor(Y₀, λ, ν, κ, ϵ, ξ₀, Λ::CGLParams; degree = 20)

Compute the solution to the ODE on the interval ``[0, ξ₀]``. Returns a
vector with four complex values, the first two are the values at `ξ₀`
and the second two are their derivatives.

The solution is computed using a Taylor expansion at `ξ = 0`. This
only works well for small values of `ξ₀` and is intended to be used
for handling the removable singularity at `ξ = 0`.
"""
function Y_zero_taylor(
    Y₀::SVector{2,Acb},
    λ::Acb,
    ν::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₀::Arb,
    Λ::CGLParams{Arb};
    degree = 20,
)
    # Compute expansion
    Y1, Y2 = cgl_linearization_equation_taylor(
        SVector{2,NTuple{2,Acb}}((Y₀[1], 0), (Y₀[2], 0)),
        λ,
        ν,
        κ,
        ϵ,
        Λ;
        degree,
    )

    remainder, remainder_derivative = _Y_zero_taylor_remainder(Y1, Y2, λ, ν, κ, ϵ, ξ₀, Λ)

    Y10, Y11 = Arblib.evaluate2(Y1, ξ₀)
    Y20, Y21 = Arblib.evaluate2(Y2, ξ₀)

    Y10 += remainder
    Y11 += remainder_derivative
    Y20 += remainder
    Y21 += remainder_derivative

    return SVector(Y10, Y20, Y11, Y21)
end
