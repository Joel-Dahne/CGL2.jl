"""
    cgl_linearization_equation(YZ, λ, κ, ϵ, ξ, Q_hat, Λ)
    cgl_linearization_equation(YZ, (λ, κ, ϵ, Q_hat, Λ), ξ)

Evaluate the right hand side of the forward ODE when written as a four
dimensional complex system. It is evaluated at the point `Y` and time
`ξ`.

For `Λ.d != 1` there is a removable singularity at `ξ = 0`. To return
a finite value we in this case required that `Y[3] = Y[4] = 0`.
"""
function cgl_linearization_equation(YZ, λ, κ, ϵ, ξ, Q_hat, Λ::CGLParams)
    (; d, ω, σ, δ) = Λ
    Y = SVector(YZ[1], YZ[2])
    Z = SVector(YZ[3], YZ[4])

    A_inv = @SMatrix[ϵ 1; -1 ϵ] / (1 + ϵ^2)
    B₁ = @SMatrix[κ 0; 0 κ]
    #B₂ = (d - 1) * @SMatrix[ϵ -1; 1 ϵ]

    C = @SMatrix[κ/σ -ω; ω κ/σ]

    J_N = let
        a, b = Q_hat(ξ)
        N₁_a = -(a^2 + b^2)^(σ - 1) * (δ * (1 + 2σ) * a^2 + 2σ * a * b + δ * b^2)
        N₁_b = -(a^2 + b^2)^(σ - 1) * (a^2 + 2δ * σ * a * b + (1 + 2σ) * b^2)
        N₂_a = (a^2 + b^2)^(σ - 1) * ((1 + 2σ) * a^2 - 2δ * σ * a * b + b^2)
        N₂_b = -(a^2 + b^2)^(σ - 1) * (δ * a^2 - 2σ * a * b + δ * (1 + 2σ) * b^2)

        @SMatrix[N₁_a N₁_b; N₂_a N₂_b]
    end

    dY = Z
    dZ = -A_inv * (C + J_N - λ * I) * Y
    if iszero(ξ)
        @assert iszero(Z)
    else
        dZ -= (A_inv * B₁ * ξ + (d - 1) * I / ξ) * Z
    end

    return vcat(dY, dZ)
end

# For use with ODEProblem
cgl_linearization_equation(u, (λ, κ, ϵ, Q_hat, Λ), ξ) =
    cgl_linearization_equation(u, λ, κ, ϵ, ξ, Q_hat, Λ)

"""
    _cgl_linearization_equation_taylor_J_N_taylor(ν, κ, ϵ, ξ₀, Λ; degree)

Compute a Taylor expansion of the term `J_N` in the equation for the
linear equation. It is done by computing the Taylor expansions of the
forward solution `Q_hat` and putting together `J_N` from this.
"""
function _cgl_linearization_equation_taylor_J_N_taylor(
    ν::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₀::Arb,
    Λ::CGLParams{Arb};
    degree::Integer = 5,
)
    (; d, ω, σ, δ) = Λ

    # Compute expansion for Q_hat.
    # We use the fact that the equation is identical to the one for Q,
    # except for the change in sign for κ and ω.
    a_hat, b_hat = cgl_equation_real_taylor(
        SVector{2,NTuple{2,Arb}}((real(ν), 0), (imag(ν), 0)),
        -κ,
        ϵ,
        zero(ξ₀),
        CGLParams(Λ, ω = -Λ.ω);
        degree,
    )

    J_N_11 =
        -(a_hat^2 + b_hat^2)^(σ - 1) *
        (δ * (1 + 2σ) * a_hat^2 + 2σ * a_hat * b_hat + δ * b_hat^2)
    J_N_12 =
        -(a_hat^2 + b_hat^2)^(σ - 1) *
        (a_hat^2 + 2δ * σ * a_hat * b_hat + (1 + 2σ) * b_hat^2)
    J_N_21 =
        (a_hat^2 + b_hat^2)^(σ - 1) *
        ((1 + 2σ) * a_hat^2 - 2δ * σ * a_hat * b_hat + b_hat^2)
    J_N_22 =
        -(a_hat^2 + b_hat^2)^(σ - 1) *
        (δ * a_hat^2 - 2σ * a_hat * b_hat + δ * (1 + 2σ) * b_hat^2)

    return @SMatrix[J_N_11 J_N_12; J_N_21 J_N_22]
end

"""
    cgl_linearization_equation_taylor(Y_0, λ, ν, κ, ϵ, Λ; degree)

Compute the Taylor expansions of `Y`. The expansion is centered at the
point `ξ = 0`, with the leading coefficients given by `Y_0` and the
linear are zero.
"""
function cgl_linearization_equation_taylor(
    Y_0::SVector{2,NTuple{2,Acb}},
    λ::Acb,
    ν::Acb,
    κ::Arb,
    ϵ::Arb,
    Λ::CGLParams{Arb};
    degree::Integer = 5,
)
    (; d, ω, σ, δ) = Λ

    Y1 = AcbSeries(Y_0[1]; degree)
    Y2 = AcbSeries(Y_0[2]; degree)

    J_N = _cgl_linearization_equation_taylor_J_N_taylor(ν, κ, ϵ, Arb(0), Λ; degree)

    for n = 0:2:(degree-2)
        v = J_N * SVector(AcbSeries(Y1, degree = n + 1), AcbSeries(Y2, degree = n + 1))

        # Explicit inverse of (n + 2) * ((n + 1) * A + B₂)
        inv_rhs = @SMatrix[ϵ 1; -1 ϵ] / ((n + 2) * (n + d) * (1 + ϵ^2))
        # Explicit value for above inverse multiplied with n * B₁ + C - λ * I
        M =
            @SMatrix[
                (ϵ * (n * κ + κ / σ - λ) + ω) (n * κ + κ / σ - λ - ϵ * ω);
                (-n * κ - κ / σ + λ + ϵ * ω) (ϵ * (n * κ + κ / σ - λ) + ω)
            ] / ((n + 2) * (n + d) * (1 + ϵ^2))

        Y1[n+2], Y2[n+2] = -M * SVector(Y1[n], Y2[n]) - inv_rhs * SVector(v[1][n], v[2][n])
    end

    return SVector(Y1, Y2)
end
