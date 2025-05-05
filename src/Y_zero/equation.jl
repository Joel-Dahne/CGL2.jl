"""
    cgl_linearization_equation(YZ, lambda, κ, ϵ, ξ, Q_hat, λ)
    cgl_linearization_equation(YZ, (lambda, κ, ϵ, Q_hat, λ), ξ)

Evaluate the right hand side of the forward ODE when written as a four
dimensional complex system. It is evaluated at the point `Y` and time
`ξ`.

For `λ.d != 1` there is a removable singularity at `ξ = 0`. To return
a finite value we in this case required that `Y[3] = Y[4] = 0`.
"""
function cgl_linearization_equation(YZ, lambda, κ, ϵ, ξ, Q_hat, λ::CGLParams)
    (; d, ω, σ, δ) = λ
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
    dZ = -A_inv * (C + J_N - lambda * I) * Y
    if iszero(ξ)
        @assert iszero(Z)
    else
        dZ -= (A_inv * B₁ * ξ + (d - 1) * I / ξ) * Z
    end

    return vcat(dY, dZ)
end

# For use with ODEProblem
cgl_linearization_equation(u, (lambda, κ, ϵ, Q_hat, λ), ξ) =
    cgl_linearization_equation(u, lambda, κ, ϵ, ξ, Q_hat, λ)

function cgl_linearization_equation_real(
    YZ_reim,
    lambda_real,
    lambda_imag,
    κ,
    ϵ,
    ξ,
    Q_hat,
    λ::CGLParams,
)
    (; d, ω, σ, δ) = λ
    Y_r = SVector(YZ_reim[1], YZ_reim[2])
    Y_i = SVector(YZ_reim[3], YZ_reim[4])
    Z_r = SVector(YZ_reim[5], YZ_reim[6])
    Z_i = SVector(YZ_reim[7], YZ_reim[8])

    a, b = Q_hat(ξ)
    N₁_a = -(a^2 + b^2)^(σ - 1) * (δ * (1 + 2σ) * a^2 + 2σ * a * b + δ * b^2)
    N₁_b = -(a^2 + b^2)^(σ - 1) * (a^2 + 2δ * σ * a * b + (1 + 2σ) * b^2)
    N₂_a = (a^2 + b^2)^(σ - 1) * ((1 + 2σ) * a^2 - 2δ * σ * a * b + b^2)
    N₂_b = -(a^2 + b^2)^(σ - 1) * (δ * a^2 - 2σ * a * b + δ * (1 + 2σ) * b^2)

    # Compute the matrices M1 and M2 following same approach as in CAPD
    begin
        N1_a_real = real(N₁_a)
        N1_a_imag = imag(N₁_a)
        N1_b_real = real(N₁_b)
        N1_b_imag = imag(N₁_b)
        N2_a_real = real(N₂_a)
        N2_a_imag = imag(N₂_a)
        N2_b_real = real(N₂_b)
        N2_b_imag = imag(N₂_b)

        M1_11_real = (-ϵ * (κ / σ + N1_a_real - lambda_real) - N2_a_real - ω) / (1 + ϵ^2)
        M1_11_imag = (-ϵ * (N1_a_imag - lambda_imag) - N2_a_imag) / (1 + ϵ^2)
        M1_12_real = (-κ / σ - ϵ * (N1_b_real - ω) - N2_b_real + lambda_real) / (1 + ϵ^2)
        M1_12_imag = (-ϵ * N1_b_imag - N2_b_imag + lambda_imag) / (1 + ϵ^2)
        M1_21_real = (κ / σ - ϵ * (N2_a_real + ω) + N1_a_real - lambda_real) / (1 + ϵ^2)
        M1_21_imag = (ϵ * N2_a_imag + N1_a_imag - lambda_imag) / (1 + ϵ^2)
        M1_22_real = (-ϵ * (κ / σ + N2_b_real - lambda_real) + N1_b_real - ω) / (1 + ϵ^2)
        M1_22_imag = (-ϵ * (N2_b_imag - lambda_imag) - N1_b_imag) / (1 + ϵ^2)

        M2_11 = κ / (1 + ϵ^2) * (-ϵ) * ξ - (d - 1) / ξ
        M2_12 = κ / (1 + ϵ^2) * (-1) * ξ
        M2_21 = -M2_12
        M2_22 = M2_11

        M1_real_cpp = @SMatrix[M1_11_real M1_12_real; M1_21_real M1_22_real]
        M1_imag_cpp = @SMatrix[M1_11_imag M1_12_imag; M1_21_imag M1_22_imag]
        M2_cpp = @SMatrix[M2_11 M2_12; M2_21 M2_22]
    end

    # Compute them using complex matrix arithmetic
    begin
        lambda = complex(lambda_real, lambda_imag)

        A_inv = @SMatrix[ϵ 1; -1 ϵ] / (1 + ϵ^2)
        B₁ = @SMatrix[κ 0; 0 κ]
        #B₂ = (d - 1) * @SMatrix[ϵ -1; 1 ϵ]

        C = @SMatrix[κ/σ -ω; ω κ/σ]

        J_N = @SMatrix[N₁_a N₁_b; N₂_a N₂_b]

        M1 = -A_inv * (C + J_N - lambda * I)
        M2 = -(A_inv * B₁ * ξ + (d - 1) * I / ξ)
    end
    # Check that they agree
    # TODO: Only compute these things once
    @assert real(M1) ≈ M1_real_cpp
    @assert imag(M1) ≈ M1_imag_cpp
    if !iszero(ξ)
        @assert M2 ≈ M2_cpp
    end

    # Compute components directly with same approach as in CAPD
    begin
        Y_r_1, Y_r_2 = Y_r
        Y_i_1, Y_i_2 = Y_i
        Z_r_1, Z_r_2 = Z_r
        Z_i_1, Z_i_2 = Z_i
        dZ_r_1 =
            M1_11_real * Y_r_1 + M1_12_real * Y_r_2 -
            (M1_11_imag * Y_i_1 + M1_12_imag * Y_i_2)
        dZ_r_2 =
            M1_21_real * Y_r_1 + M1_22_real * Y_r_2 -
            (M1_21_imag * Y_i_1 + M1_22_imag * Y_i_2)
        dZ_i_1 =
            M1_11_imag * Y_r_1 +
            M1_12_imag * Y_r_2 +
            M1_11_real * Y_i_1 +
            M1_12_real * Y_i_2
        dZ_i_2 =
            M1_21_imag * Y_r_1 +
            M1_22_imag * Y_r_2 +
            M1_21_real * Y_i_1 +
            M1_22_real * Y_i_2
        if iszero(ξ)
            @assert iszero(Z_r_1) && iszero(Z_r_2) && iszero(Z_i_1) && iszero(Z_i_2)
        else
            dZ_r_1 += M2_11 * Z_r_1 + M2_12 * Z_r_2
            dZ_r_2 += M2_21 * Z_r_1 + M2_22 * Z_r_2
            dZ_i_1 += M2_11 * Z_i_1 + M2_12 * Z_i_2
            dZ_i_2 += M2_21 * Z_i_1 + M2_22 * Z_i_2
        end
    end

    # Compute using matrix arithmetic
    dY_r = Z_r
    dY_i = Z_i
    dZ_r = real(M1) * Y_r - imag(M1) * Y_i
    dZ_i = imag(M1) * Y_r + real(M1) * Y_i
    if iszero(ξ)
        @assert iszero(Z_r) && iszero(Z_i)
    else
        dZ_r += M2 * Z_r
        dZ_i += M2 * Z_i
    end

    # Check that they agree
    @assert dZ_r ≈ SVector(dZ_r_1, dZ_r_2)
    @assert dZ_i ≈ SVector(dZ_i_1, dZ_i_2)

    return vcat(dY_r, dY_i, dZ_r, dZ_i)
end

# For use with ODEProblem
cgl_linearization_equation_real(u, (lambda_real, lambda_imag, κ, ϵ, Q_hat, λ), ξ) =
    cgl_linearization_equation_real(u, lambda_real, lambda_imag, κ, ϵ, ξ, Q_hat, λ)

function cgl_linearization_equation_taylor(
    Y_ξ₀::SVector{2,NTuple{2,Acb}},
    lambda::Acb,
    ν::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₀::Arb,
    λ::CGLParams{Arb};
    degree::Integer = 5,
)
    (; d, ω, σ, δ) = λ

    Y1 = AcbSeries(Y_ξ₀[1]; degree)
    Y2 = AcbSeries(Y_ξ₀[2]; degree)

    if iszero(ξ₀) && !isone(d) && !iszero(Y1[1]) && !iszero(Y2[1])
        return SVector(indeterminate(Y1), indeterminate(Y2))
    end

    A = SMatrix{2,2}(ϵ, 1, -1, ϵ)
    B₁ = SMatrix{2,2}(κ, 0, 0, κ)
    B₂ = (d - 1) * A
    C = SMatrix{2,2}(κ / σ, ω, -ω, κ / σ)

    for n = 0:2:(degree-2)
        if iszero(ξ₀)
            M =
                @SMatrix[
                    (ϵ * (n * κ + κ / σ - lambda) + ω) (n * κ + κ / σ - lambda - ϵ * ω);
                    (-n * κ - κ / σ + lambda + ϵ * ω) (ϵ * (n * κ + κ / σ - lambda) + ω)
                ] / ((n + 2) * (n + d) * (1 + ϵ^2))

            Y1[n+2], Y2[n+2] = -M * SVector(Y1[n], Y2[n])
        else
            error("case ξ₀ != 0 not implemented")
        end
    end

    return SVector(Y1, Y2)
end

function cgl_linearization_equation_dλ_taylor(
    Y_ξ₀_dλ::SVector{2,NTuple{2,Acb}},
    Y1::AcbSeries,
    Y2::AcbSeries,
    lambda::Acb,
    ν::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₀::Arb,
    λ::CGLParams{Arb};
    degree::Integer = 5,
)
    @assert Arblib.degree(Y1) == Arblib.degree(Y2) == degree

    (; d, ω, σ, δ) = λ

    Y1_dλ = AcbSeries(Y_ξ₀_dλ[1]; degree)
    Y2_dλ = AcbSeries(Y_ξ₀_dλ[2]; degree)

    if iszero(ξ₀) && !isone(d) && !iszero(Y1_dλ[1]) && !iszero(Y2_dλ[1])
        return SVector(indeterminate(Y1), indeterminate(Y2))
    end

    A = SMatrix{2,2}(ϵ, 1, -1, ϵ)
    B₁ = SMatrix{2,2}(κ, 0, 0, κ)
    B₂ = (d - 1) * A
    C = SMatrix{2,2}(κ / σ, ω, -ω, κ / σ)

    for n = 0:2:(degree-2)
        if iszero(ξ₀)
            inv_rhs = @SMatrix[ϵ 1; -1 ϵ] / ((n + 2) * (n + d) * (1 + ϵ^2))
            M =
                @SMatrix[
                    (ϵ * (n * κ + κ / σ - lambda) + ω) (n * κ + κ / σ - lambda - ϵ * ω);
                    (-n * κ - κ / σ + lambda + ϵ * ω) (ϵ * (n * κ + κ / σ - lambda) + ω)
                ] / ((n + 2) * (n + d) * (1 + ϵ^2))

            Y1_dλ[n+2], Y2_dλ[n+2] =
                -M * SVector(Y1_dλ[n], Y2_dλ[n]) + inv_rhs * SVector(Y1[n], Y2[n])
        else
            error("case ξ₀ != 0 not implemented")
        end
    end

    return SVector(Y1_dλ, Y2_dλ)
end
