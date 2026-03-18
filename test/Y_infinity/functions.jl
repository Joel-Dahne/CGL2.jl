@testset "functions" begin
    # IMPROVE: Test for more parameters
    lambda = Acb(0.19028080950252219, 2.769747313863597)
    γ₁ = Acb(0.1979686615344728, 0.15613711743772288)
    γ₂ = Acb(-116.03421868436834, 101.21494871936197)
    ξ = Arb(30)
    κ = Arb(0.8073018593981386)
    ϵ = Arb(0.15002213424487343)
    λ = CGLParams{Arb}(3, 1.0, 1.0, 0.0)
    (; d, ω, σ) = λ

    ξF64 = Float64(ξ)
    κF64 = Float64(κ)
    ϵF64 = Float64(ϵ)
    lambdaF64 = ComplexF64(lambda)
    λF64 = CGLParams{Float64}(λ)

    F_Y = CGL2.FunctionEnclosures_Y(lambda, γ₁, γ₂, κ, ϵ, ξ, λ)

    # Function for computing derivative using finite differences.
    fdm = central_fdm(5, 1)
    fdm2 = central_fdm(5, 2)
    fdm3 = central_fdm(5, 3)

    @testset "P_$j" for (j, P_j, P_j_dξ, P_j_dλ, P_j_dλ_dξ) in [
        (1, CGL2.P_1, CGL2.P_1_dξ, CGL2.P_1_dλ, CGL2.P_1_dλ_dξ),
        (2, CGL2.P_2, CGL2.P_2_dξ, CGL2.P_2_dλ, CGL2.P_2_dλ_dξ),
    ]
        Pj_series = P_j(ArbSeries((ξ, 1, 0)), lambda, κ, ϵ, λ)
        Pj = Pj_series[0]
        Pj_dξ = Pj_series[1]
        Pj_dξ_dξ = 2Pj_series[2]

        # Test that it solves the ODE
        if j == 1
            @test Arblib.contains_zero(
                (ϵ + im) * Pj_dξ_dξ +
                (κ * ξ + (d - 1) * (ϵ + im) / ξ) * Pj_dξ +
                (κ / σ + ω * im - lambda) * Pj,
            )
        else
            @test Arblib.contains_zero(
                (ϵ - im) * Pj_dξ_dξ +
                (κ * ξ + (d - 1) * (ϵ - im) / ξ) * Pj_dξ +
                (κ / σ - ω * im - lambda) * Pj,
            )
        end

        @test Arblib.overlaps(
            P_j(ArbSeries((ξ, 1)), lambda, κ, ϵ, λ)[1],
            P_j_dξ(ξ, lambda, κ, ϵ, λ),
        )

        @test Arblib.overlaps(
            P_j(ξ, AcbSeries((lambda, 1)), κ, ϵ, λ)[1],
            P_j_dλ(ξ, lambda, κ, ϵ, λ),
        )

        @test P_j_dξ(ξ, lambda, κ, ϵ, λ) ≈
              fdm(ξ -> P_j(ξ, lambdaF64, κF64, ϵF64, λF64), ξF64) rtol = 1e-12

        @test P_j_dλ(ξ, lambda, κ, ϵ, λ) ≈ fdm(
            lambda_real ->
                P_j(ξF64, complex(lambda_real, imag(lambdaF64)), κF64, ϵF64, λF64),
            real(lambdaF64),
        ) rtol = 1e-8

        @test P_j_dλ_dξ(ξ, lambda, κ, ϵ, λ) ≈ fdm(
            lambda_real ->
                P_j_dξ(ξF64, complex(lambda_real, imag(lambdaF64)), κF64, ϵF64, λF64),
            real(lambdaF64),
        ) rtol = 1e-9
    end

    @testset "E_$j" for (j, E_j, E_j_dξ, E_j_dλ, E_j_dλ_dξ) in [
        (1, CGL2.E_1, CGL2.E_1_dξ, CGL2.E_1_dλ, CGL2.E_1_dλ_dξ),
        (2, CGL2.E_2, CGL2.E_2_dξ, CGL2.E_2_dλ, CGL2.E_2_dλ_dξ),
    ]
        Ej_series = E_j(ArbSeries((ξ, 1, 0)), lambda, κ, ϵ, λ)
        Ej = Ej_series[0]
        Ej_dξ = Ej_series[1]
        Ej_dξ_dξ = 2Ej_series[2]

        # Test that it solves the ODE
        if j == 1
            @test Arblib.contains_zero(
                (ϵ + im) * Ej_dξ_dξ +
                (κ * ξ + (d - 1) * (ϵ + im) / ξ) * Ej_dξ +
                (κ / σ + ω * im - lambda) * Ej,
            )
        else
            @test Arblib.contains_zero(
                (ϵ - im) * Ej_dξ_dξ +
                (κ * ξ + (d - 1) * (ϵ - im) / ξ) * Ej_dξ +
                (κ / σ - ω * im - lambda) * Ej,
            )
        end

        @test Arblib.overlaps(
            E_j(ArbSeries((ξ, 1)), lambda, κ, ϵ, λ)[1],
            E_j_dξ(ξ, lambda, κ, ϵ, λ),
        )

        @test Arblib.overlaps(
            E_j(ξ, AcbSeries((lambda, 1)), κ, ϵ, λ)[1],
            E_j_dλ(ξ, lambda, κ, ϵ, λ),
        )

        @test E_j_dξ(ξ, lambda, κ, ϵ, λ) ≈
              fdm(ξ -> E_j(ξ, lambdaF64, κF64, ϵF64, λF64), ξF64) rtol = 1e-10

        @test E_j_dλ(ξ, lambda, κ, ϵ, λ) ≈ fdm(
            lambda_real ->
                E_j(ξF64, complex(lambda_real, imag(lambdaF64)), κF64, ϵF64, λF64),
            real(lambdaF64),
        ) rtol = 1e-8

        @test E_j_dλ_dξ(ξ, lambda, κ, ϵ, λ) ≈ fdm(
            lambda_real ->
                E_j_dξ(ξF64, complex(lambda_real, imag(lambdaF64)), κF64, ϵF64, λF64),
            real(lambdaF64),
        ) rtol = 1e-8
    end

    @testset "W_$j" for (
        j,
        W_j,
        P_j,
        P_j_dξ,
        P_j_dλ,
        P_j_dλ_dξ,
        E_j,
        E_j_dξ,
        E_j_dλ,
        E_j_dλ_dξ,
    ) in [
        (
            1,
            CGL2.W_1,
            CGL2.P_1,
            CGL2.P_1_dξ,
            CGL2.P_1_dλ,
            CGL2.P_1_dλ_dξ,
            CGL2.E_1,
            CGL2.E_1_dξ,
            CGL2.E_1_dλ,
            CGL2.E_1_dλ_dξ,
        ),
        (
            2,
            CGL2.W_2,
            CGL2.P_2,
            CGL2.P_2_dξ,
            CGL2.P_2_dλ,
            CGL2.P_2_dλ_dξ,
            CGL2.E_2,
            CGL2.E_2_dξ,
            CGL2.E_2_dλ,
            CGL2.E_2_dλ_dξ,
        ),
    ]
        @test Arblib.overlaps(
            W_j(ξ, lambda, κ, ϵ, λ),
            P_j(ξ, lambda, κ, ϵ, λ) * E_j_dξ(ξ, lambda, κ, ϵ, λ) -
            P_j_dξ(ξ, lambda, κ, ϵ, λ) * E_j(ξ, lambda, κ, ϵ, λ),
        )
    end

    @testset "J_P_$j" for (j, J_P_j, J_P_j_dξ, J_P_j_dλ, P_j, P_j_dλ, W_j) in [
        (1, CGL2.J_P_1, CGL2.J_P_1_dξ, CGL2.J_P_1_dλ, CGL2.P_1, CGL2.P_1_dλ, CGL2.W_1),
        (2, CGL2.J_P_2, CGL2.J_P_2_dξ, CGL2.J_P_2_dλ, CGL2.P_2, CGL2.P_2_dλ, CGL2.W_2),
    ]
        @test Arblib.overlaps(
            J_P_j(ξ, lambda, κ, ϵ, λ),
            P_j(ξ, lambda, κ, ϵ, λ) / W_j(ξ, lambda, κ, ϵ, λ),
        )

        @test Arblib.overlaps(
            J_P_j_dξ(ξ, lambda, κ, ϵ, λ),
            (P_j(
                AcbSeries((ξ, 1)),
                lambda,
                κ,
                ϵ,
                λ,
            )/W_j(AcbSeries((ξ, 1)), lambda, κ, ϵ, λ))[1],
        )

        @test Arblib.overlaps(
            J_P_j_dξ(ξ, lambda, κ, ϵ, λ),
            J_P_j(AcbSeries((ξ, 1)), lambda, κ, ϵ, λ)[1],
        )

        @test Arblib.overlaps(
            J_P_j_dλ(ξ, lambda, κ, ϵ, λ),
            (P_j(
                ξ,
                AcbSeries((lambda, 1)),
                κ,
                ϵ,
                λ,
            )/W_j(ξ, AcbSeries((lambda, 1)), κ, ϵ, λ))[1],
        )

        @test Arblib.overlaps(
            J_P_j_dλ(ξ, lambda, κ, ϵ, λ),
            J_P_j(ξ, AcbSeries((lambda, 1)), κ, ϵ, λ)[1],
        )
    end

    @testset "J_E_$j" for (j, J_E_j, J_E_j_dξ, J_E_j_dλ, E_j, E_j_dλ, W_j) in [
        (1, CGL2.J_E_1, CGL2.J_E_1_dξ, CGL2.J_E_1_dλ, CGL2.E_1, CGL2.E_1_dλ, CGL2.W_1),
        (2, CGL2.J_E_2, CGL2.J_E_2_dξ, CGL2.J_E_2_dλ, CGL2.E_2, CGL2.E_2_dλ, CGL2.W_2),
    ]
        @test Arblib.overlaps(
            J_E_j(ξ, lambda, κ, ϵ, λ),
            E_j(ξ, lambda, κ, ϵ, λ) / W_j(ξ, lambda, κ, ϵ, λ),
        )

        @test Arblib.overlaps(
            J_E_j_dξ(ξ, lambda, κ, ϵ, λ),
            (E_j(
                AcbSeries((ξ, 1)),
                lambda,
                κ,
                ϵ,
                λ,
            )/W_j(AcbSeries((ξ, 1)), lambda, κ, ϵ, λ))[1],
        )

        @test Arblib.overlaps(
            J_E_j_dξ(ξ, lambda, κ, ϵ, λ),
            J_E_j(AcbSeries((ξ, 1)), lambda, κ, ϵ, λ)[1],
        )

        @test Arblib.overlaps(
            J_E_j_dλ(ξ, lambda, κ, ϵ, λ),
            (E_j(
                ξ,
                AcbSeries((lambda, 1)),
                κ,
                ϵ,
                λ,
            )/W_j(ξ, AcbSeries((lambda, 1)), κ, ϵ, λ))[1],
        )

        @test Arblib.overlaps(
            J_E_j_dλ(ξ, lambda, κ, ϵ, λ),
            J_E_j(ξ, AcbSeries((lambda, 1)), κ, ϵ, λ)[1],
        )
    end

    @testset "K_1 and K_2" begin
        # K1 and K2 are suppose to give solutions to the linear system
        # Ψ * v = [[0, 0]; A \ F].

        A = SMatrix{2,2}(ϵ, 1, -1, ϵ)
        P = SMatrix{2,2,Acb}(im, 1, -im, 1)

        Ψ = [F_Y.E_12 F_Y.P_12; F_Y.E_12_dξ F_Y.P_12_dξ]
        F = eltype(Ψ)[0.1, 0.25]

        v = [-F_Y.K_1 * F; -F_Y.K_2 * F]

        @test all(Arblib.overlaps.(Ψ * v, [[0, 0]; inv(P) * A * P \ F]))

        K1_dξ = CGL2.K_1_dξ(ξ, lambda, κ, ϵ, λ)
        K2_dξ = CGL2.K_2_dξ(ξ, lambda, κ, ϵ, λ)

        K1_dξ_series =
            Diagonal(getindex.(diag(CGL2.K_1(AcbSeries((ξ, 1)), lambda, κ, ϵ, λ)), 1))
        K2_dξ_series =
            Diagonal(getindex.(diag(CGL2.K_2(AcbSeries((ξ, 1)), lambda, κ, ϵ, λ)), 1))

        K1_dξ_fdm = fdm(
            ξ_real ->
                CGL2.K_1(complex(ξ_real, imag(ξF64)), lambdaF64, κF64, ϵF64, λF64),
            real(ξF64),
        )

        K2_dξ_fdm = fdm(
            ξ_real ->
                CGL2.K_2(complex(ξ_real, imag(ξF64)), lambdaF64, κF64, ϵF64, λF64),
            real(ξF64),
        )

        @test all(Arblib.overlaps.(F_Y.K_1_dξ, K1_dξ_series))
        @test all(Arblib.overlaps.(F_Y.K_2_dξ, K2_dξ_series))
        @test ComplexF64.(F_Y.K_1_dξ) ≈ K1_dξ_fdm rtol = 1e-10
        @test ComplexF64.(F_Y.K_2_dξ) ≈ K2_dξ_fdm rtol = 1e-12

        K1_dλ = CGL2.K_1_dλ(ξ, lambda, κ, ϵ, λ)
        K2_dλ = CGL2.K_2_dλ(ξ, lambda, κ, ϵ, λ)

        K1_dλ_series =
            Diagonal(getindex.(diag(CGL2.K_1(ξ, AcbSeries((lambda, 1)), κ, ϵ, λ)), 1))
        K2_dλ_series =
            Diagonal(getindex.(diag(CGL2.K_2(ξ, AcbSeries((lambda, 1)), κ, ϵ, λ)), 1))

        K1_dλ_fdm = fdm(
            lambda_real ->
                CGL2.K_1(ξF64, complex(lambda_real, imag(lambdaF64)), κF64, ϵF64, λF64),
            real(lambdaF64),
        )

        K2_dλ_fdm = fdm(
            lambda_real ->
                CGL2.K_2(ξF64, complex(lambda_real, imag(lambdaF64)), κF64, ϵF64, λF64),
            real(lambdaF64),
        )

        @test all(Arblib.overlaps.(F_Y.K_1_dλ, K1_dλ_series))
        @test all(Arblib.overlaps.(F_Y.K_2_dλ, K2_dλ_series))
        @test ComplexF64.(F_Y.K_1_dλ) ≈ K1_dλ_fdm rtol = 1e-12
        @test ComplexF64.(F_Y.K_2_dλ) ≈ K2_dλ_fdm rtol = 1e-12
    end

    @testset "I_N" begin
        P = SMatrix{2,2,Acb}(im, 1, -im, 1)

        # The precise value for a and b should not play any role in
        # the correctness, we just compute some approximation here.
        νF64 = 1.9261384880241954 + 3.0638598354170337im
        a, b, a_dξ, b_dξ =
            Arb.(CGL2.Q_hat_zero_float(real(νF64), imag(νF64), κF64, ϵF64, ξF64, λF64))

        # Compute Jacobian by going through ArbSeries
        N = (a, b) -> (a^2 + b^2)^λ.σ * SVector(-λ.δ * a - b, a - λ.δ * b)
        N_a = getindex.(N(ArbSeries((a, 1)), b), 1)
        N_b = getindex.(N(a, ArbSeries((b, 1))), 1)
        JN_direct = [N_a N_b]
        IN_direct = inv(P) * JN_direct * P

        @test all(Arblib.overlaps.(IN_direct, CGL2.I_N(Acb(a, b), λ)))

        # Compute derivative w.r.t. ξ using formula for Jacobian plus ArbSeries
        @assert isone(λ.σ)
        @assert iszero(λ.δ)
        IN(a, b) = SMatrix{2,2}(
            2im * (a^2 + b^2),
            im * (a - im * b)^2,
            -im * (a + im * b)^2,
            -2im * (a^2 + b^2),
        )

        # Check that the above implementation agrees with previous one
        @test all(Arblib.overlaps.(IN_direct, IN(a, b)))

        # Compute derivative w.r.t. ξ
        IN_dξ = getindex.(IN(ArbSeries((a, a_dξ)), ArbSeries((b, b_dξ))), 1)

        @test all(Arblib.overlaps.(IN_dξ, CGL2.I_N_dξ(Acb(a, b), Acb(a_dξ, b_dξ), λ)))
    end
end
