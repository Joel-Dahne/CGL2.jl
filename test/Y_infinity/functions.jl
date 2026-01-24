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

    (; A, B₁, B₂, C, λI) = CGL2.coeff_matrices(lambda, κ, ϵ, λ)

    F_Y = CGL2.FunctionEnclosures_Y_new(lambda, γ₁, γ₂, κ, ϵ, ξ, λ)

    # Function for computing derivative using finite differences.
    fdm = central_fdm(5, 1)
    fdm2 = central_fdm(5, 2)
    fdm3 = central_fdm(5, 3)

    Y_1, Y_1_dξ, Y_1_dξ_dξ = CGL2.Y_1, CGL2.Y_1_dξ, CGL2.Y_1_dξ_dξ
    Y_2, Y_2_dξ, Y_2_dξ_dξ = CGL2.Y_2, CGL2.Y_2_dξ, CGL2.Y_2_dξ_dξ
    Y_3, Y_3_dξ, Y_3_dξ_dξ = CGL2.Y_3, CGL2.Y_3_dξ, CGL2.Y_3_dξ_dξ
    Y_4, Y_4_dξ, Y_4_dξ_dξ = CGL2.Y_4, CGL2.Y_4_dξ, CGL2.Y_4_dξ_dξ

    Y_1_dλ, Y_1_dλ_dξ = CGL2.Y_1_dλ, CGL2.Y_1_dλ_dξ
    Y_2_dλ, Y_2_dλ_dξ = CGL2.Y_2_dλ, CGL2.Y_2_dλ_dξ
    Y_3_dλ, Y_3_dλ_dξ = CGL2.Y_3_dλ, CGL2.Y_3_dλ_dξ
    Y_4_dλ, Y_4_dλ_dξ = CGL2.Y_4_dλ, CGL2.Y_4_dλ_dξ

    @testset "Y1" begin
        @test all(
            Arblib.overlaps.(
                getindex.(Y_1(ArbSeries((ξ, 1)), lambda, κ, ϵ, λ), 1),
                Y_1_dξ(ξ, lambda, κ, ϵ, λ),
            ),
        )
        @test all(
            Arblib.overlaps.(
                2getindex.(Y_1(ArbSeries((ξ, 1, 0)), lambda, κ, ϵ, λ), 2),
                Y_1_dξ_dξ(ξ, lambda, κ, ϵ, λ),
            ),
        )

        @test Y_1_dξ(ξ, lambda, κ, ϵ, λ) ≈
              fdm(ξ -> Y_1(ξ, lambdaF64, κF64, ϵF64, λF64), ξF64) rtol = 1e-10
        @test Y_1_dξ_dξ(ξ, lambda, κ, ϵ, λ) ≈
              fdm2(ξ -> Y_1(ξ, lambdaF64, κF64, ϵF64, λF64), ξF64) rtol = 1e-7

        # Test that is solves equation
        @test all(
            Arblib.contains_zero.(
                A * Y_1_dξ_dξ(ξ, lambda, κ, ϵ, λ) +
                (B₁ * ξ + B₂ / ξ) * Y_1_dξ(ξ, lambda, κ, ϵ, λ) +
                (C - λI) * Y_1(ξ, lambda, κ, ϵ, λ),
            ),
        )
        # Check that the (normalized) error is small
        @test maximum(
            abs,
            exp(-CGL2.a1(κ, ϵ) * ξ^2) * (
                A * Y_1_dξ_dξ(ξ, lambda, κ, ϵ, λ) +
                (B₁ * ξ + B₂ / ξ) * Y_1_dξ(ξ, lambda, κ, ϵ, λ) +
                (C - λI) * Y_1(ξ, lambda, κ, ϵ, λ)
            ),
        ) < 1e-16
        let ξ = 2ξ
            @test all(
                Arblib.contains_zero.(
                    A * Y_1_dξ_dξ(ξ, lambda, κ, ϵ, λ) +
                    (B₁ * ξ + B₂ / ξ) * Y_1_dξ(ξ, lambda, κ, ϵ, λ) +
                    (C - λI) * Y_1(ξ, lambda, κ, ϵ, λ),
                ),
            )
            @test maximum(
                abs,
                exp(-CGL2.a1(κ, ϵ) * ξ^2) * (
                    A * Y_1_dξ_dξ(ξ, lambda, κ, ϵ, λ) +
                    (B₁ * ξ + B₂ / ξ) * Y_1_dξ(ξ, lambda, κ, ϵ, λ) +
                    (C - λI) * Y_1(ξ, lambda, κ, ϵ, λ)
                ),
            ) < 1e-22
        end

        # Check derivatives w.r.t. λ
        @test Y_1_dλ(ξ, lambda, κ, ϵ, λ) ≈ fdm(
            lambda_real ->
                Y_1(ξF64, complex(lambda_real, imag(lambdaF64)), κF64, ϵF64, λF64),
            real(lambdaF64),
        ) rtol = 1e-10

        @test Y_1_dλ_dξ(ξ, lambda, κ, ϵ, λ) ≈ fdm(
            lambda_real ->
                Y_1_dξ(ξF64, complex(lambda_real, imag(lambdaF64)), κF64, ϵF64, λF64),
            real(lambdaF64),
        ) rtol = 1e-10

        @test Y_1_dλ_dξ(ξ, lambda, κ, ϵ, λ) ≈
              fdm(ξ -> Y_1_dλ(ξ, lambdaF64, κF64, ϵF64, λF64), ξF64) rtol = 1e-10
    end

    @testset "Y2" begin
        @test all(
            Arblib.overlaps.(
                getindex.(Y_2(ArbSeries((ξ, 1)), lambda, κ, ϵ, λ), 1),
                Y_2_dξ(ξ, lambda, κ, ϵ, λ),
            ),
        )
        @test all(
            Arblib.overlaps.(
                2getindex.(Y_2(ArbSeries((ξ, 1, 0)), lambda, κ, ϵ, λ), 2),
                Y_2_dξ_dξ(ξ, lambda, κ, ϵ, λ),
            ),
        )

        @test Y_2_dξ(ξ, lambda, κ, ϵ, λ) ≈
              fdm(ξ -> Y_2(ξ, lambdaF64, κF64, ϵF64, λF64), ξF64) rtol = 1e-8

        @test Y_2_dξ_dξ(ξ, lambda, κ, ϵ, λ) ≈
              fdm2(ξ -> Y_2(ξ, lambdaF64, κF64, ϵF64, λF64), ξF64) rtol = 1e-7

        # Test that is solves equation
        @test all(
            Arblib.contains_zero.(
                A * Y_2_dξ_dξ(ξ, lambda, κ, ϵ, λ) +
                (B₁ * ξ + B₂ / ξ) * Y_2_dξ(ξ, lambda, κ, ϵ, λ) +
                (C - λI) * Y_2(ξ, lambda, κ, ϵ, λ),
            ),
        )
        # Check that the (normalized) error is small
        @test maximum(
            abs,
            exp(-CGL2.a2(κ, ϵ) * ξ^2) * (
                A * Y_2_dξ_dξ(ξ, lambda, κ, ϵ, λ) +
                (B₁ * ξ + B₂ / ξ) * Y_2_dξ(ξ, lambda, κ, ϵ, λ) +
                (C - λI) * Y_2(ξ, lambda, κ, ϵ, λ)
            ),
        ) < 1e-16
        let ξ = 2ξ
            @test all(
                Arblib.contains_zero.(
                    A * Y_2_dξ_dξ(ξ, lambda, κ, ϵ, λ) +
                    (B₁ * ξ + B₂ / ξ) * Y_2_dξ(ξ, lambda, κ, ϵ, λ) +
                    (C - λI) * Y_2(ξ, lambda, κ, ϵ, λ),
                ),
            )
            @test maximum(
                abs,
                exp(-CGL2.a2(κ, ϵ) * ξ^2) * (
                    A * Y_2_dξ_dξ(ξ, lambda, κ, ϵ, λ) +
                    (B₁ * ξ + B₂ / ξ) * Y_2_dξ(ξ, lambda, κ, ϵ, λ) +
                    (C - λI) * Y_2(ξ, lambda, κ, ϵ, λ)
                ),
            ) < 1e-22
        end

        # Check derivatives w.r.t. λ
        @test Y_2_dλ(ξ, lambda, κ, ϵ, λ) ≈ fdm(
            lambda_real ->
                Y_2(ξF64, complex(lambda_real, imag(lambdaF64)), κF64, ϵF64, λF64),
            real(lambdaF64),
        ) rtol = 1e-10

        @test Y_2_dλ_dξ(ξ, lambda, κ, ϵ, λ) ≈ fdm(
            lambda_real ->
                Y_2_dξ(ξF64, complex(lambda_real, imag(lambdaF64)), κF64, ϵF64, λF64),
            real(lambdaF64),
        ) rtol = 1e-10

        @test Y_2_dλ_dξ(ξ, lambda, κ, ϵ, λ) ≈
              fdm(ξ -> Y_2_dλ(ξ, lambdaF64, κF64, ϵF64, λF64), ξF64) rtol = 1e-10
    end

    @testset "Y3" begin
        @test all(
            Arblib.overlaps.(
                getindex.(Y_3(ArbSeries((ξ, 1)), lambda, κ, ϵ, λ), 1),
                Y_3_dξ(ξ, lambda, κ, ϵ, λ),
            ),
        )
        @test all(
            Arblib.overlaps.(
                2getindex.(Y_3(ArbSeries((ξ, 1, 0)), lambda, κ, ϵ, λ), 2),
                Y_3_dξ_dξ(ξ, lambda, κ, ϵ, λ),
            ),
        )

        @test Y_3_dξ(ξ, lambda, κ, ϵ, λ) ≈
              fdm(ξ -> Y_3(ξ, lambdaF64, κF64, ϵF64, λF64), ξF64) rtol = 1e-11

        @test Y_3_dξ_dξ(ξ, lambda, κ, ϵ, λ) ≈
              fdm2(ξ -> Y_3(ξ, lambdaF64, κF64, ϵF64, λF64), ξF64) rtol = 1e-8

        # Test that is solves equation
        @test all(
            Arblib.contains_zero.(
                A * Y_3_dξ_dξ(ξ, lambda, κ, ϵ, λ) +
                (B₁ * ξ + B₂ / ξ) * Y_3_dξ(ξ, lambda, κ, ϵ, λ) +
                (C - λI) * Y_3(ξ, lambda, κ, ϵ, λ),
            ),
        )
        # Check that the (normalized) error is small
        @test maximum(
            abs,
            A * Y_3_dξ_dξ(ξ, lambda, κ, ϵ, λ) +
            (B₁ * ξ + B₂ / ξ) * Y_3_dξ(ξ, lambda, κ, ϵ, λ) +
            (C - λI) * Y_3(ξ, lambda, κ, ϵ, λ),
        ) < 1e-16
        let ξ = 2ξ
            @test all(
                Arblib.contains_zero.(
                    A * Y_3_dξ_dξ(ξ, lambda, κ, ϵ, λ) +
                    (B₁ * ξ + B₂ / ξ) * Y_3_dξ(ξ, lambda, κ, ϵ, λ) +
                    (C - λI) * Y_3(ξ, lambda, κ, ϵ, λ),
                ),
            )
            @test maximum(
                abs,
                A * Y_3_dξ_dξ(ξ, lambda, κ, ϵ, λ) +
                (B₁ * ξ + B₂ / ξ) * Y_3_dξ(ξ, lambda, κ, ϵ, λ) +
                (C - λI) * Y_3(ξ, lambda, κ, ϵ, λ),
            ) < 1e-22
        end

        # Check derivatives w.r.t. λ
        @test Y_3_dλ(ξ, lambda, κ, ϵ, λ) ≈ fdm(
            lambda_real ->
                Y_3(ξF64, complex(lambda_real, imag(lambdaF64)), κF64, ϵF64, λF64),
            real(lambdaF64),
        ) rtol = 1e-10

        @test Y_3_dλ_dξ(ξ, lambda, κ, ϵ, λ) ≈ fdm(
            lambda_real ->
                Y_3_dξ(ξF64, complex(lambda_real, imag(lambdaF64)), κF64, ϵF64, λF64),
            real(lambdaF64),
        ) rtol = 1e-10

        @test Y_3_dλ_dξ(ξ, lambda, κ, ϵ, λ) ≈
              fdm(ξ -> Y_3_dλ(ξ, lambdaF64, κF64, ϵF64, λF64), ξF64) rtol = 1e-10
    end

    @testset "Y4" begin
        @test all(
            Arblib.overlaps.(
                getindex.(Y_4(ArbSeries((ξ, 1)), lambda, κ, ϵ, λ), 1),
                Y_4_dξ(ξ, lambda, κ, ϵ, λ),
            ),
        )
        @test all(
            Arblib.overlaps.(
                2getindex.(Y_4(ArbSeries((ξ, 1, 0)), lambda, κ, ϵ, λ), 2),
                Y_4_dξ_dξ(ξ, lambda, κ, ϵ, λ),
            ),
        )

        @test Y_4_dξ(ξ, lambda, κ, ϵ, λ) ≈
              fdm(ξ -> Y_4(ξ, lambdaF64, κF64, ϵF64, λF64), ξF64) rtol = 1e-11

        @test Y_4_dξ_dξ(ξ, lambda, κ, ϵ, λ) ≈
              fdm2(ξ -> Y_4(ξ, lambdaF64, κF64, ϵF64, λF64), ξF64) rtol = 1e-8

        # Test that is solves equation
        @test all(
            Arblib.contains_zero.(
                A * Y_4_dξ_dξ(ξ, lambda, κ, ϵ, λ) +
                (B₁ * ξ + B₂ / ξ) * Y_4_dξ(ξ, lambda, κ, ϵ, λ) +
                (C - λI) * Y_4(ξ, lambda, κ, ϵ, λ),
            ),
        )
        # Check that the (normalized) error is small
        @test maximum(
            abs,
            A * Y_4_dξ_dξ(ξ, lambda, κ, ϵ, λ) +
            (B₁ * ξ + B₂ / ξ) * Y_4_dξ(ξ, lambda, κ, ϵ, λ) +
            (C - λI) * Y_4(ξ, lambda, κ, ϵ, λ),
        ) < 1e-16
        let ξ = 2ξ
            @test all(
                Arblib.contains_zero.(
                    A * Y_4_dξ_dξ(ξ, lambda, κ, ϵ, λ) +
                    (B₁ * ξ + B₂ / ξ) * Y_4_dξ(ξ, lambda, κ, ϵ, λ) +
                    (C - λI) * Y_4(ξ, lambda, κ, ϵ, λ),
                ),
            )
            @test maximum(
                abs,
                A * Y_4_dξ_dξ(ξ, lambda, κ, ϵ, λ) +
                (B₁ * ξ + B₂ / ξ) * Y_4_dξ(ξ, lambda, κ, ϵ, λ) +
                (C - λI) * Y_4(ξ, lambda, κ, ϵ, λ),
            ) < 1e-22
        end

        # Check derivatives w.r.t. λ
        @test Y_4_dλ(ξ, lambda, κ, ϵ, λ) ≈ fdm(
            lambda_real ->
                Y_4(ξF64, complex(lambda_real, imag(lambdaF64)), κF64, ϵF64, λF64),
            real(lambdaF64),
        ) rtol = 1e-10

        @test Y_4_dλ_dξ(ξ, lambda, κ, ϵ, λ) ≈ fdm(
            lambda_real ->
                Y_4_dξ(ξF64, complex(lambda_real, imag(lambdaF64)), κF64, ϵF64, λF64),
            real(lambdaF64),
        ) rtol = 1e-10

        @test Y_4_dλ_dξ(ξ, lambda, κ, ϵ, λ) ≈
              fdm(ξ -> Y_4_dλ(ξ, lambdaF64, κF64, ϵF64, λF64), ξF64) rtol = 1e-10
    end

    @testset "K1 and K2" begin
        # K1 and K2 are suppose to give solutions to the linear system
        # Ψ * v = [[0, 0]; A \ F].

        Y12 = hcat(Y_1(ξ, lambda, κ, ϵ, λ), Y_2(ξ, lambda, κ, ϵ, λ))
        Y12_dξ = hcat(Y_1_dξ(ξ, lambda, κ, ϵ, λ), Y_2_dξ(ξ, lambda, κ, ϵ, λ))
        Y34 = hcat(Y_3(ξ, lambda, κ, ϵ, λ), Y_4(ξ, lambda, κ, ϵ, λ))
        Y34_dξ = hcat(Y_3_dξ(ξ, lambda, κ, ϵ, λ), Y_4_dξ(ξ, lambda, κ, ϵ, λ))

        Ψ = [Y12 Y34; Y12_dξ Y34_dξ]
        F = eltype(Ψ)[0.1, 0.25]

        K1, K2 = CGL2.K_1_2(Y12, Y12_dξ, Y34, Y34_dξ, A)

        v12 = -K1 * F
        v34 = -K2 * F
        v = [v12; v34]

        @test all(Arblib.overlaps.(Ψ * v, [[0, 0]; A \ F]))

        # Test derivatives w.r.t. λ
        Y12_dλ = hcat(Y_1_dλ(ξ, lambda, κ, ϵ, λ), Y_2_dλ(ξ, lambda, κ, ϵ, λ))
        Y12_dλ_dξ = hcat(Y_1_dλ_dξ(ξ, lambda, κ, ϵ, λ), Y_2_dλ_dξ(ξ, lambda, κ, ϵ, λ))
        Y34_dλ = hcat(Y_3_dλ(ξ, lambda, κ, ϵ, λ), Y_4_dλ(ξ, lambda, κ, ϵ, λ))
        Y34_dλ_dξ = hcat(Y_3_dλ_dξ(ξ, lambda, κ, ϵ, λ), Y_4_dλ_dξ(ξ, lambda, κ, ϵ, λ))

        K1_dλ, K2_dλ =
            CGL2.K_1_2_dλ(Y12, Y12_dξ, Y34, Y34_dξ, Y12_dλ, Y12_dλ_dξ, Y34_dλ, Y34_dλ_dξ, A)

        K1_dλ_fdm = fdm(
            lambda_real -> CGL2.K_1_2(
                ξF64,
                complex(lambda_real, imag(lambdaF64)),
                κF64,
                ϵF64,
                λF64,
            )[1],
            real(lambdaF64),
        )

        K2_dλ_fdm = fdm(
            lambda_real -> CGL2.K_1_2(
                ξF64,
                complex(lambda_real, imag(lambdaF64)),
                κF64,
                ϵF64,
                λF64,
            )[2],
            real(lambdaF64),
        )

        @test K1_dλ ≈ K1_dλ_fdm rtol = 1e-12
        @test K2_dλ ≈ K2_dλ_fdm rtol = 1e-12
    end

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

    @testset "J_P_$j" for (j, J_P_j, J_P_j_dλ, P_j, P_j_dλ, W_j) in [
        (1, CGL2.J_P_1, CGL2.J_P_1_dλ, CGL2.P_1, CGL2.P_1_dλ, CGL2.W_1),
        (2, CGL2.J_P_2, CGL2.J_P_2_dλ, CGL2.P_2, CGL2.P_2_dλ, CGL2.W_2),
    ]
        @test Arblib.overlaps(
            J_P_j(ξ, lambda, κ, ϵ, λ),
            P_j(ξ, lambda, κ, ϵ, λ) / W_j(ξ, lambda, κ, ϵ, λ),
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

    @testset "J_E_$j" for (j, J_E_j, J_E_j_dλ, E_j, E_j_dλ, W_j) in [
        (1, CGL2.J_E_1, CGL2.J_E_1_dλ, CGL2.E_1, CGL2.E_1_dλ, CGL2.W_1),
        (2, CGL2.J_E_2, CGL2.J_E_2_dλ, CGL2.E_2, CGL2.E_2_dλ, CGL2.W_2),
    ]
        @test Arblib.overlaps(
            J_E_j(ξ, lambda, κ, ϵ, λ),
            E_j(ξ, lambda, κ, ϵ, λ) / W_j(ξ, lambda, κ, ϵ, λ),
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

    @testset "K1_new and K2_new" begin
        # K1 and K2 are suppose to give solutions to the linear system
        # Ψ * v = [[0, 0]; A \ F].

        M = SMatrix{2,2}(im, 1, -im, 1)
        Ψ = [M * F_Y.E_12 M * F_Y.P_12; M * F_Y.E_12_dξ M * F_Y.P_12_dξ]
        F = eltype(Ψ)[0.1, 0.25]

        v = [-F_Y.K_1 * F; -F_Y.K_2 * F]

        @test all(Arblib.overlaps.(Ψ * v, [[0, 0]; A \ F]))

        K1_dλ, K2_dλ = CGL2.K_1_2_dλ_new(ξ, lambda, κ, ϵ, λ)

        K1_dλ_series = getindex.(CGL2.K_1_2_new(ξ, AcbSeries((lambda, 1)), κ, ϵ, λ)[1], 1)
        K2_dλ_series = getindex.(CGL2.K_1_2_new(ξ, AcbSeries((lambda, 1)), κ, ϵ, λ)[2], 1)

        K1_dλ_fdm = fdm(
            lambda_real -> CGL2.K_1_2_new(
                ξF64,
                complex(lambda_real, imag(lambdaF64)),
                κF64,
                ϵF64,
                λF64,
            )[1],
            real(lambdaF64),
        )

        K2_dλ_fdm = fdm(
            lambda_real -> CGL2.K_1_2_new(
                ξF64,
                complex(lambda_real, imag(lambdaF64)),
                κF64,
                ϵF64,
                λF64,
            )[2],
            real(lambdaF64),
        )

        @test all(Arblib.overlaps.(F_Y.K_1_dλ, K1_dλ_series))
        @test all(Arblib.overlaps.(F_Y.K_2_dλ, K2_dλ_series))
        @test F_Y.K_1_dλ ≈ K1_dλ_fdm rtol = 1e-12
        @test F_Y.K_2_dλ ≈ K2_dλ_fdm rtol = 1e-12
    end

    @testset "JN" begin
        # The precise value for a and b should not play any role in
        # the correctness, we just compute some approximation here.
        νF64 = 1.9261384880241954 + 3.0638598354170337im
        a, b =
            Arb.(CGL2.Q_hat_zero_float(real(νF64), imag(νF64), κF64, ϵF64, ξF64, λF64)[1:2])

        # Compute Jacobian by going through ArbSeries
        N = (a, b) -> (a^2 + b^2)^λ.σ * SVector(-λ.δ * a - b, a - λ.δ * b)
        N_a = getindex.(N(ArbSeries((a, 1)), b), 1)
        N_b = getindex.(N(a, ArbSeries((b, 1))), 1)
        JN_direct = [N_a N_b]

        @test all(Arblib.overlaps.(JN_direct, CGL2.J_N(Acb(a, b), λ)))
    end
end
