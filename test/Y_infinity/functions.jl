@testset "functions" begin
    # IMPROVE: Test for more parameters
    ξ = Arb(30)
    κ = Arb(0.8073018593981386)
    ϵ = Arb(0.15002213424487343)
    lambda = Acb(0.2 + 1.3im)
    λ = CGLParams{Arb}(3, 1.0, 1.0, 0.0)

    ξF64 = Float64(ξ)
    κF64 = Float64(κ)
    ϵF64 = Float64(ϵ)
    lambdaF64 = ComplexF64(lambda)
    λF64 = CGLParams{Float64}(λ)

    A = SMatrix{2,2}(ϵ, 1, -1, ϵ)
    B₁ = SMatrix{2,2}(κ, 0, 0, κ)
    B₂ = (λ.d - 1) * A
    C = SMatrix{2,2}(κ / λ.σ, λ.ω, -λ.ω, κ / λ.σ)
    λI = SMatrix{2,2}(lambda, 0, 0, lambda)

    # Function for computing derivative using finite differences.
    fdm = central_fdm(5, 1)
    fdm2 = central_fdm(5, 2)
    fdm3 = central_fdm(5, 3)

    Y_1, Y_1_dξ, Y_1_dξ_dξ = CGL2.Y_1, CGL2.Y_1_dξ, CGL2.Y_1_dξ_dξ
    Y_2, Y_2_dξ, Y_2_dξ_dξ = CGL2.Y_2, CGL2.Y_2_dξ, CGL2.Y_2_dξ_dξ
    Y_3, Y_3_dξ, Y_3_dξ_dξ = CGL2.Y_3, CGL2.Y_3_dξ, CGL2.Y_3_dξ_dξ
    Y_4, Y_4_dξ, Y_4_dξ_dξ = CGL2.Y_4, CGL2.Y_4_dξ, CGL2.Y_4_dξ_dξ

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

        @test real(Y_1_dξ(ξ, lambda, κ, ϵ, λ)) ≈
              fdm(ξ -> real(Y_1(ξ, lambdaF64, κF64, ϵF64, λF64)), ξF64) rtol = 1e-8

        @test imag(Y_1_dξ(ξ, lambda, κ, ϵ, λ)) ≈
              fdm(ξ -> imag(Y_1(ξ, lambdaF64, κF64, ϵF64, λF64)), ξF64) rtol = 1e-8

        @test real(Y_1_dξ_dξ(ξ, lambda, κ, ϵ, λ)) ≈
              fdm2(ξ -> real(Y_1(ξ, lambdaF64, κF64, ϵF64, λF64)), ξF64) rtol = 1e-7

        @test imag(Y_1_dξ_dξ(ξ, lambda, κ, ϵ, λ)) ≈
              fdm2(ξ -> imag(Y_1(ξ, lambdaF64, κF64, ϵF64, λF64)), ξF64) rtol = 1e-7

        # Test that is solves equation, at least approximately

        @test exp(-CGL2.real_a2(κ, ϵ) * ξ^2) * norm(
            A * Y_1_dξ_dξ(ξ, lambda, κ, ϵ, λ) +
            (B₁ * ξ + B₂ / ξ) * Y_1_dξ(ξ, lambda, κ, ϵ, λ) +
            (C - λI) * Y_1(ξ, lambda, κ, ϵ, λ),
        ) < 1e-13
        let ξ = 2ξ
            @test exp(-CGL2.real_a2(κ, ϵ) * ξ^2) * norm(
                A * Y_1_dξ_dξ(ξ, lambda, κ, ϵ, λ) +
                (B₁ * ξ + B₂ / ξ) * Y_1_dξ(ξ, lambda, κ, ϵ, λ) +
                (C - λI) * Y_1(ξ, lambda, κ, ϵ, λ),
            ) < 1e-15
        end
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

        @test real(Y_2_dξ(ξ, lambda, κ, ϵ, λ)) ≈
              fdm(ξ -> real(Y_2(ξ, lambdaF64, κF64, ϵF64, λF64)), ξF64) rtol = 1e-8

        @test imag(Y_2_dξ(ξ, lambda, κ, ϵ, λ)) ≈
              fdm(ξ -> imag(Y_2(ξ, lambdaF64, κF64, ϵF64, λF64)), ξF64) rtol = 1e-8

        @test real(Y_2_dξ_dξ(ξ, lambda, κ, ϵ, λ)) ≈
              fdm2(ξ -> real(Y_2(ξ, lambdaF64, κF64, ϵF64, λF64)), ξF64) rtol = 1e-7

        @test imag(Y_2_dξ_dξ(ξ, lambda, κ, ϵ, λ)) ≈
              fdm2(ξ -> imag(Y_2(ξ, lambdaF64, κF64, ϵF64, λF64)), ξF64) rtol = 1e-7

        @test exp(-CGL2.real_a2(κ, ϵ) * ξ^2) * norm(
            A * Y_2_dξ_dξ(ξ, lambda, κ, ϵ, λ) +
            (B₁ * ξ + B₂ / ξ) * Y_2_dξ(ξ, lambda, κ, ϵ, λ) +
            (C - λI) * Y_2(ξ, lambda, κ, ϵ, λ),
        ) < 1e-14
        let ξ = 2ξ
            @test exp(-CGL2.real_a2(κ, ϵ) * ξ^2) * norm(
                A * Y_2_dξ_dξ(ξ, lambda, κ, ϵ, λ) +
                (B₁ * ξ + B₂ / ξ) * Y_2_dξ(ξ, lambda, κ, ϵ, λ) +
                (C - λI) * Y_2(ξ, lambda, κ, ϵ, λ),
            ) < 1e-19
        end
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

        @test real(Y_3_dξ(ξ, lambda, κ, ϵ, λ)) ≈
              fdm(ξ -> real(Y_3(ξ, lambdaF64, κF64, ϵF64, λF64)), ξF64) rtol = 1e-11

        @test imag(Y_3_dξ(ξ, lambda, κ, ϵ, λ)) ≈
              fdm(ξ -> imag(Y_3(ξ, lambdaF64, κF64, ϵF64, λF64)), ξF64) rtol = 1e-11

        @test real(Y_3_dξ_dξ(ξ, lambda, κ, ϵ, λ)) ≈
              fdm2(ξ -> real(Y_3(ξ, lambdaF64, κF64, ϵF64, λF64)), ξF64) rtol = 1e-8

        @test imag(Y_3_dξ_dξ(ξ, lambda, κ, ϵ, λ)) ≈
              fdm2(ξ -> imag(Y_3(ξ, lambdaF64, κF64, ϵF64, λF64)), ξF64) rtol = 1e-8

        # Test that is solves equation, at least approximately
        @test norm(
            A * Y_3_dξ_dξ(ξ, lambda, κ, ϵ, λ) +
            (B₁ * ξ + B₂ / ξ) * Y_3_dξ(ξ, lambda, κ, ϵ, λ) +
            (C - λI) * Y_3(ξ, lambda, κ, ϵ, λ),
        ) < 1e-12
        let ξ = 2ξ
            @test norm(
                A * Y_3_dξ_dξ(ξ, lambda, κ, ϵ, λ) +
                (B₁ * ξ + B₂ / ξ) * Y_3_dξ(ξ, lambda, κ, ϵ, λ) +
                (C - λI) * Y_3(ξ, lambda, κ, ϵ, λ),
            ) < 1e-15
        end
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

        @test real(Y_4_dξ(ξ, lambda, κ, ϵ, λ)) ≈
              fdm(ξ -> real(Y_4(ξ, lambdaF64, κF64, ϵF64, λF64)), ξF64) rtol = 1e-11

        @test imag(Y_4_dξ(ξ, lambda, κ, ϵ, λ)) ≈
              fdm(ξ -> imag(Y_4(ξ, lambdaF64, κF64, ϵF64, λF64)), ξF64) rtol = 1e-11

        @test real(Y_4_dξ_dξ(ξ, lambda, κ, ϵ, λ)) ≈
              fdm2(ξ -> real(Y_4(ξ, lambdaF64, κF64, ϵF64, λF64)), ξF64) rtol = 1e-8

        @test imag(Y_4_dξ_dξ(ξ, lambda, κ, ϵ, λ)) ≈
              fdm2(ξ -> imag(Y_4(ξ, lambdaF64, κF64, ϵF64, λF64)), ξF64) rtol = 1e-8

        # Test that is solves equation, at least approximately
        @test norm(
            A * Y_4_dξ_dξ(ξ, lambda, κ, ϵ, λ) +
            (B₁ * ξ + B₂ / ξ) * Y_4_dξ(ξ, lambda, κ, ϵ, λ) +
            (C - λI) * Y_4(ξ, lambda, κ, ϵ, λ),
        ) < 1e-12
        let ξ = 2ξ
            @test norm(
                A * Y_4_dξ_dξ(ξ, lambda, κ, ϵ, λ) +
                (B₁ * ξ + B₂ / ξ) * Y_4_dξ(ξ, lambda, κ, ϵ, λ) +
                (C - λI) * Y_4(ξ, lambda, κ, ϵ, λ),
            ) < 1e-15
        end
    end

    @testset "K1 and K2" begin
        # K1 and K2 are suppose to give solutions to the linear system
        # Ψ * v = [[0, 0]; A \ F].

        Y12 = hcat(Y_1(ξ₁, lambda, κ, ϵ, λ), Y_2(ξ₁, lambda, κ, ϵ, λ))
        Y12_dξ = hcat(Y_1_dξ(ξ₁, lambda, κ, ϵ, λ), Y_2_dξ(ξ₁, lambda, κ, ϵ, λ))
        Y34 = hcat(Y_3(ξ₁, lambda, κ, ϵ, λ), Y_4(ξ₁, lambda, κ, ϵ, λ))
        Y34_dξ = hcat(Y_3_dξ(ξ₁, lambda, κ, ϵ, λ), Y_4_dξ(ξ₁, lambda, κ, ϵ, λ))

        Ψ = [Y12 Y34; Y12_dξ Y34_dξ]
        F = eltype(Ψ)[0.1, 0.25]

        K1, K2 = CGL2.K_1_2(Y12, Y12_dξ, Y34, Y34_dξ, A)

        v12 = -K1 * F
        v34 = -K2 * F
        v = [v12; v34]

        @test all(Arblib.overlaps.(Ψ * v, [[0, 0]; A \ F]))
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
