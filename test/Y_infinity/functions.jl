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

    Y_3, Y_3_dξ, Y_3_dξ_dξ = CGL2.Y_3, CGL2.Y_3_dξ, CGL2.Y_3_dξ_dξ
    Y_4, Y_4_dξ, Y_4_dξ_dξ = CGL2.Y_4, CGL2.Y_4_dξ, CGL2.Y_4_dξ_dξ

    @testset "Y3" begin
        @test all(
            Arblib.overlaps.(
                getindex.(Y_3(ArbSeries((ξ, 1)), lambda, κ, ϵ, λ), 1),
                Y_3_dξ(ξ, lambda, κ, ϵ, λ),
            ),
        )

        #@test Arblib.overlaps(Y_3(ξ, ArbSeries((κ, 1)), ϵ, λ)[1], Y_3_dκ(ξ, κ, ϵ, λ))

        @test real(Y_3_dξ(ξ, lambda, κ, ϵ, λ)) ≈
              fdm(ξ -> real(Y_3(ξ, lambdaF64, κF64, ϵF64, λF64)), ξF64) rtol = 1e-11

        @test imag(Y_3_dξ(ξ, lambda, κ, ϵ, λ)) ≈
              fdm(ξ -> imag(Y_3(ξ, lambdaF64, κF64, ϵF64, λF64)), ξF64) rtol = 1e-11

        @test real(Y_3_dξ_dξ(ξ, lambda, κ, ϵ, λ)) ≈
              fdm2(ξ -> real(Y_3(ξ, lambdaF64, κF64, ϵF64, λF64)), ξF64) rtol = 1e-8

        @test imag(Y_3_dξ_dξ(ξ, lambda, κ, ϵ, λ)) ≈
              fdm2(ξ -> imag(Y_3(ξ, lambdaF64, κF64, ϵF64, λF64)), ξF64) rtol = 1e-8

        #@test Y_3_dξ_dξ_dξ(ξ, κ, ϵ, λ) ≈ fdm(ξ -> Y_3_dξ_dξ(ξ, κF64, ϵF64, λF64), ξF64) rtol =
        #    1e-12
        #@test Y_3_dξ_dξ_dξ(ξ, κ, ϵ, λ) ≈ fdm2(ξ -> Y_3_dξ(ξ, κF64, ϵF64, λF64), ξF64) rtol =
        #    1e-8
        #@test Y_3_dξ_dξ_dξ(ξ, κ, ϵ, λ) ≈ fdm3(ξ -> Y_3(ξ, κF64, ϵF64, λF64), ξF64) rtol = 1e-6

        #@test Y_3_dκ(ξ, κ, ϵ, λ) ≈ fdm(κ -> Y_3(ξF64, κ, ϵF64, λF64), κF64) rtol = 1e-12

        #@test Y_3_dξ_dκ(ξ, κ, ϵ, λ) ≈ fdm(κ -> Y_3_dξ(ξF64, κ, ϵF64, λF64), κF64) rtol = 1e-12

        #@test Y_3_dξ_dξ_dκ(ξ, κ, ϵ, λ) ≈ fdm(κ -> Y_3_dξ_dξ(ξF64, κ, ϵF64, λF64), κF64) rtol =
        #    1e-12

        #@test Y_3_dϵ(ξ, κ, ϵ, λ) ≈ fdm(ϵ -> Y_3(ξF64, κF64, ϵ, λF64), ϵF64) rtol = 1e-12

        #@test Y_3_dξ_dϵ(ξ, κ, ϵ, λ) ≈ fdm(ϵ -> Y_3_dξ(ξF64, κF64, ϵ, λF64), ϵF64) rtol = 1e-12

        #@test Y_3_dξ_dξ_dϵ(ξ, κ, ϵ, λ) ≈ fdm(ϵ -> Y_3_dξ_dξ(ξF64, κF64, ϵ, λF64), ϵF64) rtol =
        #    1e-12

        # Test that is solves equation, at least approximately
        @test norm(
            A * Y_3_dξ_dξ(ξ, lambda, κ, ϵ, λ) +
            (B₁ * ξ + B₂ / ξ) * Y_3_dξ(ξ, lambda, κ, ϵ, λ) +
            (C - λI) * Y_3(ξ, lambda, κ, ϵ, λ),
        ) < 1e-14
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

        #@test Arblib.overlaps(Y_4(ξ, ArbSeries((κ, 1)), ϵ, λ)[1], Y_4_dκ(ξ, κ, ϵ, λ))

        @test real(Y_4_dξ(ξ, lambda, κ, ϵ, λ)) ≈
              fdm(ξ -> real(Y_4(ξ, lambdaF64, κF64, ϵF64, λF64)), ξF64) rtol = 1e-12

        @test imag(Y_4_dξ(ξ, lambda, κ, ϵ, λ)) ≈
              fdm(ξ -> imag(Y_4(ξ, lambdaF64, κF64, ϵF64, λF64)), ξF64) rtol = 1e-12

        @test real(Y_4_dξ_dξ(ξ, lambda, κ, ϵ, λ)) ≈
              fdm2(ξ -> real(Y_4(ξ, lambdaF64, κF64, ϵF64, λF64)), ξF64) rtol = 1e-8

        @test imag(Y_4_dξ_dξ(ξ, lambda, κ, ϵ, λ)) ≈
              fdm2(ξ -> imag(Y_4(ξ, lambdaF64, κF64, ϵF64, λF64)), ξF64) rtol = 1e-8

        #@test Y_4_dξ_dξ_dξ(ξ, κ, ϵ, λ) ≈ fdm(ξ -> Y_4_dξ_dξ(ξ, κF64, ϵF64, λF64), ξF64) rtol =
        #    1e-12
        #@test Y_4_dξ_dξ_dξ(ξ, κ, ϵ, λ) ≈ fdm2(ξ -> Y_4_dξ(ξ, κF64, ϵF64, λF64), ξF64) rtol =
        #    1e-8
        #@test Y_4_dξ_dξ_dξ(ξ, κ, ϵ, λ) ≈ fdm4(ξ -> Y_4(ξ, κF64, ϵF64, λF64), ξF64) rtol = 1e-6

        #@test Y_4_dκ(ξ, κ, ϵ, λ) ≈ fdm(κ -> Y_4(ξF64, κ, ϵF64, λF64), κF64) rtol = 1e-12

        #@test Y_4_dξ_dκ(ξ, κ, ϵ, λ) ≈ fdm(κ -> Y_4_dξ(ξF64, κ, ϵF64, λF64), κF64) rtol = 1e-12

        #@test Y_4_dξ_dξ_dκ(ξ, κ, ϵ, λ) ≈ fdm(κ -> Y_4_dξ_dξ(ξF64, κ, ϵF64, λF64), κF64) rtol =
        #    1e-12

        #@test Y_4_dϵ(ξ, κ, ϵ, λ) ≈ fdm(ϵ -> Y_4(ξF64, κF64, ϵ, λF64), ϵF64) rtol = 1e-12

        #@test Y_4_dξ_dϵ(ξ, κ, ϵ, λ) ≈ fdm(ϵ -> Y_4_dξ(ξF64, κF64, ϵ, λF64), ϵF64) rtol = 1e-12

        #@test Y_4_dξ_dξ_dϵ(ξ, κ, ϵ, λ) ≈ fdm(ϵ -> Y_4_dξ_dξ(ξF64, κF64, ϵ, λF64), ϵF64) rtol =
        #    1e-12

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

end
