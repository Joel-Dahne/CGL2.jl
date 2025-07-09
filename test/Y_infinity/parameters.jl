@testset "Y_infinity/parameters" begin
    # IMPROVE: Test for more parameters
    κ = Arb(0.8073018593981386)
    ϵ = Arb(0.15002213424487343)
    lambda = Acb(0.2 + 1.3im)
    λ = CGLParams{Arb}(3, 1.0, 1.0, 0.0)

    (; A, B₁, B₂, C, λI) = CGL2.coeff_matrices(lambda, κ, ϵ, λ)

    @testset "a₂ != 0" begin
        # Check that s and A0 are eigenvalues of the appropriate
        # matrix
        @test all(
            Arblib.overlaps.(
                -0.5(A \ B₁) * CGL2.A1(typeof(κ)),
                CGL2.a1(κ, ϵ) * CGL2.A1(typeof(κ)),
            ),
        )
        @test all(
            Arblib.overlaps.(
                -0.5(A \ B₁) * CGL2.A2(typeof(κ)),
                CGL2.a2(κ, ϵ) * CGL2.A2(typeof(κ)),
            ),
        )

        # Check that s and A2 satisfy the required condition
        #let a₂ = CGL2.a21(κ, ϵ), s = CGL2.s11(lambda, κ, λ)
        #    @test all(
        #        Arblib.overlaps.(
        #            (4a₂ * A + 2B₁) * CGL2.A211(typeof(κ)),
        #            -inv(a₂) *
        #            (2a₂ * A + 2a₂ * B₂ + C - λI - s * (4a₂ * A + B₁)) *
        #            CGL2.A011(typeof(κ)),
        #        ),
        #    )
        #end
        #let a₂ = CGL2.a22(κ, ϵ), s = CGL2.s12(lambda, κ, λ)
        #    @test all(
        #        Arblib.overlaps.(
        #            (4a₂ * A + 2B₁) * CGL2.A212(typeof(κ)),
        #            -inv(a₂) *
        #            (2a₂ * A + 2a₂ * B₂ + C - λI - s * (4a₂ * A + B₁)) *
        #            CGL2.A012(typeof(κ)),
        #        ),
        #    )
        #end

        # Check that An satisfy the expected recurrence relation
        N = 10

        let a₂ = CGL2.a21(κ, ϵ), s = CGL2.s11(lambda, κ, λ)
            #As = CGL2.As_11(N, lambda, κ, ϵ, λ)
            cs, As = CGL2.As_11(N, lambda, κ, ϵ, λ)

            @test all(iszero, As[2:2:end]) # Every other is zero
            for n = 2:2:N
                @test all(
                    Arblib.contains_zero.(
                        (2a₂ * (A + B₂) + C - λI - (n + s) * (4a₂ * A + B₁)) * As[n+1] +
                        ((n + s - 2) * (n + s + 1 - 2) * A - B₂ * (n + s - 2)) *
                        As[(n-2)+1],
                    ),
                )

                @test all(
                    Arblib.contains_zero.(
                        n * κ * As[n+1] +
                        ((n + s - 2)^2 - (d - 2) * (n + s - 2)) * A * As[(n-2)+1],
                    ),
                )

                @test all(
                    Arblib.contains_zero.(
                        (
                            cs[n+1] * n * κ * I +
                            cs[(n-2)+1] * ((n + s - 2)^2 - (d - 2) * (n + s - 2)) * A
                        ) * As[1],
                    ),
                )
            end
        end

        let s = CGL2.s12(lambda, κ, λ), As = CGL2.As_12(N, lambda, κ, ϵ, λ)
            #@test all(iszero, As[2:2:end]) # Every other is zero
            for n = 2:2:N
                #@test all(
                #    Arblib.overlaps.(
                #        (C - λI - (n + s) * B₁) * As[n+1],
                #        -((n + s - 2) * (n + s - 1) * A - B₂ * (s + n - 2)) * As[(n-2)+1],
                #    ),
                #)
            end
        end

        @test isequal(CGL2.a21(κ, ϵ), conj(CGL2.a22(κ, ϵ)))
        @test isequal(CGL2.real_a2(κ, ϵ), real(CGL2.a21(κ, ϵ)))
        @test isequal(CGL2.real_a2(κ, ϵ), real(CGL2.a22(κ, ϵ)))

        #@test CGL2.s11(lambda, κ, ϵ, λ) isa Complex # TODO: Test something
        #@test CGL2.s12(lambda, κ, ϵ, λ) isa Complex # TODO: Test something
    end

    @testset "a₂ == 0" begin
        # Check that s and A0 are eigenvalues of the appropriate
        # matrix
        @test all(
            Arblib.overlaps.(
                B₁ \ (C - λI) * CGL2.A021(typeof(κ)),
                CGL2.s21(lambda, κ, λ) * CGL2.A021(typeof(κ)),
            ),
        )
        @test all(
            Arblib.overlaps.(
                B₁ \ (C - λI) * CGL2.A022(typeof(κ)),
                CGL2.s22(lambda, κ, λ) * CGL2.A022(typeof(κ)),
            ),
        )

        # Check that they have the expected correspondence
        @test Arblib.overlaps(
            CGL2.s21(lambda, κ, λ),
            CGL2.s22(lambda, κ, λ) + 2λ.ω / κ * im,
        )
        @test isequal(CGL2.real_s2(lambda, κ, λ), real(CGL2.s21(lambda, κ, λ)))
        @test isequal(CGL2.real_s2(lambda, κ, λ), real(CGL2.s22(lambda, κ, λ)))

        # Check that An satisfy the expected recurrence relation
        N = 10

        let s = CGL2.s21(lambda, κ, λ), As = CGL2.As_21(N, lambda, κ, ϵ, λ)
            @test all(iszero, As[2:2:end]) # Every other is zero
            for n = 2:2:N
                @test all(
                    Arblib.overlaps.(
                        (C - λI - (n + s) * B₁) * As[n+1],
                        -((n + s - 2) * (n + s - 1) * A - B₂ * (s + n - 2)) * As[(n-2)+1],
                    ),
                )
            end
        end

        let s = CGL2.s22(lambda, κ, λ), As = CGL2.As_22(N, lambda, κ, ϵ, λ)
            @test all(iszero, As[2:2:end]) # Every other is zero
            for n = 2:2:N
                @test all(
                    Arblib.overlaps.(
                        (C - λI - (n + s) * B₁) * As[n+1],
                        -((n + s - 2) * (n + s - 1) * A - B₂ * (s + n - 2)) * As[(n-2)+1],
                    ),
                )
            end
        end
    end
end
