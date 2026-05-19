@testset "parameter" begin
    ξ = Arb(30)
    κ = Arb(0.493223)
    ϵ = Arb(0.1)
    Λ = CGLParams{Arb}(1, 1.0, 2.3, 0.2)

    ξF64 = Float64(ξ)
    κF64 = Float64(κ)
    ϵF64 = Float64(ϵ)
    ΛF64 = CGLParams{Float64}(Λ)

    # Function for computing derivative using finite differences.
    fdm = central_fdm(5, 1)

    @testset "a" begin
        @test CGL2._a(κ, ϵ, Λ) ≈ CGL2._a(κF64, ϵF64, ΛF64)

        @test CGL2._a_dκ(κ, ϵ, Λ) ≈ CGL2._a_dκ(κF64, ϵF64, ΛF64)
        @test Arblib.overlaps(CGL2._a_dκ(κ, ϵ, Λ), CGL2._a(ArbSeries((κ, 1)), ϵ, Λ)[1])
        @test CGL2._a_dκ(κF64, ϵF64, ΛF64) ≈ fdm(κ -> CGL2._a(κ, ϵF64, ΛF64), κF64)
    end

    @testset "b" begin
        @test CGL2._b(κ, ϵ, Λ) ≈ CGL2._b(κF64, ϵF64, ΛF64)
    end

    @testset "c" begin
        @test CGL2._c(κ, ϵ, Λ) ≈ CGL2._c(κF64, ϵF64, ΛF64)

        @test CGL2._c_dκ(κ, ϵ, Λ) ≈ CGL2._c_dκ(κF64, ϵF64, ΛF64)
        @test Arblib.overlaps(CGL2._c_dκ(κ, ϵ, Λ), CGL2._c(ArbSeries((κ, 1)), ϵ, Λ)[1])
        @test CGL2._c_dκ(κF64, ϵF64, ΛF64) ≈ fdm(κ -> CGL2._c(κ, ϵF64, ΛF64), κF64)

        @test CGL2._c_dϵ(κ, ϵ, Λ) ≈ CGL2._c_dϵ(κF64, ϵF64, ΛF64)
        @test Arblib.overlaps(CGL2._c_dϵ(κ, ϵ, Λ), CGL2._c(κ, ArbSeries((ϵ, 1)), Λ)[1])
        @test CGL2._c_dϵ(κF64, ϵF64, ΛF64) ≈ fdm(ϵ -> CGL2._c(κF64, ϵ, ΛF64), ϵF64)
    end

    @testset "B_W" begin
        @test CGL2.B_W(κ, ϵ, Λ) ≈ CGL2.B_W(κF64, ϵF64, ΛF64)

        @test CGL2.B_W_dκ(κ, ϵ, Λ) ≈ CGL2.B_W_dκ(κF64, ϵF64, ΛF64)
        @test Arblib.overlaps(CGL2.B_W_dκ(κ, ϵ, Λ), CGL2.B_W(ArbSeries((κ, 1)), ϵ, Λ)[1])
        @test CGL2.B_W_dκ(κF64, ϵF64, ΛF64) ≈ fdm(κ -> CGL2.B_W(κ, ϵF64, ΛF64), κF64)

        @test CGL2.B_W_dϵ(κ, ϵ, Λ) ≈ CGL2.B_W_dϵ(κF64, ϵF64, ΛF64)
        @test Arblib.overlaps(CGL2.B_W_dϵ(κ, ϵ, Λ), CGL2.B_W(κ, ArbSeries((ϵ, 1)), Λ)[1])
        @test CGL2.B_W_dϵ(κF64, ϵF64, ΛF64) ≈ fdm(ϵ -> CGL2.B_W(κF64, ϵ, ΛF64), ϵF64)
    end
end
