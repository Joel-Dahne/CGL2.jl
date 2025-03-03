@testset "Y_infinity/parameters" begin
    # IMPROVE: Test for more parameters
    κ = 0.8073018593981386
    ϵ = 0.15002213424487343
    lambda = 0.2 + 1.3im
    λ = CGLParams(3, 1.0, 1.0, 0.0)

    @test isequal(CGL2.a21(κ, ϵ), conj(CGL2.a22(κ, ϵ)))
    @test isequal(CGL2.real_a2(κ, ϵ), real(CGL2.a21(κ, ϵ)))
    @test isequal(CGL2.real_a2(κ, ϵ), real(CGL2.a22(κ, ϵ)))

    @test CGL2.s11(lambda, κ, ϵ, λ) isa Complex # TODO: Test something
    @test CGL2.s12(lambda, κ, ϵ, λ) isa Complex # TODO: Test something

    @test isequal(CGL2.s21(lambda, κ, λ), CGL2.s22(lambda, κ, λ) + 2λ.ω / κ * im)
    @test isequal(CGL2.real_s2(lambda, κ, λ), real(CGL2.s21(lambda, κ, λ)))
    @test isequal(CGL2.real_s2(lambda, κ, λ), real(CGL2.s22(lambda, κ, λ)))
end
