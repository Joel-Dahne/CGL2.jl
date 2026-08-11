@testset "Y_zero" begin
    Y₀ = SVector(Acb(0.1, 0.2), Acb(0.3, 0.4))
    λ = Acb(0.19028080950252219, 2.769747313863597)
    ν = Acb(2.1134119144964454, 2.5294091893480206)
    κ = Arb(0.8073018593981386)
    ϵ = Arb(0.15002213424487343)
    ξ₁ = Arb(30)
    Λ = CGLParams{Arb}(3, 1.0, 1.0, 0.0)

    Y₀F64 = ComplexF64.(Y₀)
    λF64 = ComplexF64(λ)
    νF64 = ComplexF64(ν)
    κF64 = Float64(κ)
    ϵF64 = Float64(ϵ)
    ξ₁F64 = Float64(ξ₁)
    ΛF64 = CGLParams{Float64}(Λ)

    res = CGL2.Y_zero(Y₀, λ, ν, κ, ϵ, ξ₁, Λ)

    # Test with different settings for Taylor expansion
    res2 = CGL2.Y_zero_capd(Y₀, λ, ν, κ, ϵ, ξ₁, Λ, ξ₀ = Arb(1e-1), degree = 6)
    res3 = CGL2.Y_zero_capd(Y₀, λ, ν, κ, ϵ, ξ₁, Λ, ξ₀ = Arb(1e-3))
    @test all(Arblib.overlaps.(res, res2))
    @test all(Arblib.overlaps.(res, res3))

    # Compare with Float64 version
    resF64 = CGL2.Y_zero(Y₀F64, λF64, νF64, κF64, ϵF64, ξ₁F64, ΛF64)
    @test ComplexF64.(res) ≈ resF64 rtol = 1e-10
end
