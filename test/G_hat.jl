@testset "G_hat" begin
    ν = Acb(2.1134119144964454, 2.5294091893480206)
    γ₁ = Acb(0.19239527572579992, 0.14584905671816073)
    γ₂ = Acb(-1121.5787719034975, -987.6373883459842)
    κ = Arb(0.8073018593981386)
    ϵ = Arb(0.15002213424487343)
    ξ₁ = Arb(15) # A smaller value makes the tests more effective
    Λ = CGLParams{Arb}(3, 1.0, 1.0, 0.0)

    νF64 = ComplexF64(ν)
    γ₁F64 = ComplexF64(γ₁)
    γ₂F64 = ComplexF64(γ₂)
    κF64 = Float64(κ)
    ϵF64 = Float64(ϵ)
    ξ₁F64 = Float64(ξ₁)
    ΛF64 = CGLParams{Float64}(Λ)

    G_hat_jacobian = Float64.(CGL2.G_hat_jacobian(ν, γ₁, γ₂, κ, ϵ, ξ₁, Λ))
    G_hat_jacobian_fdm = jacobian(
        central_fdm(5, 1, factor = 1e5),
        ((ν_real, ν_imag, γ₂_real, γ₂_imag),) -> CGL2.G_hat(
            complex(ν_real, ν_imag),
            γ₁F64,
            complex(γ₂_real, γ₂_imag),
            κF64,
            ϵF64,
            ξ₁F64,
            ΛF64,
        ),
        [real(νF64), imag(νF64), real(γ₂F64), imag(γ₂F64)],
    )[1]

    @test G_hat_jacobian ≈ G_hat_jacobian_fdm rtol = 1e-9
    @test all(isapprox.(G_hat_jacobian, G_hat_jacobian_fdm, rtol = 1e-8))
end
