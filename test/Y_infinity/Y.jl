@testset "Y_infinity" begin
    c_0 = SVector(Acb(1), Acb(0))
    lambda = Acb(0.19028080950252219, 2.769747313863597)
    ν = Acb(2.1134119144964454, 2.5294091893480206)
    γ₁ = Acb(0.1979686615344728, 0.15613711743772288)
    γ₂ = Acb(-116.03421868436834, 101.21494871936197)
    κ = Arb(0.8073018593981386)
    ϵ = Arb(0.15002213424487343)
    ξ₁ = Arb(30)
    λ = CGLParams{Arb}(3, 1.0, 1.0, 0.0)

    c_0F64 = ComplexF64.(c_0)
    lambdaF64 = ComplexF64(lambda)
    νF64 = ComplexF64(ν)
    γ₁F64 = ComplexF64(γ₁)
    γ₂F64 = ComplexF64(γ₂)
    κF64 = Float64(κ)
    ϵF64 = Float64(ϵ)
    ξ₁F64 = Float64(ξ₁)
    λF64 = CGLParams{Float64}(λ)

    # Function for computing derivative using finite differences.
    fdm = central_fdm(5, 1)

    res = CGL2.Y_infinity(c_0, lambda, γ₁, γ₂, κ, ϵ, ξ₁, λ)

    # Compare with Float64 version
    resF64 = CGL2.Y_infinity(c_0F64, lambdaF64, γ₁F64, γ₂F64, κF64, ϵF64, ξ₁F64, λF64)
    @test ComplexF64.(res) ≈ resF64 rtol = 1e-13

    res_J = CGL2.Y_infinity_derivative(c_0, lambda, γ₁, γ₂, κ, ϵ, ξ₁, λ)

    # Compare with Float64 version
    resF64_J =
        CGL2.Y_infinity_derivative(c_0F64, lambdaF64, γ₁F64, γ₂F64, κF64, ϵF64, ξ₁F64, λF64)
    @test ComplexF64.(res_J) ≈ resF64_J rtol = 1e-13

    # Compare with fdm method
    resF64_J_fdm = fdm(
        lambda_real -> CGL2.Y_infinity(
            c_0F64,
            complex(lambda_real, imag(lambdaF64)),
            γ₁F64,
            γ₂F64,
            κF64,
            ϵF64,
            ξ₁F64,
            λF64,
        ),
        real(lambdaF64),
    )

    # Since we use a very rough enclosure for I_K_2_dλ the res_J
    # version is not particularly precise. We therefore only check
    # with rtol = 1e-4.
    @test ComplexF64.(res_J) ≈ resF64_J_fdm rtol = 1e-4
    # We also check that res_J contains the approximate version, due
    # to the finite difference is of course not guaranteed. But in
    # practice the errors from the finite difference should be much
    # smaller than the error bounds for res_J
    @test all(Arblib.contains.(res_J, Acb.(resF64_J_fdm)))
end
