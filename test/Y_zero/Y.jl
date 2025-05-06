@testset "Y_zero" begin
    x = Acb(-0.17330123921225876, 0.1393584447527974)
    lambda = Acb(0.19028080950252219, 2.769747313863597)
    ν = Acb(2.1134119144964454, 2.5294091893480206)
    κ = Arb(0.8073018593981386)
    ϵ = Arb(0.15002213424487343)
    ξ₁ = Arb(30)
    λ = CGLParams{Arb}(3, 1.0, 1.0, 0.0)

    xF64 = ComplexF64(x)
    lambdaF64 = ComplexF64(lambda)
    νF64 = ComplexF64(ν)
    κF64 = Float64(κ)
    ϵF64 = Float64(ϵ)
    ξ₁F64 = Float64(ξ₁)
    λF64 = CGLParams{Float64}(λ)

    res = CGL2.Y_zero(x, lambda, ν, κ, ϵ, ξ₁, λ)

    # Larger interval and lower degree for Taylor expansion
    # TODO: Need remainder terms!
    res2 = CGL2.Y_zero_capd(x, lambda, ν, κ, ϵ, Arb(1e-1), ξ₁, λ, degree = 5)
    @test all(Arblib.overlaps.(res, res2)) broken = true

    # Compare with Float64 version
    resF64 = CGL2.Y_zero(xF64, lambdaF64, νF64, κF64, ϵF64, ξ₁F64, λF64)
    @test ComplexF64.(res) ≈ resF64 rtol = 1e-10

    res_J = CGL2.Y_zero_jacobian(x, lambda, ν, κ, ϵ, ξ₁, λ)

    # Larger interval and lower degree for Taylor expansion
    # TODO: Need remainder terms!
    res2_J = CGL2.Y_zero_jacobian_capd(x, lambda, ν, κ, ϵ, Arb(1e-1), ξ₁, λ, degree = 5)
    @test all(Arblib.overlaps.(res_J, res2_J)) broken = true

    # Compare with Float64 version
    resF64_J = CGL2.Y_zero_jacobian(xF64, lambdaF64, νF64, κF64, ϵF64, ξ₁F64, λF64)
    @test ComplexF64.(res_J) ≈ resF64_J rtol = 1e-11

    # Check that Y_zero_float and Y_zero_float_real agree
    let Q_hat = CGL2.Q_hat_zero_float_curve(real(νF64), imag(νF64), κF64, ϵF64, ξ₁F64, λF64)
        r1 = CGL2.Y_zero_float(xF64, lambdaF64, κF64, ϵF64, ξ₁F64, Q_hat, λF64)
        r2 = CGL2.Y_zero_float_real(
            real(xF64),
            imag(xF64),
            real(lambdaF64),
            imag(lambdaF64),
            κF64,
            ϵF64,
            ξ₁F64,
            Q_hat,
            λF64,
        )
        r2 = complex.(
            SVector(r2[1], r2[2], r2[5], r2[6]),
            SVector(r2[3], r2[4], r2[7], r2[8]),
        )
        @test r1 ≈ r2
    end
end
