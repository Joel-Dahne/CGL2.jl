@testset "Q_hat_zero" begin
    ν = Acb(2.1134119144964454, 2.5294091893480206)
    κ = Arb(0.8073018593981386)
    ϵ = Arb(0.15002213424487343)
    ξ₁ = Arb(30)
    Λ = CGLParams{Arb}(3, 1.0, 1.0, 0.0)

    νF64 = ComplexF64(ν)
    κF64 = Float64(κ)
    ϵF64 = Float64(ϵ)
    ξ₁F64 = Float64(ξ₁)
    ΛF64 = CGLParams{Float64}(Λ)

    res = CGL2.Q_hat_zero(ν, κ, ϵ, ξ₁, Λ)
    resF64 = CGL2.Q_hat_zero(νF64, κF64, ϵF64, ξ₁F64, ΛF64)
    @test res ≈ resF64 rtol = 1e-11

    res_J = CGL2.Q_hat_zero_jacobian(ν, κ, ϵ, ξ₁, Λ)
    resF64_J = CGL2.Q_hat_zero_jacobian(νF64, κF64, ϵF64, ξ₁F64, ΛF64)
    resF64_J_fdm = jacobian(
        central_fdm(5, 1, factor = 1e5),
        ((νF64_real, νF64_imag),) ->
            CGL2.Q_hat_zero(complex(νF64_real, νF64_imag), κF64, ϵF64, ξ₁F64, ΛF64),
        [real(νF64), imag(νF64)],
    )[1]

    @test res_J ≈ resF64_J rtol = 1e-9
    @test real(resF64_J) ≈ resF64_J_fdm[1:2:3, :] rtol = 1e-8
    @test imag(resF64_J) ≈ resF64_J_fdm[2:2:4, :] rtol = 1e-8
end
