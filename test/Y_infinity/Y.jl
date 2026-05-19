@testset "Y_infinity" begin
    c_0 = SVector(Acb(1), Acb(0))
    lambda = Acb(0.19028080950252219, 2.769747313863597)
    ν = Acb(2.1134119144964454, 2.5294091893480206)
    γ₁ = Acb(0.1979686615344728, 0.15613711743772288)
    γ₂ = Acb(-116.03421868436834, 101.21494871936197)
    κ = Arb(0.8073018593981386)
    ϵ = Arb(0.15002213424487343)
    ξ₁ = Arb(30)
    Λ = CGLParams{Arb}(3, 1, 1, 0)

    # Function for computing derivative using finite differences.
    fdm = central_fdm(5, 1)

    res = CGL2.Y_infinity(c_0, lambda, γ₁, γ₂, κ, ϵ, ξ₁, Λ)

    # Compute derivative w.r.t. ξ using finite differences in Float64
    resF64_dξ = fdm(
        ξ -> ComplexF64.(CGL2.Y_infinity(c_0, lambda, γ₁, γ₂, κ, ϵ, Arb(ξ), Λ)),
        Float64(ξ₁),
    )

    # Check that the derivative computed with finite differences
    # agree. Note that the precision is quite low, so rtol is
    # relatively large.
    @test ComplexF64.(res[3:4]) ≈ resF64_dξ[1:2] rtol = 1e-3

    # Check if the finite difference value approximately satisfies the
    # equation

    (; d, ω, σ) = Λ
    I = SMatrix{2,2}(1, 0, 0, 1)
    J = SMatrix{2,2}(0, 1, -1, 0)
    A = ϵ * I + J
    B_1 = κ * I
    B_2 = (d - 1) * ϵ * I + (d - 1) * J
    C = κ / σ * I + ω * J
    a, b = CGL2.Q_hat_infinity(γ₁, γ₂, κ, ϵ, ξ₁, Λ)
    J_N = SMatrix{2,2}(-2a * b, 3a^2 + b^2, -(a^2 + 3b^2), 2a * b)

    Y = ComplexF64.(res[1:2])
    Y_dξ = ComplexF64.(res[3:4])
    Y_dξ_dξ = resF64_dξ[3:4]

    # The of the error should be compared to the norm of the inputs.
    # We check that it is substantially smaller than the largest norm
    # of the input.
    @test norm(
        ComplexF64.(
            A * Y_dξ_dξ + (B_1 * ξ₁ + B_2 * ξ₁^-1) * Y_dξ + (C + J_N - lambda * I) * Y,
        ),
    ) < 1e-5max(norm(Y), norm(Y_dξ), norm(Y_dξ_dξ))
end
