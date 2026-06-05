@testset "Y_infinity" begin
    c_0 = SVector(Acb(1), Acb(2))
    λ = Acb(0.19, 2.76)
    ν = Acb(2.11, 2.52)
    γ₁ = Acb(0.19, 0.15)
    γ₂ = Acb(-116.03, 101.21)
    κ = Arb(0.80)
    ϵ = Arb(0.15)
    # We use a relatively large value for ξ₁ to reduce the errors when
    # we check if it solves the equation.
    ξ₁ = Arb(50)
    Λ = CGLParams{Arb}(3, 1, 1, 0)

    # Function for computing derivative using finite differences.
    fdm = central_fdm(5, 1)

    res = CGL2.Y_infinity(c_0, λ, γ₁, γ₂, κ, ϵ, ξ₁, Λ)

    # Compute derivative w.r.t. ξ using finite differences in Float64
    resF64_dξ =
        fdm(ξ -> ComplexF64.(CGL2.Y_infinity(c_0, λ, γ₁, γ₂, κ, ϵ, Arb(ξ), Λ)), Float64(ξ₁))

    # Check that the derivative computed with finite differences
    # agree. Note that the precision is quite low, so rtol is
    # relatively large.
    @test ComplexF64.(res[3:4]) ≈ resF64_dξ[1:2] rtol = 1e-5

    # Check if the finite difference value approximately satisfies the
    # equation

    (; d, ω, σ) = Λ
    I = SMatrix{2,2}(1, 0, 0, 1)
    J = SMatrix{2,2}(0, 1, -1, 0)
    A = ϵ * I + J
    B_1 = κ * I
    B_2 = (d - 1) * ϵ * I + (d - 1) * J
    C = κ / σ * I + ω * J
    a_hat, b_hat = CGL2.Q_hat_infinity(γ₁, γ₂, κ, ϵ, ξ₁, Λ)
    J_N = SMatrix{2,2}(
        -2a_hat * b_hat,
        3a_hat^2 + b_hat^2,
        -(a_hat^2 + 3b_hat^2),
        2a_hat * b_hat,
    )

    Y = ComplexF64.(res[1:2])
    Y_dξ = ComplexF64.(res[3:4])
    Y_dξ_dξ = resF64_dξ[3:4]

    # The error should be compared to the norm of the inputs. We check
    # that it is substantially smaller than the largest norm of the
    # input.
    @test norm(
        ComplexF64.(A * Y_dξ_dξ + (B_1 * ξ₁ + B_2 * ξ₁^-1) * Y_dξ + (C + J_N - λ * I) * Y),
    ) < 2e-6max(norm(Y), norm(Y_dξ), norm(Y_dξ_dξ))
end
