@testset "functions" begin
    ξ = Arb(30)
    κ = Arb(0.8073018593981386)
    ϵ = Arb(0.15002213424487343)
    λ = CGLParams{Arb}(3, 1, 1, 0.0)
    (; d, ω, σ, δ) = λ

    c = CGL2._c(κ, ϵ, λ)

    ξF64 = Float64(ξ)
    κF64 = Float64(κ)
    ϵF64 = Float64(ϵ)
    λF64 = CGLParams{Float64}(λ)

    # Function for computing derivative using finite differences.
    fdm = central_fdm(5, 1)

    P_hat, P_hat_dξ = CGL2.P_hat, CGL2.P_hat_dξ
    E_hat, E_hat_dξ = CGL2.E_hat, CGL2.E_hat_dξ
    W_hat = CGL2.W_hat
    B_W_hat = CGL2.B_W_hat
    J_E_hat = CGL2.J_E_hat
    J_P_hat = CGL2.J_P_hat

    @testset "P_hat" begin
        @test Arblib.overlaps(P_hat(ArbSeries((ξ, 1)), κ, ϵ, λ)[1], P_hat_dξ(ξ, κ, ϵ, λ))

        @test P_hat_dξ(ξ, κ, ϵ, λ) ≈ fdm(ξ -> P_hat(ξ, κF64, ϵF64, λF64), ξF64) rtol = 1e-12
    end

    @testset "E_hat" begin
        @test Arblib.overlaps(E_hat(ArbSeries((ξ, 1)), κ, ϵ, λ)[1], E_hat_dξ(ξ, κ, ϵ, λ))

        @test E_hat_dξ(ξ, κ, ϵ, λ) ≈ fdm(ξ -> E_hat(ξ, κF64, ϵF64, λF64), ξF64) rtol = 1e-10
    end

    @testset "W_hat" begin
        @test Arblib.overlaps(
            W_hat(ξ, κ, ϵ, λ),
            P_hat(ξ, κ, ϵ, λ) * E_hat_dξ(ξ, κ, ϵ, λ) -
            P_hat_dξ(ξ, κ, ϵ, λ) * E_hat(ξ, κ, ϵ, λ),
        )
    end

    @testset "J_E_hat" begin
        @test Arblib.overlaps(
            J_E_hat(ξ, κ, ϵ, λ),
            -(1 + im * δ) / (1 - im * ϵ) * E_hat(ξ, κ, ϵ, λ) / W_hat(ξ, κ, ϵ, λ),
        )

        @test Arblib.overlaps(
            J_E_hat(ξ, κ, ϵ, λ),
            B_W_hat(κ, ϵ, λ) * E_hat(ξ, κ, ϵ, λ) * exp(c * ξ^2) * ξ^(λ.d - 1),
        )
    end

    @testset "J_P_hat" begin
        @test Arblib.overlaps(
            J_P_hat(ξ, κ, ϵ, λ),
            -(1 + im * δ) / (1 - im * ϵ) * P_hat(ξ, κ, ϵ, λ) / W_hat(ξ, κ, ϵ, λ),
        )

        @test Arblib.overlaps(
            J_P_hat(ξ, κ, ϵ, λ),
            B_W_hat(κ, ϵ, λ) * P_hat(ξ, κ, ϵ, λ) * exp(c * ξ^2) * ξ^(λ.d - 1),
        )
    end
end
