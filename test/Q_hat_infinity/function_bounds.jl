@testset "Q_infinity_function_bounds" begin
    # IMPROVE: Test for more parameters
    ξ = Arb(30)
    κ = Arb(0.8073018593981386)
    ϵ = Arb(0.15002213424487343)
    ξ₁ = Arb(15)
    λ = CGLParams{Arb}(3, 1.0, 1.0, 0.0)
    (; d, ω, σ, δ) = λ

    C = CGL2.FunctionBounds_hat(κ, ϵ, ξ₁, λ)

    P_hat, P_hat_dξ = CGL2.P_hat, CGL2.P_hat_dξ
    E_hat, E_hat_dξ = CGL2.E_hat, CGL2.E_hat_dξ
    J_P_hat = CGL2.J_P_hat
    J_E_hat = CGL2.J_E_hat

    for ξ in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64] .* ξ₁
        ####
        ## P_hat
        ####
        @test abs(P_hat(ξ, κ, ϵ, λ)) <= C.P_hat * ξ^(-1 / σ)
        @test abs(P_hat(ξ, κ, ϵ, λ)) >= 0.99C.P_hat * ξ^(-1 / σ)

        ####
        ## E_hat
        ####
        @test abs(E_hat(ξ, κ, ϵ, λ)) <= C.E_hat * exp(-real(c) * ξ^2) * ξ^(1 / σ - d)
        @test abs(E_hat(ξ, κ, ϵ, λ)) >= 0.95C.E_hat * exp(-real(c) * ξ^2) * ξ^(1 / σ - d)

        ######
        ## J_P_hat
        ######
        @test abs(J_P_hat(ξ, κ, ϵ, λ)) <=
              C.J_P_hat * exp(real(c) * ξ^2) * ξ^(-1 / σ + d - 1)
        @test abs(J_P_hat(ξ, κ, ϵ, λ)) >=
              0.99C.J_P_hat * exp(real(c) * ξ^2) * ξ^(-1 / σ + d - 1)

        ######
        ## J_E_hat
        ######
        @test abs(J_E_hat(ξ, κ, ϵ, λ)) <= C.J_E_hat * ξ^(1 / σ - 1)
        @test abs(J_E_hat(ξ, κ, ϵ, λ)) >= 0.95C.J_E_hat * ξ^(1 / σ - 1)
    end
end
