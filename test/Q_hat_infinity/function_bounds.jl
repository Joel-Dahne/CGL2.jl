@testset "Q_infinity_function_bounds" begin
    # IMPROVE: Test for more parameters
    ξ = Arb(30)
    κ = Arb(0.8073018593981386)
    ϵ = Arb(0.15002213424487343)
    ξ₁ = Arb(15)
    Λ = CGLParams{Arb}(3, 1.0, 1.0, 0.0)
    (; d, ω, σ, δ) = Λ
    a, b, c = CGL2._abc(κ, ϵ, Λ)

    C = CGL2.FunctionBounds_hat(κ, ϵ, ξ₁, Λ)

    P_hat, P_hat_dξ = CGL2.P_hat, CGL2.P_hat_dξ
    E_hat, E_hat_dξ = CGL2.E_hat, CGL2.E_hat_dξ
    J_P_hat = CGL2.J_P_hat
    J_E_hat = CGL2.J_E_hat

    # Test bounds from Lemma REF(lemma:P_hat-E_hat-bounds)
    for ξ in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64] .* ξ₁
        ####
        ## P_hat
        ####
        @test abs(P_hat(ξ, κ, ϵ, Λ)) <= C.P_hat * ξ^(-1 / σ)
        @test abs(P_hat(ξ, κ, ϵ, Λ)) >= 0.99C.P_hat * ξ^(-1 / σ)

        @test abs(P_hat_dξ(ξ, κ, ϵ, Λ)) <= C.P_hat_dξ * ξ^(-1 / σ - 1)
        @test abs(P_hat_dξ(ξ, κ, ϵ, Λ)) >= 0.95C.P_hat_dξ * ξ^(-1 / σ - 1)

        ####
        ## E_hat
        ####
        @test abs(E_hat(ξ, κ, ϵ, Λ)) <= C.E_hat * exp(-real(c) * ξ^2) * ξ^(1 / σ - d)
        @test abs(E_hat(ξ, κ, ϵ, Λ)) >= 0.95C.E_hat * exp(-real(c) * ξ^2) * ξ^(1 / σ - d)

        @test abs(E_hat_dξ(ξ, κ, ϵ, Λ)) <=
              C.E_hat_dξ * exp(-real(c) * ξ^2) * ξ^(1 / σ - d + 1)
        @test abs(E_hat_dξ(ξ, κ, ϵ, Λ)) >=
              0.95C.E_hat_dξ * exp(-real(c) * ξ^2) * ξ^(1 / σ - d + 1)

        ######
        ## J_P_hat
        ######
        @test abs(J_P_hat(ξ, κ, ϵ, Λ)) <=
              C.J_P_hat * exp(real(c) * ξ^2) * ξ^(-1 / σ + d - 1)
        @test abs(J_P_hat(ξ, κ, ϵ, Λ)) >=
              0.99C.J_P_hat * exp(real(c) * ξ^2) * ξ^(-1 / σ + d - 1)

        ######
        ## J_E_hat
        ######
        @test abs(J_E_hat(ξ, κ, ϵ, Λ)) <= C.J_E_hat * ξ^(1 / σ - 1)
        @test abs(J_E_hat(ξ, κ, ϵ, Λ)) >= 0.95C.J_E_hat * ξ^(1 / σ - 1)
    end

    # Test bounds from Lemma REF(lemma:I_E_hat-I_P_hat-bounds)
    for ξ in [1, 2, 4] .* ξ₁
        # We compute an approximation of the integral using Q_hat =
        # ξ^(-1 / σ)
        Q_hat = ξ -> ξ^(-1 / σ)

        ######
        ## I_E_hat
        ######
        # We approximate infinity by 5ξ
        I_E_hat_approx = Arblib.integrate(ξ, 5ξ, rtol = 1e-5) do ξ
            J_E_hat(ξ, κ, ϵ, Λ) * Q_hat(ξ)^(2σ + 1)
        end

        # Note that the norm of Q_hat is 1 by construction
        @test abs(I_E_hat_approx) <= C.I_E_hat * ξ^-2
        @test abs(I_E_hat_approx) >= 0.8C.I_E_hat * ξ^-2

        ######
        ## I_P_hat
        ######
        I_P_hat_approx = Arblib.integrate(ξ₁, ξ, rtol = 1e-5) do ξ
            J_P_hat(ξ, κ, ϵ, Λ) * Q_hat(ξ)^(2σ + 1)
        end

        # We only check the upper bound here. For ξ = ξ₁ the integral
        # is zero and in general the bound is very rough.
        # Note that the norm of Q_hat is 1 by construction
        @test abs(I_P_hat_approx) <= C.I_P_hat * exp(real(c) * ξ^2) * ξ^(-2 / σ + d - 2)
    end
end
