@testset "Y_infinity_function_bounds" begin
    λ = Acb(0.19, 2.76)
    γ₁ = Acb(0.19, 0.15)
    γ₂ = Acb(-116.03, 101.21)
    κ = Arb(0.80)
    ϵ = Arb(0.15)
    ξ₁ = Arb(15)
    Λ = CGLParams{Arb}(3, 1.0, 1.0, 0.0)
    (; d, σ) = Λ
    c = CGL2._c(κ, ϵ, Λ)
    A = SMatrix{2,2}(ϵ, 1, -1, ϵ)
    V = SMatrix{2,2}(im, 1, -im, 1)

    C_Y = CGL2.FunctionBounds_Y(λ, γ₁, γ₂, κ, ϵ, ξ₁, Λ)

    norm_inf = CGL2.norm_inf

    E_1 = CGL2.E_1
    E_2 = CGL2.E_2
    P_1 = CGL2.P_1
    P_2 = CGL2.P_2

    J_E_1 = CGL2.J_E_1
    J_E_2 = CGL2.J_E_2
    J_P_1 = CGL2.J_P_1
    J_P_2 = CGL2.J_P_2

    K_1, K_2 = CGL2.K_1, CGL2.K_2

    for ξ in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64] .* ξ₁
        ####
        ## E_1 and E_2
        ####

        @test abs(E_1(ξ, λ, κ, ϵ, Λ)) <=
              C_Y.E_1 * exp(-real(c) * ξ^2) * ξ^(1 / σ - d - real(λ) / κ)
        @test abs(E_2(ξ, λ, κ, ϵ, Λ)) <=
              C_Y.E_2 * exp(-real(c) * ξ^2) * ξ^(1 / σ - d - real(λ) / κ)
        @test abs(E_1(ξ, λ, κ, ϵ, Λ)) >=
              0.95C_Y.E_1 * exp(-real(c) * ξ^2) * ξ^(1 / σ - d - real(λ) / κ)
        @test abs(E_2(ξ, λ, κ, ϵ, Λ)) >=
              0.85C_Y.E_2 * exp(-real(c) * ξ^2) * ξ^(1 / σ - d - real(λ) / κ)

        ####
        ## P_1 and P_2
        ####

        @test abs(P_1(ξ, λ, κ, ϵ, Λ)) <= C_Y.P_1 * ξ^(-1 / σ + real(λ) / κ)
        @test abs(P_2(ξ, λ, κ, ϵ, Λ)) <= C_Y.P_2 * ξ^(-1 / σ + real(λ) / κ)
        @test abs(P_1(ξ, λ, κ, ϵ, Λ)) >= 0.95C_Y.P_1 * ξ^(-1 / σ + real(λ) / κ)
        @test abs(P_2(ξ, λ, κ, ϵ, Λ)) >= 0.90C_Y.P_2 * ξ^(-1 / σ + real(λ) / κ)

        ######
        ## J_E_1 and J_E_2
        ######

        @test abs(J_E_1(ξ, λ, κ, ϵ, Λ)) <= C_Y.J_E_1 * ξ^(1 / σ - real(λ) / κ - 1)
        @test abs(J_E_2(ξ, λ, κ, ϵ, Λ)) <= C_Y.J_E_2 * ξ^(1 / σ - real(λ) / κ - 1)
        @test abs(J_E_1(ξ, λ, κ, ϵ, Λ)) >= 0.95C_Y.J_E_1 * ξ^(1 / σ - real(λ) / κ - 1)
        @test abs(J_E_2(ξ, λ, κ, ϵ, Λ)) >= 0.85C_Y.J_E_2 * ξ^(1 / σ - real(λ) / κ - 1)

        ######
        ## J_P_1 and J_P_2
        ######

        @test abs(J_P_1(ξ, λ, κ, ϵ, Λ)) <=
              C_Y.J_P_1 * exp(real(c) * ξ^2) * ξ^(- 1 / σ + d + real(λ) / κ - 1)
        @test abs(J_P_2(ξ, λ, κ, ϵ, Λ)) <=
              C_Y.J_P_2 * exp(real(c) * ξ^2) * ξ^(- 1 / σ + d + real(λ) / κ - 1)
        @test abs(J_P_1(ξ, λ, κ, ϵ, Λ)) >=
              0.95C_Y.J_P_1 * exp(real(c) * ξ^2) * ξ^(- 1 / σ + d + real(λ) / κ - 1)
        @test abs(J_P_2(ξ, λ, κ, ϵ, Λ)) >=
              0.9C_Y.J_P_2 * exp(real(c) * ξ^2) * ξ^(- 1 / σ + d + real(λ) / κ - 1)

        ######
        ## K_1, K_2
        ######

        @test norm_inf(K_1(ξ, λ, κ, ϵ, Λ), 1) <=
              C_Y.K_1_1 * exp(real(c) * ξ^2) * ξ^(- 1 / σ + d + real(λ) / κ - 1)
        @test norm_inf(K_1(ξ, λ, κ, ϵ, Λ), 2) <=
              C_Y.K_1_2 * exp(real(c) * ξ^2) * ξ^(- 1 / σ + d + real(λ) / κ - 1)
        @test norm_inf(K_1(ξ, λ, κ, ϵ, Λ), 1) >=
              0.95C_Y.K_1_1 * exp(real(c) * ξ^2) * ξ^(- 1 / σ + d + real(λ) / κ - 1)
        @test norm_inf(K_1(ξ, λ, κ, ϵ, Λ), 2) >=
              0.9C_Y.K_1_2 * exp(real(c) * ξ^2) * ξ^(- 1 / σ + d + real(λ) / κ - 1)

        @test norm_inf(K_2(ξ, λ, κ, ϵ, Λ), 1) <= C_Y.K_2_1 * ξ^(1 / σ - real(λ) / κ - 1)
        @test norm_inf(K_2(ξ, λ, κ, ϵ, Λ), 2) <= C_Y.K_2_2 * ξ^(1 / σ - real(λ) / κ - 1)
        @test norm_inf(K_2(ξ, λ, κ, ϵ, Λ), 1) >= 0.9C_Y.K_2_1 * ξ^(1 / σ - real(λ) / κ - 1)
        @test norm_inf(K_2(ξ, λ, κ, ϵ, Λ), 2) >= 0.85C_Y.K_2_2 * ξ^(1 / σ - real(λ) / κ - 1)

        ####
        ## I_N
        ####

        # NOTE: This is not quite correct since γ₁ and γ₂ only
        # parametrize Q_hat_infinity at ξ₁ and not at ξ > ξ₁. The
        # difference should however be of lower order and it still
        # gives some information as a test.
        Q_hat, _ = CGL2.Q_hat_infinity(γ₁, γ₂, κ, ϵ, ξ, Λ)

        @test norm_inf(CGL2.I_N(Q_hat, Λ)) <= C_Y.I_N * ξ^-2
        @test norm_inf(CGL2.I_N(Q_hat, Λ)) >= 0.8C_Y.I_N * ξ^-2

        # NOTE: There is no clear way for how to directly test
        # I_K_1_1, I_K_1_2, I_K_2_1, I_K_2_2 and T_12. There are
        # therefore no tests for them here.
    end
end
