@testset "Y_infinity_function_bounds" begin
    lambda = Acb(0.19028080950252219, 2.769747313863597)
    γ₁ = Acb(0.1979686615344728, 0.15613711743772288)
    γ₂ = Acb(-116.03421868436834, 101.21494871936197)
    κ = Arb(0.8073018593981386)
    ϵ = Arb(0.15002213424487343)
    ξ₁ = Arb(15)
    λ = CGLParams{Arb}(3, 1.0, 1.0, 0.0)
    v = Arb("0.1")
    (; d, σ) = λ
    _, _, c = CGL2._abc(κ, ϵ, λ)
    A = SMatrix{2,2}(ϵ, 1, -1, ϵ)
    M = SMatrix{2,2}(im, 1, -im, 1)

    C_Y = CGL2.FunctionBounds_Y(lambda, γ₁, γ₂, κ, ϵ, ξ₁, v, λ)

    norm_inf = CGL2.norm_inf

    E_1, E_1_dξ = CGL2.E_1, CGL2.E_1_dξ
    E_2, E_2_dξ = CGL2.E_2, CGL2.E_2_dξ
    P_1, P_1_dξ = CGL2.P_1, CGL2.P_1_dξ
    P_2, P_2_dξ = CGL2.P_2, CGL2.P_2_dξ

    # Compute derivative of exp(c * ξ^2) * E_1(ξ)
    function exp_E_1_dξ(ξ, lambda, κ, ϵ, λ)
        a, b, c = CGL2._abc(κ, ϵ, λ)
        return exp(c * ξ^2) *
               (2c * ξ * E_1(ξ, lambda, κ, ϵ, λ) + E_1_dξ(ξ, lambda, κ, ϵ, λ))
    end
    # Compute derivative of exp(conj(c) * ξ^2) * E_2(ξ)
    function exp_E_2_dξ(ξ, lambda, κ, ϵ, λ)
        a, b, c = CGL2._abc(κ, ϵ, λ)
        return exp(conj(c) * ξ^2) *
               (2conj(c) * ξ * E_2(ξ, lambda, κ, ϵ, λ) + E_2_dξ(ξ, lambda, κ, ϵ, λ))
    end

    function exp_P_1_dξ(ξ, lambda, κ, ϵ, λ)
        a, b, c = CGL2._abc(κ, ϵ, λ)
        return exp(c * ξ^2) *
               (2c * ξ * P_1(ξ, lambda, κ, ϵ, λ) + P_1_dξ(ξ, lambda, κ, ϵ, λ))
    end
    # Compute derivative of exp(conj(c) * ξ^2) * P_2(ξ)
    function exp_P_2_dξ(ξ, lambda, κ, ϵ, λ)
        a, b, c = CGL2._abc(κ, ϵ, λ)
        return exp(conj(c) * ξ^2) *
               (2conj(c) * ξ * P_2(ξ, lambda, κ, ϵ, λ) + P_2_dξ(ξ, lambda, κ, ϵ, λ))
    end

    J_E_1 = CGL2.J_E_1
    J_E_2 = CGL2.J_E_2
    J_P_1 = CGL2.J_P_1
    J_P_2 = CGL2.J_P_2

    J_E_1_dξ = CGL2.J_E_1_dξ
    J_E_2_dξ = CGL2.J_E_2_dξ

    K_1, K_2 = CGL2.K_1, CGL2.K_2
    K_1_dξ, K_2_dξ = CGL2.K_1_dξ, CGL2.K_2_dξ

    for ξ in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64] .* ξ₁
        ####
        ## I_N
        ####

        # TODO: This is not quite correct since γ₁ and γ₂ only
        # parametrize Q_hat_infinity at ξ₁ and not at ξ > ξ₁. The
        # difference should however be of lower order and it still
        # gives some information as a test.
        Q_hat, Q_hat_dξ = CGL2.Q_hat_infinity(γ₁, γ₂, κ, ϵ, ξ, λ)

        @test norm_inf(CGL2.I_N(Q_hat, λ)) <= C_Y.I_N * ξ^-2
        @test norm_inf(CGL2.I_N(Q_hat, λ)) >= 0.8C_Y.I_N * ξ^-2
        # IMPROVE: We could maybe improve on this bound?
        @test norm_inf(CGL2.I_N_dξ(Q_hat, Q_hat_dξ, λ)) <= C_Y.I_N_dξ * ξ^-3
        @test norm_inf(CGL2.I_N_dξ(Q_hat, Q_hat_dξ, λ)) >= 0.6C_Y.I_N_dξ * ξ^-3

        ####
        ## E_1 and E_2
        ####

        @test abs(E_1(ξ, lambda, κ, ϵ, λ)) <=
              C_Y.E_1 * exp(-real(c) * ξ^2) * ξ^(1 / σ - d - real(lambda) / κ)
        @test abs(E_2(ξ, lambda, κ, ϵ, λ)) <=
              C_Y.E_2 * exp(-real(c) * ξ^2) * ξ^(1 / σ - d - real(lambda) / κ)
        @test abs(E_1(ξ, lambda, κ, ϵ, λ)) >=
              0.95C_Y.E_1 * exp(-real(c) * ξ^2) * ξ^(1 / σ - d - real(lambda) / κ)
        @test abs(E_2(ξ, lambda, κ, ϵ, λ)) >=
              0.85C_Y.E_2 * exp(-real(c) * ξ^2) * ξ^(1 / σ - d - real(lambda) / κ)

        @test abs(E_1_dξ(ξ, lambda, κ, ϵ, λ)) <=
              C_Y.E_1_dξ * exp(-real(c) * ξ^2) * ξ^(1 / σ - d - real(lambda) / κ + 1)
        @test abs(E_2_dξ(ξ, lambda, κ, ϵ, λ)) <=
              C_Y.E_2_dξ * exp(-real(c) * ξ^2) * ξ^(1 / σ - d - real(lambda) / κ + 1)
        @test abs(E_1_dξ(ξ, lambda, κ, ϵ, λ)) >=
              0.95C_Y.E_1_dξ * exp(-real(c) * ξ^2) * ξ^(1 / σ - d - real(lambda) / κ + 1)
        @test abs(E_2_dξ(ξ, lambda, κ, ϵ, λ)) >=
              0.85C_Y.E_2_dξ * exp(-real(c) * ξ^2) * ξ^(1 / σ - d - real(lambda) / κ + 1)

        @test abs(exp_E_1_dξ(ξ, lambda, κ, ϵ, λ)) <=
              C_Y.exp_E_1_dξ * ξ^(1 / σ - d - real(lambda) / κ - 1)
        @test abs(exp_E_2_dξ(ξ, lambda, κ, ϵ, λ)) <=
              C_Y.exp_E_2_dξ * ξ^(1 / σ - d - real(lambda) / κ - 1)
        @test abs(exp_E_1_dξ(ξ, lambda, κ, ϵ, λ)) >=
              0.95C_Y.exp_E_1_dξ * ξ^(1 / σ - d - real(lambda) / κ - 1)
        @test abs(exp_E_2_dξ(ξ, lambda, κ, ϵ, λ)) >=
              0.85C_Y.exp_E_2_dξ * ξ^(1 / σ - d - real(lambda) / κ - 1)

        ####
        ## P_1 and P_2
        ####

        @test abs(P_1(ξ, lambda, κ, ϵ, λ)) <= C_Y.P_1 * ξ^(-1 / σ + real(lambda) / κ)
        @test abs(P_2(ξ, lambda, κ, ϵ, λ)) <= C_Y.P_2 * ξ^(-1 / σ + real(lambda) / κ)
        @test abs(P_1(ξ, lambda, κ, ϵ, λ)) >= 0.95C_Y.P_1 * ξ^(-1 / σ + real(lambda) / κ)
        @test abs(P_2(ξ, lambda, κ, ϵ, λ)) >= 0.90C_Y.P_2 * ξ^(-1 / σ + real(lambda) / κ)

        @test abs(P_1_dξ(ξ, lambda, κ, ϵ, λ)) <=
              C_Y.P_1_dξ * ξ^(-1 / σ + real(lambda) / κ - 1)
        @test abs(P_2_dξ(ξ, lambda, κ, ϵ, λ)) <=
              C_Y.P_2_dξ * ξ^(-1 / σ + real(lambda) / κ - 1)
        @test abs(P_1_dξ(ξ, lambda, κ, ϵ, λ)) >=
              0.95C_Y.P_1_dξ * ξ^(-1 / σ + real(lambda) / κ - 1)
        @test abs(P_2_dξ(ξ, lambda, κ, ϵ, λ)) >=
              0.85C_Y.P_2_dξ * ξ^(-1 / σ + real(lambda) / κ - 1)

        @test abs(exp_P_1_dξ(ξ, lambda, κ, ϵ, λ)) <=
              C_Y.exp_P_1_dξ * exp(real(c) * ξ^2) * ξ^(-1 / σ + real(lambda) / κ + 1)
        @test abs(exp_P_2_dξ(ξ, lambda, κ, ϵ, λ)) <=
              C_Y.exp_P_2_dξ * exp(real(c) * ξ^2) * ξ^(-1 / σ + real(lambda) / κ + 1)
        @test abs(exp_P_1_dξ(ξ, lambda, κ, ϵ, λ)) >=
              0.95C_Y.exp_P_1_dξ * exp(real(c) * ξ^2) * ξ^(-1 / σ + real(lambda) / κ + 1)
        @test abs(exp_P_2_dξ(ξ, lambda, κ, ϵ, λ)) >=
              0.9C_Y.exp_P_2_dξ * exp(real(c) * ξ^2) * ξ^(-1 / σ + real(lambda) / κ + 1)

        ######
        ## J_E_1 and J_E_2
        ######

        @test abs(J_E_1(ξ, lambda, κ, ϵ, λ)) <= C_Y.J_E_1 * ξ^(1 / σ - real(lambda) / κ - 1)
        @test abs(J_E_1(ξ, lambda, κ, ϵ, λ)) >=
              0.95C_Y.J_E_1 * ξ^(1 / σ - real(lambda) / κ - 1)

        @test abs(J_E_2(ξ, lambda, κ, ϵ, λ)) <= C_Y.J_E_2 * ξ^(1 / σ - real(lambda) / κ - 1)
        @test abs(J_E_2(ξ, lambda, κ, ϵ, λ)) >=
              0.85C_Y.J_E_2 * ξ^(1 / σ - real(lambda) / κ - 1)

        @test abs(J_E_1_dξ(ξ, lambda, κ, ϵ, λ)) <=
              C_Y.J_E_1_dξ * ξ^(1 / σ - real(lambda) / κ - 2)
        @test abs(J_E_1_dξ(ξ, lambda, κ, ϵ, λ)) >=
              0.95C_Y.J_E_1_dξ * ξ^(1 / σ - real(lambda) / κ - 2)

        @test abs(J_E_2_dξ(ξ, lambda, κ, ϵ, λ)) <=
              C_Y.J_E_2_dξ * ξ^(1 / σ - real(lambda) / κ - 2)
        @test abs(J_E_2_dξ(ξ, lambda, κ, ϵ, λ)) >=
              0.85C_Y.J_E_2_dξ * ξ^(1 / σ - real(lambda) / κ - 2)

        ######
        ## J_P_1 and J_P_2
        ######

        @test abs(J_P_1(ξ, lambda, κ, ϵ, λ)) <=
              C_Y.J_P_1 * exp(real(c) * ξ^2) * ξ^(- 1 / σ + d + real(lambda) / κ - 1)
        @test abs(J_P_1(ξ, lambda, κ, ϵ, λ)) >=
              0.95C_Y.J_P_1 * exp(real(c) * ξ^2) * ξ^(- 1 / σ + d + real(lambda) / κ - 1)

        @test abs(J_P_2(ξ, lambda, κ, ϵ, λ)) <=
              C_Y.J_P_2 * exp(real(c) * ξ^2) * ξ^(- 1 / σ + d + real(lambda) / κ - 1)
        @test abs(J_P_2(ξ, lambda, κ, ϵ, λ)) >=
              0.9C_Y.J_P_2 * exp(real(c) * ξ^2) * ξ^(- 1 / σ + d + real(lambda) / κ - 1)

        ######
        ## K_1, K_2
        ######

        @test norm_inf(K_1(ξ, lambda, κ, ϵ, λ), 1) <=
              C_Y.K_1_1 * exp(real(c) * ξ^2) * ξ^(- 1 / σ + d + real(lambda) / κ - 1)
        @test norm_inf(K_1(ξ, lambda, κ, ϵ, λ), 1) >=
              0.95C_Y.K_1_1 * exp(real(c) * ξ^2) * ξ^(- 1 / σ + d + real(lambda) / κ - 1)
        @test norm_inf(K_1(ξ, lambda, κ, ϵ, λ), 2) <=
              C_Y.K_1_2 * exp(real(c) * ξ^2) * ξ^(- 1 / σ + d + real(lambda) / κ - 1)
        @test norm_inf(K_1(ξ, lambda, κ, ϵ, λ), 2) >=
              0.9C_Y.K_1_2 * exp(real(c) * ξ^2) * ξ^(- 1 / σ + d + real(lambda) / κ - 1)

        @test norm_inf(K_2(ξ, lambda, κ, ϵ, λ), 1) <=
              C_Y.K_2_1 * ξ^(1 / σ - real(lambda) / κ - 1)
        @test norm_inf(K_2(ξ, lambda, κ, ϵ, λ), 1) >=
              0.9C_Y.K_2_1 * ξ^(1 / σ - real(lambda) / κ - 1)
        @test norm_inf(K_2(ξ, lambda, κ, ϵ, λ), 2) <=
              C_Y.K_2_2 * ξ^(1 / σ - real(lambda) / κ - 1)
        @test norm_inf(K_2(ξ, lambda, κ, ϵ, λ), 2) >=
              0.85C_Y.K_2_2 * ξ^(1 / σ - real(lambda) / κ - 1)

        @test norm_inf(K_2_dξ(ξ, lambda, κ, ϵ, λ), 1) <=
              C_Y.K_2_dξ_1 * ξ^(1 / σ - real(lambda) / κ - 2)
        @test norm_inf(K_2_dξ(ξ, lambda, κ, ϵ, λ), 1) >=
              0.4C_Y.K_2_dξ_1 * ξ^(1 / σ - real(lambda) / κ - 2)
        @test norm_inf(K_2_dξ(ξ, lambda, κ, ϵ, λ), 2) <=
              C_Y.K_2_dξ_2 * ξ^(1 / σ - real(lambda) / κ - 2)
        @test norm_inf(K_2_dξ(ξ, lambda, κ, ϵ, λ), 2) >=
              0.5C_Y.K_2_dξ_2 * ξ^(1 / σ - real(lambda) / κ - 2)

        #####
        ## H_OJ
        #####

        Q_hat, Q_hat_dξ = CGL2.Q_hat_infinity(γ₁, γ₂, κ, ϵ, ξ, λ)
        H = K_2(ξ, lambda, κ, ϵ, λ) * CGL2.I_N(Q_hat, λ)
        H_dξ =
            K_2_dξ(ξ, lambda, κ, ϵ, λ) * CGL2.I_N(Q_hat, λ) +
            K_2(ξ, lambda, κ, ϵ, λ) * CGL2.I_N_dξ(Q_hat, Q_hat_dξ, λ)

        @test abs(H[1, 1]) <= C_Y.H_11 * ξ^(1 / σ - real(lambda) / κ - 3)
        @test abs(H[1, 2]) <= C_Y.H_12 * ξ^(1 / σ - real(lambda) / κ - 3)
        @test abs(H[2, 1]) <= C_Y.H_21 * ξ^(1 / σ - real(lambda) / κ - 3)
        @test abs(H[2, 2]) <= C_Y.H_22 * ξ^(1 / σ - real(lambda) / κ - 3)
        @test abs(H[1, 1]) >= 0.8C_Y.H_11 * ξ^(1 / σ - real(lambda) / κ - 3)
        @test abs(H[1, 2]) >= 0.8C_Y.H_12 * ξ^(1 / σ - real(lambda) / κ - 3)
        @test abs(H[2, 1]) >= 0.7C_Y.H_21 * ξ^(1 / σ - real(lambda) / κ - 3)
        @test abs(H[2, 2]) >= 0.7C_Y.H_22 * ξ^(1 / σ - real(lambda) / κ - 3)

        @test abs(H_dξ[1, 1]) <= C_Y.H_11_dξ * ξ^(1 / σ - real(lambda) / κ - 4)
        @test abs(H_dξ[1, 2]) <= C_Y.H_12_dξ * ξ^(1 / σ - real(lambda) / κ - 4)
        @test abs(H_dξ[2, 1]) <= C_Y.H_21_dξ * ξ^(1 / σ - real(lambda) / κ - 4)
        @test abs(H_dξ[2, 2]) <= C_Y.H_22_dξ * ξ^(1 / σ - real(lambda) / κ - 4)
        # IMPROVE: These bounds are pretty bad
        @test abs(H_dξ[1, 1]) >= 0.3C_Y.H_11_dξ * ξ^(1 / σ - real(lambda) / κ - 4)
        @test abs(H_dξ[1, 2]) >= 0.5C_Y.H_12_dξ * ξ^(1 / σ - real(lambda) / κ - 4)
        @test abs(H_dξ[2, 1]) >= 0.2C_Y.H_21_dξ * ξ^(1 / σ - real(lambda) / κ - 4)
        @test abs(H_dξ[2, 2]) >= 0.3C_Y.H_22_dξ * ξ^(1 / σ - real(lambda) / κ - 4)
    end

    # TODO: Add tests for remaining constants
end
