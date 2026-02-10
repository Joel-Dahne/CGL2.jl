@testset "Y_infinity_function_bounds" begin
    lambda = Acb(0.19028080950252219, 2.769747313863597)
    γ₁ = Acb(0.1979686615344728, 0.15613711743772288)
    γ₂ = Acb(-116.03421868436834, 101.21494871936197)
    κ = Arb(0.8073018593981386)
    ϵ = Arb(0.15002213424487343)
    ξ₁ = Arb(15)
    λ = CGLParams{Arb}(3, 1.0, 1.0, 0.0)
    (; d, σ) = λ
    _, _, c = CGL2._abc(κ, ϵ, λ)
    A = SMatrix{2,2}(ϵ, 1, -1, ϵ)
    M = SMatrix{2,2}(im, 1, -im, 1)

    C_Y = CGL2.FunctionBounds_Y(lambda, γ₁, γ₂, κ, ϵ, ξ₁, λ, include_dλ = true)

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

    E_1_dλ, E_1_dλ_dξ = CGL2.E_1_dλ, CGL2.E_1_dλ_dξ
    E_2_dλ, E_2_dλ_dξ = CGL2.E_2_dλ, CGL2.E_2_dλ_dξ
    P_1_dλ, P_1_dλ_dξ = CGL2.P_1_dλ, CGL2.P_1_dλ_dξ
    P_2_dλ, P_2_dλ_dξ = CGL2.P_2_dλ, CGL2.P_2_dλ_dξ

    J_E_1 = CGL2.J_E_1
    J_E_2 = CGL2.J_E_2
    J_P_1 = CGL2.J_P_1
    J_P_2 = CGL2.J_P_2

    J_E_1_dξ = CGL2.J_E_1_dξ
    J_E_2_dξ = CGL2.J_E_2_dξ

    J_E_1_dλ = CGL2.J_E_1_dλ
    J_E_2_dλ = CGL2.J_E_2_dλ
    J_P_1_dλ = CGL2.J_P_1_dλ
    J_P_2_dλ = CGL2.J_P_2_dλ

    K_1_2 = CGL2.K_1_2
    K_1_2_dξ = CGL2.K_1_2_dξ
    K_1_2_dλ = CGL2.K_1_2_dλ

    for ξ in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64] .* ξ₁
        ####
        ## J_N
        ####

        # TODO: This is not quite correct since γ₁ and γ₂ only
        # parametrize Q_hat_infinity at ξ₁ and not at ξ > ξ₁. The
        # difference should however be of lower order and it still
        # gives some information as a test.
        Q_hat, Q_hat_dξ = CGL2.Q_hat_infinity(γ₁, γ₂, κ, ϵ, ξ, λ)

        # IMPROVE: We could maybe improve on this bound?
        @test norm_inf(CGL2.J_N(Q_hat, λ)) <= C_Y.J_N * ξ^-2
        @test norm_inf(CGL2.J_N(Q_hat, λ)) >= 0.6C_Y.J_N * ξ^-2
        @test norm_inf(CGL2.J_N_dξ(Q_hat, Q_hat_dξ, λ)) <= C_Y.J_N_dξ * ξ^-3
        @test norm_inf(CGL2.J_N_dξ(Q_hat, Q_hat_dξ, λ)) >= 0.6C_Y.J_N_dξ * ξ^-3

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

        @test abs(E_1_dλ(ξ, lambda, κ, ϵ, λ)) <=
              C_Y.E_1_dλ * exp(-real(c) * ξ^2) * log(ξ) * ξ^(1 / σ - d - real(lambda) / κ)
        @test abs(E_2_dλ(ξ, lambda, κ, ϵ, λ)) <=
              C_Y.E_2_dλ * exp(-real(c) * ξ^2) * log(ξ) * ξ^(1 / σ - d - real(lambda) / κ)
        @test abs(E_1_dλ(ξ, lambda, κ, ϵ, λ)) >=
              0.6C_Y.E_1_dλ *
              exp(-real(c) * ξ^2) *
              log(ξ) *
              ξ^(1 / σ - d - real(lambda) / κ)
        @test abs(E_2_dλ(ξ, lambda, κ, ϵ, λ)) >=
              0.55C_Y.E_2_dλ *
              exp(-real(c) * ξ^2) *
              log(ξ) *
              ξ^(1 / σ - d - real(lambda) / κ)

        @test abs(E_1_dλ_dξ(ξ, lambda, κ, ϵ, λ)) <=
              C_Y.E_1_dλ_dξ *
              exp(-real(c) * ξ^2) *
              log(ξ) *
              ξ^(1 / σ - d - real(lambda) / κ + 1)
        @test abs(E_2_dλ_dξ(ξ, lambda, κ, ϵ, λ)) <=
              C_Y.E_2_dλ_dξ *
              exp(-real(c) * ξ^2) *
              log(ξ) *
              ξ^(1 / σ - d - real(lambda) / κ + 1)
        @test abs(E_1_dλ_dξ(ξ, lambda, κ, ϵ, λ)) >=
              0.6C_Y.E_1_dλ_dξ *
              exp(-real(c) * ξ^2) *
              log(ξ) *
              ξ^(1 / σ - d - real(lambda) / κ + 1)
        @test abs(E_2_dλ_dξ(ξ, lambda, κ, ϵ, λ)) >=
              0.5C_Y.E_2_dλ_dξ *
              exp(-real(c) * ξ^2) *
              log(ξ) *
              ξ^(1 / σ - d - real(lambda) / κ + 1)

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

        @test abs(P_1_dλ(ξ, lambda, κ, ϵ, λ)) <=
              C_Y.P_1_dλ * log(ξ) * ξ^(-1 / σ + real(lambda) / κ)
        @test abs(P_2_dλ(ξ, lambda, κ, ϵ, λ)) <=
              C_Y.P_2_dλ * log(ξ) * ξ^(-1 / σ + real(lambda) / κ)
        @test abs(P_1_dλ(ξ, lambda, κ, ϵ, λ)) >=
              0.6C_Y.P_1_dλ * log(ξ) * ξ^(-1 / σ + real(lambda) / κ)
        @test abs(P_2_dλ(ξ, lambda, κ, ϵ, λ)) >=
              0.6C_Y.P_2_dλ * log(ξ) * ξ^(-1 / σ + real(lambda) / κ)

        @test abs(P_1_dλ_dξ(ξ, lambda, κ, ϵ, λ)) <=
              C_Y.P_1_dλ_dξ * log(ξ) * ξ^(-1 / σ + real(lambda) / κ - 1)
        @test abs(P_2_dλ_dξ(ξ, lambda, κ, ϵ, λ)) <=
              C_Y.P_2_dλ_dξ * log(ξ) * ξ^(-1 / σ + real(lambda) / κ - 1)
        @test abs(P_1_dλ_dξ(ξ, lambda, κ, ϵ, λ)) >=
              0.4C_Y.P_1_dλ_dξ * log(ξ) * ξ^(-1 / σ + real(lambda) / κ - 1)
        @test abs(P_2_dλ_dξ(ξ, lambda, κ, ϵ, λ)) >=
              0.4C_Y.P_2_dλ_dξ * log(ξ) * ξ^(-1 / σ + real(lambda) / κ - 1)

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
              0.4C_Y.J_E_1_dξ * ξ^(1 / σ - real(lambda) / κ - 2)

        @test abs(J_E_2_dξ(ξ, lambda, κ, ϵ, λ)) <=
              C_Y.J_E_2_dξ * ξ^(1 / σ - real(lambda) / κ - 2)
        @test abs(J_E_2_dξ(ξ, lambda, κ, ϵ, λ)) >=
              0.5C_Y.J_E_2_dξ * ξ^(1 / σ - real(lambda) / κ - 2)

        @test abs(J_E_1_dλ(ξ, lambda, κ, ϵ, λ)) <=
              C_Y.J_E_1_dλ * log(ξ) * ξ^(1 / σ - real(lambda) / κ - 1)
        @test abs(J_E_1_dλ(ξ, lambda, κ, ϵ, λ)) >=
              0.4C_Y.J_E_1_dλ * log(ξ) * ξ^(1 / σ - real(lambda) / κ - 1)

        @test abs(J_E_2_dλ(ξ, lambda, κ, ϵ, λ)) <=
              C_Y.J_E_2_dλ * log(ξ) * ξ^(1 / σ - real(lambda) / κ - 1)
        @test abs(J_E_2_dλ(ξ, lambda, κ, ϵ, λ)) >=
              0.4C_Y.J_E_2_dλ * log(ξ) * ξ^(1 / σ - real(lambda) / κ - 1)

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

        @test abs(J_P_1_dλ(ξ, lambda, κ, ϵ, λ)) <=
              C_Y.J_P_1_dλ *
              exp(real(c) * ξ^2) *
              log(ξ) *
              ξ^(- 1 / σ + d + real(lambda) / κ - 1)
        @test abs(J_P_1_dλ(ξ, lambda, κ, ϵ, λ)) >=
              0.4C_Y.J_P_1_dλ *
              exp(real(c) * ξ^2) *
              log(ξ) *
              ξ^(- 1 / σ + d + real(lambda) / κ - 1)

        @test abs(J_P_2_dλ(ξ, lambda, κ, ϵ, λ)) <=
              C_Y.J_P_2_dλ *
              exp(real(c) * ξ^2) *
              log(ξ) *
              ξ^(- 1 / σ + d + real(lambda) / κ - 1)
        @test abs(J_P_2_dλ(ξ, lambda, κ, ϵ, λ)) >=
              0.4C_Y.J_P_2_dλ *
              exp(real(c) * ξ^2) *
              log(ξ) *
              ξ^(- 1 / σ + d + real(lambda) / κ - 1)

        ######
        ## K_1, K_2
        ######

        K_1, K_2 = K_1_2(ξ, lambda, κ, ϵ, λ)
        @test norm_inf(K_1, 1) <=
              C_Y.K_1_1 * exp(real(c) * ξ^2) * ξ^(- 1 / σ + d + real(lambda) / κ - 1)
        @test norm_inf(K_1, 1) >=
              0.95C_Y.K_1_1 * exp(real(c) * ξ^2) * ξ^(- 1 / σ + d + real(lambda) / κ - 1)
        @test norm_inf(K_1, 2) <=
              C_Y.K_1_2 * exp(real(c) * ξ^2) * ξ^(- 1 / σ + d + real(lambda) / κ - 1)
        @test norm_inf(K_1, 2) >=
              0.9C_Y.K_1_2 * exp(real(c) * ξ^2) * ξ^(- 1 / σ + d + real(lambda) / κ - 1)

        @test norm_inf(K_2, 1) <= C_Y.K_2_1 * ξ^(1 / σ - real(lambda) / κ - 1)
        @test norm_inf(K_2, 1) >= 0.9C_Y.K_2_1 * ξ^(1 / σ - real(lambda) / κ - 1)
        @test norm_inf(K_2, 2) <= C_Y.K_2_2 * ξ^(1 / σ - real(lambda) / κ - 1)
        @test norm_inf(K_2, 2) >= 0.85C_Y.K_2_2 * ξ^(1 / σ - real(lambda) / κ - 1)

        K_1_dλ, K_2_dλ = K_1_2_dλ(ξ, lambda, κ, ϵ, λ)
        @test norm_inf(K_1_dλ, 1) <=
              C_Y.K_1_dλ_1 *
              exp(real(c) * ξ^2) *
              log(ξ) *
              ξ^(- 1 / σ + d + real(lambda) / κ - 1)
        @test norm_inf(K_1_dλ, 1) >=
              0.4C_Y.K_1_dλ_1 *
              exp(real(c) * ξ^2) *
              log(ξ) *
              ξ^(- 1 / σ + d + real(lambda) / κ - 1)
        @test norm_inf(K_1_dλ, 2) <=
              C_Y.K_1_dλ_2 *
              exp(real(c) * ξ^2) *
              log(ξ) *
              ξ^(- 1 / σ + d + real(lambda) / κ - 1)
        @test norm_inf(K_1_dλ, 2) >=
              0.4C_Y.K_1_dλ_2 *
              exp(real(c) * ξ^2) *
              log(ξ) *
              ξ^(- 1 / σ + d + real(lambda) / κ - 1)

        @test norm_inf(K_2_dλ, 1) <=
              C_Y.K_2_dλ_1 * log(ξ) * ξ^(1 / σ - real(lambda) / κ - 1)
        @test norm_inf(K_2_dλ, 1) >=
              0.4C_Y.K_2_dλ_1 * log(ξ) * ξ^(1 / σ - real(lambda) / κ - 1)
        @test norm_inf(K_2_dλ, 2) <=
              C_Y.K_2_dλ_2 * log(ξ) * ξ^(1 / σ - real(lambda) / κ - 1)
        @test norm_inf(K_2_dλ, 2) >=
              0.4C_Y.K_2_dλ_2 * log(ξ) * ξ^(1 / σ - real(lambda) / κ - 1)

        #####
        ## H_OJ
        #####

        Q_hat, Q_hat_dξ = CGL2.Q_hat_infinity(γ₁, γ₂, κ, ϵ, ξ, λ)
        H = K_1_2(ξ, lambda, κ, ϵ, λ)[2] * CGL2.J_N(Q_hat, λ) * M
        H_dξ =
            K_1_2_dξ(ξ, lambda, κ, ϵ, λ)[2] * CGL2.J_N(Q_hat, λ) * M +
            K_1_2(ξ, lambda, κ, ϵ, λ)[2] * CGL2.J_N_dξ(Q_hat, Q_hat_dξ, λ) * M

        @test abs(H[1, 1]) <= C_Y.H_1j * ξ^(1 / σ - real(lambda) / κ - 3)
        @test abs(H[1, 2]) <= C_Y.H_1j * ξ^(1 / σ - real(lambda) / κ - 3)
        @test abs(H[2, 1]) <= C_Y.H_2j * ξ^(1 / σ - real(lambda) / κ - 3)
        @test abs(H[2, 2]) <= C_Y.H_2j * ξ^(1 / σ - real(lambda) / κ - 3)
        # IMPROVE: These bounds are very bad... It would be nice to improve them
        @test abs(H[1, 1]) >= 0.2C_Y.H_1j * ξ^(1 / σ - real(lambda) / κ - 3)
        @test abs(H[1, 2]) >= 0.1C_Y.H_1j * ξ^(1 / σ - real(lambda) / κ - 3)
        @test abs(H[2, 1]) >= 0.075C_Y.H_2j * ξ^(1 / σ - real(lambda) / κ - 3)
        @test abs(H[2, 2]) >= 0.1C_Y.H_2j * ξ^(1 / σ - real(lambda) / κ - 3)

        @test abs(H_dξ[1, 1]) <= C_Y.H_1j_dξ * ξ^(1 / σ - real(lambda) / κ - 4)
        @test abs(H_dξ[1, 2]) <= C_Y.H_1j_dξ * ξ^(1 / σ - real(lambda) / κ - 4)
        @test abs(H_dξ[2, 1]) <= C_Y.H_2j_dξ * ξ^(1 / σ - real(lambda) / κ - 4)
        @test abs(H_dξ[2, 2]) <= C_Y.H_2j_dξ * ξ^(1 / σ - real(lambda) / κ - 4)
        # IMPROVE: These bounds are very bad... It would be nice to improve them
        @test abs(H_dξ[1, 1]) >= 0.075C_Y.H_1j_dξ * ξ^(1 / σ - real(lambda) / κ - 4)
        @test abs(H_dξ[1, 2]) >= 0.05C_Y.H_1j_dξ * ξ^(1 / σ - real(lambda) / κ - 4)
        @test abs(H_dξ[2, 1]) >= 0.02C_Y.H_2j_dξ * ξ^(1 / σ - real(lambda) / κ - 4)
        @test abs(H_dξ[2, 2]) >= 0.075C_Y.H_2j_dξ * ξ^(1 / σ - real(lambda) / κ - 4)
    end
end
