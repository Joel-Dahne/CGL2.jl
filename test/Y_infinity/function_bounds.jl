@testset "Y_infinity_function_bounds" begin
    lambda = Acb(0.19028080950252219, 2.769747313863597)
    ν = Acb(2.1134119144964454, 2.5294091893480206)
    γ₁ = Acb(0.1979686615344728, 0.15613711743772288)
    γ₂ = Acb(-116.03421868436834, 101.21494871936197)
    κ = Arb(0.8073018593981386)
    ϵ = Arb(0.15002213424487343)
    ξ₁ = Arb(15)
    λ = CGLParams{Arb}(3, 1.0, 1.0, 0.0)
    (; d, σ) = λ
    _, _, c = CGL2._abc(κ, ϵ, λ)

    real_a1 = CGL2.real_a12(κ, ϵ)
    real_s1 = CGL2.real_s12(lambda, κ, λ)
    real_s3 = CGL2.real_s34(lambda, κ, λ)

    C_Y = CGL2.FunctionBounds_Y(lambda, γ₁, γ₂, κ, ϵ, ξ₁, λ, include_dλ = true)
    C_Y_new = CGL2.FunctionBounds_Y_new(lambda, γ₁, γ₂, κ, ϵ, ξ₁, λ, include_dλ = true)

    norm_inf = CGL2.norm_inf

    Y_1, Y_1_dξ = CGL2.Y_1, CGL2.Y_1_dξ
    Y_2, Y_2_dξ = CGL2.Y_2, CGL2.Y_2_dξ
    Y_3, Y_3_dξ = CGL2.Y_3, CGL2.Y_3_dξ
    Y_4, Y_4_dξ = CGL2.Y_4, CGL2.Y_4_dξ

    Y_1_dλ, Y_1_dλ_dξ = CGL2.Y_1_dλ, CGL2.Y_1_dλ_dξ
    Y_2_dλ, Y_2_dλ_dξ = CGL2.Y_2_dλ, CGL2.Y_2_dλ_dξ
    Y_3_dλ, Y_3_dλ_dξ = CGL2.Y_3_dλ, CGL2.Y_3_dλ_dξ
    Y_4_dλ, Y_4_dλ_dξ = CGL2.Y_4_dλ, CGL2.Y_4_dλ_dξ

    E_1, E_1_dξ = CGL2.E_1, CGL2.E_1_dξ
    E_2, E_2_dξ = CGL2.E_2, CGL2.E_2_dξ
    P_1, P_1_dξ = CGL2.P_1, CGL2.P_1_dξ
    P_2, P_2_dξ = CGL2.P_2, CGL2.P_2_dξ

    E_1_dλ, E_1_dλ_dξ = CGL2.E_1_dλ, CGL2.E_1_dλ_dξ
    E_2_dλ, E_2_dλ_dξ = CGL2.E_2_dλ, CGL2.E_2_dλ_dξ
    P_1_dλ, P_1_dλ_dξ = CGL2.P_1_dλ, CGL2.P_1_dλ_dξ
    P_2_dλ, P_2_dλ_dξ = CGL2.P_2_dλ, CGL2.P_2_dλ_dξ

    J_E_1 = CGL2.J_E_1
    J_E_2 = CGL2.J_E_2
    J_P_1 = CGL2.J_P_1
    J_P_2 = CGL2.J_P_2

    J_E_1_dλ = CGL2.J_E_1_dλ
    J_E_2_dλ = CGL2.J_E_2_dλ
    J_P_1_dλ = CGL2.J_P_1_dλ
    J_P_2_dλ = CGL2.J_P_2_dλ

    K_1_2 = CGL2.K_1_2
    K_1_2_dλ = CGL2.K_1_2_dλ
    K_1_2_new = CGL2.K_1_2_new
    K_1_2_dλ_new = CGL2.K_1_2_dλ_new

    for ξ in [1, 1.01, 1.1, 2, 4, 8, 16, 32, 64] .* ξ₁
        ####
        ## J_N
        ####

        # TODO: This is not quite correct since γ₁ and γ₂ only
        # parametrize Q_hat_infinity at ξ₁ and not at ξ > ξ₁. The
        # difference should however be of lower order and it still
        # gives some information as a test.
        Q_hat, _ = CGL2.Q_hat_infinity(γ₁, γ₂, κ, ϵ, ξ, λ)
        @test norm_inf(CGL2.J_N(Q_hat, λ)) <= C_Y.J_N * ξ^-2
        @test norm_inf(CGL2.J_N(Q_hat, λ)) >= 0.5C_Y.J_N * ξ^-2

        ####
        ## Y_12
        ####
        @test norm_inf(Y_1(ξ, lambda, κ, ϵ, λ)) <=
              C_Y.Y_12 * exp(real_a1 * ξ^2) * ξ^-real_s1
        @test norm_inf(Y_2(ξ, lambda, κ, ϵ, λ)) <=
              C_Y.Y_12 * exp(real_a1 * ξ^2) * ξ^-real_s1
        @test norm_inf(Y_1(ξ, lambda, κ, ϵ, λ)) >=
              0.5C_Y.Y_12 * exp(real_a1 * ξ^2) * ξ^-real_s1
        @test norm_inf(Y_2(ξ, lambda, κ, ϵ, λ)) >=
              0.5C_Y.Y_12 * exp(real_a1 * ξ^2) * ξ^-real_s1

        @test norm_inf(Y_1_dξ(ξ, lambda, κ, ϵ, λ)) <=
              C_Y.Y_12_dξ * exp(real_a1 * ξ^2) * ξ^(-real_s1 + 1)
        @test norm_inf(Y_2_dξ(ξ, lambda, κ, ϵ, λ)) <=
              C_Y.Y_12_dξ * exp(real_a1 * ξ^2) * ξ^(-real_s1 + 1)
        @test norm_inf(Y_1_dξ(ξ, lambda, κ, ϵ, λ)) >=
              0.5C_Y.Y_12_dξ * exp(real_a1 * ξ^2) * ξ^(-real_s1 + 1)
        @test norm_inf(Y_2_dξ(ξ, lambda, κ, ϵ, λ)) >=
              0.5C_Y.Y_12_dξ * exp(real_a1 * ξ^2) * ξ^(-real_s1 + 1)

        @test norm_inf(Y_1_dλ(ξ, lambda, κ, ϵ, λ)) <=
              C_Y.Y_12_dλ * exp(real_a1 * ξ^2) * log(ξ) * ξ^-real_s1
        @test norm_inf(Y_2_dλ(ξ, lambda, κ, ϵ, λ)) <=
              C_Y.Y_12_dλ * exp(real_a1 * ξ^2) * log(ξ) * ξ^-real_s1
        @test norm_inf(Y_1_dλ(ξ, lambda, κ, ϵ, λ)) >=
              0.5C_Y.Y_12_dλ * exp(real_a1 * ξ^2) * log(ξ) * ξ^-real_s1
        @test norm_inf(Y_2_dλ(ξ, lambda, κ, ϵ, λ)) >=
              0.5C_Y.Y_12_dλ * exp(real_a1 * ξ^2) * log(ξ) * ξ^-real_s1

        @test norm_inf(Y_1_dλ_dξ(ξ, lambda, κ, ϵ, λ)) <=
              C_Y.Y_12_dλ_dξ * exp(real_a1 * ξ^2) * log(ξ) * ξ^(-real_s1 + 1)
        @test norm_inf(Y_2_dλ_dξ(ξ, lambda, κ, ϵ, λ)) <=
              C_Y.Y_12_dλ_dξ * exp(real_a1 * ξ^2) * log(ξ) * ξ^(-real_s1 + 1)
        @test norm_inf(Y_1_dλ_dξ(ξ, lambda, κ, ϵ, λ)) >=
              0.5C_Y.Y_12_dλ_dξ * exp(real_a1 * ξ^2) * log(ξ) * ξ^(-real_s1 + 1)
        @test norm_inf(Y_2_dλ_dξ(ξ, lambda, κ, ϵ, λ)) >=
              0.5C_Y.Y_12_dλ_dξ * exp(real_a1 * ξ^2) * log(ξ) * ξ^(-real_s1 + 1)

        ####
        ## Y_34
        ####
        @test norm_inf(Y_3(ξ, lambda, κ, ϵ, λ)) <= C_Y.Y_34 * ξ^-real_s3
        @test norm_inf(Y_4(ξ, lambda, κ, ϵ, λ)) <= C_Y.Y_34 * ξ^-real_s3
        @test norm_inf(Y_3(ξ, lambda, κ, ϵ, λ)) >= 0.5C_Y.Y_34 * ξ^-real_s3
        @test norm_inf(Y_4(ξ, lambda, κ, ϵ, λ)) >= 0.5C_Y.Y_34 * ξ^-real_s3

        @test norm_inf(Y_3_dξ(ξ, lambda, κ, ϵ, λ)) <= C_Y.Y_34_dξ * ξ^(-real_s3 - 1)
        @test norm_inf(Y_4_dξ(ξ, lambda, κ, ϵ, λ)) <= C_Y.Y_34_dξ * ξ^(-real_s3 - 1)
        @test norm_inf(Y_3_dξ(ξ, lambda, κ, ϵ, λ)) >= 0.5C_Y.Y_34_dξ * ξ^(-real_s3 - 1)
        # IMPROVE: This one gives slightly worse bounds
        @test norm_inf(Y_4_dξ(ξ, lambda, κ, ϵ, λ)) >= 0.25C_Y.Y_34_dξ * ξ^(-real_s3 - 1)

        @test norm_inf(Y_3_dλ(ξ, lambda, κ, ϵ, λ)) <= C_Y.Y_34_dλ * log(ξ) * ξ^-real_s3
        @test norm_inf(Y_4_dλ(ξ, lambda, κ, ϵ, λ)) <= C_Y.Y_34_dλ * log(ξ) * ξ^-real_s3
        @test norm_inf(Y_3_dλ(ξ, lambda, κ, ϵ, λ)) >= 0.5C_Y.Y_34_dλ * log(ξ) * ξ^-real_s3
        @test norm_inf(Y_4_dλ(ξ, lambda, κ, ϵ, λ)) >= 0.5C_Y.Y_34_dλ * log(ξ) * ξ^-real_s3

        @test norm_inf(Y_3_dλ_dξ(ξ, lambda, κ, ϵ, λ)) <=
              C_Y.Y_34_dλ_dξ * log(ξ) * ξ^(-real_s3 - 1)
        @test norm_inf(Y_4_dλ_dξ(ξ, lambda, κ, ϵ, λ)) <=
              C_Y.Y_34_dλ_dξ * log(ξ) * ξ^(-real_s3 - 1)
        @test norm_inf(Y_3_dλ_dξ(ξ, lambda, κ, ϵ, λ)) >=
              0.5C_Y.Y_34_dλ_dξ * log(ξ) * ξ^(-real_s3 - 1)
        # IMPROVE: This one gives slightly worse bounds
        @test norm_inf(Y_4_dλ_dξ(ξ, lambda, κ, ϵ, λ)) >=
              0.25C_Y.Y_34_dλ_dξ * log(ξ) * ξ^(-real_s3 - 1)

        ####
        ## E_1 and E_2
        ####

        @test abs(E_1(ξ, lambda, κ, ϵ, λ)) <=
              C_Y_new.E_1 * exp(-real(c) * ξ^2) * ξ^(1 / σ - d - real(lambda) / κ)
        @test abs(E_2(ξ, lambda, κ, ϵ, λ)) <=
              C_Y_new.E_2 * exp(-real(c) * ξ^2) * ξ^(1 / σ - d - real(lambda) / κ)
        @test abs(E_1(ξ, lambda, κ, ϵ, λ)) >=
              0.95C_Y_new.E_1 * exp(-real(c) * ξ^2) * ξ^(1 / σ - d - real(lambda) / κ)
        @test abs(E_2(ξ, lambda, κ, ϵ, λ)) >=
              0.85C_Y_new.E_2 * exp(-real(c) * ξ^2) * ξ^(1 / σ - d - real(lambda) / κ)

        @test abs(E_1_dξ(ξ, lambda, κ, ϵ, λ)) <=
              C_Y_new.E_1_dξ * exp(-real(c) * ξ^2) * ξ^(1 / σ - d - real(lambda) / κ + 1)
        @test abs(E_2_dξ(ξ, lambda, κ, ϵ, λ)) <=
              C_Y_new.E_2_dξ * exp(-real(c) * ξ^2) * ξ^(1 / σ - d - real(lambda) / κ + 1)
        @test abs(E_1_dξ(ξ, lambda, κ, ϵ, λ)) >=
              0.95C_Y_new.E_1_dξ *
              exp(-real(c) * ξ^2) *
              ξ^(1 / σ - d - real(lambda) / κ + 1)
        @test abs(E_2_dξ(ξ, lambda, κ, ϵ, λ)) >=
              0.85C_Y_new.E_2_dξ *
              exp(-real(c) * ξ^2) *
              ξ^(1 / σ - d - real(lambda) / κ + 1)

        @test abs(E_1_dλ(ξ, lambda, κ, ϵ, λ)) <=
              C_Y_new.E_1_dλ *
              exp(-real(c) * ξ^2) *
              log(ξ) *
              ξ^(1 / σ - d - real(lambda) / κ)
        @test abs(E_2_dλ(ξ, lambda, κ, ϵ, λ)) <=
              C_Y_new.E_2_dλ *
              exp(-real(c) * ξ^2) *
              log(ξ) *
              ξ^(1 / σ - d - real(lambda) / κ)
        @test abs(E_1_dλ(ξ, lambda, κ, ϵ, λ)) >=
              0.6C_Y_new.E_1_dλ *
              exp(-real(c) * ξ^2) *
              log(ξ) *
              ξ^(1 / σ - d - real(lambda) / κ)
        @test abs(E_2_dλ(ξ, lambda, κ, ϵ, λ)) >=
              0.55C_Y_new.E_2_dλ *
              exp(-real(c) * ξ^2) *
              log(ξ) *
              ξ^(1 / σ - d - real(lambda) / κ)

        @test abs(E_1_dλ_dξ(ξ, lambda, κ, ϵ, λ)) <=
              C_Y_new.E_1_dλ_dξ *
              exp(-real(c) * ξ^2) *
              log(ξ) *
              ξ^(1 / σ - d - real(lambda) / κ + 1)
        @test abs(E_2_dλ_dξ(ξ, lambda, κ, ϵ, λ)) <=
              C_Y_new.E_2_dλ_dξ *
              exp(-real(c) * ξ^2) *
              log(ξ) *
              ξ^(1 / σ - d - real(lambda) / κ + 1)
        @test abs(E_1_dλ_dξ(ξ, lambda, κ, ϵ, λ)) >=
              0.6C_Y_new.E_1_dλ_dξ *
              exp(-real(c) * ξ^2) *
              log(ξ) *
              ξ^(1 / σ - d - real(lambda) / κ + 1)
        @test abs(E_2_dλ_dξ(ξ, lambda, κ, ϵ, λ)) >=
              0.5C_Y_new.E_2_dλ_dξ *
              exp(-real(c) * ξ^2) *
              log(ξ) *
              ξ^(1 / σ - d - real(lambda) / κ + 1)

        ####
        ## E_12
        ####

        @test norm_inf(hcat(E_1(ξ, lambda, κ, ϵ, λ), E_2(ξ, lambda, κ, ϵ, λ))) <=
              C_Y_new.E_12 * exp(-real(c) * ξ^2) * ξ^(1 / σ - d - real(lambda) / κ)
        @test norm_inf(hcat(E_1(ξ, lambda, κ, ϵ, λ), E_2(ξ, lambda, κ, ϵ, λ))) >=
              0.9C_Y_new.E_12 * exp(-real(c) * ξ^2) * ξ^(1 / σ - d - real(lambda) / κ)

        @test norm_inf(hcat(E_1_dξ(ξ, lambda, κ, ϵ, λ), E_2_dξ(ξ, lambda, κ, ϵ, λ))) <=
              C_Y_new.E_12_dξ * exp(-real(c) * ξ^2) * ξ^(1 / σ - d - real(lambda) / κ + 1)
        @test norm_inf(hcat(E_1_dξ(ξ, lambda, κ, ϵ, λ), E_2_dξ(ξ, lambda, κ, ϵ, λ))) >=
              0.8C_Y_new.E_12_dξ *
              exp(-real(c) * ξ^2) *
              ξ^(1 / σ - d - real(lambda) / κ + 1)

        @test norm_inf(hcat(E_1_dλ(ξ, lambda, κ, ϵ, λ), E_2_dλ(ξ, lambda, κ, ϵ, λ))) <=
              C_Y_new.E_12_dλ *
              exp(-real(c) * ξ^2) *
              log(ξ) *
              ξ^(1 / σ - d - real(lambda) / κ)
        @test norm_inf(hcat(E_1_dλ(ξ, lambda, κ, ϵ, λ), E_2_dλ(ξ, lambda, κ, ϵ, λ))) >=
              0.5C_Y_new.E_12_dλ *
              exp(-real(c) * ξ^2) *
              log(ξ) *
              ξ^(1 / σ - d - real(lambda) / κ)

        @test norm_inf(
            hcat(E_1_dλ_dξ(ξ, lambda, κ, ϵ, λ), E_2_dλ_dξ(ξ, lambda, κ, ϵ, λ)),
        ) <=
              C_Y_new.E_12_dλ_dξ *
              exp(-real(c) * ξ^2) *
              log(ξ) *
              ξ^(1 / σ - d - real(lambda) / κ + 1)
        @test norm_inf(
            hcat(E_1_dλ_dξ(ξ, lambda, κ, ϵ, λ), E_2_dλ_dξ(ξ, lambda, κ, ϵ, λ)),
        ) >=
              0.6C_Y_new.E_12_dλ_dξ *
              exp(-real(c) * ξ^2) *
              log(ξ) *
              ξ^(1 / σ - d - real(lambda) / κ + 1)

        ####
        ## P_1 and P_2
        ####

        @test abs(P_1(ξ, lambda, κ, ϵ, λ)) <= C_Y_new.P_1 * ξ^(-1 / σ + real(lambda) / κ)
        @test abs(P_2(ξ, lambda, κ, ϵ, λ)) <= C_Y_new.P_2 * ξ^(-1 / σ + real(lambda) / κ)
        @test abs(P_1(ξ, lambda, κ, ϵ, λ)) >=
              0.95C_Y_new.P_1 * ξ^(-1 / σ + real(lambda) / κ)
        @test abs(P_2(ξ, lambda, κ, ϵ, λ)) >=
              0.90C_Y_new.P_2 * ξ^(-1 / σ + real(lambda) / κ)

        @test abs(P_1_dξ(ξ, lambda, κ, ϵ, λ)) <=
              C_Y_new.P_1_dξ * ξ^(-1 / σ + real(lambda) / κ - 1)
        @test abs(P_2_dξ(ξ, lambda, κ, ϵ, λ)) <=
              C_Y_new.P_2_dξ * ξ^(-1 / σ + real(lambda) / κ - 1)
        @test abs(P_1_dξ(ξ, lambda, κ, ϵ, λ)) >=
              0.95C_Y_new.P_1_dξ * ξ^(-1 / σ + real(lambda) / κ - 1)
        @test abs(P_2_dξ(ξ, lambda, κ, ϵ, λ)) >=
              0.85C_Y_new.P_2_dξ * ξ^(-1 / σ + real(lambda) / κ - 1)

        @test abs(P_1_dλ(ξ, lambda, κ, ϵ, λ)) <=
              C_Y_new.P_1_dλ * log(ξ) * ξ^(-1 / σ + real(lambda) / κ)
        @test abs(P_2_dλ(ξ, lambda, κ, ϵ, λ)) <=
              C_Y_new.P_2_dλ * log(ξ) * ξ^(-1 / σ + real(lambda) / κ)
        @test abs(P_1_dλ(ξ, lambda, κ, ϵ, λ)) >=
              0.6C_Y_new.P_1_dλ * log(ξ) * ξ^(-1 / σ + real(lambda) / κ)
        @test abs(P_2_dλ(ξ, lambda, κ, ϵ, λ)) >=
              0.6C_Y_new.P_2_dλ * log(ξ) * ξ^(-1 / σ + real(lambda) / κ)

        @test abs(P_1_dλ_dξ(ξ, lambda, κ, ϵ, λ)) <=
              C_Y_new.P_1_dλ_dξ * log(ξ) * ξ^(-1 / σ + real(lambda) / κ - 1)
        @test abs(P_2_dλ_dξ(ξ, lambda, κ, ϵ, λ)) <=
              C_Y_new.P_2_dλ_dξ * log(ξ) * ξ^(-1 / σ + real(lambda) / κ - 1)
        @test abs(P_1_dλ_dξ(ξ, lambda, κ, ϵ, λ)) >=
              0.4C_Y_new.P_1_dλ_dξ * log(ξ) * ξ^(-1 / σ + real(lambda) / κ - 1)
        @test abs(P_2_dλ_dξ(ξ, lambda, κ, ϵ, λ)) >=
              0.4C_Y_new.P_2_dλ_dξ * log(ξ) * ξ^(-1 / σ + real(lambda) / κ - 1)

        ####
        ## P_12
        ####

        @test norm_inf(hcat(P_1(ξ, lambda, κ, ϵ, λ), P_2(ξ, lambda, κ, ϵ, λ))) <=
              C_Y_new.P_12 * ξ^(-1 / σ + real(lambda) / κ)
        @test norm_inf(hcat(P_1(ξ, lambda, κ, ϵ, λ), P_2(ξ, lambda, κ, ϵ, λ))) >=
              0.9C_Y_new.P_12 * ξ^(-1 / σ + real(lambda) / κ)

        @test norm_inf(hcat(P_1_dξ(ξ, lambda, κ, ϵ, λ), P_2_dξ(ξ, lambda, κ, ϵ, λ))) <=
              C_Y_new.P_12_dξ * ξ^(-1 / σ + real(lambda) / κ - 1)
        @test norm_inf(hcat(P_1_dξ(ξ, lambda, κ, ϵ, λ), P_2_dξ(ξ, lambda, κ, ϵ, λ))) >=
              0.8C_Y_new.P_12_dξ * ξ^(-1 / σ + real(lambda) / κ - 1)

        @test norm_inf(hcat(P_1_dλ(ξ, lambda, κ, ϵ, λ), P_2_dλ(ξ, lambda, κ, ϵ, λ))) <=
              C_Y_new.P_12_dλ * log(ξ) * ξ^(-1 / σ + real(lambda) / κ)
        @test norm_inf(hcat(P_1_dλ(ξ, lambda, κ, ϵ, λ), P_2_dλ(ξ, lambda, κ, ϵ, λ))) >=
              0.6C_Y_new.P_12_dλ * log(ξ) * ξ^(-1 / σ + real(lambda) / κ)

        @test norm_inf(
            hcat(P_1_dλ_dξ(ξ, lambda, κ, ϵ, λ), P_2_dλ_dξ(ξ, lambda, κ, ϵ, λ)),
        ) <= C_Y_new.P_12_dλ_dξ * log(ξ) * ξ^(-1 / σ + real(lambda) / κ - 1)
        @test norm_inf(
            hcat(P_1_dλ_dξ(ξ, lambda, κ, ϵ, λ), P_2_dλ_dξ(ξ, lambda, κ, ϵ, λ)),
        ) >= 0.5C_Y_new.P_12_dλ_dξ * log(ξ) * ξ^(-1 / σ + real(lambda) / κ - 1)

        ######
        ## J_E_1 and J_E_2
        ######

        @test abs(J_E_1(ξ, lambda, κ, ϵ, λ)) <=
              C_Y_new.J_E_1 * ξ^(1 / σ - real(lambda) / κ - 1)
        @test abs(J_E_1(ξ, lambda, κ, ϵ, λ)) >=
              0.95C_Y_new.J_E_1 * ξ^(1 / σ - real(lambda) / κ - 1)

        @test abs(J_E_2(ξ, lambda, κ, ϵ, λ)) <=
              C_Y_new.J_E_2 * ξ^(1 / σ - real(lambda) / κ - 1)
        @test abs(J_E_2(ξ, lambda, κ, ϵ, λ)) >=
              0.85C_Y_new.J_E_2 * ξ^(1 / σ - real(lambda) / κ - 1)

        @test abs(J_E_1_dλ(ξ, lambda, κ, ϵ, λ)) <=
              C_Y_new.J_E_1_dλ * log(ξ) * ξ^(1 / σ - real(lambda) / κ - 1)
        @test abs(J_E_1_dλ(ξ, lambda, κ, ϵ, λ)) >=
              0.4C_Y_new.J_E_1_dλ * log(ξ) * ξ^(1 / σ - real(lambda) / κ - 1)

        @test abs(J_E_2_dλ(ξ, lambda, κ, ϵ, λ)) <=
              C_Y_new.J_E_2_dλ * log(ξ) * ξ^(1 / σ - real(lambda) / κ - 1)
        @test abs(J_E_2_dλ(ξ, lambda, κ, ϵ, λ)) >=
              0.4C_Y_new.J_E_2_dλ * log(ξ) * ξ^(1 / σ - real(lambda) / κ - 1)

        ######
        ## J_P_1 and J_P_2
        ######

        @test abs(J_P_1(ξ, lambda, κ, ϵ, λ)) <=
              C_Y_new.J_P_1 * exp(real(c) * ξ^2) * ξ^(- 1 / σ + d + real(lambda) / κ - 1)
        @test abs(J_P_1(ξ, lambda, κ, ϵ, λ)) >=
              0.95C_Y_new.J_P_1 *
              exp(real(c) * ξ^2) *
              ξ^(- 1 / σ + d + real(lambda) / κ - 1)

        @test abs(J_P_2(ξ, lambda, κ, ϵ, λ)) <=
              C_Y_new.J_P_2 * exp(real(c) * ξ^2) * ξ^(- 1 / σ + d + real(lambda) / κ - 1)
        @test abs(J_P_2(ξ, lambda, κ, ϵ, λ)) >=
              0.9C_Y_new.J_P_2 * exp(real(c) * ξ^2) * ξ^(- 1 / σ + d + real(lambda) / κ - 1)

        @test abs(J_P_1_dλ(ξ, lambda, κ, ϵ, λ)) <=
              C_Y_new.J_P_1_dλ *
              exp(real(c) * ξ^2) *
              log(ξ) *
              ξ^(- 1 / σ + d + real(lambda) / κ - 1)
        @test abs(J_P_1_dλ(ξ, lambda, κ, ϵ, λ)) >=
              0.4C_Y_new.J_P_1_dλ *
              exp(real(c) * ξ^2) *
              log(ξ) *
              ξ^(- 1 / σ + d + real(lambda) / κ - 1)

        @test abs(J_P_2_dλ(ξ, lambda, κ, ϵ, λ)) <=
              C_Y_new.J_P_2_dλ *
              exp(real(c) * ξ^2) *
              log(ξ) *
              ξ^(- 1 / σ + d + real(lambda) / κ - 1)
        @test abs(J_P_2_dλ(ξ, lambda, κ, ϵ, λ)) >=
              0.4C_Y_new.J_P_2_dλ *
              exp(real(c) * ξ^2) *
              log(ξ) *
              ξ^(- 1 / σ + d + real(lambda) / κ - 1)

        ######
        ## K_1, K_2
        ######

        K_1, K_2 = K_1_2_new(ξ, lambda, κ, ϵ, λ)
        @test norm_inf(K_1) <=
              C_Y_new.K_1 * exp(real(c) * ξ^2) * ξ^(- 1 / σ + d + real(lambda) / κ - 1)
        @test norm_inf(K_1) >=
              0.95C_Y_new.K_1 * exp(real(c) * ξ^2) * ξ^(- 1 / σ + d + real(lambda) / κ - 1)

        @test norm_inf(K_2) <= C_Y_new.K_2 * ξ^(1 / σ - real(lambda) / κ - 1)
        @test norm_inf(K_2) >= 0.95C_Y_new.K_2 * ξ^(1 / σ - real(lambda) / κ - 1)

        K_1_dλ, K_2_dλ = K_1_2_dλ_new(ξ, lambda, κ, ϵ, λ)
        @test norm_inf(K_1_dλ) <=
              C_Y_new.K_1_dλ *
              exp(real(c) * ξ^2) *
              log(ξ) *
              ξ^(- 1 / σ + d + real(lambda) / κ - 1)
        @test norm_inf(K_1_dλ) >=
              0.4C_Y_new.K_1_dλ *
              exp(real(c) * ξ^2) *
              log(ξ) *
              ξ^(- 1 / σ + d + real(lambda) / κ - 1)

        @test norm_inf(K_2_dλ) <= C_Y_new.K_2_dλ * log(ξ) * ξ^(1 / σ - real(lambda) / κ - 1)
        @test norm_inf(K_2_dλ) >=
              0.4C_Y_new.K_2_dλ * log(ξ) * ξ^(1 / σ - real(lambda) / κ - 1)

        K_1, K_2 = K_1_2(ξ, lambda, κ, ϵ, λ)
        @test norm_inf(K_1) <=
              C_Y.K_1 * exp(real(c) * ξ^2) * ξ^(- 1 / σ + d + real(lambda) / κ - 1)
        @test norm_inf(K_1) >=
              0.5C_Y.K_1 * exp(real(c) * ξ^2) * ξ^-(- 1 / σ + d + real(lambda) / κ - 1)

        @test norm_inf(K_2) <= C_Y.K_2 * ξ^(1 / σ - real(lambda) / κ - 1)
        @test norm_inf(K_2) >= 0.5C_Y.K_2 * ξ^(1 / σ - real(lambda) / κ - 1)

        K_1_dλ, K_2_dλ = K_1_2_dλ(ξ, lambda, κ, ϵ, λ)
        @test norm_inf(K_1_dλ) <=
              C_Y.K_1_dλ *
              exp(real(c) * ξ^2) *
              log(ξ) *
              ξ^(- 1 / σ + d + real(lambda) / κ - 1)
        @test norm_inf(K_1_dλ) >=
              0.5C_Y.K_1_dλ *
              exp(real(c) * ξ^2) *
              log(ξ) *
              ξ^(- 1 / σ + d + real(lambda) / κ - 1)

        @test norm_inf(K_2_dλ) <= C_Y.K_1_dλ * log(ξ) * ξ^(1 / σ - real(lambda) / κ - 1)
        @test norm_inf(K_2_dλ) >= 0.5C_Y.K_1_dλ * log(ξ) * ξ^(1 / σ - real(lambda) / κ - 1)
    end
end
