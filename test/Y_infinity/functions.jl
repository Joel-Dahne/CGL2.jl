@testset "functions" begin
    # IMPROVE: Test for more parameters
    λ = Acb(0.19028080950252219, 2.769747313863597)
    γ₁ = Acb(0.1979686615344728, 0.15613711743772288)
    γ₂ = Acb(-116.03421868436834, 101.21494871936197)
    ξ = Arb(30)
    κ = Arb(0.8073018593981386)
    ϵ = Arb(0.15002213424487343)
    Λ = CGLParams{Arb}(3, 1.0, 1.0, 0.0)
    (; d, ω, σ) = Λ

    ξF64 = Float64(ξ)
    κF64 = Float64(κ)
    ϵF64 = Float64(ϵ)
    λF64 = ComplexF64(λ)
    ΛF64 = CGLParams{Float64}(Λ)

    F_Y = CGL2.FunctionEnclosures_Y(λ, γ₁, γ₂, κ, ϵ, ξ, Λ)

    # Function for computing derivative using finite differences.
    fdm = central_fdm(5, 1)
    fdm2 = central_fdm(5, 2)
    fdm3 = central_fdm(5, 3)

    @testset "P_$j" for (j, P_j, P_j_dξ) in
                        [(1, CGL2.P_1, CGL2.P_1_dξ), (2, CGL2.P_2, CGL2.P_2_dξ)]
        Pj_series = P_j(ArbSeries((ξ, 1, 0)), λ, κ, ϵ, Λ)
        Pj = Pj_series[0]
        Pj_dξ = Pj_series[1]
        Pj_dξ_dξ = 2Pj_series[2]

        # Test that it solves the ODE
        if j == 1
            @test Arblib.contains_zero(
                (ϵ + im) * Pj_dξ_dξ +
                (κ * ξ + (d - 1) * (ϵ + im) / ξ) * Pj_dξ +
                (κ / σ + ω * im - λ) * Pj,
            )
        else
            @test Arblib.contains_zero(
                (ϵ - im) * Pj_dξ_dξ +
                (κ * ξ + (d - 1) * (ϵ - im) / ξ) * Pj_dξ +
                (κ / σ - ω * im - λ) * Pj,
            )
        end

        @test Arblib.overlaps(P_j(ArbSeries((ξ, 1)), λ, κ, ϵ, Λ)[1], P_j_dξ(ξ, λ, κ, ϵ, Λ))

        @test P_j_dξ(ξ, λ, κ, ϵ, Λ) ≈ fdm(ξ -> P_j(ξ, λF64, κF64, ϵF64, ΛF64), ξF64) rtol =
            1e-12
    end

    @testset "E_$j" for (j, E_j, E_j_dξ) in
                        [(1, CGL2.E_1, CGL2.E_1_dξ), (2, CGL2.E_2, CGL2.E_2_dξ)]
        Ej_series = E_j(ArbSeries((ξ, 1, 0)), λ, κ, ϵ, Λ)
        Ej = Ej_series[0]
        Ej_dξ = Ej_series[1]
        Ej_dξ_dξ = 2Ej_series[2]

        # Test that it solves the ODE
        if j == 1
            @test Arblib.contains_zero(
                (ϵ + im) * Ej_dξ_dξ +
                (κ * ξ + (d - 1) * (ϵ + im) / ξ) * Ej_dξ +
                (κ / σ + ω * im - λ) * Ej,
            )
        else
            @test Arblib.contains_zero(
                (ϵ - im) * Ej_dξ_dξ +
                (κ * ξ + (d - 1) * (ϵ - im) / ξ) * Ej_dξ +
                (κ / σ - ω * im - λ) * Ej,
            )
        end

        @test Arblib.overlaps(E_j(ArbSeries((ξ, 1)), λ, κ, ϵ, Λ)[1], E_j_dξ(ξ, λ, κ, ϵ, Λ))

        @test E_j_dξ(ξ, λ, κ, ϵ, Λ) ≈ fdm(ξ -> E_j(ξ, λF64, κF64, ϵF64, ΛF64), ξF64) rtol =
            1e-10
    end

    @testset "W_$j" for (j, W_j, P_j, P_j_dξ, E_j, E_j_dξ) in [
        (1, CGL2.W_1, CGL2.P_1, CGL2.P_1_dξ, CGL2.E_1, CGL2.E_1_dξ),
        (2, CGL2.W_2, CGL2.P_2, CGL2.P_2_dξ, CGL2.E_2, CGL2.E_2_dξ),
    ]
        @test Arblib.overlaps(
            W_j(ξ, λ, κ, ϵ, Λ),
            P_j(ξ, λ, κ, ϵ, Λ) * E_j_dξ(ξ, λ, κ, ϵ, Λ) -
            P_j_dξ(ξ, λ, κ, ϵ, Λ) * E_j(ξ, λ, κ, ϵ, Λ),
        )
    end

    @testset "J_P_$j" for (j, J_P_j, P_j, W_j) in [
        (1, CGL2.J_P_1, CGL2.P_1, CGL2.W_1),
        (2, CGL2.J_P_2, CGL2.P_2, CGL2.W_2),
    ]
        @test Arblib.overlaps(J_P_j(ξ, λ, κ, ϵ, Λ), P_j(ξ, λ, κ, ϵ, Λ) / W_j(ξ, λ, κ, ϵ, Λ))
    end

    @testset "J_E_$j" for (j, J_E_j, E_j, W_j) in [
        (1, CGL2.J_E_1, CGL2.E_1, CGL2.W_1),
        (2, CGL2.J_E_2, CGL2.E_2, CGL2.W_2),
    ]
        @test Arblib.overlaps(J_E_j(ξ, λ, κ, ϵ, Λ), E_j(ξ, λ, κ, ϵ, Λ) / W_j(ξ, λ, κ, ϵ, Λ))
    end

    @testset "K_1 and K_2" begin
        # K1 and K2 are suppose to give solutions to the linear system
        # Ψ * v = [[0, 0]; A \ F].

        A = SMatrix{2,2}(ϵ, 1, -1, ϵ)
        P = SMatrix{2,2,Acb}(im, 1, -im, 1)

        Ψ = [F_Y.E_12 F_Y.P_12; F_Y.E_12_dξ F_Y.P_12_dξ]
        F = eltype(Ψ)[0.1, 0.25]

        v = [-F_Y.K_1 * F; -F_Y.K_2 * F]

        @test all(Arblib.overlaps.(Ψ * v, [[0, 0]; inv(P) * A * P \ F]))
    end

    @testset "I_N" begin
        P = SMatrix{2,2,Acb}(im, 1, -im, 1)

        # The precise value for a and b should not play any role in
        # the correctness, we just compute some approximation here.
        νF64 = 1.9261384880241954 + 3.0638598354170337im
        a_hat, b_hat, a_hat_dξ, b_hat_dξ =
            Arb.(CGL2.Q_hat_zero_float(real(νF64), imag(νF64), κF64, ϵF64, ξF64, ΛF64))

        # Compute Jacobian by going through ArbSeries
        N = (a_hat, b_hat) -> (a_hat^2 + b_hat^2) ^ Λ.σ * SVector(-Λ.δ * a_hat - b_hat, a_hat - Λ.δ * b_hat)
        N_a = getindex.(N(ArbSeries((a_hat, 1)), b_hat), 1)
        N_b = getindex.(N(a_hat, ArbSeries((b_hat, 1))), 1)
        J_N_direct = [N_a N_b]
        I_N_direct = inv(P) * J_N_direct * P

        @test all(Arblib.overlaps.(I_N_direct, CGL2.I_N(Acb(a_hat, b_hat), Λ)))
    end
end
