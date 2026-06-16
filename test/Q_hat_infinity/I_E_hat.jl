@testset "I_E_hat" begin
    @testset "integral_J_E_hat_P_hat" begin
        # Sample values
        κ = Arb(0.8073018593981386)
        ϵ = Arb(0.15002213424487343)
        ξ₁ = Arb(15)
        Λ = CGLParams{Arb}(3, 1, 1, 0)
        a, b, c = CGL2._abc(κ, ϵ, Λ)
        (; d) = Λ

        # Compute with function
        res_1 = CGL2.integral_J_E_hat_P_hat(κ, ϵ, ξ₁, Λ)

        # Compute by integrating using Arblib.integrate
        # Note that the error bounds computed by Arblib.integrate here
        # are not rigorous, since the integrand is not analytic. Since
        # we can't integrate to infinity anyway it doesn't matter too
        # much.
        res_2 =
            CGL2.B_W_hat(κ, ϵ, Λ) * Arblib.integrate(ξ₁, 1000ξ₁, atol = 1e-12) do η
                CGL2.U(b - a, b, c * η^2) *
                CGL2.U(a, b, -c * η^2)^2 *
                conj(CGL2.U(a, b, -c * η^2)) *
                η^(d - 1)
            end

        # Check that they roughly agree. The error in res_2 from not
        # integrating all the way to infinity is somewhat large, so
        # the error will still be relatively large.
        @test isapprox(res_1, res_2, rtol = 1e-5)
    end
end
