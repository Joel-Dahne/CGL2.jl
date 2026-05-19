"""
    cgl_hat_equation_real(Q, κ, ϵ, ξ, Λ)
    cgl_hat_equation_real(Q, (κ, ϵ, Λ), ξ)

Evaluate the right hand side of the forward ODE when written as a four
dimensional real system. It is evaluated at the point
```
Q = [a, b, α, β]
```
and time `ξ`.

For `Λ.d != 1` there is a removable singularity at `ξ = 0`. To return
a finite value we in this case required that `α = β = 0`.
"""
function cgl_hat_equation_real(Q, κ, ϵ, ξ, Λ::CGLParams)
    (; d, ω, σ, δ) = Λ
    a, b, α, β = Q

    a2b2σ = (a^2 + b^2)^σ

    # Same as cgl_equation_real except changing the sign of κ and ω
    @fastmath F1 = -κ * ξ * β - κ / σ * b - ω * a - a2b2σ * a + δ * a2b2σ * b
    @fastmath F2 = κ * ξ * α + κ / σ * a - ω * b - a2b2σ * b - δ * a2b2σ * a

    if !isone(d) && !(iszero(ξ) && iszero(α) && iszero(β))
        F1 -= (d - 1) / ξ * (α + ϵ * β)
        F2 -= (d - 1) / ξ * (β - ϵ * α)
    end

    return SVector(α, β, (F1 - ϵ * F2) / (1 + ϵ^2), (ϵ * F1 + F2) / (1 + ϵ^2))
end

# For use with ODEProblem
cgl_hat_equation_real(u, (κ, ϵ, Λ), ξ) = cgl_hat_equation_real(u, κ, ϵ, ξ, Λ)

"""
    cgl_hat_equation_real_matrix_form(Q, κ, ϵ, ξ, Λ)
    cgl_hat_equation_real_matrix_form(Q, (κ, ϵ, Λ), ξ)

Same as [`cgl_hat_equation_real`](@ref), but the computations are done
by writing the system in matrix form. This method is primarily
intended to test that the system is correctly written down.
"""
function cgl_hat_equation_real_matrix_form(Q, κ, ϵ, ξ, Λ::CGLParams)
    (; d, ω, σ, δ) = Λ

    # Rename things to follow same notation as in paper
    X = Q

    A =
        1 / (1 + ϵ^2) * @SMatrix[
            0 0 (1+ϵ^2) 0
            0 0 0 (1+ϵ^2)
            (-ω-ϵ*κ/σ) (ϵ*ω-κ/σ) (-(1 + ϵ^2)*(d-1)/ξ-ϵ*κ*ξ) -κ*ξ
            -(ϵ * ω - κ / σ) (-ω-ϵ*κ/σ) κ*ξ (-(1 + ϵ^2)*(d-1)/ξ-ϵ*κ*ξ)
        ]

    a, b, _, _ = X
    N₁ = -(a^2 + b^2)^σ * a + δ * (a^2 + b^2)^σ * b
    N₂ = -(a^2 + b^2)^σ * b - δ * (a^2 + b^2)^σ * a
    N = 1 / (1 + ϵ^2) * SVector(0, 0, N₁ - ϵ * N₂, ϵ * N₁ + N₂)

    return A * X + N
end

"""
    cgl_hat_equation_real_second_order(dQ, Q, κ, ϵ, ξ, Λ)
    cgl_hat_equation_real_second_order(dQ, Q, (κ, ϵ, Λ), ξ)

Evaluate the right hand side of the forward ODE when written as a
second order two dimensional real system.

For `Λ.d != 1` there is a removable singularity at `ξ = 0`. To return
a finite value we in this case required that `α = β = 0`.
"""
function cgl_hat_equation_real_second_order(dQ, Q, κ, ϵ, ξ, Λ::CGLParams)
    (; d, ω, σ, δ) = Λ
    a, b = Q
    α, β = dQ

    a2b2σ = (a^2 + b^2)^σ

    # Same as cgl_equation_real except changing the sign of κ and ω
    @fastmath F1 = -κ * ξ * β - κ / σ * b - ω * a - a2b2σ * a + δ * a2b2σ * b
    @fastmath F2 = κ * ξ * α + κ / σ * a - ω * b - a2b2σ * b - δ * a2b2σ * a

    if !isone(d) && !(iszero(ξ) && iszero(α) && iszero(β))
        F1 -= (d - 1) / ξ * (α + ϵ * β)
        F2 -= (d - 1) / ξ * (β - ϵ * α)
    end

    res = [(F1 - ϵ * F2) / (1 + ϵ^2), (ϵ * F1 + F2) / (1 + ϵ^2)]

    return res
end

# For use with ODEProblem
cgl_hat_equation_real_second_order(du, u, (κ, ϵ, Λ), ξ) =
    cgl_hat_equation_real_second_order(du, u, κ, ϵ, ξ, Λ)
