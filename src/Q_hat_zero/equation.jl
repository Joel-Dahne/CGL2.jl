"""
    cgl_hat_equation_real(Q, κ, ϵ, ξ, λ)
    cgl_hat_equation_real(Q, (κ, ϵ, λ), ξ)

Evaluate the right hand side of the forward ODE when written as a four
dimensional real system. It is evaluated at the point
```
Q = [a, b, α, β]
```
and time `ξ`.

For `λ.d != 1` there is a removable singularity at `ξ = 0`. To return
a finite value we in this case required that `α = β = 0`.
"""
function cgl_hat_equation_real(Q, κ, ϵ, ξ, λ::CGLParams)
    (; d, ω, σ, δ) = λ
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
cgl_hat_equation_real(u, (κ, ϵ, λ), ξ) = cgl_hat_equation_real(u, κ, ϵ, ξ, λ)
