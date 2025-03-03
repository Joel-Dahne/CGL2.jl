"""
    cgl_linearization_equation_real(Q, Q_hat, lambda, κ, ϵ, ξ, λ)
    cgl_linearization_equation_real(Q, (Q_hat, lambda, κ, ϵ, λ), ξ)

Evaluate the right hand side of the forward ODE when written as a four
dimensional real system. It is evaluated at the point
```
Q = [a, b, α, β]
```
and time `ξ`.

For `λ.d != 1` there is a removable singularity at `ξ = 0`. To return
a finite value we in this case required that `α = β = 0`.
"""
function cgl_linearization_equation_real(Y, Q_hat, lambda, κ, ϵ, ξ, λ::CGLParams)
    (; d, ω, σ, δ) = λ

    A_inv = @SMatrix[ϵ 1; -1 ϵ] / (1 + ϵ^2)
    B₁ = @SMatrix[κ 0; 0 κ]
    B₂ = (d - 1) * @SMatrix[ϵ -1; 1 ϵ]
    B = if iszero(ξ)
        @assert iszero(Y[3]) && iszero(Y[4])
        B₁ * ξ
    else
        B₁ * ξ + B₂ / ξ
    end

    C = @SMatrix[κ/σ -ω; ω κ/σ]

    J_N = let
        a, b = Q_hat(ξ)
        N₁_a = -(a^2 + b^2)^(σ - 1) * (δ * (1 + 2σ) * a^2 + 2σ * a * b + δ * b^2)
        N₁_b = -(a^2 + b^2)^(σ - 1) * (a^2 + 2δ * σ * a * b + (1 + 2σ) * b^2)
        N₂_a = (a^2 + b^2)^(σ - 1) * ((1 + 2σ) * a^2 - 2δ * σ * a * b + b^2)
        N₂_b = -(a^2 + b^2)^(σ - 1) * (δ * a^2 - 2σ * a * b + δ * (1 + 2σ) * b^2)

        @SMatrix[N₁_a N₁_b; N₂_a N₂_b]
    end

    Y₁ = SVector(Y[1], Y[2])
    Y₂ = SVector(Y[3], Y[4])

    dY₁ = Y₂
    dY₂ = A_inv * (-(C + J_N - lambda * I) * Y₁ - B * Y₂)
    if ξ == 50
        @show J_N
        display(vcat(dY₁, dY₂))
    end
    return vcat(dY₁, dY₂)
end

# For use with ODEProblem
cgl_linearization_equation_real(u, (Q_hat, lambda, κ, ϵ, λ), ξ) =
    cgl_linearization_equation_real(u, Q_hat, lambda, κ, ϵ, ξ, λ)
