"""
    linearization_real_matrix(Q₀, κ, ϵ, ξ, λ::CGLParams)

Computes the matrix corresponding to the linearization around the
solution `Q₀`. The linearized equation can be written on the form
```
Y' = B(ξ) * Y
```
for some matrix `B` depending on `Q₀`. The original ODE is of the form
```
Q` = A(ξ) * Q + N(Q)
```
where `A` is the linear part and `N` the non-linear part. The matrix
`B` is given by
```
B(ξ) = A(ξ) + J(N)(Q₀)
```
where `J(N)(Q₀)` denotes the Jacobian of `N` evaluated at `Q₀`.

This should return the same result as
```
ForwardDiff.jacobian(Q₀) do Q
    CGL2.cgl_hat_equation_real_matrix_form(Q, κ, ϵ, ξ₁, λ)
end
```
"""
function linearization_real_matrix(Q₀, κ, ϵ, ξ, λ::CGLParams)
    (; d, ω, σ, δ) = λ

    A =
        1 / (1 + ϵ^2) * @SMatrix[
            0 0 (1+ϵ^2) 0
            0 0 0 (1+ϵ^2)
            (-ω-ϵ*κ/σ) (ϵ*ω-κ/σ) (-(1 + ϵ^2)*(d-1)/ξ-ϵ*κ*ξ) -κ*ξ
            -(ϵ * ω - κ / σ) (-ω-ϵ*κ/σ) κ*ξ (-(1 + ϵ^2)*(d-1)/ξ-ϵ*κ*ξ)
        ]

    a, b, _, _ = Q₀

    N₁_a = -(a^2 + b^2)^(σ - 1) * ((1 + 2σ) * a^2 - 2δ * σ * a * b + b^2)
    N₁_b = (a^2 + b^2)^(σ - 1) * (δ * a^2 - 2σ * a * b + δ * (1 + 2σ) * b^2)
    N₂_a = -(a^2 + b^2)^(σ - 1) * (δ * (1 + 2σ) * a^2 + 2σ * a * b + δ * b^2)
    N₂_b = -(a^2 + b^2)^(σ - 1) * (a^2 + 2δ * σ * a * b + (1 + 2σ) * b^2)

    J_N =
        1 / (1 + ϵ^2) * @SMatrix[
            0 0 0 0
            0 0 0 0
            (N₁_a-ϵ*N₂_a) (N₁_b-ϵ*N₂_b) 0 0
            (ϵ*N₁_a+N₂_a) (ϵ*N₁_b+N₂_b) 0 0
        ]

    M = A + J_N

    return M
end
