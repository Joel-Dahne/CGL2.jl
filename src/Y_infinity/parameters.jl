# Matrices
function coeff_matrices(lambda, κ, ϵ, λ::CGLParams)
    (; ω, σ, d) = λ

    A = SMatrix{2,2}(ϵ, 1, -1, ϵ)
    B₁ = SMatrix{2,2}(κ, 0, 0, κ)
    B₂ = (d - 1) * A
    C = SMatrix{2,2}(κ / σ, ω, -ω, κ / σ)
    λI = SMatrix{2,2}(lambda, 0, 0, lambda)

    return (; A, B₁, B₂, C, λI)
end

function J_N_coeff_matrices(λ::CGLParams)
    (; σ, δ) = λ

    M1 = SMatrix{2,2}(δ * (1 + 2σ), -(1 + 2σ), 1, δ)
    M2 = SMatrix{2,2}(1, δ, δ, -1)
    M3 = SMatrix{2,2}(δ, -1, 1 + 2σ, δ * (1 + 2σ))

    return (; M1, M2, M3)
end
