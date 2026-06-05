"""
    H(λ, ν, γ₁, γ₂, κ, ϵ, ξ₁, Λ::CGLParams)

Compute the determinant of the `4x4` matrix with columns given by
`Y_0_1(ξ)`, `Y_0_2(ξ)`, `Y_inf_1(ξ)` and `Y_inf_2(ξ)`.
"""
function H(
    λ::Union{Complex{T},Acb},
    ν::Union{Complex{T},Acb},
    γ₁::Union{Complex{T},Acb},
    γ₂::Union{Complex{T},Acb},
    κ::T,
    ϵ::T,
    ξ₁::T,
    Λ::CGLParams{T},
) where {T}
    complex_T = ifelse(T == Arb, Acb, Complex{T})

    Y_0_1 = Y_zero(SVector{2,complex_T}(1, 0), λ, ν, κ, ϵ, ξ₁, Λ)
    Y_0_2 = Y_zero(SVector{2,complex_T}(0, 1), λ, ν, κ, ϵ, ξ₁, Λ)
    Y_inf_1 = Y_infinity(SVector{2,complex_T}(1, 0), λ, γ₁, γ₂, κ, ϵ, ξ₁, Λ)
    Y_inf_2 = Y_infinity(SVector{2,complex_T}(0, 1), λ, γ₁, γ₂, κ, ϵ, ξ₁, Λ)

    M = hcat(Y_0_1, Y_0_2, Y_inf_1, Y_inf_2)

    return det_4x4(M)
end
