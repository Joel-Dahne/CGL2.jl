"""
H(lambda, ν, γ₁, γ₂, κ, ϵ, ξ₁, λ::CGLParams)

Compute the determinant of the `4x4` matrix with columns given by
`Y_0_1(ξ)`, `Y_0_2(ξ)`, `Y_inf_1(ξ)` and `Y_inf_2(ξ)`.
"""
function H(
    lambda::Union{Complex{T},Acb},
    ν::Union{Complex{T},Acb},
    γ₁::Union{Complex{T},Acb},
    γ₂::Union{Complex{T},Acb},
    κ::T,
    ϵ::T,
    ξ₁::T,
    λ::CGLParams{T},
) where {T}
    if iswide(ν)
        # FIXME: Don't cheat!
        ν_endpoints = [
            Acb(lbound(real(ν)), lbound(imag(ν))),
            Acb(lbound(real(ν)), ubound(imag(ν))),
            Acb(ubound(real(ν)), lbound(imag(ν))),
            Acb(ubound(real(ν)), ubound(imag(ν))),
        ]
        res = tmap(ν_endpoints) do ν
            H(lambda, ν, γ₁, γ₂, κ, ϵ, ξ₁, λ)
        end
        return foldl(Arblib.union, res)
    end
    # IMPROVE: For wide values of ν this gives large overestimations.
    # Look at computing derivatives in ν to get better enclosures.
    complex_T = ifelse(T == Arb, Acb, Complex{T})
    Y_0_1 = Y_zero(SVector{2,complex_T}(1, 0), lambda, ν, κ, ϵ, ξ₁, λ)
    Y_0_2 = Y_zero(SVector{2,complex_T}(0, 1), lambda, ν, κ, ϵ, ξ₁, λ)
    Y_inf_1 = Y_infinity(SVector{2,complex_T}(1, 0), lambda, γ₁, γ₂, κ, ϵ, ξ₁, λ)
    Y_inf_2 = Y_infinity(SVector{2,complex_T}(0, 1), lambda, γ₁, γ₂, κ, ϵ, ξ₁, λ)

    M = hcat(Y_0_1, Y_0_2, Y_inf_1, Y_inf_2)

    return det_4x4(M)
end

function H(
    lambda::AcbSeries,
    ν::Acb,
    γ₁::Acb,
    γ₂::Acb,
    κ::T,
    ϵ::T,
    ξ₁::T,
    λ::CGLParams{T},
) where {T}
    @assert Arblib.degree(lambda) == 1
    lambda₀ = lambda[0]
    # IMPROVE: Compute Y_0_1 and Y_0_1_derivative together
    Y_0_1 = Y_zero(SVector{2,Acb}(1, 0), lambda₀, ν, κ, ϵ, ξ₁, λ)
    Y_0_2 = Y_zero(SVector{2,Acb}(0, 1), lambda₀, ν, κ, ϵ, ξ₁, λ)
    Y_inf_1 = Y_infinity(SVector{2,Acb}(1, 0), lambda₀, γ₁, γ₂, κ, ϵ, ξ₁, λ)
    Y_inf_2 = Y_infinity(SVector{2,Acb}(0, 1), lambda₀, γ₁, γ₂, κ, ϵ, ξ₁, λ)

    Y_0_1_derivative = Y_zero_derivative(SVector{2,Acb}(1, 0), lambda₀, ν, κ, ϵ, ξ₁, λ)
    Y_0_2_derivative = Y_zero_derivative(SVector{2,Acb}(0, 1), lambda₀, ν, κ, ϵ, ξ₁, λ)
    Y_inf_1_derivative =
        Y_infinity_derivative(SVector{2,Acb}(1, 0), lambda₀, γ₁, γ₂, κ, ϵ, ξ₁, λ)
    Y_inf_2_derivative =
        Y_infinity_derivative(SVector{2,Acb}(0, 1), lambda₀, γ₁, γ₂, κ, ϵ, ξ₁, λ)

    M = hcat(Y_0_1, Y_0_2, Y_inf_1, Y_inf_2)

    M_derivative =
        hcat(Y_0_1_derivative, Y_0_2_derivative, Y_inf_1_derivative, Y_inf_2_derivative)

    res = det_4x4(AcbSeries.(tuple.(M, M_derivative)))

    return ArbExtras.compose_zero!(res, res, lambda)
end
