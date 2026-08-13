"""
    H(λ, ν, γ₁, γ₂, κ, ϵ, ξ₁, Λ::CGLParams; bisect_ν::Int = 2)

Compute the determinant of the `4x4` matrix with columns given by
`Y_0_1(ξ)`, `Y_0_2(ξ)`, `Y_inf_1(ξ)` and `Y_inf_2(ξ)`.

To improve the computed enclosure, the value for `ν` is split into
smaller pieces and the union over them is taken. The argument
`bisect_ν` is used to determine how many times the real and imaginary
parts of `ν` are bisected.
"""
function H(
    λ::Union{Complex{T},Acb},
    ν::Union{Complex{T},Acb},
    γ₁::Union{Complex{T},Acb},
    γ₂::Union{Complex{T},Acb},
    κ::T,
    ϵ::T,
    ξ₁::T,
    Λ::CGLParams{T};
    bisect_ν::Int = 2,
) where {T}
    complex_T = ifelse(T == Arb, Acb, Complex{T})

    if bisect_ν > 0 && iswide(ν)
        νs_real =
            Arb.(ArbExtras.bisect_interval_recursive(getinterval(real(ν))..., bisect_ν))
        νs_imag =
            Arb.(ArbExtras.bisect_interval_recursive(getinterval(imag(ν))..., bisect_ν))
        νs = [Acb(ν_real, ν_imag) for ν_real in νs_real, ν_imag in νs_imag]

        Y_0_1 = reduce(
            (v, w) -> Arblib.union.(v, w),
            [Y_zero(SVector{2,complex_T}(1, 0), λ, ν, κ, ϵ, ξ₁, Λ) for ν in νs],
        )
        Y_0_2 = reduce(
            (v, w) -> Arblib.union.(v, w),
            [Y_zero(SVector{2,complex_T}(0, 1), λ, ν, κ, ϵ, ξ₁, Λ) for ν in νs],
        )
    else
        Y_0_1 = Y_zero(SVector{2,complex_T}(1, 0), λ, ν, κ, ϵ, ξ₁, Λ)
        Y_0_2 = Y_zero(SVector{2,complex_T}(0, 1), λ, ν, κ, ϵ, ξ₁, Λ)
    end

    Y_inf_1 = Y_infinity(SVector{2,complex_T}(1, 0), λ, γ₁, γ₂, κ, ϵ, ξ₁, Λ)
    Y_inf_2 = Y_infinity(SVector{2,complex_T}(0, 1), λ, γ₁, γ₂, κ, ϵ, ξ₁, Λ)

    M = hcat(Y_0_1, Y_0_2, Y_inf_1, Y_inf_2)

    return det_4x4(M)
end
