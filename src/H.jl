"""
    H(ν_real, ν_imag, c20_hat_real, c20_hat_imag, c10_hat, κ, ϵ, ξ₁, λ::CGLParams)

Compute
```
H(ξ) = Q_hat_0(ξ) - Q_hat_inf(ξ)
```
where `Q_hat_0` is given by [`Q_hat_zero`](@ref) and `Q_hat_inf` by
[`Q_hat_infinity`](@ref). This function returns a vector with four
real values, the first two are the real and imaginary values of `H` at
`ξ₁` and the second two are their derivatives.
"""
function H(
    ν_real::T,
    ν_imag::T,
    c20_hat_real::T,
    c20_hat_imag::T,
    c10_hat::Complex{T},
    κ::T,
    ϵ::T,
    ξ₁::T,
    λ::CGLParams{T},
) where {T}
    Q_hat_0, Q_hat_0_dξ = Q_hat_zero(_complex(ν_real, ν_imag), κ, ϵ, ξ₁, λ)
    Q_hat_inf, Q_hat_inf_dξ = Q_hat_infinity(c10_hat, _complex(c20_hat_real, c20_hat_imag), κ, ϵ, ξ₁, λ)

    H1 = Q_hat_0 - Q_hat_inf
    H2 = Q_hat_0_dξ - Q_hat_inf_dξ

    SVector(real(H1), imag(H1), real(H2), imag(H2))
end
