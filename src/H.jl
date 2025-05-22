"""
H(x, c, lambda, ν, γ₁, γ₂, κ, ϵ, ξ₁, λ::CGLParams)

Compute
```
H(ξ) = Y_0(ξ) - Y_inf(ξ)
```
where `Y_0` is given by [`Y_zero`](@ref) and `Y_inf` by
[`Y_infinity`](@ref). This function returns a vector with two complex
vectors, the first is the value of `H` at `ξ₁` and the second is the
derivative.
"""
function H(
    x::Union{Complex{T},Acb},
    c::SVector{2,<:Union{Complex{T},Acb}},
    lambda::Union{Complex{T},Acb},
    ν::Union{Complex{T},Acb},
    γ₁::Union{Complex{T},Acb},
    γ₂::Union{Complex{T},Acb},
    κ::T,
    ϵ::T,
    ξ₁::T,
    λ::CGLParams{T},
) where {T}
    Y_0 = Y_zero(x, lambda, ν, κ, ϵ, ξ₁, λ)
    Y_inf = Y_infinity(c, lambda, γ₁, γ₂, κ, ϵ, ξ₁, λ)

    return Y_0 - Y_inf
end

"""
H_precomputed(x, c, ν, γ₁, γ₂, κ, ϵ, ξ₁, Q_hat, λ::CGLParams)

Same as [`H`](@ref), but instead of taking `ν`, `γ₁` and `γ₂` as
arguments it take a precomputed `ODESolution` representing the forward
solution `Q_hat`.
"""
function H_precomputed(
    x::Complex{T},
    c::SVector{2,Complex{T}},
    lambda::Complex{T},
    κ::T,
    ϵ::T,
    ξ₁::T,
    Q_hat::ODESolution{T},
    λ::CGLParams{T},
) where {T}
    Y_0 = Y_zero(x, lambda, κ, ϵ, ξ₁, Q_hat, λ)
    Y_inf = Y_infinity(c, lambda, κ, ϵ, ξ₁, complex(Q_hat(ξ₁)[1:2]...), λ)

    return Y_0 - Y_inf
end

"""
H_jacobian(x, c, lambda, ν, γ₁, γ₂, κ, ϵ, ξ₁, λ::CGLParams)

This function computes the Jacobian of [`H`](@ref) w.r.t. the
parameters `x`, `c` and `lambda`.
"""
function H_jacobian(
    x::Union{Complex{T},Acb},
    c::SVector{2,<:Union{Complex{T},Acb}},
    lambda::Union{Complex{T},Acb},
    ν::Union{Complex{T},Acb},
    γ₁::Union{Complex{T},Acb},
    γ₂::Union{Complex{T},Acb},
    κ::T,
    ϵ::T,
    ξ₁::T,
    λ::CGLParams{T},
) where {T}
    Y_0_J = Y_zero_jacobian(x, lambda, ν, κ, ϵ, ξ₁, λ)
    Y_inf_J = Y_infinity_jacobian(c, lambda, γ₁, γ₂, κ, ϵ, ξ₁, λ)

    return SMatrix{4,4}(
        # Derivative with respect to x
        Y_0_J[:, 1]...,
        # Derivative with respect to c₁
        -Y_inf_J[:, 1]...,
        # Derivative with respect to c₂
        -Y_inf_J[:, 2]...,
        # Derivative with respect to lambda
        Y_0_J[:, 2] - Y_inf_J[:, 3]...,
    )
end
