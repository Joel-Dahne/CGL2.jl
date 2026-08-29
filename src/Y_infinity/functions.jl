"""
    FunctionEnclosures_Y(λ, γ₁, γ₂, κ, ϵ, ξ₁, Λ)

Contains enclosures of the functions

- [`E_12`](@ref)
- [`E_12_dξ`](@ref)
- [`P_12`](@ref)
- [`P_12_dξ`](@ref)

when evaluated at `ξ₁`.

Even though none of the functions depend on `γ₁` and `γ₂` we keep them
as arguments to mirror [`FunctionBounds_Y`](@ref).
"""
struct FunctionEnclosures_Y
    E_12::Diagonal{Acb,SVector{2,Acb}}
    E_12_dξ::Diagonal{Acb,SVector{2,Acb}}
    P_12::Diagonal{Acb,SVector{2,Acb}}
    P_12_dξ::Diagonal{Acb,SVector{2,Acb}}
end

function FunctionEnclosures_Y(
    λ::Acb,
    γ₁::Acb,
    γ₂::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    Λ::CGLParams{Arb},
)
    return FunctionEnclosures_Y(
        Diagonal(SVector(E_1(ξ₁, λ, κ, ϵ, Λ), E_2(ξ₁, λ, κ, ϵ, Λ))),
        Diagonal(SVector(E_1_dξ(ξ₁, λ, κ, ϵ, Λ), E_2_dξ(ξ₁, λ, κ, ϵ, Λ))),
        Diagonal(SVector(P_1(ξ₁, λ, κ, ϵ, Λ), P_2(ξ₁, λ, κ, ϵ, Λ))),
        Diagonal(SVector(P_1_dξ(ξ₁, λ, κ, ϵ, Λ), P_2_dξ(ξ₁, λ, κ, ϵ, Λ))),
    )
end

function P_1(ξ, λ, κ, ϵ, Λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, Λ)
    z = -c * ξ^2
    return U(a - λ, b, z)
end

function P_1_dξ(ξ, λ, κ, ϵ, Λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, Λ)
    z = -c * ξ^2
    z_dξ = -2c * ξ
    return U_dz(a - λ, b, z) * z_dξ
end

function P_2(ξ, λ, κ, ϵ, Λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, Λ)
    z = -conj(c) * ξ^2
    return U(conj(a) - λ, b, z)
end

function P_2_dξ(ξ, λ, κ, ϵ, Λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, Λ)
    z = -conj(c) * ξ^2
    z_dξ = -2conj(c) * ξ
    return U_dz(conj(a) - λ, b, z) * z_dξ
end

function E_1(ξ, λ, κ, ϵ, Λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, Λ)
    z = -c * ξ^2
    return exp(z) * U(b - a + λ, b, -z)
end

function E_1_dξ(ξ, λ, κ, ϵ, Λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, Λ)
    z = -c * ξ^2
    z_dξ = -2c * ξ
    return exp(z) * (U(b - a + λ, b, -z) - U_dz(b - a + λ, b, -z)) * z_dξ
end

function E_2(ξ, λ, κ, ϵ, Λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, Λ)
    z = -conj(c) * ξ^2
    return exp(z) * U(b - conj(a) + λ, b, -z)
end

function E_2_dξ(ξ, λ, κ, ϵ, Λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, Λ)
    z = -conj(c) * ξ^2
    z_dξ = -2conj(c) * ξ
    return exp(z) * (U(b - conj(a) + λ, b, -z) - U_dz(b - conj(a) + λ, b, -z)) * z_dξ
end

function B_W_1(λ, κ, ϵ, Λ::CGLParams)
    a, b, c = _abc(κ, ϵ, Λ)
    sgn = c isa AcbSeries ? sign(imag(c[0])) : sign(imag(c))
    return 1 // 2 * exp(sgn * im * (b - a + λ) * π) * (-c)^(b - 1)
end

function B_W_2(λ, κ, ϵ, Λ::CGLParams)
    a, b, c = _abc(κ, ϵ, Λ)
    sgn = c isa AcbSeries ? sign(imag(c[0])) : sign(imag(c))
    return 1 // 2 * exp(-sgn * im * (b - conj(a) + λ) * π) * (-conj(c))^(b - 1)
end

function J_P_1(ξ, λ, κ, ϵ, Λ::CGLParams)
    (; d) = Λ
    c = _c(κ, ϵ, Λ)
    return B_W_1(λ, κ, ϵ, Λ) * P_1(ξ, λ, κ, ϵ, Λ) * exp(c * ξ^2) * ξ^(d - 1)
end

function J_P_2(ξ, λ, κ, ϵ, Λ::CGLParams)
    (; d) = Λ
    c = _c(κ, ϵ, Λ)
    return B_W_2(λ, κ, ϵ, Λ) * P_2(ξ, λ, κ, ϵ, Λ) * exp(conj(c) * ξ^2) * ξ^(d - 1)
end

function J_E_1(ξ, λ, κ, ϵ, Λ::CGLParams)
    (; d) = Λ
    c = _c(κ, ϵ, Λ)
    return B_W_1(λ, κ, ϵ, Λ) * E_1(ξ, λ, κ, ϵ, Λ) * exp(c * ξ^2) * ξ^(d - 1)
end

function J_E_2(ξ, λ, κ, ϵ, Λ::CGLParams)
    (; d) = Λ
    c = _c(κ, ϵ, Λ)
    return B_W_2(λ, κ, ϵ, Λ) * E_2(ξ, λ, κ, ϵ, Λ) * exp(conj(c) * ξ^2) * ξ^(d - 1)
end

# These functions are only used for testing

function W_1(ξ, λ, κ, ϵ, Λ::CGLParams)
    a, b, c = _abc(κ, ϵ, Λ)
    z = -c * ξ^2
    sgn = c isa AcbSeries ? sign(imag(c[0])) : sign(imag(c))
    return -2c * exp(-sgn * im * (b - a + λ) * π) * ξ * z^-b * exp(z)
end

function W_2(ξ, λ, κ, ϵ, Λ::CGLParams)
    a, b, c = _abc(κ, ϵ, Λ)
    z = -conj(c) * ξ^2
    sgn = c isa AcbSeries ? sign(imag(c[0])) : sign(imag(c))
    return -2conj(c) * exp(sgn * im * (b - conj(a) + λ) * π) * ξ * z^-b * exp(z)
end

function K_1(ξ, λ, κ, ϵ, Λ::CGLParams)
    return Diagonal(
        SVector((ϵ - im) * J_P_1(ξ, λ, κ, ϵ, Λ), (ϵ + im) * J_P_2(ξ, λ, κ, ϵ, Λ)),
    ) / (1 + ϵ^2)
end

function K_2(ξ, λ, κ, ϵ, Λ::CGLParams)
    return Diagonal(
        SVector((ϵ - im) * J_E_1(ξ, λ, κ, ϵ, Λ), (ϵ + im) * J_E_2(ξ, λ, κ, ϵ, Λ)),
    ) / (1 + ϵ^2)
end

function I_N(Q_hat, Λ::CGLParams)
    # The formula below assumes σ = 1 and δ = 0
    @assert isone(Λ.σ)
    @assert iszero(Λ.δ)

    return SMatrix{2,2}(
        2im * abs2(Q_hat),
        conj(-im * Q_hat^2),
        -im * Q_hat^2,
        -2im * abs2(Q_hat),
    )
end
