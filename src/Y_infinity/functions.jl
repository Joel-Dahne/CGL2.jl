"""
    FunctionEnclosures_Y(κ, ϵ, ξ₁, λ)

Contains enclosures of the functions

- [`E_12`](@ref)
- [`E_12_dξ`](@ref)
- [`P_12`](@ref)
- [`P_12_dξ`](@ref)
- [`K_1`](@ref)
- [`K_2`](@ref)
- [`K_1_dξ`](@ref)
- [`K_2_dξ`](@ref)
- [`I_N`](@ref)
- [`I_N_dξ`](@ref)

when evaluated at `ξ₁`.
"""
struct FunctionEnclosures_Y
    E_12::Diagonal{Acb,SVector{2,Acb}}
    E_12_dξ::Diagonal{Acb,SVector{2,Acb}}
    P_12::Diagonal{Acb,SVector{2,Acb}}
    P_12_dξ::Diagonal{Acb,SVector{2,Acb}}
    K_1::Diagonal{Acb,SVector{2,Acb}}
    K_2::Diagonal{Acb,SVector{2,Acb}}
    K_1_dξ::Diagonal{Acb,SVector{2,Acb}}
    K_2_dξ::Diagonal{Acb,SVector{2,Acb}}
    I_N::SMatrix{2,2,Acb}
    I_N_dξ::SMatrix{2,2,Acb}
end

function FunctionEnclosures_Y(
    lambda::Acb,
    γ₁::Acb,
    γ₂::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
)
    # Compute enclosure of forward solution
    Q_hat, Q_hat_dξ = Q_hat_infinity(γ₁, γ₂, κ, ϵ, ξ₁, λ)

    return FunctionEnclosures_Y(
        Diagonal(SVector(E_1(ξ₁, lambda, κ, ϵ, λ), E_2(ξ₁, lambda, κ, ϵ, λ))),
        Diagonal(SVector(E_1_dξ(ξ₁, lambda, κ, ϵ, λ), E_2_dξ(ξ₁, lambda, κ, ϵ, λ))),
        Diagonal(SVector(P_1(ξ₁, lambda, κ, ϵ, λ), P_2(ξ₁, lambda, κ, ϵ, λ))),
        Diagonal(SVector(P_1_dξ(ξ₁, lambda, κ, ϵ, λ), P_2_dξ(ξ₁, lambda, κ, ϵ, λ))),
        K_1(ξ₁, lambda, κ, ϵ, λ),
        K_2(ξ₁, lambda, κ, ϵ, λ),
        K_1_dξ(ξ₁, lambda, κ, ϵ, λ),
        K_2_dξ(ξ₁, lambda, κ, ϵ, λ),
        I_N(Q_hat, λ),
        I_N_dξ(Q_hat, Q_hat_dξ, λ),
    )
end

function P_1(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, λ)
    z = -c * ξ^2
    return U(a - lambda / 2κ, b, z)
end

function P_1_dξ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, λ)
    z = -c * ξ^2
    z_dξ = -2c * ξ
    return U_dz(a - lambda / 2κ, b, z) * z_dξ
end

function P_2(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, λ)
    z = -conj(c) * ξ^2
    return U(conj(a) - lambda / 2κ, b, z)
end

function P_2_dξ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, λ)
    z = -conj(c) * ξ^2
    z_dξ = -2conj(c) * ξ
    return U_dz(conj(a) - lambda / 2κ, b, z) * z_dξ
end

function E_1(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, λ)
    z = -c * ξ^2
    return exp(z) * U(b - a + lambda / 2κ, b, -z)
end

function E_1_dξ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, λ)
    z = -c * ξ^2
    z_dξ = -2c * ξ
    return exp(z) *
           (U(b - a + lambda / 2κ, b, -z) - U_dz(b - a + lambda / 2κ, b, -z)) *
           z_dξ
end

function E_2(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, λ)
    z = -conj(c) * ξ^2
    return exp(z) * U(b - conj(a) + lambda / 2κ, b, -z)
end

function E_2_dξ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, λ)
    z = -conj(c) * ξ^2
    z_dξ = -2conj(c) * ξ
    return exp(z) *
           (U(b - conj(a) + lambda / 2κ, b, -z) - U_dz(b - conj(a) + lambda / 2κ, b, -z)) *
           z_dξ
end

function K_1(ξ, lambda, κ, ϵ, λ::CGLParams)
    return -Diagonal(
        SVector((ϵ - im) * J_P_1(ξ, lambda, κ, ϵ, λ), (ϵ + im) * J_P_2(ξ, lambda, κ, ϵ, λ)),
    ) / (1 + ϵ^2)
end

function K_2(ξ, lambda, κ, ϵ, λ::CGLParams)
    return Diagonal(
        SVector((ϵ - im) * J_E_1(ξ, lambda, κ, ϵ, λ), (ϵ + im) * J_E_2(ξ, lambda, κ, ϵ, λ)),
    ) / (1 + ϵ^2)
end

function K_1_dξ(ξ, lambda, κ, ϵ, λ::CGLParams)
    return -Diagonal(
        SVector(
            (ϵ - im) * J_P_1_dξ(ξ, lambda, κ, ϵ, λ),
            (ϵ + im) * J_P_2_dξ(ξ, lambda, κ, ϵ, λ),
        ),
    ) / (1 + ϵ^2)
end

function K_2_dξ(ξ, lambda, κ, ϵ, λ::CGLParams)
    return Diagonal(
        SVector(
            (ϵ - im) * J_E_1_dξ(ξ, lambda, κ, ϵ, λ),
            (ϵ + im) * J_E_2_dξ(ξ, lambda, κ, ϵ, λ),
        ),
    ) / (1 + ϵ^2)
end

function I_N(Q_hat, λ::CGLParams{T}) where {T}
    @assert isone(λ.σ)
    @assert iszero(λ.δ)

    return SMatrix{2,2}(
        2im * abs2(Q_hat),
        conj(-im * Q_hat^2),
        -im * Q_hat^2,
        -2im * abs2(Q_hat),
    )
end

function I_N_dξ(Q_hat, Q_hat_dξ, λ::CGLParams{T}) where {T}
    @assert isone(λ.σ)
    @assert iszero(λ.δ)

    return SMatrix{2,2}(
        4im * (real(Q_hat) * real(Q_hat_dξ) + imag(Q_hat) * imag(Q_hat_dξ)),
        conj(-2im * Q_hat * Q_hat_dξ),
        -2im * Q_hat * Q_hat_dξ,
        -4im * (real(Q_hat) * real(Q_hat_dξ) + imag(Q_hat) * imag(Q_hat_dξ)),
    )
end

# The rest of these functions are only used for testing

function W_1(ξ, lambda, κ, ϵ, λ::CGLParams)
    a, b, c = _abc(κ, ϵ, λ)
    z = -c * ξ^2
    sgn = c isa AcbSeries ? sign(imag(c[0])) : sign(imag(c))
    return -2c * exp(-sgn * im * (b - a + lambda / 2κ) * π) * ξ * z^-b * exp(z)
end

function W_2(ξ, lambda, κ, ϵ, λ::CGLParams)
    a, b, c = _abc(κ, ϵ, λ)
    z = -conj(c) * ξ^2
    sgn = c isa AcbSeries ? sign(imag(c[0])) : sign(imag(c))
    return -2conj(c) * exp(sgn * im * (b - conj(a) + lambda / 2κ) * π) * ξ * z^-b * exp(z)
end

function B_W_1(lambda, κ, ϵ, λ::CGLParams)
    a, b, c = _abc(κ, ϵ, λ)
    sgn = c isa AcbSeries ? sign(imag(c[0])) : sign(imag(c))
    return 1 // 2 * exp(sgn * im * (b - a + lambda / 2κ) * π) * (-c)^(b - 1)
end

function B_W_2(lambda, κ, ϵ, λ::CGLParams)
    a, b, c = _abc(κ, ϵ, λ)
    sgn = c isa AcbSeries ? sign(imag(c[0])) : sign(imag(c))
    return 1 // 2 * exp(-sgn * im * (b - conj(a) + lambda / 2κ) * π) * (-conj(c))^(b - 1)
end

function J_P_1(ξ, lambda, κ, ϵ, λ::CGLParams)
    (; d) = λ
    c = _c(κ, ϵ, λ)
    return B_W_1(lambda, κ, ϵ, λ) * P_1(ξ, lambda, κ, ϵ, λ) * exp(c * ξ^2) * ξ^(d - 1)
end

function J_P_2(ξ, lambda, κ, ϵ, λ::CGLParams)
    (; d) = λ
    c = _c(κ, ϵ, λ)
    return B_W_2(lambda, κ, ϵ, λ) * P_2(ξ, lambda, κ, ϵ, λ) * exp(conj(c) * ξ^2) * ξ^(d - 1)
end

function J_E_1(ξ, lambda, κ, ϵ, λ::CGLParams)
    (; d) = λ
    c = _c(κ, ϵ, λ)
    return B_W_1(lambda, κ, ϵ, λ) * E_1(ξ, lambda, κ, ϵ, λ) * exp(c * ξ^2) * ξ^(d - 1)
end

function J_E_2(ξ, lambda, κ, ϵ, λ::CGLParams)
    (; d) = λ
    c = _c(κ, ϵ, λ)
    return B_W_2(lambda, κ, ϵ, λ) * E_2(ξ, lambda, κ, ϵ, λ) * exp(conj(c) * ξ^2) * ξ^(d - 1)
end

function J_P_1_dξ(ξ, lambda, κ, ϵ, λ::CGLParams)
    (; d) = λ
    c = _c(κ, ϵ, λ)
    return B_W_1(lambda, κ, ϵ, λ) *
           (
               P_1_dξ(ξ, lambda, κ, ϵ, λ) * ξ^-1 +
               2c * P_1(ξ, lambda, κ, ϵ, λ) +
               (d - 1) * P_1(ξ, lambda, κ, ϵ, λ) * ξ^-2
           ) *
           exp(c * ξ^2) *
           ξ^d
end

function J_P_2_dξ(ξ, lambda, κ, ϵ, λ::CGLParams)
    (; d) = λ
    c = _c(κ, ϵ, λ)
    return B_W_2(lambda, κ, ϵ, λ) *
           (
               P_2_dξ(ξ, lambda, κ, ϵ, λ) * ξ^-1 +
               2conj(c) * P_2(ξ, lambda, κ, ϵ, λ) +
               (d - 1) * P_2(ξ, lambda, κ, ϵ, λ) * ξ^-2
           ) *
           exp(conj(c) * ξ^2) *
           ξ^d
end

function J_E_1_dξ(ξ, lambda, κ, ϵ, λ::CGLParams)
    (; d) = λ
    c = _c(κ, ϵ, λ)
    return B_W_1(lambda, κ, ϵ, λ) *
           (
               E_1_dξ(ξ, lambda, κ, ϵ, λ) * ξ^-1 +
               2c * E_1(ξ, lambda, κ, ϵ, λ) +
               (d - 1) * E_1(ξ, lambda, κ, ϵ, λ) * ξ^-2
           ) *
           exp(c * ξ^2) *
           ξ^d
end

function J_E_2_dξ(ξ, lambda, κ, ϵ, λ::CGLParams)
    (; d) = λ
    c = _c(κ, ϵ, λ)
    return B_W_2(lambda, κ, ϵ, λ) *
           (
               E_2_dξ(ξ, lambda, κ, ϵ, λ) * ξ^-1 +
               2conj(c) * E_2(ξ, lambda, κ, ϵ, λ) +
               (d - 1) * E_2(ξ, lambda, κ, ϵ, λ) * ξ^-2
           ) *
           exp(conj(c) * ξ^2) *
           ξ^d
end
