"""
    FunctionEnclosures_Y(κ, ϵ, ξ₁, Λ)

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
    Λ::CGLParams{Arb},
)
    # Compute enclosure of forward solution
    Q_hat, Q_hat_dξ = Q_hat_infinity(γ₁, γ₂, κ, ϵ, ξ₁, Λ)

    return FunctionEnclosures_Y(
        Diagonal(SVector(E_1(ξ₁, lambda, κ, ϵ, Λ), E_2(ξ₁, lambda, κ, ϵ, Λ))),
        Diagonal(SVector(E_1_dξ(ξ₁, lambda, κ, ϵ, Λ), E_2_dξ(ξ₁, lambda, κ, ϵ, Λ))),
        Diagonal(SVector(P_1(ξ₁, lambda, κ, ϵ, Λ), P_2(ξ₁, lambda, κ, ϵ, Λ))),
        Diagonal(SVector(P_1_dξ(ξ₁, lambda, κ, ϵ, Λ), P_2_dξ(ξ₁, lambda, κ, ϵ, Λ))),
        K_1(ξ₁, lambda, κ, ϵ, Λ),
        K_2(ξ₁, lambda, κ, ϵ, Λ),
        K_1_dξ(ξ₁, lambda, κ, ϵ, Λ),
        K_2_dξ(ξ₁, lambda, κ, ϵ, Λ),
        I_N(Q_hat, Λ),
        I_N_dξ(Q_hat, Q_hat_dξ, Λ),
    )
end

function P_1(ξ, lambda, κ, ϵ, Λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, Λ)
    z = -c * ξ^2
    return U(a - lambda / 2κ, b, z)
end

function P_1_dξ(ξ, lambda, κ, ϵ, Λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, Λ)
    z = -c * ξ^2
    z_dξ = -2c * ξ
    return U_dz(a - lambda / 2κ, b, z) * z_dξ
end

function P_2(ξ, lambda, κ, ϵ, Λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, Λ)
    z = -conj(c) * ξ^2
    return U(conj(a) - lambda / 2κ, b, z)
end

function P_2_dξ(ξ, lambda, κ, ϵ, Λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, Λ)
    z = -conj(c) * ξ^2
    z_dξ = -2conj(c) * ξ
    return U_dz(conj(a) - lambda / 2κ, b, z) * z_dξ
end

function E_1(ξ, lambda, κ, ϵ, Λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, Λ)
    z = -c * ξ^2
    return exp(z) * U(b - a + lambda / 2κ, b, -z)
end

function E_1_dξ(ξ, lambda, κ, ϵ, Λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, Λ)
    z = -c * ξ^2
    z_dξ = -2c * ξ
    return exp(z) *
           (U(b - a + lambda / 2κ, b, -z) - U_dz(b - a + lambda / 2κ, b, -z)) *
           z_dξ
end

function E_2(ξ, lambda, κ, ϵ, Λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, Λ)
    z = -conj(c) * ξ^2
    return exp(z) * U(b - conj(a) + lambda / 2κ, b, -z)
end

function E_2_dξ(ξ, lambda, κ, ϵ, Λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, Λ)
    z = -conj(c) * ξ^2
    z_dξ = -2conj(c) * ξ
    return exp(z) *
           (U(b - conj(a) + lambda / 2κ, b, -z) - U_dz(b - conj(a) + lambda / 2κ, b, -z)) *
           z_dξ
end

function K_1(ξ, lambda, κ, ϵ, Λ::CGLParams)
    return -Diagonal(
        SVector((ϵ - im) * J_P_1(ξ, lambda, κ, ϵ, Λ), (ϵ + im) * J_P_2(ξ, lambda, κ, ϵ, Λ)),
    ) / (1 + ϵ^2)
end

function K_2(ξ, lambda, κ, ϵ, Λ::CGLParams)
    return Diagonal(
        SVector((ϵ - im) * J_E_1(ξ, lambda, κ, ϵ, Λ), (ϵ + im) * J_E_2(ξ, lambda, κ, ϵ, Λ)),
    ) / (1 + ϵ^2)
end

function K_1_dξ(ξ, lambda, κ, ϵ, Λ::CGLParams)
    return -Diagonal(
        SVector(
            (ϵ - im) * J_P_1_dξ(ξ, lambda, κ, ϵ, Λ),
            (ϵ + im) * J_P_2_dξ(ξ, lambda, κ, ϵ, Λ),
        ),
    ) / (1 + ϵ^2)
end

function K_2_dξ(ξ, lambda, κ, ϵ, Λ::CGLParams)
    return Diagonal(
        SVector(
            (ϵ - im) * J_E_1_dξ(ξ, lambda, κ, ϵ, Λ),
            (ϵ + im) * J_E_2_dξ(ξ, lambda, κ, ϵ, Λ),
        ),
    ) / (1 + ϵ^2)
end

function I_N(Q_hat, Λ::CGLParams{T}) where {T}
    @assert isone(Λ.σ)
    @assert iszero(Λ.δ)

    return SMatrix{2,2}(
        2im * abs2(Q_hat),
        conj(-im * Q_hat^2),
        -im * Q_hat^2,
        -2im * abs2(Q_hat),
    )
end

function I_N_dξ(Q_hat, Q_hat_dξ, Λ::CGLParams{T}) where {T}
    @assert isone(Λ.σ)
    @assert iszero(Λ.δ)

    return SMatrix{2,2}(
        4im * (real(Q_hat) * real(Q_hat_dξ) + imag(Q_hat) * imag(Q_hat_dξ)),
        conj(-2im * Q_hat * Q_hat_dξ),
        -2im * Q_hat * Q_hat_dξ,
        -4im * (real(Q_hat) * real(Q_hat_dξ) + imag(Q_hat) * imag(Q_hat_dξ)),
    )
end

# The rest of these functions are only used for testing

function W_1(ξ, lambda, κ, ϵ, Λ::CGLParams)
    a, b, c = _abc(κ, ϵ, Λ)
    z = -c * ξ^2
    sgn = c isa AcbSeries ? sign(imag(c[0])) : sign(imag(c))
    return -2c * exp(-sgn * im * (b - a + lambda / 2κ) * π) * ξ * z^-b * exp(z)
end

function W_2(ξ, lambda, κ, ϵ, Λ::CGLParams)
    a, b, c = _abc(κ, ϵ, Λ)
    z = -conj(c) * ξ^2
    sgn = c isa AcbSeries ? sign(imag(c[0])) : sign(imag(c))
    return -2conj(c) * exp(sgn * im * (b - conj(a) + lambda / 2κ) * π) * ξ * z^-b * exp(z)
end

function B_W_1(lambda, κ, ϵ, Λ::CGLParams)
    a, b, c = _abc(κ, ϵ, Λ)
    sgn = c isa AcbSeries ? sign(imag(c[0])) : sign(imag(c))
    return 1 // 2 * exp(sgn * im * (b - a + lambda / 2κ) * π) * (-c)^(b - 1)
end

function B_W_2(lambda, κ, ϵ, Λ::CGLParams)
    a, b, c = _abc(κ, ϵ, Λ)
    sgn = c isa AcbSeries ? sign(imag(c[0])) : sign(imag(c))
    return 1 // 2 * exp(-sgn * im * (b - conj(a) + lambda / 2κ) * π) * (-conj(c))^(b - 1)
end

function J_P_1(ξ, lambda, κ, ϵ, Λ::CGLParams)
    (; d) = Λ
    c = _c(κ, ϵ, Λ)
    return B_W_1(lambda, κ, ϵ, Λ) * P_1(ξ, lambda, κ, ϵ, Λ) * exp(c * ξ^2) * ξ^(d - 1)
end

function J_P_2(ξ, lambda, κ, ϵ, Λ::CGLParams)
    (; d) = Λ
    c = _c(κ, ϵ, Λ)
    return B_W_2(lambda, κ, ϵ, Λ) * P_2(ξ, lambda, κ, ϵ, Λ) * exp(conj(c) * ξ^2) * ξ^(d - 1)
end

function J_E_1(ξ, lambda, κ, ϵ, Λ::CGLParams)
    (; d) = Λ
    c = _c(κ, ϵ, Λ)
    return B_W_1(lambda, κ, ϵ, Λ) * E_1(ξ, lambda, κ, ϵ, Λ) * exp(c * ξ^2) * ξ^(d - 1)
end

function J_E_2(ξ, lambda, κ, ϵ, Λ::CGLParams)
    (; d) = Λ
    c = _c(κ, ϵ, Λ)
    return B_W_2(lambda, κ, ϵ, Λ) * E_2(ξ, lambda, κ, ϵ, Λ) * exp(conj(c) * ξ^2) * ξ^(d - 1)
end

function J_P_1_dξ(ξ, lambda, κ, ϵ, Λ::CGLParams)
    (; d) = Λ
    c = _c(κ, ϵ, Λ)
    return B_W_1(lambda, κ, ϵ, Λ) *
           (
               P_1_dξ(ξ, lambda, κ, ϵ, Λ) * ξ^-1 +
               2c * P_1(ξ, lambda, κ, ϵ, Λ) +
               (d - 1) * P_1(ξ, lambda, κ, ϵ, Λ) * ξ^-2
           ) *
           exp(c * ξ^2) *
           ξ^d
end

function J_P_2_dξ(ξ, lambda, κ, ϵ, Λ::CGLParams)
    (; d) = Λ
    c = _c(κ, ϵ, Λ)
    return B_W_2(lambda, κ, ϵ, Λ) *
           (
               P_2_dξ(ξ, lambda, κ, ϵ, Λ) * ξ^-1 +
               2conj(c) * P_2(ξ, lambda, κ, ϵ, Λ) +
               (d - 1) * P_2(ξ, lambda, κ, ϵ, Λ) * ξ^-2
           ) *
           exp(conj(c) * ξ^2) *
           ξ^d
end

function J_E_1_dξ(ξ, lambda, κ, ϵ, Λ::CGLParams)
    (; d) = Λ
    c = _c(κ, ϵ, Λ)
    return B_W_1(lambda, κ, ϵ, Λ) *
           (
               E_1_dξ(ξ, lambda, κ, ϵ, Λ) * ξ^-1 +
               2c * E_1(ξ, lambda, κ, ϵ, Λ) +
               (d - 1) * E_1(ξ, lambda, κ, ϵ, Λ) * ξ^-2
           ) *
           exp(c * ξ^2) *
           ξ^d
end

function J_E_2_dξ(ξ, lambda, κ, ϵ, Λ::CGLParams)
    (; d) = Λ
    c = _c(κ, ϵ, Λ)
    return B_W_2(lambda, κ, ϵ, Λ) *
           (
               E_2_dξ(ξ, lambda, κ, ϵ, Λ) * ξ^-1 +
               2conj(c) * E_2(ξ, lambda, κ, ϵ, Λ) +
               (d - 1) * E_2(ξ, lambda, κ, ϵ, Λ) * ξ^-2
           ) *
           exp(conj(c) * ξ^2) *
           ξ^d
end
