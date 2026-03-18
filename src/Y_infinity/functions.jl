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

    E_12_dλ::Diagonal{Acb,SVector{2,Acb}}
    E_12_dλ_dξ::Diagonal{Acb,SVector{2,Acb}}
    P_12_dλ::Diagonal{Acb,SVector{2,Acb}}
    P_12_dλ_dξ::Diagonal{Acb,SVector{2,Acb}}
    K_1_dλ::Diagonal{Acb,SVector{2,Acb}}
    K_2_dλ::Diagonal{Acb,SVector{2,Acb}}

    function FunctionEnclosures_Y(
        lambda::Acb,
        γ₁::Acb,
        γ₂::Acb,
        κ::Arb,
        ϵ::Arb,
        ξ₁::Arb,
        λ::CGLParams{Arb};
        include_dλ::Bool = false,
    )
        K1, K2 = K_1_2(ξ₁, lambda, κ, ϵ, λ)
        K1_dξ, K2_dξ = K_1_2_dξ(ξ₁, lambda, κ, ϵ, λ)
        K1_dλ, K2_dλ = K_1_2_dλ(ξ₁, lambda, κ, ϵ, λ)

        # Compute enclosure of forward solution
        Q_hat_ξ₁, dQ_hat_ξ₁ = Q_hat_infinity(γ₁, γ₂, κ, ϵ, ξ₁, λ)
        IN = I_N(Q_hat_ξ₁, λ)
        IN_dξ = I_N_dξ(Q_hat_ξ₁, dQ_hat_ξ₁, λ)

        F = new(
            Diagonal(SVector(E_1(ξ₁, lambda, κ, ϵ, λ), E_2(ξ₁, lambda, κ, ϵ, λ))),
            Diagonal(SVector(E_1_dξ(ξ₁, lambda, κ, ϵ, λ), E_2_dξ(ξ₁, lambda, κ, ϵ, λ))),
            Diagonal(SVector(P_1(ξ₁, lambda, κ, ϵ, λ), P_2(ξ₁, lambda, κ, ϵ, λ))),
            Diagonal(SVector(P_1_dξ(ξ₁, lambda, κ, ϵ, λ), P_2_dξ(ξ₁, lambda, κ, ϵ, λ))),
            K1,
            K2,
            K1_dξ,
            K2_dξ,
            IN,
            IN_dξ,
            Diagonal(SVector(E_1_dλ(ξ₁, lambda, κ, ϵ, λ), E_2_dλ(ξ₁, lambda, κ, ϵ, λ))),
            Diagonal(
                SVector(E_1_dλ_dξ(ξ₁, lambda, κ, ϵ, λ), E_2_dλ_dξ(ξ₁, lambda, κ, ϵ, λ)),
            ),
            Diagonal(SVector(P_1_dλ(ξ₁, lambda, κ, ϵ, λ), P_2_dλ(ξ₁, lambda, κ, ϵ, λ))),
            Diagonal(
                SVector(P_1_dλ_dξ(ξ₁, lambda, κ, ϵ, λ), P_2_dλ_dξ(ξ₁, lambda, κ, ϵ, λ)),
            ),
            K1_dλ,
            K2_dλ,
        )

        return F
    end
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

function P_1_dλ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, λ)
    z = -c * ξ^2
    return -U_da(a - lambda / 2κ, b, z) / 2κ
end

function P_1_dλ_dξ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, λ)
    z = -c * ξ^2
    z_dξ = -2c * ξ
    return -U_dzda(a - lambda / 2κ, b, z) / 2κ * z_dξ
end

function P_2_dλ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, λ)
    z = -conj(c) * ξ^2
    return -U_da(conj(a) - lambda / 2κ, b, z) / 2κ
end

function P_2_dλ_dξ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, λ)
    z = -conj(c) * ξ^2
    z_dξ = -2conj(c) * ξ
    return -U_dzda(conj(a) - lambda / 2κ, b, z) / 2κ * z_dξ
end

function E_1_dλ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, λ)
    z = -c * ξ^2
    return exp(z) * U_da(b - a + lambda / 2κ, b, -z) / 2κ
end

function E_1_dλ_dξ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, λ)
    z = -c * ξ^2
    z_dξ = -2c * ξ
    return exp(z) *
           (U_da(b - a + lambda / 2κ, b, -z) - U_dzda(b - a + lambda / 2κ, b, -z)) / 2κ *
           z_dξ
end

function E_2_dλ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, λ)
    z = -conj(c) * ξ^2
    return exp(z) * U_da(b - conj(a) + lambda / 2κ, b, -z) / 2κ
end

function E_2_dλ_dξ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, λ)
    z = -conj(c) * ξ^2
    z_dξ = -2conj(c) * ξ
    return exp(z) * (
        U_da(b - conj(a) + lambda / 2κ, b, -z) - U_dzda(b - conj(a) + lambda / 2κ, b, -z)
    ) / 2κ * z_dξ
end

function W_1(ξ, lambda, κ, ϵ, λ::CGLParams)
    a, b, c = _abc(κ, ϵ, λ)

    z = -c * ξ^2

    sgn = if c isa AcbSeries
        sign(imag(c[0]))
    else
        sign(imag(c))
    end

    return -2c * exp(-sgn * im * (b - a + lambda / 2κ) * π) * ξ * z^-b * exp(z)
end

function W_2(ξ, lambda, κ, ϵ, λ::CGLParams)
    a, b, c = _abc(κ, ϵ, λ)

    z = -conj(c) * ξ^2

    sgn = if c isa AcbSeries
        sign(imag(c[0]))
    else
        sign(imag(c))
    end

    return -2conj(c) * exp(sgn * im * (b - conj(a) + lambda / 2κ) * π) * ξ * z^-b * exp(z)
end

function B_W_1(lambda, κ, ϵ, λ::CGLParams)
    a, b, c = _abc(κ, ϵ, λ)

    sgn = if c isa AcbSeries
        sign(imag(c[0]))
    else
        sign(imag(c))
    end

    return 1 // 2 * exp(sgn * im * (b - a + lambda / 2κ) * π) * (-c)^(b - 1)
end

function B_W_2(lambda, κ, ϵ, λ::CGLParams)
    a, b, c = _abc(κ, ϵ, λ)

    sgn = if c isa AcbSeries
        sign(imag(c[0]))
    else
        sign(imag(c))
    end

    return 1 // 2 * exp(-sgn * im * (b - conj(a) + lambda / 2κ) * π) * (-conj(c))^(b - 1)
end

function B_W_1_dλ(lambda, κ, ϵ, λ::CGLParams)
    a, b, c = _abc(κ, ϵ, λ)

    sgn = if c isa AcbSeries
        sign(imag(c[0]))
    else
        sign(imag(c))
    end

    return 1 // 2 *
           (sgn * im / 2κ * π) *
           exp(sgn * im * (b - a + lambda / 2κ) * π) *
           (-c)^(b - 1)
end

function B_W_2_dλ(lambda, κ, ϵ, λ::CGLParams)
    a, b, c = _abc(κ, ϵ, λ)

    sgn = if c isa AcbSeries
        sign(imag(c[0]))
    else
        sign(imag(c))
    end

    return 1 // 2 *
           (-sgn * im / 2κ * π) *
           exp(-sgn * im * (b - conj(a) + lambda / 2κ) * π) *
           (-conj(c))^(b - 1)
end

function J_P_1(ξ, lambda, κ, ϵ, λ::CGLParams)
    (; d) = λ
    _, _, c = _abc(κ, ϵ, λ)
    return B_W_1(lambda, κ, ϵ, λ) * P_1(ξ, lambda, κ, ϵ, λ) * exp(c * ξ^2) * ξ^(d - 1)
end

function J_P_2(ξ, lambda, κ, ϵ, λ::CGLParams)
    (; d) = λ
    _, _, c = _abc(κ, ϵ, λ)
    return B_W_2(lambda, κ, ϵ, λ) * P_2(ξ, lambda, κ, ϵ, λ) * exp(conj(c) * ξ^2) * ξ^(d - 1)
end

function J_E_1(ξ, lambda, κ, ϵ, λ::CGLParams)
    (; d) = λ
    _, _, c = _abc(κ, ϵ, λ)
    return B_W_1(lambda, κ, ϵ, λ) * E_1(ξ, lambda, κ, ϵ, λ) * exp(c * ξ^2) * ξ^(d - 1)
end

function J_E_2(ξ, lambda, κ, ϵ, λ::CGLParams)
    (; d) = λ
    _, _, c = _abc(κ, ϵ, λ)
    return B_W_2(lambda, κ, ϵ, λ) * E_2(ξ, lambda, κ, ϵ, λ) * exp(conj(c) * ξ^2) * ξ^(d - 1)
end

function J_P_1_dξ(ξ, lambda, κ, ϵ, λ::CGLParams)
    (; d) = λ
    _, _, c = _abc(κ, ϵ, λ)
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
    _, _, c = _abc(κ, ϵ, λ)
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
    _, _, c = _abc(κ, ϵ, λ)
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
    _, _, c = _abc(κ, ϵ, λ)
    return B_W_2(lambda, κ, ϵ, λ) *
           (
               E_2_dξ(ξ, lambda, κ, ϵ, λ) * ξ^-1 +
               2conj(c) * E_2(ξ, lambda, κ, ϵ, λ) +
               (d - 1) * E_2(ξ, lambda, κ, ϵ, λ) * ξ^-2
           ) *
           exp(conj(c) * ξ^2) *
           ξ^d
end

function J_P_1_dλ(ξ, lambda, κ, ϵ, λ::CGLParams)
    (; d) = λ
    _, _, c = _abc(κ, ϵ, λ)
    return (
               B_W_1_dλ(lambda, κ, ϵ, λ) * P_1(ξ, lambda, κ, ϵ, λ) +
               B_W_1(lambda, κ, ϵ, λ) * P_1_dλ(ξ, lambda, κ, ϵ, λ)
           ) *
           exp(c * ξ^2) *
           ξ^(d - 1)
end

function J_P_2_dλ(ξ, lambda, κ, ϵ, λ::CGLParams)
    (; d) = λ
    _, _, c = _abc(κ, ϵ, λ)
    return (
               B_W_2_dλ(lambda, κ, ϵ, λ) * P_2(ξ, lambda, κ, ϵ, λ) +
               B_W_2(lambda, κ, ϵ, λ) * P_2_dλ(ξ, lambda, κ, ϵ, λ)
           ) *
           exp(conj(c) * ξ^2) *
           ξ^(d - 1)
end

function J_E_1_dλ(ξ, lambda, κ, ϵ, λ::CGLParams)
    (; d) = λ
    _, _, c = _abc(κ, ϵ, λ)
    return (
               B_W_1_dλ(lambda, κ, ϵ, λ) * E_1(ξ, lambda, κ, ϵ, λ) +
               B_W_1(lambda, κ, ϵ, λ) * E_1_dλ(ξ, lambda, κ, ϵ, λ)
           ) *
           exp(c * ξ^2) *
           ξ^(d - 1)
end

function J_E_2_dλ(ξ, lambda, κ, ϵ, λ::CGLParams)
    (; d) = λ
    _, _, c = _abc(κ, ϵ, λ)
    return (
               B_W_2_dλ(lambda, κ, ϵ, λ) * E_2(ξ, lambda, κ, ϵ, λ) +
               B_W_2(lambda, κ, ϵ, λ) * E_2_dλ(ξ, lambda, κ, ϵ, λ)
           ) *
           exp(conj(c) * ξ^2) *
           ξ^(d - 1)
end

function I_N(Q_hat_ξ, λ::CGLParams{T}) where {T}
    @assert isone(λ.σ)
    @assert iszero(λ.δ)

    return SMatrix{2,2}(
        2im * abs2(Q_hat_ξ),
        conj(-im * Q_hat_ξ^2),
        -im * Q_hat_ξ^2,
        -2im * abs2(Q_hat_ξ),
    )
end

function I_N_dξ(Q_hat_ξ, dQ_hat_ξ, λ::CGLParams{T}) where {T}
    @assert isone(λ.σ)
    @assert iszero(λ.δ)

    return SMatrix{2,2}(
        4im * (real(Q_hat_ξ) * real(dQ_hat_ξ) + imag(Q_hat_ξ) * imag(dQ_hat_ξ)),
        conj(-2im * Q_hat_ξ * dQ_hat_ξ),
        -2im * Q_hat_ξ * dQ_hat_ξ,
        -4im * (real(Q_hat_ξ) * real(dQ_hat_ξ) + imag(Q_hat_ξ) * imag(dQ_hat_ξ)),
    )
end

function K_1_2(ξ, lambda, κ, ϵ, λ::CGLParams)
    K_1 =
        -Diagonal(
            SVector(
                (ϵ - im) * J_P_1(ξ, lambda, κ, ϵ, λ),
                (ϵ + im) * J_P_2(ξ, lambda, κ, ϵ, λ),
            ),
        ) / (1 + ϵ^2)
    K_2 =
        Diagonal(
            SVector(
                (ϵ - im) * J_E_1(ξ, lambda, κ, ϵ, λ),
                (ϵ + im) * J_E_2(ξ, lambda, κ, ϵ, λ),
            ),
        ) / (1 + ϵ^2)
    return K_1, K_2
end

function K_1_2_dξ(ξ, lambda, κ, ϵ, λ::CGLParams)
    K_1_dξ =
        -Diagonal(
            SVector(
                (ϵ - im) * J_P_1_dξ(ξ, lambda, κ, ϵ, λ),
                (ϵ + im) * J_P_2_dξ(ξ, lambda, κ, ϵ, λ),
            ),
        ) / (1 + ϵ^2)
    K_2_dξ =
        Diagonal(
            SVector(
                (ϵ - im) * J_E_1_dξ(ξ, lambda, κ, ϵ, λ),
                (ϵ + im) * J_E_2_dξ(ξ, lambda, κ, ϵ, λ),
            ),
        ) / (1 + ϵ^2)
    return K_1_dξ, K_2_dξ
end

function K_1_2_dλ(ξ, lambda, κ, ϵ, λ::CGLParams)
    K_1_dλ =
        -Diagonal(
            SVector(
                (ϵ - im) * J_P_1_dλ(ξ, lambda, κ, ϵ, λ),
                (ϵ + im) * J_P_2_dλ(ξ, lambda, κ, ϵ, λ),
            ),
        ) / (1 + ϵ^2)
    K_2_dλ =
        Diagonal(
            SVector(
                (ϵ - im) * J_E_1_dλ(ξ, lambda, κ, ϵ, λ),
                (ϵ + im) * J_E_2_dλ(ξ, lambda, κ, ϵ, λ),
            ),
        ) / (1 + ϵ^2)
    return K_1_dλ, K_2_dλ
end
