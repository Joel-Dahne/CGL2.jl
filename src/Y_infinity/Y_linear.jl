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
