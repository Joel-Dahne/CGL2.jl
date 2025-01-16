function P_hat(ξ, κ, ϵ, λ::CGLParams)
    a, b, c = _abc(κ, ϵ, λ)

    z = -c * ξ^2

    return U(a, b, z)
end

function P_hat_dξ(ξ, κ, ϵ, λ::CGLParams)
    a, b, c = _abc(κ, ϵ, λ)

    z = -c * ξ^2
    z_dξ = -2c * ξ

    return U_dz(a, b, z) * z_dξ
end

function E_hat(ξ, κ, ϵ, λ::CGLParams)
    a, b, c = _abc(κ, ϵ, λ)

    z = -c * ξ^2

    return exp(z) * U(b - a, b, -z)
end

function E_hat_dξ(ξ, κ, ϵ, λ::CGLParams)
    a, b, c = _abc(κ, ϵ, λ)

    z = -c * ξ^2
    z_dξ = -2c * ξ

    return exp(z) * (U(b - a, b, -z) - U_dz(b - a, b, -z)) * z_dξ
end
