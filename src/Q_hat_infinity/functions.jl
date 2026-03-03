struct FunctionEnclosures_hat
    P_hat::Acb
    P_hat_dξ::Acb
    E_hat::Acb
    E_hat_dξ::Acb
    J_P_hat::Acb
    J_E_hat::Acb

    FunctionEnclosures_hat() = new(
        indeterminate(Acb),
        indeterminate(Acb),
        indeterminate(Acb),
        indeterminate(Acb),
        indeterminate(Acb),
        indeterminate(Acb),
    )
end

function FunctionEnclosures_hat(κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    F_hat = FunctionEnclosures_hat()

    F_hat.P_hat[] = P_hat(ξ₁, κ, ϵ, λ)
    F_hat.P_hat_dξ[] = P_hat_dξ(ξ₁, κ, ϵ, λ)

    F_hat.E_hat[] = E_hat(ξ₁, κ, ϵ, λ)
    F_hat.E_hat_dξ[] = E_hat_dξ(ξ₁, κ, ϵ, λ)

    F_hat.J_P_hat[] = J_P_hat(ξ₁, κ, ϵ, λ)
    F_hat.J_E_hat[] = J_E_hat(ξ₁, κ, ϵ, λ)

    return F_hat
end

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

function W_hat(ξ, κ, ϵ, λ::CGLParams)
    a, b, c = _abc(κ, ϵ, λ)

    z = -c * ξ^2

    sgn = if c isa AcbSeries
        sign(imag(c[0]))
    else
        sign(imag(c))
    end

    return -2c * exp(-sgn * im * (b - a) * π) * ξ * z^-b * exp(z)
end

function B_W_hat(κ, ϵ, λ::CGLParams)
    (; δ) = λ
    a, b, c = _abc(κ, ϵ, λ)

    sgn = if c isa Acb
        sign(Arblib.imagref(c))
    else
        sign(imag(c))
    end

    return _complex(δ, -1) / κ * exp(sgn * im * (b - a) * π) * (-c)^b
end

function J_P_hat(ξ, κ, ϵ, λ::CGLParams)
    c = _c(κ, ϵ, λ)

    return B_W_hat(κ, ϵ, λ) * P_hat(ξ, κ, ϵ, λ) * exp(c * ξ^2) * ξ^(λ.d - 1)
end

function J_E_hat(ξ, κ, ϵ, λ::CGLParams)
    c = _c(κ, ϵ, λ)

    return B_W_hat(κ, ϵ, λ) * E_hat(ξ, κ, ϵ, λ) * exp(c * ξ^2) * ξ^(λ.d - 1)
end
