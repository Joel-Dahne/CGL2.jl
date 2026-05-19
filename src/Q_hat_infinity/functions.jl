"""
    FunctionEnclosures_hat(κ, ϵ, ξ₁, Λ)

Contains enclosures of the functions

- [`P_hat`](@ref)
- [`P_hat_dξ`](@ref)
- [`E_hat`](@ref)
- [`E_hat_dξ`](@ref)
- [`J_E_hat`](@ref)
- [`J_P_hat`](@ref)

when evaluated at `ξ₁`.
"""
struct FunctionEnclosures_hat
    P_hat::Acb
    P_hat_dξ::Acb
    E_hat::Acb
    E_hat_dξ::Acb
    J_E_hat::Acb
    J_P_hat::Acb

    FunctionEnclosures_hat() = new(
        indeterminate(Acb),
        indeterminate(Acb),
        indeterminate(Acb),
        indeterminate(Acb),
        indeterminate(Acb),
        indeterminate(Acb),
    )
end

function FunctionEnclosures_hat(κ::Arb, ϵ::Arb, ξ₁::Arb, Λ::CGLParams{Arb})
    F_hat = FunctionEnclosures_hat()

    F_hat.P_hat[] = P_hat(ξ₁, κ, ϵ, Λ)
    F_hat.P_hat_dξ[] = P_hat_dξ(ξ₁, κ, ϵ, Λ)

    F_hat.E_hat[] = E_hat(ξ₁, κ, ϵ, Λ)
    F_hat.E_hat_dξ[] = E_hat_dξ(ξ₁, κ, ϵ, Λ)

    F_hat.J_E_hat[] = J_E_hat(ξ₁, κ, ϵ, Λ)
    F_hat.J_P_hat[] = J_P_hat(ξ₁, κ, ϵ, Λ)

    return F_hat
end

function P_hat(ξ, κ, ϵ, Λ::CGLParams)
    a, b, c = _abc(κ, ϵ, Λ)
    z = -c * ξ^2
    return U(a, b, z)
end

function P_hat_dξ(ξ, κ, ϵ, Λ::CGLParams)
    a, b, c = _abc(κ, ϵ, Λ)
    z = -c * ξ^2
    z_dξ = -2c * ξ
    return U_dz(a, b, z) * z_dξ
end

function E_hat(ξ, κ, ϵ, Λ::CGLParams)
    a, b, c = _abc(κ, ϵ, Λ)
    z = -c * ξ^2
    return exp(z) * U(b - a, b, -z)
end

function E_hat_dξ(ξ, κ, ϵ, Λ::CGLParams)
    a, b, c = _abc(κ, ϵ, Λ)
    z = -c * ξ^2
    z_dξ = -2c * ξ
    return exp(z) * (U(b - a, b, -z) - U_dz(b - a, b, -z)) * z_dξ
end

# This is only used for testing.
function W_hat(ξ, κ, ϵ, Λ::CGLParams)
    a, b, c = _abc(κ, ϵ, Λ)
    z = -c * ξ^2
    sgn = c isa AcbSeries ? sign(imag(c[0])) : sign(imag(c))
    return -2c * exp(-sgn * im * (b - a) * π) * ξ * z^-b * exp(z)
end

function B_W_hat(κ, ϵ, Λ::CGLParams)
    (; δ) = Λ
    a, b, c = _abc(κ, ϵ, Λ)
    sgn = c isa AcbSeries ? sign(imag(c[0])) : sign(imag(c))
    return -_complex(δ, -1) / κ * exp(sgn * im * (b - a) * π) * (-c)^b
end

function J_E_hat(ξ, κ, ϵ, Λ::CGLParams)
    c = _c(κ, ϵ, Λ)
    return B_W_hat(κ, ϵ, Λ) * E_hat(ξ, κ, ϵ, Λ) * exp(c * ξ^2) * ξ^(Λ.d - 1)
end

function J_P_hat(ξ, κ, ϵ, Λ::CGLParams)
    c = _c(κ, ϵ, Λ)
    return B_W_hat(κ, ϵ, Λ) * P_hat(ξ, κ, ϵ, Λ) * exp(c * ξ^2) * ξ^(Λ.d - 1)
end
