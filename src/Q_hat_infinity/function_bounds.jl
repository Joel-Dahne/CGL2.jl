function B_W_hat(κ, ϵ, λ::CGLParams)
    (; δ) = λ
    a, b, c = _abc(κ, ϵ, λ)

    sgn = if c isa AcbSeries
        sign(Arblib.imagref(Arblib.ref(c, 0)))
    elseif c isa Acb
        sign(Arblib.imagref(c))
    else
        sign(imag(c))
    end

    return _complex(-δ, 1) / (-κ) * exp(sgn * im * (b - a) * π) * (-c)^b
end

function C_P_hat(κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    a, b, c = _abc(κ, ϵ, λ)

    return C_U(a, b, -c * ξ₁^2) * abs((-c)^-a)
end

function C_E_hat(κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    a, b, c = _abc(κ, ϵ, λ)

    return C_U(b - a, b, c * ξ₁^2) * abs(c^(a - b))
end

C_J_P_hat(κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}) =
    abs(B_W_hat(κ, ϵ, λ)) * C_P_hat(κ, ϵ, ξ₁, λ)

C_J_E_hat(κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}) =
    abs(B_W_hat(κ, ϵ, λ)) * C_E_hat(κ, ϵ, ξ₁, λ)
