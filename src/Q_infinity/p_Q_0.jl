# TODO: Add tests and documentation for this function
function integral_J_E_P(κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    a, b, c = _abc(κ, ϵ, λ)
    (; d, σ) = λ

    I_U = zero(Acb)

    n = 4

    p_U = map(0:(n-1)) do k
        rising(a, k) * rising(a - b + 1, k) / (factorial(k) * (-c)^k)
    end
    p_U_bma = map(0:(n-1)) do k
        rising(b - a, k) * rising(-a + 1, k) / (factorial(k) * c^k)
    end

    C_R_U = CGL2.C_R_U(n, a, b, c * ξ₁^2)
    C_R_U_bma = CGL2.C_R_U(n, b - a, b, -c * ξ₁^2)

    for i_1 = 0:n
        for i_2 = 0:n
            for i_3 = 0:n
                for i_4 = 0:n
                    # Exponent after integration
                    exponent = -2 / σ - 2(i_1 + i_2 + i_3 + i_4)

                    if i_1 == n || i_2 == n || i_3 == n || i_4 == n
                        # At least one of the terms is a remainder
                        # term. We only bound the absolute value of
                        # the integral.
                        coefficient = Arb(1)
                        coefficient *= i_1 == n ? C_R_U_bma : abs(p_U_bma[i_1+1])
                        coefficient *= i_2 == n ? C_R_U : abs(p_U[i_2+1])
                        coefficient *= i_3 == n ? C_R_U : abs(p_U[i_3+1])
                        coefficient *= i_4 == n ? C_R_U : abs(p_U[i_4+1])

                        term_bound = coefficient * ξ₁^exponent / abs(exponent)

                        term = add_error(Acb(0), term_bound)

                        I_U += term
                    else
                        # None of the terms are remainder terms, we
                        # compute an enclosure of the integral.
                        I_U +=
                            p_U_bma[i_1+1] *
                            p_U[i_2+1] *
                            p_U[i_3+1] *
                            conj(p_U[i_4+1]) *
                            ξ₁^exponent / abs(exponent)
                    end
                end
            end
        end
    end

    return B_W(κ, ϵ, λ) * (-c)^(a - b) * abs(c^-a)^2 * c^-a * I_U
end

function I_E_infty_enclosure(γ::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    a, b, c = _abc(κ, ϵ, λ)
    v = Arb("0.001")
    (; σ) = λ

    CU = UBounds(a, b, c, ξ₁)
    C = FunctionBounds(κ, ϵ, ξ₁, λ, CU)
    norms = NormBounds(γ, κ, ϵ, ξ₁, v, λ, C)

    I_E_main = abs(γ)^2 * γ * integral_J_E_P(κ, ϵ, ξ₁, λ)

    C_R_Q =
        (C.P * C_I_E(κ, ϵ, ξ₁, v, λ, C) + C.E * C_I_P(κ, ϵ, ξ₁, v, λ, C)) *
        norms.Q^(2σ + 1) *
        ξ₁^((2σ + 1) * v - 2)
    R_I_E_bound =
        C.J_E * (3abs(γ)^2 * C.P^2 * C_R_Q + 3abs(γ) * C.P * C_R_Q^2 + C_R_Q^3) / (2 / σ) *
        ξ₁^(-2 / σ)
    R_I_E = add_error(zero(Acb), R_I_E_bound)

    return I_E_main + R_I_E
end

function p_Q_0(γ::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    a, b, c = _abc(κ, ϵ, λ)
    return c^-a * (γ + I_E_infty_enclosure(γ, κ, ϵ, ξ₁, λ))
end

# IMPROVE: Do we need this function?
p_Q_0(γ::ComplexF64, κ::Float64, ϵ::Float64, ξ₁::Float64, λ::CGLParams{Float64}) =
    Float64(p_Q_0(Acb(γ), Arb(κ), Arb(ϵ), Arb(ξ₁), CGLParams{Arb}(λ)))
