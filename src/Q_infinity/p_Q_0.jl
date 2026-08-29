"""
    integral_J_E_P(κ::Arb, ϵ::Arb, ξ₁::Arb, Λ::CGLParams{Arb})

Compute an enclosure of the integral in `η` from `ξ₁` to infinity of
the function

```
J_E(η) * abs(P(η))^2 * P(η)
```

The approach is based on Lemma REF(lemma:I_E_infty-integral). The lemma
reduces it to the integral with the integrand given by a product of
asymptotic series with remainder terms. To enclose the integral we
compute the coefficients in the series and bounds for the remainder
terms.

The multiplication of the series is handled using nested for loops.
For terms in the product where at least one of the factors is a
remainder term we compute a bound for the absolute value. For terms
were all factors are from the series we enclose the integral by
integrating it explicitly.
"""
function integral_J_E_P(κ::Arb, ϵ::Arb, ξ₁::Arb, Λ::CGLParams{Arb})
    a, b, c = _abc(κ, ϵ, Λ)
    (; d, σ) = Λ

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

    return B_W(κ, ϵ, Λ) * (-c)^(a - b) * abs(c^-a)^2 * c^-a * I_U
end

"""
    I_E_infty_enclosure(γ, κ, ϵ, ξ₁, Λ)

Compute an enclosure of ``I_{E,∞}`` from Lemma REF(lemma:p_Q_0). The
enclosure is based on Lemma REF(lemma:I_E_infty).
"""
function I_E_infty_enclosure(γ::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, Λ::CGLParams{Arb})
    a, b, c = _abc(κ, ϵ, Λ)
    (; σ) = Λ

    v = Arb(0) # v is zero in the entire section
    CU = UBounds(a, b, c, ξ₁)
    C = FunctionBounds(κ, ϵ, ξ₁, Λ, CU)
    CI = IBounds(κ, ϵ, ξ₁, v, Λ, C)
    norms = NormBounds(γ, κ, ϵ, ξ₁, v, Λ, C, CI)

    I_E_main = abs(γ)^2 * γ * integral_J_E_P(κ, ϵ, ξ₁, Λ)

    C_R_1 = C.P * C.J_E / 2 * norms.Q^(2σ + 1)
    C_R_2 = C.E * CI.I_P * norms.Q^(2σ + 1)
    A_max = max(C_R_1, C_R_2) * ξ₁^-2
    R_I_E_bound =
        C.J_E * (C_R_1 + C_R_2) / 4 *
        (3abs(γ)^2 * C.P^2 + 3abs(γ) * C.P * A_max + A_max^2) *
        ξ₁^-4
    R_I_E = add_error(zero(Acb), R_I_E_bound)

    return I_E_main + R_I_E
end

"""
    p_Q_0(γ, κ, ϵ, ξ₁, Λ)

Compute the coefficient ``p_{Q,0}`` from Lemma REF(lemma:p_Q_0). Up to
a factor `c^-a`, it is the coefficient for the leading order term in
the asymptotic expansion of `Q`.
"""
function p_Q_0(γ::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, Λ::CGLParams{Arb})
    return γ + I_E_infty_enclosure(γ, κ, ϵ, ξ₁, Λ)
end
