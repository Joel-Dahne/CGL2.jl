# TODO: Add tests and documentation for this function
function integral_J_E_hat_P_hat(κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    a, b, c = _abc(κ, ϵ, λ)
    (; d, σ) = λ

    I_U = zero(Acb)

    n = 4

    p_U = map(0:(n-1)) do k
        rising(a, k) * rising(a - b + 1, k) / (factorial(k) * c^k)
    end
    p_U_bma = map(0:(n-1)) do k
        rising(b - a, k) * rising(-a + 1, k) / (factorial(k) * (-c)^k)
    end

    C_R_U = CGL2.C_R_U(n, a, b, -c * ξ₁^2)
    C_R_U_bma = CGL2.C_R_U(n, b - a, b, c * ξ₁^2)

    for i_1 = 0:n
        for i_2 = 0:n
            for i_3 = 0:n
                for i_4 = 0:n
                    # Exponent after integration
                    exponent = -2 / σ - 2(i_1 + i_2 + i_3 + i_4)

                    if i_1 == n || i_2 == n || i_3 == n || i_4 == n
                        # One of the terms is a remainder term. We
                        # only bound the absolute value of the
                        # integral.
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

    #I_W = Arblib.integrate(ξ₁, 20ξ₁) do η
    #    U(b - a, b, c * η^2) * U(a, b, -c * η^2)^2 * conj(U(a, b, -c * η^2)) * η^(d - 1)
    #end / (c^(a - b) * ((-c)^(-a))^2 * conj((-c)^(-a)))

    return B_W_hat(κ, ϵ, λ) * c^(a - b) * abs((-c)^-a)^2 * (-c)^-a * I_U
end

function I_E_hat_enclosure(
    γ₁::Acb,
    γ₂::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    F::FunctionEnclosures_hat,
    C::FunctionBounds_hat,
    norms::NormBounds_hat,
)
    a, b, c = _abc(κ, ϵ, λ)
    (; d, σ) = λ

    # Requirements of Lemma REF(lemma:I_E_hat-enclosure)
    # The requirements from Lemma REF(lemma:Q-hat-leading-term) are
    # checked internally by the function C_R_Q_hat.
    @assert isone(σ)

    I_E_hat_main = abs(γ₁)^2 * γ₁ * integral_J_E_hat_P_hat(κ, ϵ, ξ₁, λ)

    C_R_Q_hat_1 = CGL2.C_R_Q_hat_1(γ₁, γ₂, κ, ϵ, ξ₁, v, λ, F, C, norms)
    R_I_E_hat_bound =
        C.J_E_hat *
        (
            3abs(γ₁)^2 * C.P_hat^2 * C_R_Q_hat_1 / abs(-2 / σ + (2σ + 1) * v - 2) +
            3abs(γ₁) * C.P_hat * C_R_Q_hat_1^2 / abs(-2 / σ + 2(2σ + 1) * v - 4) *
            ξ₁^((2σ + 1) * v - 2) +
            C_R_Q_hat_1^3 / abs(-2 / σ + 3(2σ + 1) * v - 6) * ξ₁^(2(2σ + 1) * v - 4)
        ) *
        ξ₁^(-2 / σ + (2σ + 1) * v - 2)
    R_I_E_hat = add_error(zero(Acb), R_I_E_hat_bound)

    return I_E_hat_main + R_I_E_hat
end

function I_E_hat_dγ₂_enclosure(
    γ₁::Acb,
    γ₂::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    F::FunctionEnclosures_hat,
    C::FunctionBounds_hat,
    norms::NormBounds_hat,
)
    (; σ) = λ

    # Requirements of Lemma REF(lemma:I_E_hat-I_P_hat-dgamma-bounds)
    # are checked in the computation of `C`, where the associated
    # constants are computed.

    I_E_hat_dγ₂_bound =
        (2σ + 1) * C.I_E_hat * ξ₁^((2σ + 1) * v - 2) * norms.Q_hat^2σ * norms.Q_hat_dγ₂

    return add_error(zero(Acb), I_E_hat_dγ₂_bound)
end
