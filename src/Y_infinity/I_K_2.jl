"""
    I_K_2_enclosure(c_0, lambda, κ, ϵ, ξ₁, v, λ, Z, F_Z, C_Z, norms_Z)

Compute an enclosure of ``I_{K_2}(ξ)`` at the point `ξ = ξ₁`.

It uses two different approaches depending on the accuracy of the
enclosure for `Z`.

For wide enclosures of `Z`, it bounds the value using Lemma
REF(lemma:I_K_1-I_K_2-bounds). In this case we classify wide as the
enclosure overlapping zero.

For tighter enclosures of `Z` the bound is based on Lemma REF(TODO).
TODO: Add more documentation once the Lemma is finalized.
"""
function I_K_2_enclosure(
    c_0::SVector{2,Acb},
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    Z::SVector{2,Acb},
    F_Z::FunctionEnclosures_Y,
    C_Z::FunctionBounds_Y,
    norms_Z::NormBounds_Y,
)
    (; d, σ) = λ
    c = _c(κ, ϵ, λ)

    if all(Arblib.contains_zero, Z)
        # This is the wide case. The bound is computed using
        # REF(lemma:I_K_1-I_K_2-bounds)

        exponent = 2 / σ - d - 2real(lambda) / κ + v - 4

        # These are the requirements of Lemma REF(lemma:I_K_1-I_K_2-bounds)
        @assert v > 0
        @assert real(c) > 0
        @assert exponent < 0

        I_K_2_1_bound = C_Z.I_K_2_1 * exp(-real(c) * ξ₁^2) * ξ₁^exponent * norms_Z.Z
        I_K_2_2_bound = C_Z.I_K_2_2 * exp(-real(c) * ξ₁^2) * ξ₁^exponent * norms_Z.Z

        return add_error.(zero(c_0), SVector(I_K_2_1_bound, I_K_2_2_bound))
    end

    # This is the tight case. The bound is based on Lemma REF(TODO).

    # The requirements of Lemma REF(TODO) related to Lemma REF(TODO)
    # are checked in the computation of C_Z. The other requirements
    # are:
    @assert v > 0
    @assert real(c) > 0

    H = F_Z.K_2 * F_Z.I_N
    main = inv(2c) * H * Z

    # Bound remainder term
    C_Z_1 =
        C_Z.E_1 * abs(c_0[1]) +
        (C_Z.E_1 * C_Z.I_K_1_1 + C_Z.P_1 * C_Z.I_K_2_1 * ξ₁^-2) * ξ₁^(v - 2) * norms_Z.Z
    C_Z_2 =
        C_Z.E_2 * abs(c_0[2]) +
        (C_Z.E_2 * C_Z.I_K_1_2 + C_Z.P_2 * C_Z.I_K_2_2 * ξ₁^-2) * ξ₁^(v - 2) * norms_Z.Z

    C_exp_Z_1_dξ =
        C_Z.exp_E_1_dξ * abs(c_0[1]) +
        (
            C_Z.exp_E_1_dξ * C_Z.I_K_1_1 +
            C_Z.E_1 * C_Z.K_1_1 * C_Z.I_N +
            C_Z.exp_P_1_dξ * C_Z.I_K_2_1 +
            C_Z.P_1 * C_Z.K_2_1 * C_Z.I_N
        ) * ξ₁^(v - 2)
    C_exp_Z_2_dξ =
        C_Z.exp_E_2_dξ * abs(c_0[2]) +
        (
            C_Z.exp_E_2_dξ * C_Z.I_K_1_2 +
            C_Z.E_2 * C_Z.K_1_2 * C_Z.I_N +
            C_Z.exp_P_2_dξ * C_Z.I_K_2_2 +
            C_Z.P_2 * C_Z.K_2_2 * C_Z.I_N
        ) * ξ₁^(v - 2)

    C_D_11 = C_Z.H_11 * C_Z_1
    C_D_12 = C_Z.H_12 * C_Z_2
    C_D_21 = C_Z.H_21 * C_Z_1
    C_D_22 = C_Z.H_22 * C_Z_2

    C_D_11_dξ = C_Z.H_11_dξ * C_Z_1 + C_Z.H_11 * C_exp_Z_1_dξ
    C_D_12_dξ = C_Z.H_12_dξ * C_Z_2 + C_Z.H_12 * C_exp_Z_2_dξ
    C_D_21_dξ = C_Z.H_21_dξ * C_Z_1 + C_Z.H_21 * C_exp_Z_1_dξ
    C_D_22_dξ = C_Z.H_22_dξ * C_Z_2 + C_Z.H_22 * C_exp_Z_2_dξ

    β = 2 / σ - d - 2real(lambda) / κ - 5
    # Enclosure of integral of exp(-real(c)η^2) * η^β from ξ₁ to infinity
    integral_exponential = inv(2real(c)^((β + 1) / 2)) * gamma((β + 1) / 2, real(c) * ξ₁^2)
    I_K_2_11_bound = (C_D_11 + C_D_11_dξ) * integral_exponential
    I_K_2_12_bound = (C_D_12 + C_D_12_dξ) * integral_exponential
    I_K_2_21_bound = (C_D_21 + C_D_21_dξ) * integral_exponential
    I_K_2_22_bound = (C_D_22 + C_D_22_dξ) * integral_exponential

    remainder_bound =
        SVector(I_K_2_11_bound + I_K_2_12_bound, I_K_2_21_bound + I_K_2_22_bound)

    return add_error.(main, remainder_bound)
end
