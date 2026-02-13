function I_K_2_enclosure(
    c_0::SVector{2,Acb},
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    Z::SVector{2,Acb},
    dZ::SVector{2,Acb},
    F_Z::FunctionEnclosures_Y,
    C_Z::FunctionBounds_Y,
    C_I_K_j::I_K_j_Bounds,
    norms_Z::NormBounds_Y,
)
    (; d, σ) = λ
    _, _, c = _abc(κ, ϵ, λ)

    @assert v > 0
    @assert real(c) > 0

    # In the first iteration dY is non-finite. We then compute a
    # zeroth order bound.
    if !all(isfinite, dZ)
        exponent = 2 / σ - d - 2real(lambda) / κ + v - 4
        @assert exponent < 0

        I_K_2_1_bound = C_I_K_j.C_I_K_2_1 * exp(-real(c) * ξ₁^2) * ξ₁^exponent * norms_Z.Z
        I_K_2_2_bound = C_I_K_j.C_I_K_2_2 * exp(-real(c) * ξ₁^2) * ξ₁^exponent * norms_Z.Z

        return add_error.(zero(c_0), SVector(I_K_2_1_bound, I_K_2_2_bound))
    end

    # In later iterations we have finite enclosures of both Y and dY.
    # We then compute an enclosure coming from integration by parts
    # two times. See Lemma REF(XXX).

    H = F_Z.K_2 * F_Z.I_N

    D_11 = H[1, 1] * exp(c * ξ₁^2) * Z[1]
    D_21 = H[2, 1] * exp(c * ξ₁^2) * Z[2]
    D_12 = H[1, 2] * exp(conj(c) * ξ₁^2) * Z[1]
    D_22 = H[2, 2] * exp(conj(c) * ξ₁^2) * Z[2]

    # IMPROVE: We could explicitly cancel the exponentials here.
    I_K_2_11_main = inv(2c) * exp(-c * ξ₁^2) * D_11
    I_K_2_21_main = inv(2c) * exp(-c * ξ₁^2) * D_21
    I_K_2_12_main = inv(2c) * exp(-conj(c) * ξ₁^2) * D_12
    I_K_2_22_main = inv(2c) * exp(-conj(c) * ξ₁^2) * D_22

    # Bound remainder term
    C_Z_1 =
        C_Z.E_1 * abs(c_0[1]) +
        (C_Z.E_1 * C_I_K_j.C_I_K_1_1 + C_Z.P_1 * C_I_K_j.C_I_K_2_1 * ξ₁^-2) *
        ξ₁^(v - 2) *
        norms_Z.Z
    C_Z_2 =
        C_Z.E_2 * abs(c_0[2]) +
        (C_Z.E_2 * C_I_K_j.C_I_K_1_2 + C_Z.P_2 * C_I_K_j.C_I_K_2_2 * ξ₁^-2) *
        ξ₁^(v - 2) *
        norms_Z.Z

    C_exp_Z_1_dξ =
        C_Z.exp_E_1_dξ * abs(c_0[1]) +
        (
            C_Z.exp_E_1_dξ * C_I_K_j.C_I_K_1_1 +
            C_Z.E_1 * C_Z.K_1_1 * C_Z.I_N +
            C_Z.exp_P_1_dξ * C_I_K_j.C_I_K_2_1 +
            C_Z.P_1 * C_Z.K_2_1 * C_Z.I_N
        ) * ξ₁^(v - 2)
    C_exp_Z_2_dξ =
        C_Z.exp_E_2_dξ * abs(c_0[2]) +
        (
            C_Z.exp_E_2_dξ * C_I_K_j.C_I_K_1_2 +
            C_Z.E_2 * C_Z.K_1_2 * C_Z.I_N +
            C_Z.exp_P_2_dξ * C_I_K_j.C_I_K_2_2 +
            C_Z.P_2 * C_Z.K_2_2 * C_Z.I_N
        ) * ξ₁^(v - 2)

    C_D_11 = C_Z.H_1j * C_Z_1
    C_D_12 = C_Z.H_1j * C_Z_2
    C_D_21 = C_Z.H_2j * C_Z_1
    C_D_22 = C_Z.H_2j * C_Z_2

    C_D_11_dξ = C_Z.H_1j_dξ * C_Z_1 + C_Z.H_1j * C_exp_Z_1_dξ
    C_D_12_dξ = C_Z.H_1j_dξ * C_Z_2 + C_Z.H_1j * C_exp_Z_2_dξ
    C_D_21_dξ = C_Z.H_2j_dξ * C_Z_1 + C_Z.H_2j * C_exp_Z_1_dξ
    C_D_22_dξ = C_Z.H_2j_dξ * C_Z_2 + C_Z.H_2j * C_exp_Z_2_dξ

    exponent = 2 / σ - d - 2real(lambda) / κ - 6
    @assert exponent < 0
    I_K_2_11_bound = (C_D_11 + C_D_11_dξ) / 2real(c) * exp(-real(c) * ξ₁^2) * ξ₁^exponent
    I_K_2_12_bound = (C_D_12 + C_D_12_dξ) / 2real(c) * exp(-real(c) * ξ₁^2) * ξ₁^exponent
    I_K_2_21_bound = (C_D_21 + C_D_21_dξ) / 2real(c) * exp(-real(c) * ξ₁^2) * ξ₁^exponent
    I_K_2_22_bound = (C_D_22 + C_D_22_dξ) / 2real(c) * exp(-real(c) * ξ₁^2) * ξ₁^exponent

    main = SVector(I_K_2_11_main + I_K_2_12_main, I_K_2_21_main + I_K_2_22_main)
    bound = SVector(I_K_2_11_bound + I_K_2_12_bound, I_K_2_21_bound + I_K_2_22_bound)

    return add_error.(main, bound)
end

function I_K_2_dλ_enclosure(
    c_0::SVector{2,Acb},
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C_Z::FunctionBounds_Y,
    C_I_K_j::I_K_j_Bounds,
    norms_Z::NormBounds_Y,
)
    return I_K_2_dλ_1_enclosure(c_0, lambda, κ, ϵ, ξ₁, v, λ, C_Z, C_I_K_j, norms_Z) +
           I_K_2_dλ_2_enclosure(c_0, lambda, κ, ϵ, ξ₁, v, λ, C_Z, C_I_K_j, norms_Z)
end

function I_K_2_dλ_1_enclosure(
    c_0::SVector{2,Acb},
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C_Z::FunctionBounds_Y,
    C_I_K_j::I_K_j_Bounds,
    norms_Z::NormBounds_Y,
)
    (; d, σ) = λ
    _, _, c = _abc(κ, ϵ, λ)
    exponent = 2 / σ - d - 2real(lambda) / κ + v - 4

    @assert v > 0
    @assert -real(c) < 0
    @assert exponent < 0

    I_K_2_dλ_1_1_bound =
        C_I_K_j.C_I_K_2_dλ_1_1 * exp(-real(c) * ξ₁^2) * log(ξ₁) * ξ₁^exponent * norms_Z.Z
    I_K_2_dλ_1_2_bound =
        C_I_K_j.C_I_K_2_dλ_1_2 * exp(-real(c) * ξ₁^2) * log(ξ₁) * ξ₁^exponent * norms_Z.Z

    return add_error.(zero(c_0), SVector(I_K_2_dλ_1_1_bound, I_K_2_dλ_1_2_bound))
end

function I_K_2_dλ_2_enclosure(
    c_0::SVector{2,Acb},
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C_Z::FunctionBounds_Y,
    C_I_K_j::I_K_j_Bounds,
    norms_Z::NormBounds_Y,
)
    (; d, σ) = λ
    _, _, c = _abc(κ, ϵ, λ)
    exponent = 2 / σ - d - 2real(lambda) / κ + v - 4

    @assert v > 0
    @assert -real(c) < 0
    @assert exponent < 0

    I_K_2_dλ_2_1_bound =
        C_I_K_j.C_I_K_2_dλ_2_1 * exp(-real(c) * ξ₁^2) * ξ₁^exponent * norms_Z.Z_dλ
    I_K_2_dλ_2_2_bound =
        C_I_K_j.C_I_K_2_dλ_2_2 * exp(-real(c) * ξ₁^2) * ξ₁^exponent * norms_Z.Z_dλ

    return add_error.(zero(c_0), SVector(I_K_2_dλ_2_1_bound, I_K_2_dλ_2_2_bound))
end
