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
    _, _, c = _abc(κ, ϵ, λ) # TODO: We use the name c for two things...
    M = SMatrix{2,2}(im, 1, -im, 1)
    exponent = 2 / σ - d - 2real(lambda) / κ + v - 4

    @assert v > 0
    @assert -real(c) < 0
    @assert exponent < 0

    # In the first iteration dY is non-finite. We then compute a
    # zeroth order bound.

    I_K_2_1_bound = C_I_K_j.C_I_K_2_1 * exp(-real(c) * ξ₁^2) * ξ₁^exponent * norms_Z.Z
    I_K_2_2_bound = C_I_K_j.C_I_K_2_2 * exp(-real(c) * ξ₁^2) * ξ₁^exponent * norms_Z.Z

    return add_error.(zero(c_0), SVector(I_K_2_1_bound, I_K_2_2_bound))

    # TODO: Implement the below

    if !all(isfinite, dY)
        return add_error.(zero(c_0), SVector(I_K_2_1_bound, I_K_2_2_bound))
    end

    # In later iterations we have finite enclosures of both Y and dY.
    # We then compute an enclosure coming from integration by parts
    # two times. See Lemma REF(XXX).

    H = F_Y.K_2 * F_Y.J_N * M

    D_11 = H[1, 1] * (exp(c*ξ₁^2)*inv(M)*Y)[1]
    D_21 = H[2, 1] * (exp(c*ξ₁^2)*inv(M)*Y)[1]
    D_12 = H[1, 2] * (exp(conj(c)*ξ₁^2)*inv(M)*Y)[2]
    D_22 = H[2, 2] * (exp(conj(c)*ξ₁^2)*inv(M)*Y)[2]

    # IMPROVE: We could explicitly cancel the exponentials here.
    I_K_2_11_main = inv(2c) * exp(-c * ξ₁^2) * D_11
    I_K_2_21_main = inv(2c) * exp(-c * ξ₁^2) * D_21
    I_K_2_12_main = inv(2c) * exp(-conj(c) * ξ₁^2) * D_12
    I_K_2_22_main = inv(2c) * exp(-conj(c) * ξ₁^2) * D_22

    # Bound remainder term

    # TODO: Check that these are correct
    C_D_11 = C_Y.J_E_1 * C_Y.J_N / sqrt(1 + ϵ^2)
    C_D_12 = C_Y.J_E_1 * C_Y.J_N / sqrt(1 + ϵ^2)
    C_D_21 = C_Y.J_E_2 * C_Y.J_N / sqrt(1 + ϵ^2)
    C_D_22 = C_Y.J_E_2 * C_Y.J_N / sqrt(1 + ϵ^2)

    # FIXME
    C_D_11_dξ = Arb(1)
    C_D_12_dξ = Arb(1)
    C_D_21_dξ = Arb(1)
    C_D_22_dξ = Arb(1)

    # FIXME: Implement this and probably move it elsewhere?
    exp_pow_integral_bound(c, b) = inv(2c) * exp(-c * ξ₁^2) * ξ₁^(b - 1)

    # FIXME
    I_K_2_11_bound =
        inv(2abs(c)) * (
            C_D_11_dξ * exp_pow_integral_bound(real(c), exponent) +
            C_D_11 * exp_pow_integral_bound(real(c), exponent)
        )
    I_K_2_21_bound =
        inv(2abs(c)) * (
            C_D_21_dξ * exp_pow_integral_bound(real(c), exponent) +
            C_D_21 * exp_pow_integral_bound(real(c), exponent)
        )
    I_K_2_12_bound =
        inv(2abs(c)) * (
            C_D_12_dξ * exp_pow_integral_bound(real(c), exponent) +
            C_D_12 * exp_pow_integral_bound(real(c), exponent)
        )
    I_K_2_22_bound =
        inv(2abs(c)) * (
            C_D_22_dξ * exp_pow_integral_bound(real(c), exponent) +
            C_D_22 * exp_pow_integral_bound(real(c), exponent)
        )

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
