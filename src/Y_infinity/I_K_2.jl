function I_K_2_enclosure(
    c::SVector{2,Acb},
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C_Y::FunctionBounds_Y,
    norms_Y::NormBounds_Y,
)
    (; d, σ) = λ
    exponent = 2 / σ - d - 2real(lambda) / κ + v - 2
    real_a = real_a12(κ, ϵ)

    @assert v > 0
    @assert real_a < 0
    @assert exponent < 0

    C_I_K_2 = C_Y.K_1 * C_Y.J_N / abs(exponent)

    I_K_2_bound = C_I_K_2 * exp(real_a12(κ, ϵ) * ξ₁^2) * ξ₁^(exponent) * norms_Y.Y

    return add_error.(zero(c), I_K_2_bound)
end

function I_K_2_dλ_enclosure(
    c::SVector{2,Acb},
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C_Y::FunctionBounds_Y,
    norms_Y::NormBounds_Y,
)
    return I_K_2_dλ_1_enclosure(c, lambda, κ, ϵ, ξ₁, v, λ, C_Y, norms_Y) +
           I_K_2_dλ_2_enclosure(c, lambda, κ, ϵ, ξ₁, v, λ, C_Y, norms_Y)
end

function I_K_2_dλ_1_enclosure(
    c::SVector{2,Acb},
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C_Y::FunctionBounds_Y,
    norms_Y::NormBounds_Y,
)
    (; d, σ) = λ
    exponent = 2 / σ - d - 2real(lambda) / κ + v - 2
    real_a = real_a12(κ, ϵ)

    @assert v > 0
    @assert real_a < 0
    @assert exponent < 0

    C_I_K_2_dλ_1_1 = C_Y.K_1_dλ * C_Y.J_N / (exponent)^2
    C_I_K_2_dλ_1_2 = C_Y.K_1_dλ * C_Y.J_N / abs(exponent)

    I_K_2_dλ_bound =
        (C_I_K_2_dλ_1_1 * log(ξ₁) + C_I_K_2_dλ_1_2) *
        exp(real_a12(κ, ϵ) * ξ₁^2) *
        ξ₁^(exponent) *
        norms_Y.Y

    return add_error.(zero(c), I_K_2_dλ_bound)
end

function I_K_2_dλ_2_enclosure(
    c::SVector{2,Acb},
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C_Y::FunctionBounds_Y,
    norms_Y::NormBounds_Y,
)
    (; d, σ) = λ
    exponent = 2 / σ - d - 2real(lambda) / κ + v - 2
    real_a = real_a12(κ, ϵ)

    @assert v > 0
    @assert real_a < 0
    @assert exponent < 0

    C_I_K_dλ_2 = C_Y.K_1 * C_Y.J_N / abs(exponent)

    I_K_2_dλ_bound = C_I_K_dλ_2 * exp(real_a12(κ, ϵ) * ξ₁^2) * ξ₁^(exponent) * norms_Y.Y

    return add_error.(zero(c), I_K_2_dλ_bound)
end
