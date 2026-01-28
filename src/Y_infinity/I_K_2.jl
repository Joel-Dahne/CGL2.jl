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
    _, _, c_ = _abc(κ, ϵ, λ) # TODO: We use the name c for two things...
    exponent = 2 / σ - d - 2real(lambda) / κ + v - 4

    @assert v > 0
    @assert -real(c_) < 0
    @assert exponent < 0

    I_K_2_2_bound =
        C_I_K_2_2(lambda, κ, ϵ, v, λ, C_Y) * exp(-real(c_) * ξ₁^2) * ξ₁^exponent * norms_Y.Y
    I_K_2_1_bound =
        C_I_K_2_1(lambda, κ, ϵ, v, λ, C_Y) * exp(-real(c_) * ξ₁^2) * ξ₁^exponent * norms_Y.Y

    return add_error.(zero(c), SVector(I_K_2_1_bound, I_K_2_2_bound))
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
    _, _, c_ = _abc(κ, ϵ, λ) # TODO: We use the name c for two things...
    exponent = 2 / σ - d - 2real(lambda) / κ + v - 4

    @assert v > 0
    @assert -real(c_) < 0
    @assert exponent < 0

    I_K_2_dλ_1_1_bound =
        C_I_K_2_dλ_1_1(lambda, κ, ϵ, ξ₁, v, λ, C_Y) *
        exp(-real(c_) * ξ₁^2) *
        log(ξ₁) *
        ξ₁^exponent *
        norms_Y.Y
    I_K_2_dλ_1_2_bound =
        C_I_K_2_dλ_1_2(lambda, κ, ϵ, ξ₁, v, λ, C_Y) *
        exp(-real(c_) * ξ₁^2) *
        log(ξ₁) *
        ξ₁^exponent *
        norms_Y.Y

    return add_error.(zero(c), SVector(I_K_2_dλ_1_1_bound, I_K_2_dλ_1_2_bound))
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
    _, _, c_ = _abc(κ, ϵ, λ) # TODO: We use the name c for two things...
    exponent = 2 / σ - d - 2real(lambda) / κ + v - 4

    @assert v > 0
    @assert -real(c_) < 0
    @assert exponent < 0

    I_K_2_dλ_2_1_bound =
        C_I_K_2_dλ_2_1(lambda, κ, ϵ, v, λ, C_Y) *
        exp(-real(c_) * ξ₁^2) *
        ξ₁^exponent *
        norms_Y.Y_dλ
    I_K_2_dλ_2_2_bound =
        C_I_K_2_dλ_2_2(lambda, κ, ϵ, v, λ, C_Y) *
        exp(-real(c_) * ξ₁^2) *
        ξ₁^exponent *
        norms_Y.Y_dλ

    return add_error.(zero(c), SVector(I_K_2_dλ_2_1_bound, I_K_2_dλ_2_2_bound))
end
