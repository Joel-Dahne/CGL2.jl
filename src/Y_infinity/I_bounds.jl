# Lemma I_K_1-I_K_2-bounds

function C_I_K_1_1(v::Arb, C_Y::FunctionBounds_Y)
    @assert v - 2 < 0
    return C_Y.K_1_1 * C_Y.J_N / abs(v - 2)
end

function C_I_K_1_2(v::Arb, C_Y::FunctionBounds_Y)
    @assert v - 2 < 0
    return C_Y.K_1_2 * C_Y.J_N / abs(v - 2)
end

function C_I_K_2_1(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C_Y::FunctionBounds_Y,
)
    _, _, c = _abc(κ, ϵ, λ)
    return C_Y.K_2_1 * C_Y.J_N / (2real(c))
end

function C_I_K_2_2(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C_Y::FunctionBounds_Y,
)
    _, _, c = _abc(κ, ϵ, λ)
    return C_Y.K_2_2 * C_Y.J_N / (2real(c))
end

# Lemma Y-I_K_1_I_K_2-lambda-1-bounds

function C_I_K_1_dλ_1_1(ξ₁::Arb, v::Arb, C_Y::FunctionBounds_Y)
    @assert v - 2 < 0
    return C_Y.K_1_dλ_1 * C_Y.J_N * (abs(v - 2) + inv(log(ξ₁))) / (v - 2)^2
end

function C_I_K_1_dλ_1_2(ξ₁::Arb, v::Arb, C_Y::FunctionBounds_Y)
    @assert v - 2 < 0
    return C_Y.K_1_dλ_2 * C_Y.J_N * (abs(v - 2) + inv(log(ξ₁))) / (v - 2)^2
end


function C_I_K_2_dλ_1_1(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C_Y::FunctionBounds_Y,
)
    _, _, c = _abc(κ, ϵ, λ)
    exponent = 2 / λ.σ - λ.d - 2real(lambda) / κ + v - 4
    @assert ξ₁ > exp(-inv(exponent))
    return C_Y.K_2_dλ_1 * C_Y.J_N / (2real(c))
end

function C_I_K_2_dλ_1_2(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C_Y::FunctionBounds_Y,
)
    _, _, c = _abc(κ, ϵ, λ)
    exponent = 2 / λ.σ - λ.d - 2real(lambda) / κ + v - 4
    @assert ξ₁ > exp(-inv(exponent))
    return C_Y.K_2_dλ_2 * C_Y.J_N / (2real(c))
end

# Lemma Y-I_K_1_I_K_2-lambda-2-bounds

function C_I_K_1_dλ_2_1(v::Arb, C_Y::FunctionBounds_Y)
    return C_I_K_1_2(v, C_Y)
end

function C_I_K_1_dλ_2_2(v::Arb, C_Y::FunctionBounds_Y)
    return C_I_K_1_2(v, C_Y)
end

function C_I_K_2_dλ_2_1(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C_Y::FunctionBounds_Y,
)
    return C_I_K_2_1(lambda, κ, ϵ, v, λ, C_Y)
end

function C_I_K_2_dλ_2_2(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C_Y::FunctionBounds_Y,
)
    return C_I_K_2_2(lambda, κ, ϵ, v, λ, C_Y)
end
