function C_I_K_1(v::Arb, C_Y::Union{FunctionBounds_Y,FunctionBounds_Y_new})
    return C_Y.K_1 * C_Y.J_N / abs(v - 2)
end

function C_I_K_2(
    lambda::Acb,
    κ::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C_Y::Union{FunctionBounds_Y,FunctionBounds_Y_new},
)
    exponent = real_s34(lambda, κ, λ) - real_s12(lambda, κ, λ) + v - 2
    return C_Y.K_2 * C_Y.J_N / abs(exponent)
end

function C_I_K_1_dλ_1_1(v::Arb, C_Y::Union{FunctionBounds_Y,FunctionBounds_Y_new})
    return C_Y.K_1_dλ * C_Y.J_N / abs(v - 2)
end

function C_I_K_1_dλ_1_2(v::Arb, C_Y::Union{FunctionBounds_Y,FunctionBounds_Y_new})
    return C_Y.K_1_dλ * C_Y.J_N / (v - 2)^2
end

function C_I_K_2_dλ_1_1(
    lambda::Acb,
    κ::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C_Y::Union{FunctionBounds_Y,FunctionBounds_Y_new},
)
    exponent = real_s34(lambda, κ, λ) - real_s12(lambda, κ, λ) + v - 2
    return C_Y.K_2_dλ * C_Y.J_N / abs(exponent)
end

function C_I_K_2_dλ_1_2(
    lambda::Acb,
    κ::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C_Y::Union{FunctionBounds_Y,FunctionBounds_Y_new},
)
    exponent = real_s34(lambda, κ, λ) - real_s12(lambda, κ, λ) + v - 2
    return C_Y.K_2_dλ * C_Y.J_N / exponent^2
end

function C_I_K_1_dλ_2(v::Arb, C_Y::Union{FunctionBounds_Y,FunctionBounds_Y_new})
    return C_I_K_1(v, C_Y)
end

function C_I_K_2_dλ_2(
    lambda::Acb,
    κ::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C_Y::Union{FunctionBounds_Y,FunctionBounds_Y_new},
)
    return C_I_K_2(lambda, κ, v, λ, C_Y)
end
