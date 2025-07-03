"""
    Y_infinity(c, lambda, γ₁, γ₂, κ, ϵ, ξ₁, λ::CGLParams)

Compute the solution to the ODE on the interval ``[ξ₁, ∞)``. Returns a
vector with two complex values, where the first is the value at `ξ₁`
and the second is the derivative.
"""
function Y_infinity(
    c::SVector{2,Acb},
    lambda::Acb,
    γ₁::Acb,
    γ₂::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
)
    v = Arb("0.1")

    # Compute enclosure of forward solution
    Q_hat, dQ_hat = Q_hat_infinity(γ₁, γ₂, κ, ϵ, ξ₁, λ)

    # Precompute functions as well as function and norm bounds
    F_Y = FunctionEnclosures_Y(lambda, γ₁, γ₂, κ, ϵ, ξ₁, λ)

    C_Y = FunctionBounds_Y(lambda, γ₁, γ₂, κ, ϵ, ξ₁, λ)

    norms_Y = NormBounds_Y(c, lambda, κ, ϵ, ξ₁, v, λ, C_Y)

    # Compute zeroth order bounds
    Y = add_error.(zero.(c), norms_Y.Y * exp(real_a12(κ, ϵ) * ξ₁^2) * ξ₁^v)
    dY = add_error.(zero.(c), norms_Y.Y_dξ * exp(real_a12(κ, ϵ) * ξ₁^2) * ξ₁^v)

    # Improve bounds
    I_K_2 = zero(Y) # FIXME: Should not be zero

    Y = F_Y.Y_12 * c + F_Y.Y_34 * I_K_2

    I_K_1_dξ = F_Y.K_1 * F_Y.J_N * Y
    I_K_2_dξ = -F_Y.K_2 * F_Y.J_N * Y

    dY =
        F_Y.Y_12_dξ * c +
        F_Y.Y_12 * I_K_1_dξ +
        F_Y.Y_34_dξ * I_K_2 +
        F_Y.Y_34 * I_K_2_dξ

    return vcat(Y, dY)
end

function Y_infinity(
    c::SVector{2,ComplexF64},
    lambda::ComplexF64,
    γ₁::ComplexF64,
    γ₂::ComplexF64,
    κ::Float64,
    ϵ::Float64,
    ξ₁::Float64,
    λ::CGLParams{Float64},
)
    Q_hat_ξ₁ = Q_hat_infinity(γ₁, γ₂, κ, ϵ, ξ₁, λ)[1]

    return Y_infinity(c, lambda, κ, ϵ, ξ₁, Q_hat_ξ₁, λ)
end

function Y_infinity(
    c::SVector{2,ComplexF64},
    lambda::ComplexF64,
    κ::Float64,
    ϵ::Float64,
    ξ₁::Float64,
    Q_hat_ξ₁::ComplexF64,
    λ::CGLParams{Float64},
)
    A = SMatrix{2,2}(ϵ, 1, -1, ϵ)

    Y12 = hcat(Y_1(ξ₁, lambda, κ, ϵ, λ), Y_2(ξ₁, lambda, κ, ϵ, λ))
    Y12_dξ = hcat(Y_1_dξ(ξ₁, lambda, κ, ϵ, λ), Y_2_dξ(ξ₁, lambda, κ, ϵ, λ))
    Y34 = hcat(Y_3(ξ₁, lambda, κ, ϵ, λ), Y_4(ξ₁, lambda, κ, ϵ, λ))
    Y34_dξ = hcat(Y_3_dξ(ξ₁, lambda, κ, ϵ, λ), Y_4_dξ(ξ₁, lambda, κ, ϵ, λ))

    JN = J_N(Q_hat_ξ₁, λ)
    K1, K2 = K_1_2(Y12, Y12_dξ, Y34, Y34_dξ, A)

    # First order approximation
    Y = Y12 * c

    I_K_1 = zero(Y) # This is actually zero
    I_K_2 = zero(Y) # TODO: Compute approximation of this

    Y = Y12 * c + Y12 * I_K_1 + Y34 * I_K_2

    I_K_1_dξ = K1 * JN * Y
    I_K_2_dξ = -K2 * JN * Y

    dY = Y12_dξ * c + Y12_dξ * I_K_1 + Y12 * I_K_1_dξ + Y34_dξ * I_K_2 + Y34 * I_K_2_dξ

    return vcat(Y, dY)
end

"""
    Y_infinity_derivative(γ, κ, ϵ, ξ₁, λ::CGLParams)

This function computes the derivative of [`Y_infinity`](@ref) w.r.t.
the parameter `lambda`.
"""
function Y_infinity_derivative(
    c::SVector{2,Acb},
    lambda::Acb,
    γ₁::Acb,
    γ₂::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
)
    v = Arb("0.1")

    C = FunctionBounds_Y(lambda, γ₁, γ₂, κ, ϵ, ξ₁, λ)

    norms = NormBounds_Y(c, lambda, κ, ϵ, ξ₁, v, λ, C)

    # Compute zeroth order bounds
    Y = add_error.(zero.(c), norms.Y * exp(real_a12(κ, ϵ) * ξ₁^2) * ξ₁^v)
    dY = add_error.(zero.(c), norms.Y_dξ * exp(real_a12(κ, ϵ) * ξ₁^2) * ξ₁^v)
    # TODO: Add proper norm bounds
    Y_dλ = add_error.(zero.(c), norms.Y * exp(real_a12(κ, ϵ) * ξ₁^2) * ξ₁^v)
    dY_dλ = add_error.(zero.(c), norms.Y_dξ * exp(real_a12(κ, ϵ) * ξ₁^2) * ξ₁^v)

    for _ = 1:3
        # TODO
    end

    return SVector(Y_dλ..., dY_dλ...)
end

function Y_infinity_derivative(
    c::SVector{2,ComplexF64},
    lambda::ComplexF64,
    γ₁::ComplexF64,
    γ₂::ComplexF64,
    κ::Float64,
    ϵ::Float64,
    ξ₁::Float64,
    λ::CGLParams{Float64},
)
    Q_hat_ξ₁ = Q_hat_infinity(γ₁, γ₂, κ, ϵ, ξ₁, λ)[1]

    return Y_infinity_derivative(c, lambda, κ, ϵ, ξ₁, Q_hat_ξ₁, λ)
end

function Y_infinity_derivative(
    c::SVector{2,ComplexF64},
    lambda::ComplexF64,
    κ::Float64,
    ϵ::Float64,
    ξ₁::Float64,
    Q_hat_ξ₁::ComplexF64,
    λ::CGLParams{Float64},
)
    A = SMatrix{2,2}(ϵ, 1, -1, ϵ)

    Y1 = Y_1(ξ₁, lambda, κ, ϵ, λ)
    Y1_dξ = Y_1_dξ(ξ₁, lambda, κ, ϵ, λ)
    Y2 = Y_2(ξ₁, lambda, κ, ϵ, λ)
    Y2_dξ = Y_2_dξ(ξ₁, lambda, κ, ϵ, λ)

    Y12 = hcat(Y1, Y2)
    Y12_dξ = hcat(Y1_dξ, Y2_dξ)
    Y34 = hcat(Y_3(ξ₁, lambda, κ, ϵ, λ), Y_4(ξ₁, lambda, κ, ϵ, λ))
    Y34_dξ = hcat(Y_3_dξ(ξ₁, lambda, κ, ϵ, λ), Y_4_dξ(ξ₁, lambda, κ, ϵ, λ))

    JN = J_N(Q_hat_ξ₁, λ)
    K1, K2 = K_1_2(Y12, Y12_dξ, Y34, Y34_dξ, A)

    # First order approximation
    Y = Y12 * c

    I_K_1 = zero(Y) # This is actually zero
    I_K_2 = zero(Y) # TODO: Compute approximation of this

    Y = Y12 * c + Y12 * I_K_1 + Y34 * I_K_2

    I_K_1_dξ = K1 * JN * Y
    I_K_2_dξ = -K2 * JN * Y

    dY = Y12_dξ * c + Y12_dξ * I_K_1 + Y12 * I_K_1_dξ + Y34_dξ * I_K_2 + Y34 * I_K_2_dξ

    # FIXME: Implement these
    Y_dλ = Y
    dY_dλ = Y

    return SVector(Y_dλ..., dY_dλ...)
end
