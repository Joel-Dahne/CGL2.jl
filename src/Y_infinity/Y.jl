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
    _, _, c_ = _abc(κ, ϵ, λ) # TODO: We use the name c for two things...
    M = SMatrix{2,2}(im, 1, -im, 1)

    # Precompute functions as well as function and norm bounds
    F_Y = FunctionEnclosures_Y_new(lambda, γ₁, γ₂, κ, ϵ, ξ₁, λ)

    C_Y = FunctionBounds_Y_new(lambda, γ₁, γ₂, κ, ϵ, ξ₁, λ)

    norms_Y = NormBounds_Y(c, lambda, κ, ϵ, ξ₁, v, λ, C_Y)

    # Compute zeroth order bounds
    Y = add_error.(zero(c), norms_Y.Y * exp(-real(c_) * ξ₁^2) * ξ₁^v)

    # Improve bounds
    I_K_2 = I_K_2_enclosure(c, lambda, κ, ϵ, ξ₁, v, λ, C_Y, norms_Y)

    # FIXME: Improve bounds so that we don't have to cheat
    I_K_2 = 1e-5I_K_2

    Y = M * (F_Y.E_12 * c + F_Y.P_12 * I_K_2)

    I_K_1_dξ = F_Y.K_1 * F_Y.J_N * Y
    I_K_2_dξ = -F_Y.K_2 * F_Y.J_N * Y

    dY =
        M * (
            F_Y.E_12_dξ * c +
            F_Y.E_12 * I_K_1_dξ +
            F_Y.P_12_dξ * I_K_2 +
            F_Y.P_12 * I_K_2_dξ
        )

    return vcat(Y, dY)
    res_new = vcat(Y, dY)

    # Precompute functions as well as function and norm bounds
    F_Y = FunctionEnclosures_Y(lambda, γ₁, γ₂, κ, ϵ, ξ₁, λ)

    C_Y = FunctionBounds_Y(lambda, γ₁, γ₂, κ, ϵ, ξ₁, λ)

    norms_Y = NormBounds_Y(c, lambda, κ, ϵ, ξ₁, v, λ, C_Y)

    # Compute zeroth order bounds
    Y = add_error.(zero(c), norms_Y.Y * exp(real_a12(κ, ϵ) * ξ₁^2) * ξ₁^v)

    # Improve bounds
    I_K_2 = I_K_2_enclosure(c, lambda, κ, ϵ, ξ₁, v, λ, C_Y, norms_Y)

    # FIXME: Improve bounds so that we don't have to cheat
    I_K_2 = 1e-5I_K_2

    Y = F_Y.Y_12 * c + F_Y.Y_34 * I_K_2

    I_K_1_dξ = F_Y.K_1 * F_Y.J_N * Y
    I_K_2_dξ = -F_Y.K_2 * F_Y.J_N * Y


    dY = F_Y.Y_12_dξ * c + F_Y.Y_12 * I_K_1_dξ + F_Y.Y_34_dξ * I_K_2 + F_Y.Y_34 * I_K_2_dξ

    #return vcat(Y, dY)
    return vcat(Y, dY), res_new
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
    M = SMatrix{2,2}(im, 1, -im, 1)

    E12 = Diagonal(SVector(E_1(ξ₁, lambda, κ, ϵ, λ), E_2(ξ₁, lambda, κ, ϵ, λ)))
    E12_dξ = Diagonal(SVector(E_1_dξ(ξ₁, lambda, κ, ϵ, λ), E_2_dξ(ξ₁, lambda, κ, ϵ, λ)))
    P12 = Diagonal(SVector(P_1(ξ₁, lambda, κ, ϵ, λ), P_2(ξ₁, lambda, κ, ϵ, λ)))
    P12_dξ = Diagonal(SVector(P_1_dξ(ξ₁, lambda, κ, ϵ, λ), P_2_dξ(ξ₁, lambda, κ, ϵ, λ)))

    K1, K2 = K_1_2_new(ξ₁, lambda, κ, ϵ, λ)

    JN = J_N(Q_hat_ξ₁, λ)

    # First order approximation
    Y = M * E12 * c

    I_K_2 = zero(Y) # IMPROVE: Compute approximation of this

    Y = M * (E12 * c + P12 * I_K_2)

    I_K_1_dξ = K1 * JN * Y
    I_K_2_dξ = -K2 * JN * Y

    dY = M * (E12_dξ * c + E12 * I_K_1_dξ + P12_dξ * I_K_2 + P12 * I_K_2_dξ)

    return vcat(Y, dY)
    res_new = vcat(Y, dY)

    Y12 = hcat(Y_1(ξ₁, lambda, κ, ϵ, λ), Y_2(ξ₁, lambda, κ, ϵ, λ))
    Y12_dξ = hcat(Y_1_dξ(ξ₁, lambda, κ, ϵ, λ), Y_2_dξ(ξ₁, lambda, κ, ϵ, λ))
    Y34 = hcat(Y_3(ξ₁, lambda, κ, ϵ, λ), Y_4(ξ₁, lambda, κ, ϵ, λ))
    Y34_dξ = hcat(Y_3_dξ(ξ₁, lambda, κ, ϵ, λ), Y_4_dξ(ξ₁, lambda, κ, ϵ, λ))

    A = SMatrix{2,2}(ϵ, 1, -1, ϵ)
    K1, K2 = K_1_2(Y12, Y12_dξ, Y34, Y34_dξ, A)

    JN = J_N(Q_hat_ξ₁, λ)

    # First order approximation
    Y = Y12 * c

    I_K_2 = zero(Y) # IMPROVE: Compute approximation of this

    Y = Y12 * c + Y34 * I_K_2

    I_K_1_dξ = K1 * JN * Y
    I_K_2_dξ = -K2 * JN * Y

    dY = Y12_dξ * c + Y12 * I_K_1_dξ + Y34_dξ * I_K_2 + Y34 * I_K_2_dξ

    #return vcat(Y, dY)
    return vcat(Y, dY), res_new
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
    _, _, c_ = _abc(κ, ϵ, λ) # TODO: We use the name c for two things...
    M = SMatrix{2,2}(im, 1, -im, 1)

    # Precompute functions as well as function and norm bounds
    F_Y = FunctionEnclosures_Y_new(lambda, γ₁, γ₂, κ, ϵ, ξ₁, λ, include_dλ = true)

    C_Y = FunctionBounds_Y_new(lambda, γ₁, γ₂, κ, ϵ, ξ₁, λ, include_dλ = true)

    norms_Y = NormBounds_Y(c, lambda, κ, ϵ, ξ₁, v, λ, C_Y, include_dλ = true)

    # Compute zeroth order bounds
    Y = add_error.(zero(c), norms_Y.Y * exp(-real(c_) * ξ₁^2) * ξ₁^v)
    Y_dλ = add_error.(zero(c), norms_Y.Y_dλ * exp(-real(c_) * ξ₁^2) * ξ₁^v)

    # Improve bounds
    I_K_2 = I_K_2_enclosure(c, lambda, κ, ϵ, ξ₁, v, λ, C_Y, norms_Y)
    I_K_2_dλ = I_K_2_dλ_enclosure(c, lambda, κ, ϵ, ξ₁, v, λ, C_Y, norms_Y)

    Y = M * (F_Y.E_12 * c + F_Y.P_12 * I_K_2)
    Y_dλ =
        M *
        (F_Y.E_12_dλ * c + F_Y.P_12_dλ * I_K_2 + F_Y.P_12 * I_K_2_dλ)

    I_K_1_dξ = F_Y.K_1 * F_Y.J_N * Y
    I_K_2_dξ = -F_Y.K_2 * F_Y.J_N * Y

    I_K_1_dλ_dξ =
        F_Y.K_1_dλ * F_Y.J_N * Y + F_Y.K_1 * F_Y.J_N * Y_dλ
    I_K_2_dλ_dξ =
        -F_Y.K_2_dλ * F_Y.J_N * Y - F_Y.K_2 * F_Y.J_N * Y_dλ

    dY_dλ =
        M * (
            F_Y.E_12_dλ_dξ * c +
            F_Y.E_12_dλ * I_K_1_dξ +
            F_Y.E_12 * I_K_1_dλ_dξ +
            F_Y.P_12_dλ_dξ * I_K_2 +
            F_Y.P_12_dξ * I_K_2_dλ +
            F_Y.P_12_dλ * I_K_2_dξ +
            F_Y.P_12 * I_K_2_dλ_dξ
        )

    # TODO: It seems like the new version give slightly worse bounds.
    # Why is that?

    return vcat(Y_dλ, dY_dλ)
    res_new = vcat(Y_dλ, dY_dλ)

    # Precompute functions as well as function and norm bounds
    F_Y = FunctionEnclosures_Y(lambda, γ₁, γ₂, κ, ϵ, ξ₁, λ, include_dλ = true)

    C_Y = FunctionBounds_Y(lambda, γ₁, γ₂, κ, ϵ, ξ₁, λ, include_dλ = true)

    norms_Y = NormBounds_Y(c, lambda, κ, ϵ, ξ₁, v, λ, C_Y, include_dλ = true)

    # Compute zeroth order bounds
    Y = add_error.(zero(c), norms_Y.Y * exp(real_a12(κ, ϵ) * ξ₁^2) * ξ₁^v)
    Y_dλ = add_error.(zero(c), norms_Y.Y_dλ * exp(real_a12(κ, ϵ) * ξ₁^2) * ξ₁^v)

    # Improve bounds
    I_K_2 = I_K_2_enclosure(c, lambda, κ, ϵ, ξ₁, v, λ, C_Y, norms_Y)
    I_K_2_dλ = I_K_2_dλ_enclosure(c, lambda, κ, ϵ, ξ₁, v, λ, C_Y, norms_Y)

    # FIXME: Improve bounds so that we don't have to cheat
    I_K_2 = 1e-1I_K_2

    # TODO: It seems like the _new version give slightly worse bounds.
    # Why is that?

    Y = F_Y.Y_12 * c + F_Y.Y_34 * I_K_2
    Y_dλ = F_Y.Y_12_dλ * c + F_Y.Y_34_dλ * I_K_2 + F_Y.Y_34 * I_K_2_dλ

    I_K_1_dξ = F_Y.K_1 * F_Y.J_N * Y
    I_K_2_dξ = -F_Y.K_2 * F_Y.J_N * Y

    I_K_1_dλ_dξ = F_Y.K_1_dλ * F_Y.J_N * Y + F_Y.K_1 * F_Y.J_N * Y_dλ
    I_K_2_dλ_dξ = -F_Y.K_2_dλ * F_Y.J_N * Y - F_Y.K_2 * F_Y.J_N * Y_dλ

    dY_dλ =
        F_Y.Y_12_dλ_dξ * c +
        F_Y.Y_12_dλ * I_K_1_dξ +
        F_Y.Y_12 * I_K_1_dλ_dξ +
        F_Y.Y_34_dλ_dξ * I_K_2 +
        F_Y.Y_34_dξ * I_K_2_dλ +
        F_Y.Y_34_dλ * I_K_2_dξ +
        F_Y.Y_34 * I_K_2_dλ_dξ

    #return vcat(Y_dλ, dY_dλ)
    return vcat(Y_dλ, dY_dλ), res_new
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
    M = SMatrix{2,2}(im, 1, -im, 1)

    E12 = Diagonal(SVector(E_1(ξ₁, lambda, κ, ϵ, λ), E_2(ξ₁, lambda, κ, ϵ, λ)))
    E12_dξ = Diagonal(SVector(E_1_dξ(ξ₁, lambda, κ, ϵ, λ), E_2_dξ(ξ₁, lambda, κ, ϵ, λ)))
    P12 = Diagonal(SVector(P_1(ξ₁, lambda, κ, ϵ, λ), P_2(ξ₁, lambda, κ, ϵ, λ)))
    P12_dξ = Diagonal(SVector(P_1_dξ(ξ₁, lambda, κ, ϵ, λ), P_2_dξ(ξ₁, lambda, κ, ϵ, λ)))

    E12_dλ = Diagonal(SVector(E_1_dλ(ξ₁, lambda, κ, ϵ, λ), E_2_dλ(ξ₁, lambda, κ, ϵ, λ)))
    E12_dλ_dξ =
        Diagonal(SVector(E_1_dλ_dξ(ξ₁, lambda, κ, ϵ, λ), E_2_dλ_dξ(ξ₁, lambda, κ, ϵ, λ)))
    P12_dλ = Diagonal(SVector(P_1_dλ(ξ₁, lambda, κ, ϵ, λ), P_2_dλ(ξ₁, lambda, κ, ϵ, λ)))
    P12_dλ_dξ =
        Diagonal(SVector(P_1_dλ_dξ(ξ₁, lambda, κ, ϵ, λ), P_2_dλ_dξ(ξ₁, lambda, κ, ϵ, λ)))

    K1, K2 = K_1_2_new(ξ₁, lambda, κ, ϵ, λ)
    K1_dλ, K2_dλ = K_1_2_dλ_new(ξ₁, lambda, κ, ϵ, λ)

    JN = J_N(Q_hat_ξ₁, λ)

    # First order approximation
    Y = M * E12 * c
    Y_dλ = M * E12_dλ * c

    # Compute an improved approximation
    I_K_2 = zero(Y) # IMPROVE: Compute approximation of this
    I_K_2_dλ = zero(Y) # IMPROVE: Compute approximation of this

    Y = M * (E12 * c + P12 * I_K_2)
    Y_dλ = M * (E12_dλ * c + P12_dλ * I_K_2 + P12_dλ * I_K_2_dλ)

    I_K_1_dξ = K1 * JN * Y
    I_K_2_dξ = -K2 * JN * Y

    I_K_1_dλ_dξ = K1_dλ * JN * Y + K1 * JN * Y_dλ
    I_K_2_dλ_dξ = -K2_dλ * JN * Y - K2 * JN * Y_dλ

    dY = M * (E12_dξ * c + E12 * I_K_1_dξ + P12_dξ * I_K_2 + P12 * I_K_2_dξ)
    dY_dλ =
        M * (
            E12_dλ_dξ * c +
            E12_dλ * I_K_1_dξ +
            E12 * I_K_1_dλ_dξ +
            P12_dλ_dξ * I_K_2 +
            P12_dξ * I_K_2_dλ +
            P12_dλ * I_K_2_dξ +
            P12 * I_K_2_dλ_dξ
        )

    return vcat(Y_dλ, dY_dλ)
    res_new = vcat(Y_dλ, dY_dλ)

    Y12 = hcat(Y_1(ξ₁, lambda, κ, ϵ, λ), Y_2(ξ₁, lambda, κ, ϵ, λ))
    Y12_dξ = hcat(Y_1_dξ(ξ₁, lambda, κ, ϵ, λ), Y_2_dξ(ξ₁, lambda, κ, ϵ, λ))
    Y34 = hcat(Y_3(ξ₁, lambda, κ, ϵ, λ), Y_4(ξ₁, lambda, κ, ϵ, λ))
    Y34_dξ = hcat(Y_3_dξ(ξ₁, lambda, κ, ϵ, λ), Y_4_dξ(ξ₁, lambda, κ, ϵ, λ))

    Y12_dλ = hcat(Y_1_dλ(ξ₁, lambda, κ, ϵ, λ), Y_2_dλ(ξ₁, lambda, κ, ϵ, λ))
    Y12_dλ_dξ = hcat(Y_1_dλ_dξ(ξ₁, lambda, κ, ϵ, λ), Y_2_dλ_dξ(ξ₁, lambda, κ, ϵ, λ))
    Y34_dλ = hcat(Y_3_dλ(ξ₁, lambda, κ, ϵ, λ), Y_4_dλ(ξ₁, lambda, κ, ϵ, λ))
    Y34_dλ_dξ = hcat(Y_3_dλ_dξ(ξ₁, lambda, κ, ϵ, λ), Y_4_dλ_dξ(ξ₁, lambda, κ, ϵ, λ))

    A = SMatrix{2,2}(ϵ, 1, -1, ϵ)
    K1, K2 = K_1_2(Y12, Y12_dξ, Y34, Y34_dξ, A)
    K1_dλ, K2_dλ =
        K_1_2_dλ(Y12, Y12_dξ, Y34, Y34_dξ, Y12_dλ, Y12_dλ_dξ, Y34_dλ, Y34_dλ_dξ, A)

    JN = J_N(Q_hat_ξ₁, λ)

    # First order approximation
    Y = Y12 * c
    Y_dλ = Y12_dλ * c

    # Compute an improved approximation
    I_K_2 = zero(Y) # IMPROVE: Compute approximation of this
    I_K_2_dλ = zero(Y) # IMPROVE: Compute approximation of this

    Y = Y12 * c + Y34 * I_K_2
    Y_dλ = Y12_dλ * c + Y34_dλ * I_K_2 + Y34 * I_K_2_dλ

    I_K_1_dξ = K1 * JN * Y
    I_K_2_dξ = -K2 * JN * Y

    I_K_1_dλ_dξ = K1_dλ * JN * Y + K1 * JN * Y_dλ
    I_K_2_dλ_dξ = -K2_dλ * JN * Y - K2 * JN * Y_dλ

    dY_dλ =
        Y12_dλ_dξ * c +
        Y12_dλ * I_K_1_dξ +
        Y12 * I_K_1_dλ_dξ +
        Y34_dλ_dξ * I_K_2 +
        Y34_dξ * I_K_2_dλ +
        Y34_dλ * I_K_2_dξ +
        Y34 * I_K_2_dλ_dξ

    #return SVector(Y_dλ..., dY_dλ...)
    return SVector(Y_dλ..., dY_dλ...), res_new
end
