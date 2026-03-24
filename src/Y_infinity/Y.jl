"""
    Y_infinity(c, lambda, γ₁, γ₂, κ, ϵ, ξ₁, λ::CGLParams)

Compute the solution to the ODE on the interval ``[ξ₁, ∞)``. Returns a
vector with two complex values, where the first is the value at `ξ₁`
and the second is the derivative.
"""
function Y_infinity(
    c_0::SVector{2,Acb},
    lambda::Acb,
    γ₁::Acb,
    γ₂::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
)
    v = Arb("0.1")
    c = _c(κ, ϵ, λ)

    # Precompute functions as well as function and norm bounds
    F_Z = FunctionEnclosures_Y(lambda, γ₁, γ₂, κ, ϵ, ξ₁, λ)

    C_Z = FunctionBounds_Y(lambda, γ₁, γ₂, κ, ϵ, ξ₁, v, λ)

    norms_Z = NormBounds_Y(c_0, lambda, κ, ϵ, ξ₁, v, λ, C_Z)

    # Compute zeroth order bounds
    Z = add_error.(zero(c_0), norms_Z.Z * exp(-real(c) * ξ₁^2) * ξ₁^v)
    dZ = indeterminate.(c_0)

    # Improve bounds iteratively.
    for _ = 1:5
        I_K_2 = I_K_2_enclosure(c_0, lambda, κ, ϵ, ξ₁, v, λ, Z, F_Z, C_Z, norms_Z)

        Z = F_Z.E_12 * c_0 + F_Z.P_12 * I_K_2

        I_K_1_dξ = F_Z.K_1 * F_Z.I_N * Z
        I_K_2_dξ = -F_Z.K_2 * F_Z.I_N * Z

        dZ =
            F_Z.E_12_dξ * c_0 +
            F_Z.E_12 * I_K_1_dξ +
            F_Z.P_12_dξ * I_K_2 +
            F_Z.P_12 * I_K_2_dξ
    end

    M = SMatrix{2,2}(im, 1, -im, 1)
    Y = M * Z
    dY = M * dZ
    return vcat(Y, dY)
end

"""
    Y_infinity_derivative(γ, κ, ϵ, ξ₁, λ::CGLParams)

This function computes the derivative of [`Y_infinity`](@ref) w.r.t.
the parameter `lambda`.
"""
function Y_infinity_derivative(
    c_0::SVector{2,Acb},
    lambda::Acb,
    γ₁::Acb,
    γ₂::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
)
    v = Arb("0.1")
    c = _c(κ, ϵ, λ)

    # Precompute functions as well as function and norm bounds
    F_Z = FunctionEnclosures_Y(lambda, γ₁, γ₂, κ, ϵ, ξ₁, λ)

    C_Z = FunctionBounds_Y(lambda, γ₁, γ₂, κ, ϵ, ξ₁, v, λ)

    norms_Z = NormBounds_Y(c_0, lambda, κ, ϵ, ξ₁, v, λ, C_Z)

    # Compute zeroth order bounds
    Z = add_error.(zero(c_0), norms_Z.Z * exp(-real(c) * ξ₁^2) * ξ₁^v)
    dZ = indeterminate.(c_0)
    Z_dλ = add_error.(zero(c_0), norms_Z.Z_dλ * exp(-real(c) * ξ₁^2) * ξ₁^v)
    dZ_dλ = indeterminate.(c_0)

    # Improve bounds iteratively.
    for _ = 1:5
        I_K_2 = I_K_2_enclosure(c_0, lambda, κ, ϵ, ξ₁, v, λ, Z, F_Z, C_Z, norms_Z)

        I_K_2_dλ = I_K_2_dλ_enclosure(c_0, lambda, κ, ϵ, ξ₁, v, λ, C_Z, norms_Z)

        Z = F_Z.E_12 * c_0 + F_Z.P_12 * I_K_2
        Z_dλ = F_Z.E_12_dλ * c_0 + F_Z.P_12_dλ * I_K_2 + F_Z.P_12 * I_K_2_dλ

        I_K_1_dξ = F_Z.K_1 * F_Z.I_N * Z
        I_K_2_dξ = -F_Z.K_2 * F_Z.I_N * Z

        I_K_1_dλ_dξ = F_Z.K_1_dλ * F_Z.I_N * Z + F_Z.K_1 * F_Z.I_N * Z_dλ
        I_K_2_dλ_dξ = -F_Z.K_2_dλ * F_Z.I_N * Z - F_Z.K_2 * F_Z.I_N * Z_dλ

        dZ_dλ =
            F_Z.E_12_dλ_dξ * c_0 +
            F_Z.E_12_dλ * I_K_1_dξ +
            F_Z.E_12 * I_K_1_dλ_dξ +
            F_Z.P_12_dλ_dξ * I_K_2 +
            F_Z.P_12_dξ * I_K_2_dλ +
            F_Z.P_12_dλ * I_K_2_dξ +
            F_Z.P_12 * I_K_2_dλ_dξ
    end

    M = SMatrix{2,2}(im, 1, -im, 1)
    Y_dλ = M * Z_dλ
    dY_dλ = M * dZ_dλ
    return vcat(Y_dλ, dY_dλ)
end
