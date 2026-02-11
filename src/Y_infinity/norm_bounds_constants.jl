function C_T_12(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C_Z::FunctionBounds_Y,
    C_I_K_j::I_K_j_Bounds,
)
    (; d, σ) = λ
    _, _, c = _abc(κ, ϵ, λ)

    @assert v > 0
    @assert real(c) > 0
    return max(C_Z.E_1 * C_I_K_j.C_I_K_1_1, C_Z.E_2 * C_I_K_j.C_I_K_1_2) +
           max(C_Z.P_1 * C_I_K_j.C_I_K_2_1, C_Z.P_2 * C_I_K_j.C_I_K_2_2) * ξ₁^-2
end

function C_Z_dλ_1(c::SVector{2,Acb}, v::Arb, C_Z::FunctionBounds_Y)
    return max(C_Z.E_1_dλ * abs(c[1]), C_Z.E_2_dλ * abs(c[2])) * exp(Arb(-1)) / v
end

function C_Z_dλ_2(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C_Z::FunctionBounds_Y,
    C_I_K_j::I_K_j_Bounds,
)
    return max(C_Z.E_1_dλ * C_I_K_j.C_I_K_1_1, C_Z.E_2_dλ * C_I_K_j.C_I_K_1_2) *
           exp(Arb(-1)) / v * ξ₁^(v - 2) +
           max(C_Z.P_1_dλ * C_I_K_j.C_I_K_2_1, C_Z.P_2_dλ * C_I_K_j.C_I_K_2_2) *
           log(ξ₁) *
           ξ₁^-4 +
           max(C_Z.E_1 * C_I_K_j.C_I_K_1_dλ_1_1, C_Z.E_2 * C_I_K_j.C_I_K_1_dλ_1_2) *
           log(ξ₁) *
           ξ₁^-2 +
           max(C_Z.P_1 * C_I_K_j.C_I_K_2_dλ_1_1, C_Z.P_2 * C_I_K_j.C_I_K_2_dλ_1_2) *
           log(ξ₁) *
           ξ₁^-4
end

function C_Z_dλ_3(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C_Z::FunctionBounds_Y,
    C_I_K_j::I_K_j_Bounds,
)
    return max(C_Z.E_1 * C_I_K_j.C_I_K_1_dλ_2_1, C_Z.E_2 * C_I_K_j.C_I_K_1_dλ_2_2) * ξ₁^-2 +
           max(C_Z.P_1 * C_I_K_j.C_I_K_2_dλ_2_1, C_Z.P_2 * C_I_K_j.C_I_K_2_dλ_2_2) * ξ₁^-4
end
