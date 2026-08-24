"""
    I_K_2_enclosure(c_0, λ, κ, ϵ, ξ₁, Λ, C_Z, norms_Z)

Compute an enclosure of ``I_{K_2}(ξ)`` at the point `ξ = ξ₁`.

It uses the bound from Lemma REF(lemma:I_K_1-I_K_2-bounds).
"""
function I_K_2_enclosure(
    c_0::SVector{2,Acb},
    λ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    Λ::CGLParams{Arb},
    C_Z::FunctionBounds_Y,
    norms_Z::NormBounds_Y,
)
    (; d, σ) = Λ
    c = _c(κ, ϵ, Λ)

    exponent = 2 / σ - d - 4real(λ) - 4

    I_K_2_1_bound = C_Z.I_K_2_1 * exp(-real(c) * ξ₁^2) * ξ₁^exponent * norms_Z.Z
    I_K_2_2_bound = C_Z.I_K_2_2 * exp(-real(c) * ξ₁^2) * ξ₁^exponent * norms_Z.Z

    return add_error.(zero(c_0), SVector(I_K_2_1_bound, I_K_2_2_bound))
end
