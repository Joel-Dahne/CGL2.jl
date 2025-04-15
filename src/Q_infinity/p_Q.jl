function p_Q(γ, κ, ϵ, ξ₁, λ::CGLParams{T}) where {T}
    a, b, c = CGL2._abc(κ, ϵ, λ)
    pQ = c^-a * γ

    if T == Arb # Only add error bounds if we are working with an Arb
        CU = CGL2.UBounds(a, b, c, ξ₁)
        C = CGL2.FunctionBounds(κ, ϵ, ξ₁, λ, CU)
        norms = CGL2.NormBounds(γ, κ, ϵ, ξ₁, v, λ, C)
        C_p_Q =
            abs(c^-a) *
            CGL2.C_I_E(κ, ϵ, ξ₁, v, λ, C) *
            norms.Q^(2σ + 1) *
            ξ₁^((2σ + 1) * v - 2)
        add_error(pQ, C_p_Q)
    end

    return pQ
end
