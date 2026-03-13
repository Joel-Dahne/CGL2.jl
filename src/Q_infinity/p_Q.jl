function p_Q(γ, κ, ϵ, ξ₁, λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, λ)
    pQ = c^-a * γ

    if T == Arb # Only add error bounds if we are working with an Arb
        v = Arb("0.001") # IMPROVE: Can we take v = 0 here?
        (; σ) = λ

        # IMPROVE: We need to improve this bound, most likely by
        # improving C_I_E.
        CU = UBounds(a, b, c, ξ₁)
        C = FunctionBounds(κ, ϵ, ξ₁, λ, CU)
        norms = NormBounds(γ, κ, ϵ, ξ₁, v, λ, C)
        C_p_Q =
            abs(c^-a) * C_I_E(κ, ϵ, ξ₁, v, λ, C) * norms.Q^(2σ + 1) * ξ₁^((2σ + 1) * v - 2)
        pQ = add_error(pQ, C_p_Q)
    end

    return pQ
end
