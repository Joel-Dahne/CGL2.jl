struct NormBounds_Y
    Z::Arb
    Z_dλ::Arb

    function NormBounds_Y(
        c::SVector{2,Acb},
        lambda::Acb,
        κ::Arb,
        ϵ::Arb,
        ξ₁::Arb,
        v::Arb,
        λ::CGLParams{Arb},
        C::FunctionBounds_Y,
    )
        norms = new(norm_bound_Z(c, lambda, κ, ϵ, ξ₁, v, λ, C), indeterminate(Arb))
        norms.Z_dλ[] = norm_bound_Z_dλ(c, lambda, κ, ϵ, ξ₁, v, λ, C, norms)
        return norms
    end
end

"""
    norm_bound_Z(c, lambda, κ, ϵ, ξ₁, v, λ::CGLParams, C_Z::FunctionBounds_Y)

Compute an upper bound for the norm of `Z` using the fixed point
formulation.
"""
function norm_bound_Z(
    c::SVector{2,Acb},
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C_Z::FunctionBounds_Y,
)
    C_T = C_T_12(lambda, κ, ϵ, ξ₁, v, λ, C_Z)
    if 2C_T * ξ₁^(-2) < 1
        return inv(1 - C_T * ξ₁^-2) * max(C_Z.E_1 * abs(c[1]), C_Z.E_2 * abs(c[2])) * ξ₁^-v
    else
        #@debug "Non-finite norm" 2C_T * ξ₁^(-2)
        return indeterminate(Arb)
    end
end

function norm_bound_Z_dλ(
    c::SVector{2,Acb},
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C_Z::FunctionBounds_Y,
    norms_Z::NormBounds_Y,
)
    num = C_Z_dλ_1(c, v, C_Z) + C_Z_dλ_2(lambda, κ, ϵ, ξ₁, v, λ, C_Z) * norms_Z.Z
    den = 1 - C_Z_dλ_3(lambda, κ, ϵ, ξ₁, v, λ, C_Z)
    if den > 0
        return num / den
    else
        #@debug "Non-finite norm" den
        return indeterminate(Arb)
    end
end
