struct NormBounds_Y
    Y::Arb
    Y_dλ::Arb

    function NormBounds_Y(
        c::SVector{2,Acb},
        lambda::Acb,
        κ::Arb,
        ϵ::Arb,
        ξ₁::Arb,
        v::Arb,
        λ::CGLParams{Arb},
        C::FunctionBounds_Y,
        C_I_K_j::I_K_j_Bounds;
        include_dλ::Bool = false,
    )
        norms = new(norm_bound_Y(c, lambda, κ, ϵ, ξ₁, v, λ, C, C_I_K_j), indeterminate(κ))

        if include_dλ
            norms.Y_dλ[] = norm_bound_Y_dλ(c, lambda, κ, ϵ, ξ₁, v, λ, C, C_I_K_j, norms)
        end

        return norms
    end
end

"""
    norm_bound_Y(c, lambda, κ, ϵ, ξ₁, v, λ::CGLParams)

Compute an upper bound for the norm of `Y` using the fixed point
formulation.
"""
function norm_bound_Y(
    c::SVector{2,Acb},
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C_Y::FunctionBounds_Y,
    C_I_K_j::I_K_j_Bounds,
)
    C_T = C_T_12(lambda, κ, ϵ, ξ₁, v, λ, C_Y, C_I_K_j)
    if 2C_T * ξ₁^(-2) < 1
        return inv(1 - C_T * ξ₁^-2) * 2max(C_Y.E_1 * abs(c[1]), C_Y.E_2 * abs(c[2])) * ξ₁^-v
    else
        #@debug "Non-finite norm" 2C_T * ξ₁^(-2)
        return indeterminate(Arb)
    end
end

function norm_bound_Y_dλ(
    c::SVector{2,Acb},
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C_Y::FunctionBounds_Y,
    C_I_K_j::I_K_j_Bounds,
    norms_Y::NormBounds_Y,
)
    num = C_Y_dλ_1(c, v, C_Y) + C_Y_dλ_2(lambda, κ, ϵ, ξ₁, v, λ, C_Y, C_I_K_j) * norms_Y.Y
    den = 1 - C_Y_dλ_3(lambda, κ, ϵ, ξ₁, v, λ, C_Y, C_I_K_j)
    if den > 0
        return num / den
    else
        #@debug "Non-finite norm" den
        return indeterminate(Arb)
    end
end
