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
        C::Union{FunctionBounds_Y,FunctionBounds_Y_new};
        include_dλ::Bool = false,
    )
        norms = new(norm_bound_Y(c, lambda, κ, ϵ, ξ₁, v, λ, C), indeterminate(κ))

        if include_dλ
            norms.Y_dλ[] = norm_bound_Y_dλ(c, lambda, κ, ϵ, ξ₁, v, λ, C, norms)
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
)
    C_T = C_T_Y(lambda, κ, ϵ, ξ₁, v, λ, C_Y)

    if 2C_T * ξ₁^(-2) < 1
        return C_Y.Y_12 * norm_inf(c) * ξ₁^-v / (1 - C_T * ξ₁^-2)
    else
        #@debug "Non-finite norm" 2C_T * ξ₁^(-2)
        return indeterminate(Arb)
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
    C_Y::FunctionBounds_Y_new,
)
    C_T = C_T_12(lambda, κ, ϵ, ξ₁, v, λ, C_Y)
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
    norms_Y::NormBounds_Y,
)
    num = C_Y_dλ_1(v, C_Y) * norm_inf(c) + C_Y_dλ_2(lambda, κ, ξ₁, v, λ, C_Y) * norms_Y.Y
    den = 1 - C_Y_dλ_3(lambda, κ, ξ₁, v, λ, C_Y)

    if den > 0
        return num / den
    else
        #@debug "Non-finite norm" den
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
    C_Y::FunctionBounds_Y_new,
    norms_Y::NormBounds_Y,
)
    num = C_Y_dλ_1(c, v, C_Y) + C_Y_dλ_2(lambda, κ, ϵ, ξ₁, v, λ, C_Y) * norms_Y.Y
    den = 1 - C_Y_dλ_3(lambda, κ, ϵ, ξ₁, v, λ, C_Y)
    if den > 0
        return num / den
    else
        #@debug "Non-finite norm" den
        return indeterminate(Arb)
    end
end
