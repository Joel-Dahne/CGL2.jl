struct NormBounds_Y
    Y::Arb
    Y_dξ::Arb

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
        norms = new(norm_bound_Y(c, lambda, κ, ϵ, ξ₁, v, λ, C), indeterminate(κ))

        norms.Y_dξ[] = norm_bound_Y_dξ(c, lambda, κ, ϵ, ξ₁, v, λ, C, norms)

        return norms
    end
end

norm_bound_Y(
    c::SVector{2,Acb},
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
) = Y_infinity_fixed_point(c, lambda, κ, ϵ, ξ₁, v, λ, C)

function norm_bound_Y_dξ(
    c::SVector{2,Acb},
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
    norms::NormBounds_Y,
)
    # TODO
    return one(Arb)
end
