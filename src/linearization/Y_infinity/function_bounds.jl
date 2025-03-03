struct FunctionBounds_Y
    J_N::Arb
    Y_12::Arb
    Y_34::Arb
    K_1::Arb
    K_2::Arb

    function FunctionBounds_Y(
        C_ab::Arb,
        lambda::Acb,
        κ::Arb,
        ϵ::Arb,
        ξ₁::Arb,
        λ::CGLParams{Arb},
    )
        C = new(
            C_J_N(C_ab, lambda, κ, ϵ, ξ₁, λ),
            C_Y_12(lambda, κ, ϵ, ξ₁, λ),
            C_Y_34(lambda, κ, ϵ, ξ₁, λ),
            indeterminate(κ),
            indeterminate(κ),
        )

        C.K_1[] = C_K_1(lambda, κ, ϵ, ξ₁, λ, C)
        C.K_2[] = C_K_2(lambda, κ, ϵ, ξ₁, λ, C)

        return C
    end
end

function C_J_N(C_ab::Arb, lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    (; σ, δ) = λ
    # TODO: opnorm is most likely not rigorous. We probably want to
    # compute explicit values for them.
    p = 2
    return 2^(σ - 1) *
           C_ab^(σ + 1) *
           (
               opnorm(@SMatrix(Float64[δ*(1+2σ) 1; -(1 + 2σ) δ]), p) +
               2σ * opnorm(@SMatrix(Float64[1 δ; δ -1]), p) +
               opnorm(@SMatrix(Float64[δ (1+2σ); -1 δ*(1+2σ)]), p)
           )
end

function C_Y_12(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    return one(lambda) # FIXME
end

function C_Y_34(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    return one(lambda) # FIXME
end

function C_K_1(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, C::FunctionBounds_Y)
    return one(lambda) # FIXME
end

function C_K_2(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, C::FunctionBounds_Y)
    return one(lambda) # FIXME
end
