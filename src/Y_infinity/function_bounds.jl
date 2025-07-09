struct FunctionBounds_Y
    J_N::Arb
    Y_12::Arb
    Y_34::Arb
    Y_12_dξ::Arb
    Y_34_dξ::Arb
    K_1::Arb
    K_2::Arb
    Y_12_dλ::Arb
    Y_34_dλ::Arb
    Y_12_dλ_dξ::Arb
    Y_34_dλ_dξ::Arb
    K_1_dλ::Arb
    K_2_dλ::Arb

    function FunctionBounds_Y(
        lambda::Acb,
        γ₁::Acb,
        γ₂::Acb,
        κ::Arb,
        ϵ::Arb,
        ξ₁::Arb,
        λ::CGLParams{Arb};
        include_dλ::Bool = false,
    )
        C_ab = Arb(1) # FIXME

        C = new(
            C_J_N(C_ab, lambda, κ, ϵ, ξ₁, λ),
            C_Y_12(lambda, κ, ϵ, ξ₁, λ),
            C_Y_34(lambda, κ, ϵ, ξ₁, λ),
            C_Y_12_dξ(lambda, κ, ϵ, ξ₁, λ),
            C_Y_34_dξ(lambda, κ, ϵ, ξ₁, λ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
        )

        C.K_1[] = C_K_1(lambda, κ, ϵ, ξ₁, λ, C)
        C.K_2[] = C_K_2(lambda, κ, ϵ, ξ₁, λ, C)

        if include_dλ
            C.Y_12_dλ[] = C_Y_12_dλ(lambda, κ, ϵ, ξ₁, λ, C)
            C.Y_34_dλ[] = C_Y_34_dλ(lambda, κ, ϵ, ξ₁, λ, C)
            C.Y_12_dλ_dξ[] = C_Y_12_dλ_dξ(lambda, κ, ϵ, ξ₁, λ, C)
            C.Y_34_dλ_dξ[] = C_Y_34_dλ_dξ(lambda, κ, ϵ, ξ₁, λ, C)
            C.K_1_dλ[] = C_K_1_dλ(lambda, κ, ϵ, ξ₁, λ, C)
            C.K_2_dλ[] = C_K_2_dλ(lambda, κ, ϵ, ξ₁, λ, C)
        end

        return C
    end
end

function C_J_N(C_ab::Arb, lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    (; σ, δ) = λ

    return 2^(σ - 1) *
           C_ab^(σ + 1) *
           (
               1 +
               abs(δ) * (1 + 2σ) +
               2σ * (1 + abs(δ)) +
               max(1 + 2σ * abs(δ), 1 + abs(δ) * (1 + 2σ))
           )
end

function C_Y_12(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    return one(lambda) # FIXME
end

function C_Y_34(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    return one(lambda) # FIXME
end

function C_Y_12_dξ(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    return one(lambda) # FIXME
end

function C_Y_34_dξ(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    return one(lambda) # FIXME
end


function C_Y_12_dλ(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    return one(lambda) # FIXME
end

function C_Y_34_dλ(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    return one(lambda) # FIXME
end

function C_Y_12_dλ_dξ(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    return one(lambda) # FIXME
end

function C_Y_34_dλ_dξ(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    return one(lambda) # FIXME
end

function C_K_1(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, C::FunctionBounds_Y)
    return one(lambda) # FIXME
end

function C_K_2(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, C::FunctionBounds_Y)
    return one(lambda) # FIXME
end

function C_K_1_dλ(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    return one(lambda) # FIXME
end

function C_K_2_dλ(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    return one(lambda) # FIXME
end
