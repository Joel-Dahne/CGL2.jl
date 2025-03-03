function C_T_Y(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    real_a₂ = real_a2(κ, ϵ)
    real_s₂ = real_s2(lambda, κ, λ)

    @assert v > 0
    @assert real_a₂ < 0
    @assert real_s₂ > 0
    @assert v + real_s₂ - 2 < 0

    return C.Y_12 * C.K_1 * C.J_N / abs(v - 2) +
           C.Y_34 * C.K_2 * C.J_N / abs(v + real_a₂ - 2)
end

"""
    Y_infinity_fixed_point(c, lambda, κ, ϵ, ξ₁, v, λ::CGLParams)

Consider the fixed point problem that is used in the paper to get a
solution on ``[ξ₁, ∞)``. This function computes `ρ` such that there
exists a solution in the ball of radius `ρ`, this solution is globally
unique.

We have a unique fixed points in a ball of radius `ρ` if
```
ρ >= C_Y_12 * norm(c) * ξ₁^-v / (1 - C_T * ξ₁^-2)
```
and
```
2C_T * ξ₁^(-2) < 1
```
We therefore only have to check that the second inequality holds, and
then the lower bound for `ρ` is given by the first inequality.
"""
function Y_infinity_fixed_point(
    c::SVector{2,Acb},
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    v::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y;
    throw_on_failure::Bool = false,
)
    (; σ) = λ
    @assert v > 0 # Required for the below bounds to be valid

    C_T = C_T_Y(lambda, κ, ϵ, ξ₁, v, λ, C)

    if !(2C_T * ξ₁^(-2) < 1)
        if throw_on_failure
            throw(
                ErrorException(
                    "second inequality not satisfied, 2C_T * ξ₁^(-2) = $(2C_T * ξ₁^(-2))",
                ),
            )
        else
            return indeterminate(κ)
        end
    end

    return C.Y_12 * norm(c) * ξ₁^-v / (1 - C_T * ξ₁^-2)
end
