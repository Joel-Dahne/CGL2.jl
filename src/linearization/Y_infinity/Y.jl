"""
    Y_infinity(γ, κ, ϵ, ξ₁, λ::CGLParams)

Compute the solution to the ODE on the interval ``[ξ₁, ∞)``. Returns a
vector with two complex values, where the first is the value at `ξ₁`
and the second is the derivative.
"""
function Y_infinity(
    C_ab::Arb,
    c::SVector{2,Acb},
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
)
    v = Arb("0.1")

    C = FunctionBounds_Y(C_ab, lambda, κ, ϵ, ξ₁, λ)

    norms = NormBounds_Y(c, lambda, κ, ϵ, ξ₁, v, λ, C)

    # Compute zeroth order bounds
    Y = add_error.(zero.(c), norms.Y * exp(real_a2(κ, ϵ) * ξ₁^2) * ξ₁^v)
    dY = add_error.(zero.(c), norms.Y_dξ * exp(real_a2(κ, ϵ) * ξ₁^2) * ξ₁^v)

    for _ = 1:3
        # TODO
    end

    return vcat(Y, dY)
end
