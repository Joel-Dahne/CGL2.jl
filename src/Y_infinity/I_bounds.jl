"""
    I_K_j_Bounds(C_Y::FunctionBounds; include_dλ = false)

This contains the bounds related to the `I_K_1` and `I_K_2` functions
and related derivatives with respect to `lambda`.

More precisely it contains the bounds from

- Lemma REF(lemma:I_K_1-I_K_2-bounds)
- Lemma REF(lemma:Y-I_K_1_I_K_2-lambda-1-bounds)
- Lemma REF(lemma:Y-I_K_1_I_K_2-lambda-2-bounds)

If `include_dλ = false` it doesn't include the bounds corresponding to
derivatives in `lambda`.

It checks all the conditions on the parameters that these lemmas
assume. If any of these conditions are not satisfied it will throw an
error. When using this struct the bounds can therefore safely be
assume to hold.
"""
struct I_K_j_Bounds
    C_I_K_1_1::Arb
    C_I_K_1_2::Arb
    C_I_K_2_1::Arb
    C_I_K_2_2::Arb
    C_I_K_1_dλ_1_1::Arb
    C_I_K_1_dλ_1_2::Arb
    C_I_K_2_dλ_1_1::Arb
    C_I_K_2_dλ_1_2::Arb
    C_I_K_1_dλ_2_1::Arb
    C_I_K_1_dλ_2_2::Arb
    C_I_K_2_dλ_2_1::Arb
    C_I_K_2_dλ_2_2::Arb

    function I_K_j_Bounds(
        lambda::Acb,
        κ::Arb,
        ϵ::Arb,
        ξ₁::Arb,
        v::Arb,
        λ::CGLParams{Arb},
        C_Y::FunctionBounds_Y;
        include_dλ::Bool = false,
    )
        (; d, σ) = λ
        _, _, c = _abc(κ, ϵ, λ)
        exponent = 2 / σ - d - 2real(lambda) / κ + v - 4

        # These are requirements of all three lemmas
        real(c) > 0 || throw(ArgumentError("real(c) > 0 not satisfied"))
        v > 0 || throw(ArgumentError("v > 0 not satisfied"))
        v - 2 < 0 || throw(ArgumentError("v - 2 < 0 not satisfied"))
        exponent < 0 || throw(ArgumentError("exponent < 0 not satisfied"))

        # This is a requirement for Lemma
        # REF(lemma:Y-I_K_1_I_K_2-lambda-1-bounds)
        ξ₁ > exp(-inv(exponent)) > 1 ||
            throw(ArgumentError("ξ₁ > exp(-inv(exponent)) not satisfied"))

        return new(
            C_I_K_1_1(v, C_Y),
            C_I_K_1_2(v, C_Y),
            C_I_K_2_1(κ, ϵ, λ, C_Y),
            C_I_K_2_2(κ, ϵ, λ, C_Y),
            include_dλ ? C_I_K_1_dλ_1_1(ξ₁, v, C_Y) : indeterminate(κ),
            include_dλ ? C_I_K_1_dλ_1_2(ξ₁, v, C_Y) : indeterminate(κ),
            include_dλ ? C_I_K_2_dλ_1_1(κ, ϵ, λ, C_Y) : indeterminate(κ),
            include_dλ ? C_I_K_2_dλ_1_2(κ, ϵ, λ, C_Y) : indeterminate(κ),
            include_dλ ? C_I_K_1_dλ_2_1(v, C_Y) : indeterminate(κ),
            include_dλ ? C_I_K_1_dλ_2_2(v, C_Y) : indeterminate(κ),
            include_dλ ? C_I_K_2_dλ_2_1(κ, ϵ, λ, C_Y) : indeterminate(κ),
            include_dλ ? C_I_K_2_dλ_2_2(κ, ϵ, λ, C_Y) : indeterminate(κ),
        )
    end
end

# Lemma I_K_1-I_K_2-bounds

function C_I_K_1_1(v::Arb, C_Y::FunctionBounds_Y)
    return C_Y.K_1_1 * C_Y.J_N / abs(v - 2)
end

function C_I_K_1_2(v::Arb, C_Y::FunctionBounds_Y)
    return C_Y.K_1_2 * C_Y.J_N / abs(v - 2)
end

function C_I_K_2_1(κ::Arb, ϵ::Arb, λ::CGLParams{Arb}, C_Y::FunctionBounds_Y)
    _, _, c = _abc(κ, ϵ, λ)
    return C_Y.K_2_1 * C_Y.J_N / (2real(c))
end

function C_I_K_2_2(κ::Arb, ϵ::Arb, λ::CGLParams{Arb}, C_Y::FunctionBounds_Y)
    _, _, c = _abc(κ, ϵ, λ)
    return C_Y.K_2_2 * C_Y.J_N / (2real(c))
end

# Lemma Y-I_K_1_I_K_2-lambda-1-bounds

function C_I_K_1_dλ_1_1(ξ₁::Arb, v::Arb, C_Y::FunctionBounds_Y)
    return C_Y.K_1_dλ_1 * C_Y.J_N * (abs(v - 2) + inv(log(ξ₁))) / (v - 2)^2
end

function C_I_K_1_dλ_1_2(ξ₁::Arb, v::Arb, C_Y::FunctionBounds_Y)
    return C_Y.K_1_dλ_2 * C_Y.J_N * (abs(v - 2) + inv(log(ξ₁))) / (v - 2)^2
end


function C_I_K_2_dλ_1_1(κ::Arb, ϵ::Arb, λ::CGLParams{Arb}, C_Y::FunctionBounds_Y)
    _, _, c = _abc(κ, ϵ, λ)
    return C_Y.K_2_dλ_1 * C_Y.J_N / (2real(c))
end

function C_I_K_2_dλ_1_2(κ::Arb, ϵ::Arb, λ::CGLParams{Arb}, C_Y::FunctionBounds_Y)
    _, _, c = _abc(κ, ϵ, λ)
    return C_Y.K_2_dλ_2 * C_Y.J_N / (2real(c))
end

# Lemma Y-I_K_1_I_K_2-lambda-2-bounds

function C_I_K_1_dλ_2_1(v::Arb, C_Y::FunctionBounds_Y)
    return C_I_K_1_2(v, C_Y)
end

function C_I_K_1_dλ_2_2(v::Arb, C_Y::FunctionBounds_Y)
    return C_I_K_1_2(v, C_Y)
end

function C_I_K_2_dλ_2_1(κ::Arb, ϵ::Arb, λ::CGLParams{Arb}, C_Y::FunctionBounds_Y)
    return C_I_K_2_1(κ, ϵ, λ, C_Y)
end

function C_I_K_2_dλ_2_2(κ::Arb, ϵ::Arb, λ::CGLParams{Arb}, C_Y::FunctionBounds_Y)
    return C_I_K_2_2(κ, ϵ, λ, C_Y)
end
