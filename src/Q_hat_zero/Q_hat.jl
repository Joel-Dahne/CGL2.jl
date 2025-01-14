"""
    Q_hat_zero(ν, κ, ϵ, ξ₁, λ::CGLParams; tol::Float64 = 1e-11)

Compute the solution to the ODE on the interval ``[0, ξ₁]``. Returns a
vector with two complex values, where the first is the value at `ξ₁`
and the second is the derivative.
"""
function Q_hat_zero(
    ν::Complex{T},
    κ::T,
    ϵ::T,
    ξ₁::T,
    λ::CGLParams{T};
    tol::Float64 = 1e-11,
) where {T}
    Q = if T == Arb
        Q_hat_zero_capd(real(ν), imag(ν), κ, ϵ, ξ₁, λ; tol)
    else
        Q_hat_zero_float(real(ν), imag(ν), κ, ϵ, ξ₁, λ; tol)
    end
    return SVector(_complex(Q[1], Q[2]), _complex(Q[3], Q[4]))
end
