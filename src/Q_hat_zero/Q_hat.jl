"""
    Q_hat_zero(ν, κ, ϵ, ξ₁, Λ::CGLParams; tol::Float64 = 1e-11)

Compute the solution to the ODE on the interval ``[0, ξ₁]``. Returns a
vector with two complex values, where the first is the value at `ξ₁`
and the second is the derivative.
"""
function Q_hat_zero(
    ν::Union{Complex{T},Acb},
    κ::T,
    ϵ::T,
    ξ₁::T,
    Λ::CGLParams{T};
    tol::Float64 = 1e-11,
) where {T}
    Q_hat = if T == Arb
        Q_hat_zero_capd(real(ν), imag(ν), κ, ϵ, ξ₁, Λ; tol)
    else
        Q_hat_zero_float(real(ν), imag(ν), κ, ϵ, ξ₁, Λ; tol)
    end
    return SVector(_complex(Q_hat[1], Q_hat[2]), _complex(Q_hat[3], Q_hat[4]))
end

"""
    Q_hat_zero_jacobian(ν, κ, ϵ, ξ₁, Λ::CGLParams; tol::Float64 = 1e-11)

This function computes the Jacobian of [`Q_hat_zero`](@ref) w.r.t. the
parameter `ν.
"""
function Q_hat_zero_jacobian(
    ν::Union{Complex{T},Acb},
    κ::T,
    ϵ::T,
    ξ₁::T,
    Λ::CGLParams{T};
    tol::Float64 = 1e-11,
) where {T}
    J = if T == Arb
        Q_hat_zero_jacobian_capd(real(ν), imag(ν), κ, ϵ, ξ₁, Λ; tol)
    else
        Q_hat_zero_jacobian_float(real(ν), imag(ν), κ, ϵ, ξ₁, Λ; tol)
    end
    return SMatrix{2,1}(_complex(J[1, 1], J[2, 1]), _complex(J[3, 1], J[4, 1]))
end
