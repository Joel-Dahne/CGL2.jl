"""
    Y_zero(Y₀, lambda, ν, κ, ϵ, ξ₁, λ::CGLParams; tol::Float64 = 1e-11)
    Y_zero(Y₀, lambda, κ, ϵ, ξ₁, Q_hat, λ::CGLParams; tol::Float64 = 1e-11)

Compute the solution to the ODE on the interval ``[0, ξ₁]`` with
initial values given by `Y₀`. Returns a vector with four complex
values, where the first two are the values at `ξ₁` and the last two
are the derivatives.
"""
function Y_zero(
    Y₀::SVector{2,Acb},
    lambda::Acb,
    ν::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    tol::Float64 = 1e-11,
) where {Arb}
    return Y_zero_capd(Y₀, lambda, ν, κ, ϵ, ξ₁, λ; tol)
end

function Y_zero(
    Y₀::SVector{2,ComplexF64},
    lambda::ComplexF64,
    ν::ComplexF64,
    κ::Float64,
    ϵ::Float64,
    ξ₁::Float64,
    λ::CGLParams{Float64};
    tol::Float64 = 1e-11,
)
    Q_hat = CGL2.Q_hat_zero_float_curve(real(ν), imag(ν), κ, ϵ, ξ₁, λ)

    return Y_zero(Y₀, lambda, κ, ϵ, ξ₁, Q_hat, λ)
end

function Y_zero(
    Y₀::SVector{2,ComplexF64},
    lambda::ComplexF64,
    κ::Float64,
    ϵ::Float64,
    ξ₁::Float64,
    Q_hat,
    λ::CGLParams{Float64};
    tol::Float64 = 1e-11,
)
    return Y_zero_float(Y₀, lambda, κ, ϵ, ξ₁, Q_hat, λ; tol)
end
