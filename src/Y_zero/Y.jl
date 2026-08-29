"""
    Y_zero(Y₀, λ, ν, κ, ϵ, ξ₁, Λ::CGLParams; tol::Float64 = 1e-11)
    Y_zero(Y₀, λ, κ, ϵ, ξ₁, Q_hat, Λ::CGLParams; tol::Float64 = 1e-11)

Compute the solution to the ODE on the interval ``[0, ξ₁]`` with
initial values given by `Y₀`. Returns a vector with four complex
values, where the first two are the values at `ξ₁` and the last two
are the derivatives.
"""
function Y_zero(
    Y₀::SVector{2,Acb},
    λ::Acb,
    ν::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    Λ::CGLParams{Arb};
    tol::Float64 = 1e-11,
)
    return Y_zero_capd(Y₀, λ, ν, κ, ϵ, ξ₁, Λ; tol)
end

function Y_zero(
    Y₀::SVector{2,ComplexF64},
    λ::ComplexF64,
    ν::ComplexF64,
    κ::Float64,
    ϵ::Float64,
    ξ₁::Float64,
    Λ::CGLParams{Float64};
    tol::Float64 = 1e-11,
)
    Q_hat = Q_hat_zero_float_curve(real(ν), imag(ν), κ, ϵ, ξ₁, Λ; tol)

    return Y_zero(Y₀, λ, κ, ϵ, ξ₁, Q_hat, Λ; tol)
end

function Y_zero(
    Y₀::SVector{2,ComplexF64},
    λ::ComplexF64,
    κ::Float64,
    ϵ::Float64,
    ξ₁::Float64,
    Q_hat,
    Λ::CGLParams{Float64};
    tol::Float64 = 1e-11,
)
    return Y_zero_float(Y₀, λ, κ, ϵ, ξ₁, Q_hat, Λ; tol)
end
