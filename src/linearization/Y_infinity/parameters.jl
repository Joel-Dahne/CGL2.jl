# Matrices
function coeff_matrices(lambda, κ, ϵ, λ::CGLParams)
    (; ω, σ, d) = λ

    A = SMatrix{2,2}(ϵ, 1, -1, ϵ)
    B₁ = SMatrix{2,2}(κ, 0, 0, κ)
    B₂ = (d - 1) * A
    C = SMatrix{2,2}(κ / σ, ω, -ω, κ / σ)
    λI = SMatrix{2,2}(lambda, 0, 0, lambda)

    return (; A, B₁, B₂, C, λI)
end

function J_N_coeff_matrices(λ::CGLParams)
    (; σ, δ) = λ

    M1 = SMatrix{2,2}(δ * (1 + 2σ), -(1 + 2σ), 1, δ)
    M2 = SMatrix{2,2}(1, δ, δ, -1)
    M3 = SMatrix{2,2}(δ, -1, 1 + 2σ, δ * (1 + 2σ))

    return (; M1, M2, M3)
end

# Parameters for asymptotic expansion of linear equation

A01(T) = SVector{2,ifelse(T == Arb, Acb, Complex{T})}(1, im)
A02(T) = SVector{2,ifelse(T == Arb, Acb, Complex{T})}(1, -im)

a21(κ, ϵ) = -κ * (ϵ + im) / 2(1 + ϵ^2)
a22(κ, ϵ) = -κ * (ϵ - im) / 2(1 + ϵ^2)

real_a2(κ, ϵ) = -κ * ϵ / 2(1 + ϵ^2)

s11(lambda, κ, λ::CGLParams) = -1 / λ.σ + λ.d + lambda / κ + λ.ω / κ * im
s12(lambda, κ, λ::CGLParams) = -1 / λ.σ + λ.d + lambda / κ - λ.ω / κ * im
s21(lambda, κ, λ::CGLParams) = 1 / λ.σ - lambda / κ - λ.ω / κ * im
s22(lambda, κ, λ::CGLParams) = 1 / λ.σ - lambda / κ + λ.ω / κ * im

real_s1(lambda, κ, λ::CGLParams) = -1 / λ.σ + λ.d + real(lambda) / κ
real_s2(lambda, κ, λ::CGLParams) = 1 / λ.σ - real(lambda) / κ

function _cs_odd(N, lambda, κ, ϵ, λ::CGLParams{T}; s, sgn1, sgn2) where {T}
    (; d) = λ

    cs_odd = ifelse(T == Arb, Acb, Complex{T})[1]
    for n = 2:2:N
        cn =
            sgn1 * ((s + n - 2)^2 - (d - 2) * (s + n - 2)) * (ϵ + sgn2 * im) / (n * κ) *
            cs_odd[end]
        push!(cs_odd, cn)
    end

    return cs_odd
end

cs_odd_11(N, lambda, κ, ϵ, λ::CGLParams{T}) where {T} =
    _cs_odd(N, lambda, κ, ϵ, λ; s = s11(lambda, κ, λ), sgn1 = -1, sgn2 = -1)
cs_odd_12(N, lambda, κ, ϵ, λ::CGLParams{T}) where {T} =
    _cs_odd(N, lambda, κ, ϵ, λ; s = s12(lambda, κ, λ), sgn1 = -1, sgn2 = 1)
cs_odd_21(N, lambda, κ, ϵ, λ::CGLParams{T}) where {T} =
    _cs_odd(N, lambda, κ, ϵ, λ; s = s21(lambda, κ, λ), sgn1 = 1, sgn2 = -1)
cs_odd_22(N, lambda, κ, ϵ, λ::CGLParams{T}) where {T} =
    _cs_odd(N, lambda, κ, ϵ, λ; s = s22(lambda, κ, λ), sgn1 = 1, sgn2 = 1)
