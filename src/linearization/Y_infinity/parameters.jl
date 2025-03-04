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

# Parameters for Case I, a₂ != 0

a21(κ, ϵ) = -κ * (ϵ + im) / 2(1 + ϵ^2)
a22(κ, ϵ) = -κ * (ϵ - im) / 2(1 + ϵ^2)

real_a2(κ, ϵ) = -κ * ϵ / 2(1 + ϵ^2)

s11(lambda, κ, λ::CGLParams) = -1 / λ.σ + λ.d + lambda / κ + λ.ω / κ * im
s12(lambda, κ, λ::CGLParams) = -1 / λ.σ + λ.d + lambda / κ - λ.ω / κ * im

real_s1(lambda, κ, λ::CGLParams) = -1 / λ.σ + λ.d + real(lambda) / κ

A011(T) =
    if T == Arb
        SVector{2,Acb}(1, im)
    else
        SVector{2,Complex{T}}(1, im)
    end
A012(T) =
    if T == Arb
        SVector{2,Acb}(1, -im)
    else
        SVector{2,Complex{T}}(1, -im)
    end

# TODO: We probably want a method returning only the coefficients
function _As_1(N, lambda, κ, ϵ, λ::CGLParams{T}; s, A0, sgn) where {T}
    (; d) = λ

    (; A, B₁, C, λI) = coeff_matrices(lambda, κ, ϵ, λ)

    cs = eltype(A0)[1, 0]

    for n = 2:N
        cn =
            -((s + n - 2)^2 - (s + n - 2) * (d - 2)) * (ϵ + sgn * im) / (n * κ) *
            cs[(n-2)+1]
        push!(cs, cn)
    end

    return cs .* Ref(A0)
end

As_11(N, lambda, κ, ϵ, λ::CGLParams{T}) where {T} =
    _As_1(N, lambda, κ, ϵ, λ; s = s11(lambda, κ, λ), A0 = A011(T), sgn = -1)

As_12(N, lambda, κ, ϵ, λ::CGLParams{T}) where {T} =
    _As_1(N, lambda, κ, ϵ, λ; s = s12(lambda, κ, λ), A0 = A012(T), sgn = 1)

# Parameters for Case II, a₂ = 0

s21(lambda, κ, λ::CGLParams) = 1 / λ.σ - lambda / κ + λ.ω / κ * im
s22(lambda, κ, λ::CGLParams) = 1 / λ.σ - lambda / κ - λ.ω / κ * im

real_s2(lambda, κ, λ::CGLParams) = 1 / λ.σ - real(lambda) / κ

A021(T) =
    if T == Arb
        SVector{2,Acb}(1, -im)
    else
        SVector{2,Complex{T}}(1, -im)
    end
A022(T) =
    if T == Arb
        SVector{2,Acb}(1, im)
    else
        SVector{2,Complex{T}}(1, im)
    end

function _As_2(N, lambda, κ, ϵ, λ::CGLParams{T}; s, A0) where {T}
    (; d) = λ

    (; A, B₁, C, λI) = coeff_matrices(lambda, κ, ϵ, λ)

    As = [A0, zero(A0)]
    for n = 2:N
        An =
            -((s + n - 2)^2 - (s + n - 2) * (d - 2)) *
            inv(C - λI - (s + n) * B₁) *
            A *
            As[(n-2)+1]
        push!(As, An)
    end

    return As
end

As_21(N, lambda, κ, ϵ, λ::CGLParams{T}) where {T} =
    _As_2(N, lambda, κ, ϵ, λ; s = s21(lambda, κ, λ), A0 = A021(T))

As_22(N, lambda, κ, ϵ, λ::CGLParams{T}) where {T} =
    _As_2(N, lambda, κ, ϵ, λ; s = s22(lambda, κ, λ), A0 = A022(T))
