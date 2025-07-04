A1(T) = SVector{2,ifelse(T == Arb, Acb, Complex{T})}(1, im)
A2(T) = SVector{2,ifelse(T == Arb, Acb, Complex{T})}(1, -im)
A3(T) = SVector{2,ifelse(T == Arb, Acb, Complex{T})}(1, im)
A4(T) = SVector{2,ifelse(T == Arb, Acb, Complex{T})}(1, -im)

a1(κ, ϵ) = -κ * (ϵ + im) / 2(1 + ϵ^2)
a2(κ, ϵ) = -κ * (ϵ - im) / 2(1 + ϵ^2)

# Real part of a₁ and a₂
real_a12(κ, ϵ) = -κ * ϵ / 2(1 + ϵ^2)

s1(lambda, κ, λ::CGLParams) = -1 / λ.σ + λ.d + lambda / κ + λ.ω / κ * im
s2(lambda, κ, λ::CGLParams) = -1 / λ.σ + λ.d + lambda / κ - λ.ω / κ * im
s3(lambda, κ, λ::CGLParams) = 1 / λ.σ - lambda / κ - λ.ω / κ * im
s4(lambda, κ, λ::CGLParams) = 1 / λ.σ - lambda / κ + λ.ω / κ * im

s1_dλ(lambda, κ, λ::CGLParams) = 1 / κ
s2_dλ(lambda, κ, λ::CGLParams) = 1 / κ
s3_dλ(lambda, κ, λ::CGLParams) = -1 / κ
s4_dλ(lambda, κ, λ::CGLParams) = -1 / κ

# Real part of s₁ and s₂
real_s12(lambda, κ, λ::CGLParams) = -1 / λ.σ + λ.d + real(lambda) / κ
# Real part of s₃ and s₄
real_s34(lambda, κ, λ::CGLParams) = 1 / λ.σ - real(lambda) / κ

function _c2ns(N::Int, sgn1::Int, sgn2::Int, s, κ, ϵ, λ::CGLParams{T}) where {T}
    (; d) = λ

    c2ns = ifelse(T == Arb, Acb, Complex{T})[1]
    for n = 1:N
        c2n =
            sgn1 * ((s + 2n - 2)^2 - (d - 2) * (s + 2n - 2)) * (ϵ + sgn2 * im) / (2n * κ) *
            c2ns[end]
        push!(c2ns, c2n)
    end

    return c2ns
end

c2ns1(N, s, κ, ϵ, λ::CGLParams) = _c2ns(N, -1, -1, s, κ, ϵ, λ)
c2ns2(N, s, κ, ϵ, λ::CGLParams) = _c2ns(N, -1, 1, s, κ, ϵ, λ)
c2ns3(N, s, κ, ϵ, λ::CGLParams) = _c2ns(N, 1, -1, s, κ, ϵ, λ)
c2ns4(N, s, κ, ϵ, λ::CGLParams) = _c2ns(N, 1, 1, s, κ, ϵ, λ)

function _c2ns_dλ(N::Int, sgn1::Int, sgn2::Int, s, s_dλ, κ, ϵ, λ::CGLParams{T}) where {T}
    (; d) = λ

    c2ns = _c2ns(N, sgn1, sgn2, s, κ, ϵ, λ)

    c2ns_dλ = ifelse(T == Arb, Acb, Complex{T})[0]
    for n = 1:N
        c2n_dλ =
            sgn1 * (2s_dλ * (s + 2n - 2) - (d - 2) * s_dλ) * (ϵ + sgn2 * im) / (2n * κ) *
            c2ns[n] +
            sgn1 * ((s + 2n - 2)^2 - (d - 2) * (s + 2n - 2)) * (ϵ + sgn2 * im) / (2n * κ) *
            c2ns_dλ[end]
        push!(c2ns_dλ, c2n_dλ)
    end

    return c2ns_dλ
end

c2ns1_dλ(N, s, s_dλ, κ, ϵ, λ::CGLParams) = _c2ns_dλ(N, -1, -1, s, s_dλ, κ, ϵ, λ)
c2ns2_dλ(N, s, s_dλ, κ, ϵ, λ::CGLParams) = _c2ns_dλ(N, -1, 1, s, s_dλ, κ, ϵ, λ)
c2ns3_dλ(N, s, s_dλ, κ, ϵ, λ::CGLParams) = _c2ns_dλ(N, 1, -1, s, s_dλ, κ, ϵ, λ)
c2ns4_dλ(N, s, s_dλ, κ, ϵ, λ::CGLParams) = _c2ns_dλ(N, 1, 1, s, s_dλ, κ, ϵ, λ)

function C_R_Y_1(N::Integer, ξ, s, κ, ϵ, λ::CGLParams{T}) where {T}
    # FIXME: For now we just take 10 times the size of the last term
    # as the remainder
    return 10abs(c2ns1(N, s, κ, ϵ, λ)[end])
end

function C_R_Y_1_dξ(N::Integer, ξ, s, κ, ϵ, λ::CGLParams{T}) where {T}
    # FIXME: For now we just take 10 times the size of the last term
    # as the remainder
    return 10abs((2N + s) * c2ns1(N, s, κ, ϵ, λ)[end])
end

function C_R_Y_1_dξ_dξ(N::Integer, ξ, s, κ, ϵ, λ::CGLParams{T}) where {T}
    # FIXME: For now we just take 10 times the size of the last term
    # as the remainder
    return 10abs((2N + s) * (2N + 1 + s) * c2ns1(N, s, κ, ϵ, λ)[end])
end

function C_R_Y_2(N::Integer, ξ, s, κ, ϵ, λ::CGLParams{T}) where {T}
    # FIXME: For now we just take 10 times the size of the last term
    # as the remainder
    return 10abs(c2ns2(N, s, κ, ϵ, λ)[end])
end

function C_R_Y_2_dξ(N::Integer, ξ, s, κ, ϵ, λ::CGLParams{T}) where {T}
    # FIXME: For now we just take 10 times the size of the last term
    # as the remainder
    return 10abs((2N + s) * c2ns2(N, s, κ, ϵ, λ)[end])
end

function C_R_Y_2_dξ_dξ(N::Integer, ξ, s, κ, ϵ, λ::CGLParams{T}) where {T}
    # FIXME: For now we just take 10 times the size of the last term
    # as the remainder
    return 10abs((2N + s) * (2N + 1 + s) * c2ns2(N, s, κ, ϵ, λ)[end])
end

function C_R_Y_3(N::Integer, ξ, s, κ, ϵ, λ::CGLParams{T}) where {T}
    # FIXME: For now we just take 10 times the size of the last term
    # as the remainder
    return 10abs(c2ns3(N, s, κ, ϵ, λ)[end])
end

function C_R_Y_3_dξ(N::Integer, ξ, s, κ, ϵ, λ::CGLParams{T}) where {T}
    # FIXME: For now we just take 10 times the size of the last term
    # as the remainder
    return 10abs((2N + s) * c2ns3(N, s, κ, ϵ, λ)[end])
end

function C_R_Y_3_dξ_dξ(N::Integer, ξ, s, κ, ϵ, λ::CGLParams{T}) where {T}
    # FIXME: For now we just take 10 times the size of the last term
    # as the remainder
    return 10abs((2N + s) * (2N + 1 + s) * c2ns3(N, s, κ, ϵ, λ)[end])
end


function C_R_Y_4(N::Integer, ξ, s, κ, ϵ, λ::CGLParams{T}) where {T}
    # FIXME: For now we just take 10 times the size of the last term
    # as the remainder
    return 10abs(c2ns4(N, s, κ, ϵ, λ)[end])
end

function C_R_Y_4_dξ(N::Integer, ξ, s, κ, ϵ, λ::CGLParams{T}) where {T}
    # FIXME: For now we just take 10 times the size of the last term
    # as the remainder
    return 10abs((2N + s) * c2ns4(N, s, κ, ϵ, λ)[end])
end

function C_R_Y_4_dξ_dξ(N::Integer, ξ, s, κ, ϵ, λ::CGLParams{T}) where {T}
    # FIXME: For now we just take 10 times the size of the last term
    # as the remainder
    return 10abs((2N + s) * (2N + 1 + s) * c2ns4(N, s, κ, ϵ, λ)[end])
end
