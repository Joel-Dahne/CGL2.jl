A1(T) = SVector{2,ifelse(T == Arb, Acb, Complex{T})}(1, im)
A2(T) = SVector{2,ifelse(T == Arb, Acb, Complex{T})}(1, -im)
A3(T) = SVector{2,ifelse(T == Arb, Acb, Complex{T})}(1, im)
A4(T) = SVector{2,ifelse(T == Arb, Acb, Complex{T})}(1, -im)

a1(κ, ϵ) = -κ * (ϵ + im) / 2(1 + ϵ^2)
a2(κ, ϵ) = -κ * (ϵ - im) / 2(1 + ϵ^2)

s1(lambda, κ, λ::CGLParams) = -1 / λ.σ + λ.d + lambda / κ + λ.ω / κ * im
s2(lambda, κ, λ::CGLParams) = -1 / λ.σ + λ.d + lambda / κ - λ.ω / κ * im
s3(lambda, κ, λ::CGLParams) = 1 / λ.σ - lambda / κ - λ.ω / κ * im
s4(lambda, κ, λ::CGLParams) = 1 / λ.σ - lambda / κ + λ.ω / κ * im

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
