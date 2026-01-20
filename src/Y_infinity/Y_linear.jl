function _Y_k(ξ, Aₖ, sₖ, aₖ, c₂ₙₖs, C_R_Y_k)
    N = length(c₂ₙₖs)
    ns = 0:(N-1)

    S = sum(zip(ns, c₂ₙₖs)) do (n, c₂ₙₖ)
        c₂ₙₖ * ξ^(-2n)
    end

    if S isa Arblib.AcbOrRef # FIXME: Handle AcbSeries
        S = add_error(S, C_R_Y_k * abs(ξ^(-2N)))
    end

    if iszero(aₖ)
        return ξ^(-sₖ) * S * Aₖ
    else
        return exp(aₖ * ξ^2) * ξ^(-sₖ) * S * Aₖ
    end
end

function _Y_k_dξ(ξ, Aₖ, sₖ, aₖ, c₂ₙₖs, C_R_Y_k, C_R_Y_k_dξ)
    N = length(c₂ₙₖs)
    ns = 0:(N-1)

    # Derivative of sum w.r.t ξ
    S_dξ = -sum(zip(ns, c₂ₙₖs)) do (n, c₂ₙₖ)
        (2n + sₖ) * c₂ₙₖ * ξ^(-2n - 1)
    end

    if S_dξ isa Arblib.AcbOrRef # FIXME: Handle AcbSeries
        S_dξ = add_error(S_dξ, C_R_Y_k_dξ * abs(ξ^(-2N - 1)))
    end

    if iszero(aₖ)
        return ξ^(-sₖ) * S_dξ * Aₖ
    else
        S = sum(zip(ns, c₂ₙₖs)) do (n, c₂ₙₖ)
            c₂ₙₖ * ξ^(-2n)
        end

        if S isa Arblib.AcbOrRef # FIXME: Handle AcbSeries
            S = add_error(S, C_R_Y_k * abs(ξ^(-2N)))
        end

        return exp(aₖ * ξ^2) * ξ^(-sₖ) * (2aₖ * ξ * S + S_dξ) * Aₖ
    end
end

function _Y_k_dξ_dξ(ξ, Aₖ, sₖ, aₖ, c₂ₙₖs, C_R_Y_k, C_R_Y_k_dξ, C_R_Y_k_dξ_dξ)
    N = length(c₂ₙₖs)
    ns = 0:(N-1)

    S_dξ_dξ = sum(zip(ns, c₂ₙₖs)) do (n, c₂ₙₖ)
        (2n + sₖ) * (2n + 1 + sₖ) * c₂ₙₖ * ξ^(-2n - 2)
    end

    if S_dξ_dξ isa Arblib.AcbOrRef # FIXME: Handle AcbSeries
        S_dξ_dξ = add_error(S_dξ_dξ, C_R_Y_k_dξ_dξ * abs(ξ^(-2N - 2)))
    end

    if iszero(aₖ)
        return ξ^(-sₖ) * S_dξ_dξ * Aₖ
    else
        S = sum(zip(ns, c₂ₙₖs)) do (n, c₂ₙₖ)
            c₂ₙₖ * ξ^(-2n)
        end
        S_dξ = -sum(zip(ns, c₂ₙₖs)) do (n, c₂ₙₖ)
            (2n + sₖ) * c₂ₙₖ * ξ^(-2n - 1)
        end

        if S isa Arblib.AcbOrRef # FIXME: Handle AcbSeries
            S = add_error(S, C_R_Y_k * abs(ξ^(-2N)))
            S_dξ = add_error(S_dξ, C_R_Y_k_dξ * abs(ξ^(-2N - 1)))
        end

        return exp(aₖ * ξ^2) *
               ξ^(-sₖ) *
               (((2aₖ * ξ)^2 + 2aₖ) * S + 4aₖ * ξ * S_dξ + S_dξ_dξ) *
               Aₖ
    end
end

function _Y_k_dλ(ξ, Aₖ, sₖ, sₖ_dλ, aₖ, c₂ₙₖs, c₂ₙₖs_dλ, C_R_Y_k)
    N = length(c₂ₙₖs)
    ns = 0:(N-1)

    S = sum(zip(ns, c₂ₙₖs)) do (n, c₂ₙₖ)
        c₂ₙₖ * ξ^(-2n)
    end

    S_dλ = sum(zip(ns, c₂ₙₖs_dλ)) do (n, c₂ₙₖ_dλ)
        c₂ₙₖ_dλ * ξ^(-2n)
    end

    if S isa Arblib.AcbOrRef # FIXME: Handle AcbSeries
        S = add_error(S, C_R_Y_k * abs(ξ^(-2N)))
        S_dλ = add_error(S_dλ, C_R_Y_k * abs(ξ^(-2N))) # FIXME
    end

    if iszero(aₖ)
        return ξ^(-sₖ) * (-sₖ_dλ * log(ξ) * S + S_dλ) * Aₖ
    else
        return exp(aₖ * ξ^2) * ξ^(-sₖ) * (-sₖ_dλ * log(ξ) * S + S_dλ) * Aₖ
    end
end

function _Y_k_dλ_dξ(ξ, Aₖ, sₖ, sₖ_dλ, aₖ, c₂ₙₖs, c₂ₙₖs_dλ, C_R_Y_k, C_R_Y_k_dξ)
    N = length(c₂ₙₖs)
    ns = 0:(N-1)

    S_dξ = -sum(zip(ns, c₂ₙₖs)) do (n, c₂ₙₖ)
        (2n + sₖ) * c₂ₙₖ * ξ^(-2n - 1)
    end

    S_dλ_dξ = -sum(zip(ns, c₂ₙₖs, c₂ₙₖs_dλ)) do (n, c₂ₙₖ, c₂ₙₖ_dλ)
        (sₖ_dλ * c₂ₙₖ + (2n + sₖ) * c₂ₙₖ_dλ) * ξ^(-2n - 1)
    end

    if S_dξ isa Arblib.AcbOrRef # FIXME: Handle AcbSeries
        S_dξ = add_error(S_dξ, C_R_Y_k_dξ * abs(ξ^(-2N - 1)))
        S_dλ_dξ = add_error(S_dλ_dξ, C_R_Y_k_dξ * abs(ξ^(-2N - 1))) # FIXME
    end

    if iszero(aₖ)
        return ξ^(-sₖ) * (-sₖ_dλ * log(ξ) * S_dξ + S_dλ_dξ) * Aₖ
    else
        S = sum(zip(ns, c₂ₙₖs)) do (n, c₂ₙₖ)
            c₂ₙₖ * ξ^(-2n)
        end

        S_dλ = sum(zip(ns, c₂ₙₖs_dλ)) do (n, c₂ₙₖ_dλ)
            c₂ₙₖ_dλ * ξ^(-2n)
        end

        if S isa Arblib.AcbOrRef # FIXME: Handle AcbSeries
            S = add_error(S, C_R_Y_k * abs(ξ^(-2N)))
            S_dλ = add_error(S_dλ, C_R_Y_k * abs(ξ^(-2N))) # FIXME
        end

        return exp(aₖ * ξ^2) *
               ξ^(-sₖ) *
               (-sₖ_dλ * log(ξ) * (2aₖ * ξ * S + S_dξ) + 2aₖ * ξ * S_dλ + S_dλ_dξ) *
               Aₖ
    end
end

function Y_1(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    A₁ = A1(T)
    s₁ = s1(lambda, κ, λ)
    a₁ = a1(κ, ϵ)

    N = 10 # TODO: Choose N
    c2ns = c2ns1(N - 1, s₁, κ, ϵ, λ)

    C_R_Y_k = T == Arb ? C_R_Y_1(N, ξ, s₁, κ, ϵ, λ) : zero(T)

    return _Y_k(ξ, A₁, s₁, a₁, c2ns, C_R_Y_k)
end

function Y_1_dξ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    A₁ = A1(T)
    s₁ = s1(lambda, κ, λ)
    a₁ = a1(κ, ϵ)

    N = 10 # TODO: Choose N
    c2ns = c2ns1(N - 1, s₁, κ, ϵ, λ)

    C_R_Y_k = T == Arb ? C_R_Y_1(N, ξ, s₁, κ, ϵ, λ) : zero(T)
    C_R_Y_k_dξ = T == Arb ? C_R_Y_1_dξ(N, ξ, s₁, κ, ϵ, λ) : zero(T)

    return _Y_k_dξ(ξ, A₁, s₁, a₁, c2ns, C_R_Y_k, C_R_Y_k_dξ)
end

function Y_1_dξ_dξ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    A₁ = A1(T)
    s₁ = s1(lambda, κ, λ)
    a₁ = a1(κ, ϵ)

    N = 10 # TODO: Choose N
    c2ns = c2ns1(N - 1, s₁, κ, ϵ, λ)

    C_R_Y_k = T == Arb ? C_R_Y_1(N, ξ, s₁, κ, ϵ, λ) : zero(T)
    C_R_Y_k_dξ = T == Arb ? C_R_Y_1_dξ(N, ξ, s₁, κ, ϵ, λ) : zero(T)
    C_R_Y_k_dξ_dξ = T == Arb ? C_R_Y_1_dξ_dξ(N, ξ, s₁, κ, ϵ, λ) : zero(T)

    return _Y_k_dξ_dξ(ξ, A₁, s₁, a₁, c2ns, C_R_Y_k, C_R_Y_k_dξ, C_R_Y_k_dξ_dξ)
end

function Y_2(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    A₂ = A2(T)
    s₂ = s2(lambda, κ, λ)
    a₂ = a2(κ, ϵ)

    N = 10 # TODO: Choose N
    c2ns = c2ns2(N - 1, s₂, κ, ϵ, λ)

    C_R_Y_k = T == Arb ? C_R_Y_2(N, ξ, s₂, κ, ϵ, λ) : zero(T)

    return _Y_k(ξ, A₂, s₂, a₂, c2ns, C_R_Y_k)
end

function Y_2_dξ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    A₂ = A2(T)
    s₂ = s2(lambda, κ, λ)
    a₂ = a2(κ, ϵ)

    N = 10 # TODO: Choose N
    c2ns = c2ns2(N - 1, s₂, κ, ϵ, λ)

    C_R_Y_k = T == Arb ? C_R_Y_2(N, ξ, s₂, κ, ϵ, λ) : zero(T)
    C_R_Y_k_dξ = T == Arb ? C_R_Y_2_dξ(N, ξ, s₂, κ, ϵ, λ) : zero(T)

    return _Y_k_dξ(ξ, A₂, s₂, a₂, c2ns, C_R_Y_k, C_R_Y_k_dξ)
end

function Y_2_dξ_dξ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    A₂ = A2(T)
    s₂ = s2(lambda, κ, λ)
    a₂ = a2(κ, ϵ)

    N = 10 # TODO: Choose N
    c2ns = c2ns2(N - 1, s₂, κ, ϵ, λ)

    C_R_Y_k = T == Arb ? C_R_Y_2(N, ξ, s₂, κ, ϵ, λ) : zero(T)
    C_R_Y_k_dξ = T == Arb ? C_R_Y_2_dξ(N, ξ, s₂, κ, ϵ, λ) : zero(T)
    C_R_Y_k_dξ_dξ = T == Arb ? C_R_Y_2_dξ_dξ(N, ξ, s₂, κ, ϵ, λ) : zero(T)

    return _Y_k_dξ_dξ(ξ, A₂, s₂, a₂, c2ns, C_R_Y_k, C_R_Y_k_dξ, C_R_Y_k_dξ_dξ)
end

function Y_3(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    A₃ = A3(T)
    s₃ = s3(lambda, κ, λ)
    a₃ = zero(s₃)

    N = 10 # TODO: Choose N
    c2ns = c2ns3(N - 1, s₃, κ, ϵ, λ)

    C_R_Y_k = T == Arb ? C_R_Y_3(N, ξ, s₃, κ, ϵ, λ) : zero(T)

    return _Y_k(ξ, A₃, s₃, a₃, c2ns, C_R_Y_k)
end

function Y_3_dξ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    A₃ = A3(T)
    s₃ = s3(lambda, κ, λ)
    a₃ = zero(s₃)

    N = 10 # TODO: Choose N
    c2ns = c2ns3(N - 1, s₃, κ, ϵ, λ)

    C_R_Y_k = T == Arb ? C_R_Y_3(N, ξ, s₃, κ, ϵ, λ) : zero(T)
    C_R_Y_k_dξ = T == Arb ? C_R_Y_3_dξ(N, ξ, s₃, κ, ϵ, λ) : zero(T)

    return _Y_k_dξ(ξ, A₃, s₃, a₃, c2ns, C_R_Y_k, C_R_Y_k_dξ)
end

function Y_3_dξ_dξ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    A₃ = A3(T)
    s₃ = s3(lambda, κ, λ)
    a₃ = zero(s₃)

    N = 10 # TODO: Choose N
    c2ns = c2ns3(N - 1, s₃, κ, ϵ, λ)

    C_R_Y_k = T == Arb ? C_R_Y_3(N, ξ, s₃, κ, ϵ, λ) : zero(T)
    C_R_Y_k_dξ = T == Arb ? C_R_Y_3_dξ(N, ξ, s₃, κ, ϵ, λ) : zero(T)
    C_R_Y_k_dξ_dξ = T == Arb ? C_R_Y_3_dξ_dξ(N, ξ, s₃, κ, ϵ, λ) : zero(T)

    return _Y_k_dξ_dξ(ξ, A₃, s₃, a₃, c2ns, C_R_Y_k, C_R_Y_k_dξ, C_R_Y_k_dξ_dξ)
end

function Y_4(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    A₄ = A4(T)
    s₄ = s4(lambda, κ, λ)
    a₄ = zero(s₄)

    N = 10 # TODO: Choose N
    c2ns = c2ns4(N - 1, s₄, κ, ϵ, λ)

    C_R_Y_k = T == Arb ? C_R_Y_4(N, ξ, s₄, κ, ϵ, λ) : zero(T)

    return _Y_k(ξ, A₄, s₄, a₄, c2ns, C_R_Y_k)
end

function Y_4_dξ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    A₄ = A4(T)
    s₄ = s4(lambda, κ, λ)
    a₄ = zero(s₄)

    N = 10 # TODO: Choose N
    c2ns = c2ns4(N - 1, s₄, κ, ϵ, λ)

    C_R_Y_k = T == Arb ? C_R_Y_4(N, ξ, s₄, κ, ϵ, λ) : zero(T)
    C_R_Y_k_dξ = T == Arb ? C_R_Y_4_dξ(N, ξ, s₄, κ, ϵ, λ) : zero(T)

    return _Y_k_dξ(ξ, A₄, s₄, a₄, c2ns, C_R_Y_k, C_R_Y_k_dξ)
end

function Y_4_dξ_dξ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    A₄ = A4(T)
    s₄ = s4(lambda, κ, λ)
    a₄ = zero(s₄)

    N = 10 # TODO: Choose N
    c2ns = c2ns4(N - 1, s₄, κ, ϵ, λ)

    C_R_Y_k = T == Arb ? C_R_Y_4(N, ξ, s₄, κ, ϵ, λ) : zero(T)
    C_R_Y_k_dξ = T == Arb ? C_R_Y_4_dξ(N, ξ, s₄, κ, ϵ, λ) : zero(T)
    C_R_Y_k_dξ_dξ = T == Arb ? C_R_Y_4_dξ_dξ(N, ξ, s₄, κ, ϵ, λ) : zero(T)

    return _Y_k_dξ_dξ(ξ, A₄, s₄, a₄, c2ns, C_R_Y_k, C_R_Y_k_dξ, C_R_Y_k_dξ_dξ)
end

function Y_1_dλ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    A₁ = A1(T)
    s₁ = s1(lambda, κ, λ)
    s₁_dλ = s1_dλ(lambda, κ, λ)
    a₁ = a1(κ, ϵ)

    N = 10 # TODO: Choose N
    c2ns = c2ns1(N - 1, s₁, κ, ϵ, λ)
    c2ns_dλ = c2ns1_dλ(N - 1, s₁, s₁_dλ, κ, ϵ, λ)

    C_R_Y_k = T == Arb ? C_R_Y_1(N, ξ, s₁, κ, ϵ, λ) : zero(T)

    return _Y_k_dλ(ξ, A₁, s₁, s₁_dλ, a₁, c2ns, c2ns_dλ, C_R_Y_k)
end

function Y_2_dλ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    A₂ = A2(T)
    s₂ = s2(lambda, κ, λ)
    s₂_dλ = s2_dλ(lambda, κ, λ)
    a₂ = a2(κ, ϵ)

    N = 10 # TODO: Choose N
    c2ns = c2ns2(N - 1, s₂, κ, ϵ, λ)
    c2ns_dλ = c2ns2_dλ(N - 1, s₂, s₂_dλ, κ, ϵ, λ)

    C_R_Y_k = T == Arb ? C_R_Y_2(N, ξ, s₂, κ, ϵ, λ) : zero(T)

    return _Y_k_dλ(ξ, A₂, s₂, s₂_dλ, a₂, c2ns, c2ns_dλ, C_R_Y_k)
end

function Y_3_dλ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    A₃ = A3(T)
    s₃ = s3(lambda, κ, λ)
    s₃_dλ = s3_dλ(lambda, κ, λ)
    a₃ = zero(s₃)

    N = 10 # TODO: Choose N
    c2ns = c2ns3(N - 1, s₃, κ, ϵ, λ)
    c2ns_dλ = c2ns3_dλ(N - 1, s₃, s₃_dλ, κ, ϵ, λ)

    C_R_Y_k = T == Arb ? C_R_Y_1(N, ξ, s₃, κ, ϵ, λ) : zero(T)

    return _Y_k_dλ(ξ, A₃, s₃, s₃_dλ, a₃, c2ns, c2ns_dλ, C_R_Y_k)
end

function Y_4_dλ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    A₄ = A4(T)
    s₄ = s4(lambda, κ, λ)
    s₄_dλ = s4_dλ(lambda, κ, λ)
    a₄ = zero(s₄)

    N = 10 # TODO: Choose N
    c2ns = c2ns4(N - 1, s₄, κ, ϵ, λ)
    c2ns_dλ = c2ns4_dλ(N - 1, s₄, s₄_dλ, κ, ϵ, λ)

    C_R_Y_k = T == Arb ? C_R_Y_1(N, ξ, s₄, κ, ϵ, λ) : zero(T)

    return _Y_k_dλ(ξ, A₄, s₄, s₄_dλ, a₄, c2ns, c2ns_dλ, C_R_Y_k)
end

function Y_1_dλ_dξ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    A₁ = A1(T)
    s₁ = s1(lambda, κ, λ)
    s₁_dλ = s1_dλ(lambda, κ, λ)
    a₁ = a1(κ, ϵ)

    N = 10 # TODO: Choose N
    c2ns = c2ns1(N - 1, s₁, κ, ϵ, λ)
    c2ns_dλ = c2ns1_dλ(N - 1, s₁, s₁_dλ, κ, ϵ, λ)

    C_R_Y_k = T == Arb ? C_R_Y_1(N, ξ, s₁, κ, ϵ, λ) : zero(T)
    C_R_Y_k_dξ = T == Arb ? C_R_Y_1_dξ(N, ξ, s₁, κ, ϵ, λ) : zero(T)

    return _Y_k_dλ_dξ(ξ, A₁, s₁, s₁_dλ, a₁, c2ns, c2ns_dλ, C_R_Y_k, C_R_Y_k_dξ)
end

function Y_2_dλ_dξ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    A₂ = A2(T)
    s₂ = s2(lambda, κ, λ)
    s₂_dλ = s2_dλ(lambda, κ, λ)
    a₂ = a2(κ, ϵ)

    N = 10 # TODO: Choose N
    c2ns = c2ns2(N - 1, s₂, κ, ϵ, λ)
    c2ns_dλ = c2ns2_dλ(N - 1, s₂, s₂_dλ, κ, ϵ, λ)

    C_R_Y_k = T == Arb ? C_R_Y_2(N, ξ, s₂, κ, ϵ, λ) : zero(T)
    C_R_Y_k_dξ = T == Arb ? C_R_Y_2_dξ(N, ξ, s₂, κ, ϵ, λ) : zero(T)

    return _Y_k_dλ_dξ(ξ, A₂, s₂, s₂_dλ, a₂, c2ns, c2ns_dλ, C_R_Y_k, C_R_Y_k_dξ)
end

function Y_3_dλ_dξ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    A₃ = A3(T)
    s₃ = s3(lambda, κ, λ)
    s₃_dλ = s3_dλ(lambda, κ, λ)
    a₃ = zero(s₃)

    N = 10 # TODO: Choose N
    c2ns = c2ns3(N - 1, s₃, κ, ϵ, λ)
    c2ns_dλ = c2ns3_dλ(N - 1, s₃, s₃_dλ, κ, ϵ, λ)

    C_R_Y_k = T == Arb ? C_R_Y_1(N, ξ, s₃, κ, ϵ, λ) : zero(T)
    C_R_Y_k_dξ = T == Arb ? C_R_Y_3_dξ(N, ξ, s₃, κ, ϵ, λ) : zero(T)

    return _Y_k_dλ_dξ(ξ, A₃, s₃, s₃_dλ, a₃, c2ns, c2ns_dλ, C_R_Y_k, C_R_Y_k_dξ)
end

function Y_4_dλ_dξ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    A₄ = A4(T)
    s₄ = s4(lambda, κ, λ)
    s₄_dλ = s4_dλ(lambda, κ, λ)
    a₄ = zero(s₄)

    N = 10 # TODO: Choose N
    c2ns = c2ns4(N - 1, s₄, κ, ϵ, λ)
    c2ns_dλ = c2ns4_dλ(N - 1, s₄, s₄_dλ, κ, ϵ, λ)

    C_R_Y_k = T == Arb ? C_R_Y_1(N, ξ, s₄, κ, ϵ, λ) : zero(T)
    C_R_Y_k_dξ = T == Arb ? C_R_Y_4_dξ(N, ξ, s₄, κ, ϵ, λ) : zero(T)

    return _Y_k_dλ_dξ(ξ, A₄, s₄, s₄_dλ, a₄, c2ns, c2ns_dλ, C_R_Y_k, C_R_Y_k_dξ)
end

function P_1(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, λ)
    z = -c * ξ^2
    return U(a - lambda / 2κ, b, z)
end

function P_1_dξ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, λ)
    z = -c * ξ^2
    z_dξ = -2c * ξ
    return U_dz(a - lambda / 2κ, b, z) * z_dξ
end

function P_2(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, λ)
    z = -conj(c) * ξ^2
    return U(conj(a) - lambda / 2κ, b, z)
end

function P_2_dξ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, λ)
    z = -conj(c) * ξ^2
    z_dξ = -2conj(c) * ξ
    return U_dz(conj(a) - lambda / 2κ, b, z) * z_dξ
end

function E_1(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, λ)
    z = -c * ξ^2
    return exp(z) * U(b - a + lambda / 2κ, b, -z)
end

function E_1_dξ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, λ)
    z = -c * ξ^2
    z_dξ = -2c * ξ
    return exp(z) *
           (U(b - a + lambda / 2κ, b, -z) - U_dz(b - a + lambda / 2κ, b, -z)) *
           z_dξ
end

function E_2(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, λ)
    z = -conj(c) * ξ^2
    return exp(z) * U(b - conj(a) + lambda / 2κ, b, -z)
end

function E_2_dξ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, λ)
    z = -conj(c) * ξ^2
    z_dξ = -2conj(c) * ξ
    return exp(z) *
           (U(b - conj(a) + lambda / 2κ, b, -z) - U_dz(b - conj(a) + lambda / 2κ, b, -z)) *
           z_dξ
end

function P_1_dλ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, λ)
    z = -c * ξ^2
    return -U_da(a - lambda / 2κ, b, z) / 2κ
end

function P_1_dλ_dξ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, λ)
    z = -c * ξ^2
    z_dξ = -2c * ξ
    return -U_dzda(a - lambda / 2κ, b, z) / 2κ * z_dξ
end

function P_2_dλ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, λ)
    z = -conj(c) * ξ^2
    return -U_da(conj(a) - lambda / 2κ, b, z) / 2κ
end

function P_2_dλ_dξ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, λ)
    z = -conj(c) * ξ^2
    z_dξ = -2conj(c) * ξ
    return -U_dzda(conj(a) - lambda / 2κ, b, z) / 2κ * z_dξ
end

function E_1_dλ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, λ)
    z = -c * ξ^2
    return exp(z) * U_da(b - a + lambda / 2κ, b, -z) / 2κ
end

function E_1_dλ_dξ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, λ)
    z = -c * ξ^2
    z_dξ = -2c * ξ
    return exp(z) *
           (U_da(b - a + lambda / 2κ, b, -z) - U_dzda(b - a + lambda / 2κ, b, -z)) / 2κ *
           z_dξ
end

function E_2_dλ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, λ)
    z = -conj(c) * ξ^2
    return exp(z) * U_da(b - conj(a) + lambda / 2κ, b, -z) / 2κ
end

function E_2_dλ_dξ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    a, b, c = _abc(κ, ϵ, λ)
    z = -conj(c) * ξ^2
    z_dξ = -2conj(c) * ξ
    return exp(z) * (
        U_da(b - conj(a) + lambda / 2κ, b, -z) - U_dzda(b - conj(a) + lambda / 2κ, b, -z)
    ) / 2κ * z_dξ
end

function W_1(ξ, lambda, κ, ϵ, λ::CGLParams)
    a, b, c = _abc(κ, ϵ, λ)

    z = -c * ξ^2

    sgn = if c isa AcbSeries
        sign(imag(c[0]))
    else
        sign(imag(c))
    end

    return -2c * exp(-sgn * im * (b - a + lambda / 2κ) * π) * ξ * z^-b * exp(z)
end

function W_2(ξ, lambda, κ, ϵ, λ::CGLParams)
    a, b, c = _abc(κ, ϵ, λ)

    z = -conj(c) * ξ^2

    sgn = if c isa AcbSeries
        sign(imag(c[0]))
    else
        sign(imag(c))
    end

    return -2conj(c) * exp(sgn * im * (b - conj(a) + lambda / 2κ) * π) * ξ * z^-b * exp(z)
end

function B_W_1(lambda, κ, ϵ, λ::CGLParams)
    a, b, c = _abc(κ, ϵ, λ)

    sgn = if c isa AcbSeries
        sign(imag(c[0]))
    else
        sign(imag(c))
    end

    return 1 // 2 * exp(sgn * im * (b - a + lambda / 2κ) * π) * (-c)^(b - 1)
end

function B_W_2(lambda, κ, ϵ, λ::CGLParams)
    a, b, c = _abc(κ, ϵ, λ)

    sgn = if c isa AcbSeries
        sign(imag(c[0]))
    else
        sign(imag(c))
    end

    return 1 // 2 * exp(-sgn * im * (b - conj(a) + lambda / 2κ) * π) * (-conj(c))^(b - 1)
end

function B_W_1_dλ(lambda, κ, ϵ, λ::CGLParams)
    a, b, c = _abc(κ, ϵ, λ)

    sgn = if c isa AcbSeries
        sign(imag(c[0]))
    else
        sign(imag(c))
    end

    return 1 // 2 *
           (sgn * im / 2κ * π) *
           exp(sgn * im * (b - a + lambda / 2κ) * π) *
           (-c)^(b - 1)
end

function B_W_2_dλ(lambda, κ, ϵ, λ::CGLParams)
    a, b, c = _abc(κ, ϵ, λ)

    sgn = if c isa AcbSeries
        sign(imag(c[0]))
    else
        sign(imag(c))
    end

    return 1 // 2 *
           (-sgn * im / 2κ * π) *
           exp(-sgn * im * (b - conj(a) + lambda / 2κ) * π) *
           (-conj(c))^(b - 1)
end

function J_P_1(ξ, lambda, κ, ϵ, λ::CGLParams)
    (; d) = λ
    _, _, c = _abc(κ, ϵ, λ)
    return B_W_1(lambda, κ, ϵ, λ) * P_1(ξ, lambda, κ, ϵ, λ) * exp(c * ξ^2) * ξ^(d - 1)
end

function J_P_2(ξ, lambda, κ, ϵ, λ::CGLParams)
    (; d) = λ
    _, _, c = _abc(κ, ϵ, λ)
    return B_W_2(lambda, κ, ϵ, λ) * P_2(ξ, lambda, κ, ϵ, λ) * exp(conj(c) * ξ^2) * ξ^(d - 1)
end

function J_E_1(ξ, lambda, κ, ϵ, λ::CGLParams)
    (; d) = λ
    _, _, c = _abc(κ, ϵ, λ)
    return B_W_1(lambda, κ, ϵ, λ) * E_1(ξ, lambda, κ, ϵ, λ) * exp(c * ξ^2) * ξ^(d - 1)
end

function J_E_2(ξ, lambda, κ, ϵ, λ::CGLParams)
    (; d) = λ
    _, _, c = _abc(κ, ϵ, λ)
    return B_W_2(lambda, κ, ϵ, λ) * E_2(ξ, lambda, κ, ϵ, λ) * exp(conj(c) * ξ^2) * ξ^(d - 1)
end

function J_P_1_dλ(ξ, lambda, κ, ϵ, λ::CGLParams)
    (; d) = λ
    _, _, c = _abc(κ, ϵ, λ)
    return B_W_1_dλ(lambda, κ, ϵ, λ) * P_1(ξ, lambda, κ, ϵ, λ) * exp(c * ξ^2) * ξ^(d - 1) +
           B_W_1(lambda, κ, ϵ, λ) * P_1_dλ(ξ, lambda, κ, ϵ, λ) * exp(c * ξ^2) * ξ^(d - 1)
end

function J_P_2_dλ(ξ, lambda, κ, ϵ, λ::CGLParams)
    (; d) = λ
    _, _, c = _abc(κ, ϵ, λ)
    return B_W_2_dλ(lambda, κ, ϵ, λ) *
           P_2(ξ, lambda, κ, ϵ, λ) *
           exp(conj(c) * ξ^2) *
           ξ^(d - 1) +
           B_W_2(lambda, κ, ϵ, λ) *
           P_2_dλ(ξ, lambda, κ, ϵ, λ) *
           exp(conj(c) * ξ^2) *
           ξ^(d - 1)
end

function J_E_1_dλ(ξ, lambda, κ, ϵ, λ::CGLParams)
    (; d) = λ
    _, _, c = _abc(κ, ϵ, λ)
    return B_W_1_dλ(lambda, κ, ϵ, λ) * E_1(ξ, lambda, κ, ϵ, λ) * exp(c * ξ^2) * ξ^(d - 1) +
           B_W_1(lambda, κ, ϵ, λ) * E_1_dλ(ξ, lambda, κ, ϵ, λ) * exp(c * ξ^2) * ξ^(d - 1)
end

function J_E_2_dλ(ξ, lambda, κ, ϵ, λ::CGLParams)
    (; d) = λ
    _, _, c = _abc(κ, ϵ, λ)
    return B_W_2_dλ(lambda, κ, ϵ, λ) *
           E_2(ξ, lambda, κ, ϵ, λ) *
           exp(conj(c) * ξ^2) *
           ξ^(d - 1) +
           B_W_2(lambda, κ, ϵ, λ) *
           E_2_dλ(ξ, lambda, κ, ϵ, λ) *
           exp(conj(c) * ξ^2) *
           ξ^(d - 1)
end
