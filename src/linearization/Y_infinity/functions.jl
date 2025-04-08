function Y_1(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    a₂ = a21(κ, ϵ)
    s = s11(lambda, κ, λ)

    # TODO: Choose N and add remainder

    N = 10

    ns = 0:2:N
    cs_odd = cs_odd_11(N, lambda, κ, ϵ, λ)

    res = sum(zip(ns, cs_odd)) do (n, cn)
        cn * ξ^(-n - s)
    end

    return exp(a₂ * ξ^2) * res * A01(T)
end

function Y_1_dξ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    a₂ = a21(κ, ϵ)
    s = s11(lambda, κ, λ)

    # TODO: Choose N and add remainder

    N = 10

    ns = 0:2:N
    cs_odd = cs_odd_11(N, lambda, κ, ϵ, λ)

    res1 = sum(zip(ns, cs_odd)) do (n, cn)
        cn * ξ^(-n - s)
    end

    res2 = -sum(zip(ns, cs_odd)) do (n, cn)
        (n + s) * cn * ξ^(-n - 1 - s)
    end

    return exp(a₂ * ξ^2) * (2a₂ * ξ * res1 + res2) * A01(T)
end

function Y_1_dξ_dξ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    a₂ = a21(κ, ϵ)
    s = s11(lambda, κ, λ)

    # TODO: Choose N and add remainder

    N = 10

    ns = 0:2:N
    cs_odd = cs_odd_11(N, lambda, κ, ϵ, λ)

    res1 = sum(zip(ns, cs_odd)) do (n, cn)
        cn * ξ^(-n - s)
    end

    res2 = -sum(zip(ns, cs_odd)) do (n, cn)
        (n + s) * cn * ξ^(-n - 1 - s)
    end

    res3 = sum(zip(ns, cs_odd)) do (n, cn)
        (n + s) * (n + 1 + s) * cn * ξ^(-n - 2 - s)
    end

    return exp(a₂ * ξ^2) * (((2a₂ * ξ)^2 + 2a₂) * res1 + 4a₂ * ξ * res2 + res3) * A01(T)
end

function Y_2(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    a₂ = a22(κ, ϵ)
    s = s12(lambda, κ, λ)

    # TODO: Choose N and add remainder
    N = 10

    ns = 0:2:N
    cs_odd = cs_odd_12(N, lambda, κ, ϵ, λ)

    res = sum(zip(ns, cs_odd)) do (n, cn)
        cn * ξ^(-n - s)
    end

    return exp(a₂ * ξ^2) * res * A02(T)
end

function Y_2_dξ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    a₂ = a22(κ, ϵ)
    s = s12(lambda, κ, λ)

    # TODO: Choose N and add remainder

    N = 10

    ns = 0:2:N
    cs_odd = cs_odd_12(N, lambda, κ, ϵ, λ)

    res1 = sum(zip(ns, cs_odd)) do (n, cn)
        cn * ξ^(-n - s)
    end

    res2 = -sum(zip(ns, cs_odd)) do (n, cn)
        (n + s) * cn * ξ^(-n - 1 - s)
    end

    return exp(a₂ * ξ^2) * (2a₂ * ξ * res1 + res2) * A02(T)
end

function Y_2_dξ_dξ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    a₂ = a22(κ, ϵ)
    s = s12(lambda, κ, λ)

    # TODO: Choose N and add remainder

    N = 10

    ns = 0:2:N
    cs_odd = cs_odd_12(N, lambda, κ, ϵ, λ)

    res1 = sum(zip(ns, cs_odd)) do (n, cn)
        cn * ξ^(-n - s)
    end

    res2 = -sum(zip(ns, cs_odd)) do (n, cn)
        (n + s) * cn * ξ^(-n - 1 - s)
    end

    res3 = sum(zip(ns, cs_odd)) do (n, cn)
        (n + s) * (n + 1 + s) * cn * ξ^(-n - 2 - s)
    end

    return exp(a₂ * ξ^2) * (((2a₂ * ξ)^2 + 2a₂) * res1 + 4a₂ * ξ * res2 + res3) * A02(T)
end


function Y_3(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    s = s21(lambda, κ, λ)

    # TODO: Choose N and add remainder

    N = 10

    ns = 0:2:N
    cs_odd = cs_odd_21(N, lambda, κ, ϵ, λ)

    res = sum(zip(ns, cs_odd)) do (n, cn)
        cn * ξ^(-n - s)
    end

    return res * A01(T)
end

function Y_3_dξ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    s = s21(lambda, κ, λ)

    # TODO: Choose N and add remainder

    N = 10

    ns = 0:2:N
    cs_odd = cs_odd_21(N, lambda, κ, ϵ, λ)

    res = -sum(zip(ns, cs_odd)) do (n, cn)
        (n + s) * cn * ξ^(-n - 1 - s)
    end

    return res * A01(T)
end

function Y_3_dξ_dξ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    s = s21(lambda, κ, λ)

    # TODO: Choose N and add remainder

    N = 10

    ns = 0:2:N
    cs_odd = cs_odd_21(N, lambda, κ, ϵ, λ)

    res = sum(zip(ns, cs_odd)) do (n, cn)
        (n + s) * (n + 1 + s) * cn * ξ^(-n - 2 - s)
    end

    return res * A01(T)
end

function Y_4(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    s = s22(lambda, κ, λ)

    # TODO: Choose N and add remainder

    N = 10

    ns = 0:2:N
    cs_odd = cs_odd_22(N, lambda, κ, ϵ, λ)

    res = sum(zip(ns, cs_odd)) do (n, cn)
        cn * ξ^(-n - s)
    end

    return res * A02(T)
end

function Y_4_dξ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    s = s22(lambda, κ, λ)

    # TODO: Choose N and add remainder

    N = 10

    ns = 0:2:N
    cs_odd = cs_odd_22(N, lambda, κ, ϵ, λ)

    res = -sum(zip(ns, cs_odd)) do (n, cn)
        (n + s) * cn * ξ^(-n - 1 - s)
    end

    return res * A02(T)
end

function Y_4_dξ_dξ(ξ, lambda, κ, ϵ, λ::CGLParams{T}) where {T}
    s = s22(lambda, κ, λ)

    # TODO: Choose N and add remainder

    N = 10

    ns = 0:2:N
    cs_odd = cs_odd_22(N, lambda, κ, ϵ, λ)

    res = sum(zip(ns, cs_odd)) do (n, cn)
        (n + s) * (n + s + 1) * cn * ξ^(-n - 2 - s)
    end

    return res * A02(T)
end
