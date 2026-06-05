function branch_eigenvalue_points(
    μs::Vector{Arb},
    γs::Vector{Acb},
    κs::Vector{Arb},
    ϵs::Vector{Arb},
    ξ₁s::Vector{Arb},
    Λ::CGLParams{Arb};
    verbose = false,
)
    res = map(μs, γs, κs, ϵs, ξ₁s) do μ, γ, κ, ϵ, ξ₁
        eigenvalue_solve(μ, γ, κ, ϵ, ξ₁, Λ; verbose)
    end

    xs = getindex.(res, 1)
    cs = getindex.(res, 2)
    λs = getindex.(res, 3)

    return xs, cs, λs
end
