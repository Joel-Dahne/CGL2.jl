function H_solve(
    lambda::Acf,
    ν::Acb,
    γ₁::Acb,
    γ₂::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{T};
    expansion_rate = 0.05,
    max_iterations = 10,
    verbose = false,
) where {T}
    root, _ = verify_root_from_approximation(
        lambda -> H(lambda, ν, γ₁, γ₂, κ, ϵ, ξ₁, λ),
        Acb(lambda);
        expansion_rate,
        max_iterations,
        verbose,
    )

    return root
end
