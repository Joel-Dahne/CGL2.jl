function H_solve(
    lambda::Acf,
    ν::Acb,
    γ₁::Acb,
    γ₂::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{T};
    return_uniqueness::Union{Val{false},Val{true}} = Val{false}(),
    try_expand_uniqueness = return_uniqueness isa Val{true},
    expansion_rate = 0.05,
    max_iterations = 10,
    verbose = false,
    extra_verbose = false,
) where {T}
    H_x = lambda -> H(lambda, ν, γ₁, γ₂, κ, ϵ, ξ₁, λ)

    root, root_uniqueness = verify_root_from_approximation(
        H_x,
        Acb(lambda);
        expansion_rate,
        max_iterations,
        verbose,
    )

    if try_expand_uniqueness
        verbose && @info "Expanding region for uniqueness"
        # TODO: Extend this to Acb
        root_uniqueness = expand_uniqueness(H_x, root_uniqueness; verbose, extra_verbose)
    end

    if return_uniqueness isa Val{true}
        return root, root_uniqueness
    else
        return root
    end
end
