function G_hat_approximate(
    γ₁::Complex{T},
    κ::T,
    ϵ::T,
    ξ₁::T,
    Λ::CGLParams{T};
    return_convergence::Union{Val{false},Val{true}} = Val{false}(),
    verbose = false,
) where {T}
    F((ν, γ₂), (γ₁, κ, ϵ, ξ₁, Λ)) = G_hat(ν, γ₁, γ₂, κ, ϵ, ξ₁, Λ)
    x₀ = SVector(zero(γ₁), zero(γ₁)) # IMPROVE: Pick this in a smarter way
    prob = NonlinearProblem{false}(F, x₀, (γ₁, κ, ϵ, ξ₁, Λ))
    sol = try
        solve(
            prob,
            NewtonRaphson(autodiff = AutoFiniteDiff()),
            maxiters = 100, # Should be enough to saturate converge
            reltol = 1e-10,
            abstol = 1e-10,
        )
    catch e
        e isa SingularException || rethrow(e)
        verbose && @warn "Encountered" e ϵ

        if return_convergence isa Val{false}
            return zero(γ₁), zero(γ₁)
        else
            return false, zero(γ₁), zero(γ₁)
        end
    end

    if verbose && !(norm(sol.resid) < 1e-1)
        @warn "Very low precision when refining approximation" ϵ sol.resid
    end

    converged = NonlinearSolve.SciMLBase.successful_retcode(sol)

    ν, γ₂ = if converged
        sol[1], sol[2]
    else
        zero(γ₁), zero(γ₁)
    end
    if return_convergence isa Val{false}
        return ν, γ₂
    else
        return converged, ν, γ₂
    end
end

function G_hat_solve(
    ν::Acf,
    γ₁::Acb,
    γ₂::Acf,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    Λ::CGLParams{Arb};
    expansion_rate = 0.05,
    max_iterations = 10,
    verbose = false,
)
    G_hat_x = ((ν, γ₂),) -> G_hat(ν, γ₁, γ₂, κ, ϵ, ξ₁, Λ)
    dG_hat_x = ((ν, γ₂),) -> G_hat_jacobian(ν, γ₁, γ₂, κ, ϵ, ξ₁, Λ)

    root, _ = verify_root_from_approximation(
        G_hat_x,
        dG_hat_x,
        SVector{2,Acb}(ν, γ₂);
        expansion_rate,
        max_iterations,
        verbose,
    )

    # Manually perform a couple more Newton iterations. It might be
    # better to have verify_root_from_approximation do this by itself,
    # but this is a simple solution.
    for _ = 1:2
        root = newton_step(G_hat_x, dG_hat_x, root)
    end

    return root
end
