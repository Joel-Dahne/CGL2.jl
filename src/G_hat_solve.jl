function G_hat_approximate(
    γ₁::Complex{T},
    κ::T,
    ϵ::T,
    ξ₁::T,
    λ::CGLParams{T};
    return_convergence::Union{Val{false},Val{true}} = Val{false}(),
    verbose = false,
) where {T}
    F((ν, γ₂), (γ₁, κ, ϵ, ξ₁, λ)) = G_hat(ν, γ₁, γ₂, κ, ϵ, ξ₁, λ)
    x₀ = SVector(zero(γ₁), zero(γ₁)) # IMPROVE: Pick this in a smarter way
    prob = NonlinearProblem{false}(F, x₀, (γ₁, κ, ϵ, ξ₁, λ))
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
