function H_approximate(
    μ::T,
    γ::Complex{T},
    κ::T,
    ϵ::T,
    ξ₁::T,
    λ::CGLParams{T};
    return_convergence::Union{Val{false},Val{true}} = Val{false}(),
    verbose = false,
) where {T}
    F(ν, (μ, γ, κ, ϵ, ξ₁, λ)) = H(ν[1], ν[2], μ, γ, κ, ϵ, ξ₁, λ)

    prob = NonlinearProblem{false}(F, SVector(μ, zero(μ)), (μ, γ, κ, ϵ, ξ₁, λ))
    sol = try
        solve(
            prob,
            NewtonRaphson(autodiff = AutoFiniteDiff()),
            maxiters = 50, # Should be enough to saturate converge
            reltol = 1e-10,
            abstol = 1e-10,
        )
    catch e
        e isa SingularException || rethrow(e)
        verbose && @warn "Encountered" e ϵ

        if return_convergence isa Val{false}
            return _complex(μ, zero(μ))
        else
            return false, _complex(μ, zero(μ))
        end
    end

    converged = norm(sol.resid) < 1e-1

    if verbose && !converged
        @warn "Very low precision when refining approximation" ϵ sol.resid
    end

    if return_convergence isa Val{false}
        return _complex(sol.u...)
    else
        return converged, _complex(sol.u...)
    end
end
