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
    # Coefficient for leading term in asymptotic expansion
    c10_hat = let
        a, b, c = CGL2._abc(κ, ϵ, λ)
        γ * c^-a / (-c)^-a
    end

    F(x, (c10_hat, κ, ϵ, ξ₁, λ)) = H(x..., c10_hat, κ, ϵ, ξ₁, λ)
    x₀ = SVector(0.0, 0.0, 0.0, 0.0) # IMPROVE: Pick this in a smarter way
    prob = NonlinearProblem{false}(F, x₀, (c10_hat, κ, ϵ, ξ₁, λ))
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
            return _complex(zero(μ), zero(μ)), _complex(zero(μ), zero(μ))
        else
            return false, _complex(μ, zero(μ)), _complex(zero(μ), zero(μ))
        end
    end

    if verbose && !(norm(sol.resid) < 1e-1)
        @warn "Very low precision when refining approximation" ϵ sol.resid
    end

    converged = NonlinearSolve.SciMLBase.successful_retcode(sol)

    ν, c20_hat = if converged
        _complex(sol.u[1:2]...), _complex(sol.u[3:4]...)
    else
        _complex(zero(μ), zero(μ)), _complex(zero(μ), zero(μ))
    end
    if return_convergence isa Val{false}
        return ν, c20_hat
    else
        return converged, ν, c20_hat
    end
end
