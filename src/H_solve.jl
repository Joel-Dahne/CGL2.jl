function H_approximate(
    lambda_approx::Complex{T},
    ν::Complex{T},
    κ::T,
    ϵ::T,
    ξ₁::T,
    λ::CGLParams{T};
    return_convergence::Union{Val{false},Val{true}} = Val{false}(),
    verbose = false,
) where {T}
    Q_hat = CGL2.Q_hat_zero_float_curve(real(ν), imag(ν), κ, ϵ, ξ₁, λ)
    Q_hat_ξ₁ = complex(Q_hat(ξ₁)[1:2]...)

    F((x, c1, c2, lambda), (κ, ϵ, ξ₁, Q_hat, Q_hat_ξ₁, λ)) =
        H(x, SVector(c1, c2), lambda, κ, ϵ, ξ₁, Q_hat, Q_hat_ξ₁, λ)
    x₀ = SVector{4,Complex{T}}(0, 0, 0, lambda_approx) # IMPROVE: Pick this in a smarter way
    prob = NonlinearProblem{false}(F, x₀, (κ, ϵ, ξ₁, Q_hat, Q_hat_ξ₁, λ))
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
            return _complex(zero(μ), zero(μ)), _complex(zero(μ), zero(μ))
        else
            return false, _complex(μ, zero(μ)), _complex(zero(μ), zero(μ))
        end
    end

    if verbose && !(norm(sol.resid) < 1e-1)
        @warn "Very low precision when refining approximation" ϵ sol.resid
    end

    converged = NonlinearSolve.SciMLBase.successful_retcode(sol)

    x, c, lambda = if converged
        sol[1], SVector(sol[2], sol[3]), sol[4]
    else
        zero(Complex{T}), SVector(zero(Complex{T}), zero(Complex{T})), zero(Complex{T})
    end
    if return_convergence isa Val{false}
        return x, c, lambda
    else
        return converged, x, c, lambda
    end
end

function H_solve(
    x::Acf,
    c::SVector{2,Acf},
    lambda::Acf,
    ν::Acb,
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
    # TODO: Implement these methods
    H_x = ((x, c1, c2, lambda),) -> H(x, SVector(c1, c2), lambda, ν, κ, ϵ, ξ₁, λ)
    dH_x = ((x, c1, c2, lambda),) -> H_jacobian(x, SVector(c1, c2), lambda, ν, κ, ϵ, ξ₁, λ)

    root, root_uniqueness = verify_root_from_approximation(
        H_x,
        dH_x,
        SVector{4,Acb}(x, c[1], c[2], lambda);
        expansion_rate,
        max_iterations,
        verbose,
    )

    if try_expand_uniqueness
        verbose && @info "Expanding region for uniqueness"
        # TODO: Extend this to Acb
        root_uniqueness =
            expand_uniqueness(H_x, dH_x, root_uniqueness; verbose, extra_verbose)
    end

    if return_uniqueness isa Val{true}
        return root, root_uniqueness
    else
        return root
    end
end
