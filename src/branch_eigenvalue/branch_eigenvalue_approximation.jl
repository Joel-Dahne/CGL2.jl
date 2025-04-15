function branch_eigenvalue_approximation(
    μs::Vector{T},
    γs::Vector{Complex{T}},
    κs::Vector{T},
    ϵs::Vector{T},
    ξ₁s::Vector{T},
    λ::CGLParams{T};
    verbose = false,
) where {T}
    res = map(μs, γs, κs, ϵs, ξ₁s) do μ, γ, κ, ϵ, ξ₁
        verbose && @info "Iteration ϵ = $ϵ"

        ###
        # Step 1: Compute p_Q
        ###
        pQ = CGL2.p_Q(γ, κ, ϵ, ξ₁, λ)

        ###
        # Step 2.1: Solve for γ₁, giving asymptotic behavior of Q_hat at infinity
        ###
        a, b, c = CGL2._abc(κ, ϵ, λ)
        γ₁ = pQ / (-c)^-a

        @assert CGL2.p_Q_hat(γ₁, κ, ϵ, λ) ≈ pQ

        ###
        # Step 2.2: Solve for ν and γ₂
        ###
        ν, γ₂ = CGL2.G_hat_approximate(γ₁, κ, ϵ, ξ₁, λ)

        @assert isapprox(norm(CGL2.G_hat(ν, γ₁, γ₂, κ, ϵ, ξ₁, λ)), 0, atol = 1e-9)

        ###
        # Step 3.1: Find approximate eigenvalue from finite difference method
        ###
        λs = filter(
            lambda -> imag(lambda) > 0,
            CGL2.linearization_eigenvalues_real_1(ν, κ, ϵ, ξ₁, λ)[1],
        )

        lambda_approx = λs[findmax(real, λs)[2]]

        ###
        # Step 3: Solve for Z and λ
        ###
        Q_hat = CGL2.Q_hat_zero_float_curve(real(ν), imag(ν), κ, ϵ, ξ₁, λ)
        Q_hat_ξ₁ = CGL2._complex(Q_hat(ξ₁)[1:2]...)

        @assert Q_hat_ξ₁ ≈ CGL2.Q_hat_infinity(γ₁, γ₂, κ, ϵ, ξ₁, λ)[1]

        x, c, lambda = CGL2.H_approximate(lambda_approx, κ, ϵ, ξ₁, Q_hat, Q_hat_ξ₁, λ)

        @assert isapprox(
            norm(CGL2.H(x, c, lambda, κ, ϵ, ξ₁, Q_hat, Q_hat_ξ₁, λ)),
            0,
            atol = 1e-9,
        )
        @assert isapprox(lambda, lambda_approx, rtol = 1e-3)

        (x, c, lambda)
    end

    xs = getindex.(res, 1)
    cs = getindex.(res, 2)
    lambdas = getindex.(res, 3)

    return xs, cs, lambdas
end
