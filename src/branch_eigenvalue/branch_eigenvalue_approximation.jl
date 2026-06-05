function branch_eigenvalue_approximation(
    μs::Vector{T},
    γs::Vector{Complex{T}},
    κs::Vector{T},
    ϵs::Vector{T},
    ξ₁s::Vector{T},
    Λ::CGLParams{T};
    verbose = false,
) where {T}
    res = map(μs, γs, κs, ϵs, ξ₁s) do μ, γ, κ, ϵ, ξ₁
        verbose && @info "Iteration ϵ = $ϵ"

        ###
        # Step 1: Compute p_Q_0
        ###
        pQ0 = p_Q_0(γ, κ, ϵ, ξ₁, Λ)

        ###
        # Step 2.1: Solve for γ₁, giving asymptotic behavior of Q_hat at infinity
        ###
        a, b, c = _abc(κ, ϵ, Λ)
        γ₁ = pQ0 / (-c)^-a

        @assert p_Q_hat(γ₁, κ, ϵ, Λ) ≈ pQ0

        ###
        # Step 2.2: Solve for ν and γ₂
        ###
        ν, γ₂ = G_hat_approximate(γ₁, κ, ϵ, ξ₁, Λ)

        @assert isapprox(norm(G_hat(ν, γ₁, γ₂, κ, ϵ, ξ₁, Λ)), 0, atol = 1e-9)

        ###
        # Step 3.1: Find approximate eigenvalue from finite difference method
        ###
        Λs = filter(λ -> imag(λ) > 0, linearization_eigenvalues_real_1(ν, κ, ϵ, ξ₁, Λ)[1])

        λ_approx = Λs[findmax(real, Λs)[2]]

        ###
        # Step 3: Solve for Z and Λ
        ###
        x, c, λ = H_approximate(λ_approx, ν, κ, ϵ, ξ₁, Λ)

        @assert isapprox(norm(H(x, c, λ, ν, γ₁, γ₂, κ, ϵ, ξ₁, Λ)), 0, atol = 1e-9)
        @assert isapprox(λ, λ_approx, rtol = 1e-3)

        (x, c, λ)
    end

    xs = getindex.(res, 1)
    cs = getindex.(res, 2)
    λs = getindex.(res, 3)

    return xs, cs, λs
end
