function eigenvalue_solve(
    μ::Arb,
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    max_iterations::Int = 10,
    verbose = false,
)
    verbose && @info "Iteration ϵ = $ϵ"

    κF64 = Float64(κ)
    ϵF64 = Float64(ϵ)
    ξ₁F64 = Float64(ξ₁)
    λF64 = CGLParams{Float64}(λ)

    ###
    # Step 1: Compute p_Q_0
    ###
    pQ0 = p_Q_0(γ, κ, ϵ, ξ₁, λ)
    #return pQ0
    ###
    # Step 2.1: Solve for γ₁, giving asymptotic behavior of Q_hat at infinity
    ###
    verbose && @info "Solving for γ₁"

    a, b, c = _abc(κ, ϵ, λ)
    γ₁ = pQ0 / (-c)^-a

    @assert Arblib.overlaps(p_Q_hat(γ₁, κ, ϵ, λ), pQ0)
    #return γ₁
    verbose && @info "Got" γ₁

    ###
    # Step 2.2: Solve for ν and γ₂
    ###
    verbose && @info "Solving for ν and γ₂"

    νF64_approx, γ₂F64_approx = G_hat_approximate(ComplexF64(γ₁), κF64, ϵF64, ξ₁F64, λF64)

    ν, γ₂ = G_hat_solve(Acf(νF64_approx), γ₁, Acf(γ₂F64_approx), κ, ϵ, ξ₁, λ; verbose)

    @show ComplexF64(ν) ≈ νF64_approx
    @show ComplexF64(γ₂) ≈ γ₂F64_approx

    verbose && @info "Got" ν γ₂

    # Return ν and γ₂
    #return ν, γ₂

    ###
    # Step 3.1: Find approximate eigenvalue
    ###
    verbose && @info "Computing approximate eigenvalue"

    ###
    # Step 3.1.1: Find approximation using finite difference method
    ###
    verbose && @info "Computing approximation using finite difference"

    λsF64 = filter(
        lambda -> imag(lambda) > 0,
        linearization_eigenvalues_real_1(νF64_approx, κF64, ϵF64, ξ₁F64, λF64)[1],
    )

    lambdaF64_approx = λsF64[findmax(real, λsF64)[2]]

    verbose && @info "Got" lambdaF64_approx

    ###
    # Step 3.1.2: Refine approximation using Newton
    ###
    verbose && @info "Refining approximation using Newton"

    # IMPROVE: Consider running more iterations?
    H_approx = H(
        AcbSeries((lambdaF64_approx, 1)),
        midpoint(Acb, ν),
        midpoint(Acb, γ₁),
        midpoint(Acb, γ₂),
        midpoint(Arb, κ),
        ϵ,
        ξ₁,
        λ,
    )
    lambda_approx =
        lambdaF64_approx - midpoint(Acf, H_approx[0]) / midpoint(Acf, H_approx[1])

    verbose && @info "Got" lambda_approx

    ###
    # Step 3.2: Solve λ
    ###
    verbose && @info "Solving for λ"

    #verbose && @info "Solving for λ using midpoint of ν"
    #lambda_mid = H_solve(lambda_approx, midpoint(Acb, ν), γ₁, γ₂, κ, ϵ, ξ₁, λ; verbose)

    # Return parameters
    #return lambda_approx, ν, γ₁, γ₂

    # FIXME: Improve enclosures so that we don't need this scaling
    ν_radius_scaling = Mag(8e-3)
    verbose && @info "Solving for λ using ν with radius scaled by" ν_radius_scaling
    ν = Acb(
        setball(Arb, midpoint(real(ν)), ν_radius_scaling * Arblib.radius(real(ν))),
        setball(Arb, midpoint(imag(ν)), ν_radius_scaling * Arblib.radius(imag(ν))),
    )
    lambda = H_solve(lambda_approx, ν, γ₁, γ₂, κ, ϵ, ξ₁, λ; verbose, max_iterations)

    # Note that the finite difference approximation is quite bad
    @show isapprox(ComplexF64(lambda), lambdaF64_approx, rtol = 1e-4)

    return lambda
end
