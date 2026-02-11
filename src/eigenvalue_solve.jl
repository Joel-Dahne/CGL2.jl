function eigenvalue_solve(
    μ::Arb,
    γ::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    verbose = false,
)
    verbose && @info "Iteration ϵ = $ϵ"

    κF64 = Float64(κ)
    ϵF64 = Float64(ϵ)
    ξ₁F64 = Float64(ξ₁)
    λF64 = CGLParams{Float64}(λ)

    ###
    # Step 1: Compute p_Q
    ###
    pQ = CGL2.p_Q(γ, κ, ϵ, ξ₁, λ)

    ###
    # Step 2.1: Solve for γ₁, giving asymptotic behavior of Q_hat at infinity
    ###
    verbose && @info "Solving for γ₁"

    a, b, c = CGL2._abc(κ, ϵ, λ)
    γ₁ = pQ / (-c)^-a

    @assert Arblib.overlaps(CGL2.p_Q_hat(γ₁, κ, ϵ, λ), pQ)

    verbose && @info "Got" γ₁

    ###
    # Step 2.2: Solve for ν and γ₂
    ###
    verbose && @info "Solving for ν and γ₂"

    νF64, γ₂F64 = CGL2.G_hat_approximate(ComplexF64(γ₁), κF64, ϵF64, ξ₁F64, λF64)

    # TODO: Implement rigorous version
    ν, γ₂ = CGL2.G_hat_solve(Acf(νF64), γ₁, Acf(γ₂F64), κ, ϵ, ξ₁, λ; verbose)

    @show ComplexF64(ν) ≈ νF64
    @show ComplexF64(γ₂) ≈ γ₂F64

    verbose && @info "Got" ν γ₂

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
        CGL2.linearization_eigenvalues_real_1(νF64, κF64, ϵF64, ξ₁F64, λF64)[1],
    )

    lambdaF64_approx = λsF64[findmax(real, λsF64)[2]]

    verbose && @info "Got" lambdaF64_approx

    ###
    # Step 3.1.2: Refine approximation using Newton
    ###
    verbose && @info "Refining approximation using Newton"

    # IMPROVE: Consider running more iterations?
    H_approx = H(AcbSeries((lambdaF64_approx, 1)), ν, γ₁, γ₂, κ, ϵ, ξ₁, λ)
    lambda_approx =
        lambdaF64_approx - midpoint(Acf, H_approx[0]) / midpoint(Acf, H_approx[1])

    verbose && @info "Got" lambda_approx

    ###
    # Step 3.2: Solve λ
    ###
    verbose && @info "Solving for λ"

    verbose && @info "Solving for λ using midpoint of ν"
    lambda_mid = CGL2.H_solve(lambda_approx, midpoint(Acb, ν), γ₁, γ₂, κ, ϵ, ξ₁, λ; verbose)

    # FIXME: Improve enclosures so that we don't need this scaling
    ν_radius_scaling = Mag(1e-3)
    verbose && @info "Solving for λ using ν with radius scaled by" ν_radius_scaling
    Arblib.mul!(Arblib.radref(Arblib.realref(ν)), radius(real(ν)), ν_radius_scaling)
    Arblib.mul!(Arblib.radref(Arblib.imagref(ν)), radius(imag(ν)), ν_radius_scaling)
    lambda = CGL2.H_solve(lambda_approx, ν, γ₁, γ₂, κ, ϵ, ξ₁, λ; verbose)

    # Note that the finite difference approximation is quite bad
    @show isapprox(ComplexF64(lambda), lambdaF64_approx, rtol = 1e-4)

    return lambda
end
