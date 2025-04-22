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

    @assert ComplexF64(ν) ≈ νF64
    @assert ComplexF64(γ₂) ≈ γ₂F64

    verbose && @info "Got" ν γ₂

    ###
    # Step 3.1: Find approximate eigenvalue from finite difference method
    ###
    verbose && @info "Computing approximate eigenvalue"
    λsF64 = filter(
        lambda -> imag(lambda) > 0,
        CGL2.linearization_eigenvalues_real_1(νF64, κF64, ϵF64, ξ₁F64, λF64)[1],
    )

    lambdaF64_approx = λsF64[findmax(real, λsF64)[2]]

    verbose && @info "Got" lambdaF64_approx

    ###
    # Step 3: Solve for x and λ
    ###
    verbose && @info "Solving for x and λ"

    Q_hat_F64 = CGL2.Q_hat_zero_float_curve(real(νF64), imag(νF64), κF64, ϵF64, ξ₁F64, λF64)
    Q_hat_ξ₁_F64 = complex(Q_hat_F64(ξ₁F64)[1:2]...)

    xF64, cF64, lambdaF64 = CGL2.H_approximate(
        lambdaF64_approx,
        κF64,
        ϵF64,
        ξ₁F64,
        Q_hat_F64,
        Q_hat_ξ₁_F64,
        λF64,
    )

    # TODO: Implement rigorous version
    x, c, lambda = Acb(xF64), Acb.(cF64), Acb(lambdaF64)

    @assert ComplexF64(x) ≈ xF64
    @assert ComplexF64.(c) ≈ cF64
    @assert ComplexF64(lambda) ≈ lambdaF64

    return x, c, lambda
end
