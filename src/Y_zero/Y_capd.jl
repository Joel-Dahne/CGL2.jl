"""
    _Y_zero_capd(
        Y_ξ₀::SVector{4,Complex{Interval{Float64}}},
        lambda::Complex{Interval{Float64}},
        ν::Complex{Interval{Float64}},
        κ::Interval{Float64},
        ϵ::Interval{Float64},
        ξ₀::Interval{Float64},
        ξ₁::Interval{Float64},
        λ::CGLParams{Interval{Float64}};
        output_jacobian::Union{Val{false},Val{true}} = Val{false}(),
        tol::Float64 = 1e-11,
    )

Internal function for calling the CAPD program.

It computes the solution to the ODE on the interval ``[ξ₀, ξ₁]`` and
returns four complex values, the first two are the values at `ξ₁` and
the second two are their derivatives.

If `output_jacobian = Val{true}()` it also computes the Jacobian
w.r.t. `Y_ξ₀` and `lambda`.
"""
function _Y_zero_capd(
    Q_hat_ξ₀::SVector{4,Interval{Float64}},
    Y_ξ₀::SVector{4,Complex{Interval{Float64}}},
    lambda::Complex{Interval{Float64}},
    κ::Interval{Float64},
    ϵ::Interval{Float64},
    ξ₀::Interval{Float64},
    ξ₁::Interval{Float64},
    λ::CGLParams{Interval{Float64}};
    output_jacobian::Union{Val{false},Val{true}} = Val{false}(),
    tol::Float64 = 1e-11,
)
    # Build the input to the CAPD program
    input_Q_hat_ξ₀ = ""
    for x in Q_hat_ξ₀
        input_Q_hat_ξ₀ *= "[$(inf(x)), $(sup(x))]\n"
    end
    input_Y_ξ₀ = ""
    input_Y_ξ₀ *= "[$(inf(real(Y_ξ₀[1]))), $(sup(real(Y_ξ₀[1])))]\n"
    input_Y_ξ₀ *= "[$(inf(real(Y_ξ₀[2]))), $(sup(real(Y_ξ₀[2])))]\n"
    input_Y_ξ₀ *= "[$(inf(imag(Y_ξ₀[1]))), $(sup(imag(Y_ξ₀[1])))]\n"
    input_Y_ξ₀ *= "[$(inf(imag(Y_ξ₀[2]))), $(sup(imag(Y_ξ₀[2])))]\n"
    input_Y_ξ₀ *= "[$(inf(real(Y_ξ₀[3]))), $(sup(real(Y_ξ₀[3])))]\n"
    input_Y_ξ₀ *= "[$(inf(real(Y_ξ₀[4]))), $(sup(real(Y_ξ₀[4])))]\n"
    input_Y_ξ₀ *= "[$(inf(imag(Y_ξ₀[3]))), $(sup(imag(Y_ξ₀[3])))]\n"
    input_Y_ξ₀ *= "[$(inf(imag(Y_ξ₀[4]))), $(sup(imag(Y_ξ₀[4])))]\n"
    input_params = "$(λ.d)\n"
    for x in [real(lambda), imag(lambda), κ, ϵ, λ.ω, λ.σ, λ.δ]
        input_params *= "[$(inf(x)), $(sup(x))]\n"
    end
    input_ξspan = ""
    for x in [ξ₀, ξ₁]
        input_ξspan *= "[$(inf(x)), $(sup(x))]\n"
    end
    input_output_jacobian = ifelse(output_jacobian isa Val{true}, "1\n", "0\n")
    input_tol = "$tol\n"

    input = join([input_Q_hat_ξ₀, input_Y_ξ₀, input_params, input_ξspan, input_output_jacobian, input_tol])

    # IMPROVE: Write directly to stdout of cmd instead of using echo
    program = pkgdir(@__MODULE__, "capd", "build", "Y")
    cmd = pipeline(`echo $input`, `$program`)

    output = try
        readchomp(cmd)
    catch e
        # If NaN occurs during the computation the program aborts. We
        # catch this and handle it in the same way as if an exception
        # was thrown during the computations.
        e isa ProcessFailedException || rethrow(e)

        "Exception"
    end

    #println(output)

    if contains(output, "Exception")
        n = if output_jacobian isa Val{true}
            80
        else
            8
        end
        res = fill(IntervalArithmetic.emptyinterval(Interval{Float64}), n)
    else
        res = parse.(Interval{Float64}, split(output, "\n"))::Vector{Interval{Float64}}
    end

    if output_jacobian isa Val{true}
        # IMPROVE: In principle it should be enough to compute only
        # the derivatives with respect to the real parts since
        # everything is analytic.
        J = SMatrix{8,10}(res)

        return SMatrix{4,5,Complex{Interval{Float64}}}(
            complex(J[1, 1], J[3, 1]), # ∂Y₁ / ∂Y₁(0)
            complex(J[2, 1], J[4, 1]), # ∂Y₂ / ∂Y₁(0)
            complex(J[5, 1], J[7, 1]), # ∂Z₁ / ∂Y₁(0)
            complex(J[6, 1], J[8, 1]), # ∂Z₂ / ∂Y₁(0)
            complex(J[1, 2], J[3, 2]), # ∂Y₁ / ∂Y₂(0)
            complex(J[2, 2], J[4, 2]), # ∂Y₂ / ∂Y₂(0)
            complex(J[5, 2], J[7, 2]), # ∂Z₁ / ∂Y₂(0)
            complex(J[6, 2], J[8, 2]), # ∂Z₂ / ∂Y₂(0)
            complex(J[1, 5], J[3, 5]), # ∂Y₁ / ∂Z₁(0)
            complex(J[2, 5], J[4, 5]), # ∂Y₂ / ∂Z₁(0)
            complex(J[5, 5], J[7, 5]), # ∂Z₁ / ∂Z₁(0)
            complex(J[6, 5], J[8, 5]), # ∂Z₂ / ∂Z₁(0)
            complex(J[1, 6], J[3, 6]), # ∂Y₁ / ∂Z₂(0)
            complex(J[2, 6], J[4, 6]), # ∂Y₂ / ∂Z₂(0)
            complex(J[5, 6], J[7, 6]), # ∂Z₁ / ∂Z₂(0)
            complex(J[6, 6], J[8, 6]), # ∂Z₂ / ∂Z₂(0)
            complex(J[1, 9], J[3, 9]), # ∂Y₁ / ∂λ
            complex(J[2, 9], J[4, 9]), # ∂Y₂ / ∂λ
            complex(J[5, 9], J[7, 9]), # ∂Z₁ / ∂λ
            complex(J[6, 9], J[8, 9]), # ∂Z₂ / ∂λ
        )
    else
        Y = SVector(
            complex(res[1], res[3]),
            complex(res[2], res[4]),
            complex(res[5], res[7]),
            complex(res[6], res[8]),
        )

        return SVector{4,Complex{Interval{Float64}}}(Y)
    end
end

"""
    Y_zero_capd(x, lambda, ν, κ, ϵ, ξ₁, λ::CGLParams; tol::Float64 = 1e-11, degree = 20)
    Y_zero_capd(x, lambda, ν, κ, ϵ, ξ₀, ξ₁, λ::CGLParams; tol::Float64 = 1e-11, degree = 20)

Compute the solution to the ODE on the interval ``[0, ξ₁]``. Returns a
vector with four complex values, the first two are the values at `ξ₁`
and the second two are their derivatives.

The solution is computed using the rigorous CAPD integrator.

If `ξ₀` is given then it uses a single Taylor expansion on the
interval `[0, ξ₀]` (of degree `degree`) and CAPD on `[ξ₀, ξ₁]`. For
`λ.d != 1` this is automatically used (`ξ₀ = 1e-2` by default) to
handle the removable singularity at zero.

If the given `ξ₀` gives a non-finite enclosure on `[0, ξ₀]`, then it
tries with half that value. If it fails again it tries to halve it
once more, iterating like this for a maximum of a few times.
"""
Y_zero_capd(
    x::Acb,
    lambda::Acb,
    ν::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    tol::Float64 = 1e-11,
    degree = 20,
) = Y_zero_capd(
    x,
    lambda,
    ν,
    κ,
    ϵ,
    ifelse(isone(λ.d), zero(Arb), Arb(1e-2)),
    ξ₁,
    λ;
    tol,
    degree,
)

function Y_zero_capd(
    x::Acb,
    lambda::Acb,
    ν::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₀::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    tol::Float64 = 1e-11,
    degree = 20,
)
    S = Interval{Float64}

    Q_hat_ξ₀ = if !iszero(ξ₀) && !iszero(ν)
        @assert 0 < ξ₀ < ξ₁
        # Integrate system on [0, ξ₀] using Taylor expansion at zero
        Q_hat_ξ₀ = Q_hat_zero_taylor(real(ν), imag(ν), κ, ϵ, ξ₀, λ)
        if !all(isfinite, Q_hat_ξ₀)
            iterations = 0
            while !all(isfinite, Q_hat_ξ₀) && iterations < 5
                iterations += 1
                ξ₀ /= 2
                Q_hat_ξ₀ = Q_hat_zero_taylor(real(ν), imag(ν), κ, ϵ, ξ₀, λ)
            end
            iterations == 5 && @debug "Non-finite enclosure for smallest ξ₀" ξ₀
        end
        convert(SVector{4,S}, Q_hat_ξ₀)
    else
        SVector{4,S}(real(ν), imag(ν), interval(0.0), interval(0.0))
    end

    Y_ξ₀ = if !iszero(ξ₀)
        @assert 0 < ξ₀ < ξ₁
        # Integrate system on [0, ξ₀] using Taylor expansion at zero
        Y_ξ₀ = Y_zero_taylor(x, lambda, ν, κ, ϵ, ξ₀, λ; degree)
        if !all(isfinite, Y_ξ₀)
            iterations = 0
            while !all(isfinite, Y_ξ₀) && iterations < 5
                iterations += 1
                ξ₀ /= 2
                Y_ξ₀ = Y_zero_taylor(x, lambda, ν, κ, ϵ, ξ₀, λ; degree)
            end
            iterations == 5 && @debug "Non-finite enclosure for smallest ξ₀" ξ₀
        end
        convert(SVector{4,Complex{S}}, Y_ξ₀)
    else
        SVector{4,Complex{S}}(x, interval(1.0), interval(0.0), interval(0.0))
    end

    # Integrate system on [ξ₀, ξ₁] using capd.
    # We use the fact that the equation is identical to the one for Q,
    # except for the change in sign for κ and ω.
    Y = _Y_zero_capd(
        Q_hat_ξ₀,
        Y_ξ₀,
        convert(Complex{S}, lambda),
        convert(S, κ),
        convert(S, ϵ),
        convert(S, ξ₀),
        convert(S, ξ₁),
        CGLParams{S}(λ);
        tol,
    )

    return Acb.(Y)
end

"""
    Y_zero_jacobian_capd(x, lambda, ν, κ, ϵ, ξ₁, λ::CGLParams; tol::Float64 = 1e-11, degree = 20)
    Y_zero_jacobian_capd(x, lambda, ν, κ, ϵ, ξ₀, ξ₁, λ::CGLParams; tol::Float64 = 1e-11, degree = 20)

This function computes the Jacobian of [`Y_zero_capd`](@ref) w.r.t.
the parameters `x` and `lambda`.

Similar to [`Q_zero_capd`](@ref) the solution is computed using the
rigorous CAPD integrator.
"""
Y_zero_jacobian_capd(
    x::Acb,
    lambda::Acb,
    ν::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    tol::Float64 = 1e-11,
    degree = 20,
) = Y_zero_jacobian_capd(
    x,
    lambda,
    ν,
    κ,
    ϵ,
    ifelse(isone(λ.d), zero(Arb), Arb(1e-2)),
    ξ₁,
    λ;
    tol,
    degree,
)

function Y_zero_jacobian_capd(
    x::Acb,
    lambda::Acb,
    ν::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₀::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    tol::Float64 = 1e-11,
    degree = 20,
)
    S = Interval{Float64}

    Q_hat_ξ₀ = if !iszero(ξ₀)
        @assert 0 < ξ₀ < ξ₁
        # Integrate system on [0, ξ₀] using Taylor expansion at zero
        Q_hat_ξ₀ = Q_hat_zero_taylor(real(ν), imag(ν), κ, ϵ, ξ₀, λ)
        if !all(isfinite, Q_hat_ξ₀)
            iterations = 0
            while !all(isfinite, Q_hat_ξ₀) && iterations < 5
                iterations += 1
                ξ₀ /= 2
                Q_hat_ξ₀ = Q_hat_zero_taylor(real(ν), imag(ν), κ, ϵ, ξ₀, λ)
            end
            iterations == 5 && @debug "Non-finite enclosure for smallest ξ₀" ξ₀
        end
        convert(SVector{4,S}, Q_hat_ξ₀)
    else
        SVector{4,S}(real(ν), imag(ν), interval(0.0), interval(0.0))
    end

    Y_ξ₀, J_ξ₀ = let
        if !iszero(ξ₀)
            #@assert 0 < ξ₀ < ξ₁
            # Integrate system on [0, ξ₀] using Taylor expansion at zero
            Y_ξ₀, J_ξ₀ = Y_zero_jacobian_taylor(x, lambda, ν, κ, ϵ, ξ₀, λ; degree)
            if !all(isfinite, Y_ξ₀)
                iterations = 0
                while !all(isfinite, Y_ξ₀) && iterations < 5
                    iterations += 1
                    ξ₀ /= 2
                    Y_ξ₀, J_ξ₀ =
                        Y_zero_jacobian_taylor(x, lambda, ν, κ, ϵ, ξ₀, λ; degree)
                end
                iterations == 5 && @debug "Non-finite enclosure for smallest ξ₀" ξ₀
            end
            Y_ξ₀ = convert(SVector{4,Complex{S}}, Y_ξ₀)
            J_ξ₀ = convert(SMatrix{4,2,Complex{S}}, J_ξ₀)
        else
            Y_ξ₀ = SVector{4,Complex{S}}(x, interval(1.0), interval(0.0), interval(0.0))
            # Empty integration so the only non-zero derivative is the
            # one of Q_ξ₀[1] w.r.t. μ, which is 1.
            J_ξ₀ = SMatrix{4,2,Complex{S}}(
                interval(1.0),
                interval(0.0),
                interval(0.0),
                interval(0.0),
                interval(0.0),
                interval(0.0),
                interval(0.0),
                interval(0.0),
            )
        end

        # J_ξ₀ now contains derivatives of Q_ξ₀. We want to add a row
        # [0, 1] for the derivative of lambda.
        Y_ξ₀, vcat(J_ξ₀, SMatrix{1,2,Complex{S}}(interval(0.0), interval(1.0)))
    end

    # Integrate system on [ξ₀, ξ₁] using capd.
    # We use the fact that the equation is identical to the one for Q,
    # except for the change in sign for κ and ω.
    J_ξ₀_ξ₁ = _Y_zero_capd(
        Q_hat_ξ₀,
        Y_ξ₀,
        convert(Complex{S}, lambda),
        convert(S, κ),
        convert(S, ϵ),
        convert(S, ξ₀),
        convert(S, ξ₁),
        CGLParams{S}(λ),
        output_jacobian = Val{true}();
        tol,
    )

    # The Jacobian on the interval [0, ξ₁] is the product of the one
    # on [0, ξ₀] and the one on [ξ₀, ξ₁].
    J = J_ξ₀_ξ₁ * J_ξ₀

    return Acb.(J)
end
