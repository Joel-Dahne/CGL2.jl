"""
    _Y_zero_capd(
        Y_ξ₀::SVector{4,Acb},
        κ::Arb,
        ξ₀::Arb,
        ξ₁::Arb,
        Λ::CGLParams{Arb};
        output_jacobian::Union{Val{false},Val{true}} = Val(false),
        wrt_epsilon::Bool = false,
        include_parameter_derivatives::Union{Val{false},Val{true}} = Val(true),
        output_curve::Union{Val{false},Val{true}} = Val(true),
        tol::Float64 = 1e-11,
    )

Internal function for calling the CAPD program.

It computes the solution to the ODE on the interval ``[ξ₀, ξ₁]`` using
the CAPD C++ library. See `CAPD/README.md` for more information about
the C++ implementation. The initial value at `ξ₀` is given by
`Q_hat_ξ₀` and `Y_ξ₀`.

It returns 4 complex values, the first two are the values at `ξ₁` and
the second two are their derivatives.
"""
function _Y_zero_capd(
    Q_hat_ξ₀::SVector{4,Arb},
    Y_ξ₀::SVector{4,Acb},
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₀::Arb,
    ξ₁::Arb,
    Λ::CGLParams{Arb};
    tol::Float64 = 1e-11,
)
    # Convenience functions for computing lower and upper bounds as Float64
    _inf(x::Arb) = Float64(lbound(x), RoundDown)
    _sup(x::Arb) = Float64(ubound(x), RoundUp)

    # Run the C++ program
    exit_success, output = try
        # If the program aborts it writes to stderr. To avoid
        # cluttering the output we redirect this to devnull. (What it
        # writes is not particularly useful information.)
        open(
            pipeline(`$(pkgdir(@__MODULE__, "CAPD", "build", "Y_zero"))`, stderr = devnull),
            "w+",
        ) do io
            # Write initial value
            for x in Q_hat_ξ₀
                println(io, "[$(_inf(x)), $(_sup(x))]")
            end
            println(io, "[$(inf(real(Y_ξ₀[1]))), $(sup(real(Y_ξ₀[1])))]")
            println(io, "[$(inf(real(Y_ξ₀[2]))), $(sup(real(Y_ξ₀[2])))]")
            println(io, "[$(inf(imag(Y_ξ₀[1]))), $(sup(imag(Y_ξ₀[1])))]")
            println(io, "[$(inf(imag(Y_ξ₀[2]))), $(sup(imag(Y_ξ₀[2])))]")
            println(io, "[$(inf(real(Y_ξ₀[3]))), $(sup(real(Y_ξ₀[3])))]")
            println(io, "[$(inf(real(Y_ξ₀[4]))), $(sup(real(Y_ξ₀[4])))]")
            println(io, "[$(inf(imag(Y_ξ₀[3]))), $(sup(imag(Y_ξ₀[3])))]")
            println(io, "[$(inf(imag(Y_ξ₀[4]))), $(sup(imag(Y_ξ₀[4])))]")
            # Write parameters
            println(io, Λ.d)
            for x in [real(lambda), imag(lambda), κ, ϵ, Λ.ω, Λ.σ, Λ.δ]
                println(io, "[$(_inf(x)), $(_sup(x))]")
            end
            # Write integration interval
            println(io, "[$(_inf(ξ₀)), $(_sup(ξ₀))]")
            println(io, "[$(_inf(ξ₁)), $(_sup(ξ₁))]")
            # Write settings
            println(io, tol)
            close(io.in)

            output = readchomp(io)
            exit_success = success(io)::Bool
            exit_success, output
        end
    catch e
        e isa ProcessFailedException || rethrow(e)

        # If NaN occurs during the computation the program aborts. We
        # catch this return a failed exit code and empty output.
        false, ""
    end

    # To simplify the implementation we use parse(Interval{Float64},
    # str) for parsing the output from the program.

    if !exit_success
        Y_reim = fill(nai(Float64), 8)
    else
        Y_reim = parse.(Interval{Float64}, split(output, "\n"))::Vector{Interval{Float64}}
    end

    Y = SVector(
        complex(Y_reim[1], Y_reim[3]),
        complex(Y_reim[2], Y_reim[4]),
        complex(Y_reim[5], Y_reim[7]),
        complex(Y_reim[6], Y_reim[8]),
    )

    return SVector{4,Acb}(Y)
end

"""
    Y_zero_capd(Y₀, lambda, ν, κ, ϵ, ξ₁, λ::CGLParams; ξ₀, tol, degree = 20)

Compute the solution to the ODE on the interval ``[0, ξ₁]`` with
initial values given by `Y₀`. Returns a vector with four complex
values, where the first two are the values at `ξ₁` and the last two
are the derivatives.

The solution is computed using the rigorous CAPD integrator.

If `ξ₀` is non-zero it uses a single Taylor expansion on the interval
`[0, ξ₀]` and the CAPD integrator on `[ξ₀, ξ₁]`. This is needed to
avoid the removable singularity at `ξ = 0` which CAPD cannot handle
directly. For `Λ.d = 1` there is no removable singularity and the
default value is `ξ₀ = 0`, otherwise the default value is `ξ₀ = 1e-2`.

If the given `ξ₀` gives a non-finite enclosure on `[0, ξ₀]`, then it
tries with half that value. If it fails again it tries to halve it
once more, iterating like this for a maximum of a few times.
"""
function Y_zero_capd(
    Y₀::SVector{2,Acb},
    lambda::Acb,
    ν::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    ξ₀::Arb = ifelse(isone(λ.d), zero(Arb), Arb(1e-2)),
    tol::Float64 = 1e-11,
    degree = 20,
)
    Q_hat_ξ₀, Y_ξ₀ = if !iszero(ξ₀)
        @assert 0 < ξ₀ < ξ₁
        # Integrate system on [0, ξ₀] using Taylor expansion at zero
        Q_hat_ξ₀ = Q_hat_zero_taylor(real(ν), imag(ν), κ, ϵ, ξ₀, λ)
        Y_ξ₀ = Y_zero_taylor(Y₀, lambda, ν, κ, ϵ, ξ₀, λ; degree)
        if !(all(isfinite, Q_hat_ξ₀) && all(isfinite, Y_ξ₀))
            iterations = 0
            while !(all(isfinite, Q_hat_ξ₀) && all(isfinite, Y_ξ₀)) && iterations < 5
                iterations += 1
                ξ₀ /= 2
                Q_hat_ξ₀ = Q_hat_zero_taylor(real(ν), imag(ν), κ, ϵ, ξ₀, λ)
                Y_ξ₀ = Y_zero_taylor(Y₀, lambda, ν, κ, ϵ, ξ₀, λ; degree)
            end
            iterations == 5 && @debug "Non-finite enclosure for smallest ξ₀" ξ₀
        end
        Q_hat_ξ₀, Y_ξ₀
    else
        SVector{4,Arb}(real(ν), imag(ν), 0, 0),
        SVector{4,Acb}(interval(Y₀[1]), interval(Y₀[2]), 0, 0)
    end

    # Integrate system on [ξ₀, ξ₁] using capd.
    return _Y_zero_capd(Q_hat_ξ₀, Y_ξ₀, lambda, κ, ϵ, ξ₀, ξ₁, λ; tol)
end
