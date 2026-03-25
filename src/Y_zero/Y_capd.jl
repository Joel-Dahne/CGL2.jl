"""
    _Y_zero_capd(
        Q_hat_ξ₀::SVector{4,Interval{Float64}},
        Y_ξ₀::SVector{4,Complex{Interval{Float64}}},
        lambda::Complex{Interval{Float64}},
        ν::Complex{Interval{Float64}},
        κ::Interval{Float64},
        ϵ::Interval{Float64},
        ξ₀::Interval{Float64},
        ξ₁::Interval{Float64},
        λ::CGLParams{Interval{Float64}};
        tol::Float64 = 1e-11,
    )

Internal function for calling the CAPD program.

It computes the solution to the ODE on the interval ``[ξ₀, ξ₁]`` and
returns four complex values, the first two are the values at `ξ₁` and
the second two are their derivatives.
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
    input_tol = "$tol\n"

    input = join([input_Q_hat_ξ₀, input_Y_ξ₀, input_params, input_ξspan, input_tol])

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

    if contains(output, "Exception")
        res = fill(emptyinterval(Interval{Float64}), 8)
    else
        res = parse.(Interval{Float64}, split(output, "\n"))::Vector{Interval{Float64}}
    end

    Y = SVector(
        complex(res[1], res[3]),
        complex(res[2], res[4]),
        complex(res[5], res[7]),
        complex(res[6], res[8]),
    )

    return SVector{4,Complex{Interval{Float64}}}(Y)
end

"""
    Y_zero_capd(Y₀, lambda, ν, κ, ϵ, ξ₁, λ::CGLParams; tol::Float64 = 1e-11, degree = 20)
    Y_zero_capd(Y₀, lambda, ν, κ, ϵ, ξ₀, ξ₁, λ::CGLParams; tol::Float64 = 1e-11, degree = 20)

Compute the solution to the ODE on the interval ``[0, ξ₁]`` with
initial values given by `Y₀`. Returns a vector with four complex
values, where the first two are the values at `ξ₁` and the last two
are the derivatives.

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
    Y₀::SVector{2,Acb},
    lambda::Acb,
    ν::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb};
    tol::Float64 = 1e-11,
    degree = 20,
) = Y_zero_capd(
    Y₀,
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
    Y₀::SVector{2,Acb},
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

    Q_hat_ξ₀, Y_ξ₀ = if !iszero(ξ₀) && !iszero(ν)
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
        convert(SVector{4,S}, Q_hat_ξ₀), convert(SVector{4,Complex{S}}, Y_ξ₀)
    else
        SVector{4,S}(real(ν), imag(ν), interval(0.0), interval(0.0)),
        SVector{4,Complex{S}}(
            interval(Y₀[1]),
            interval(Y₀[2]),
            interval(0.0),
            interval(0.0),
        )
    end

    # Integrate system on [ξ₀, ξ₁] using capd.
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
