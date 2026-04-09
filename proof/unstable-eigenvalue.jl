### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ e914aa7a-31fe-11f1-b1c4-7dea491a45e5
begin
    using Pkg
    Pkg.activate(".", io = devnull)
    using ArbExtras
    using Arblib
    using CairoMakie
    using CGL2
    using LaTeXStrings
    using OhMyThreads
    using PlutoUI

    import CairoMakie.Makie.GeometryBasics: Rect2f

    # Set the bits of precision used for the computations with Arblib
    setprecision(Arb, 128)

    nothing
end

# ╔═╡ 8e9f1365-8c0e-41cd-b094-663a3760c0d1
md"""
# Existence of an unstable eigenvalue
"""

# ╔═╡ baa97573-bc3e-433b-88f1-a7c0a826d27a
TableOfContents()

# ╔═╡ abc2a340-ab72-4340-838c-bfff029ddd77
md"""
Check this box to set the code to save the figures.
- Save figures $(@bind save_figures CheckBox(default = false))
"""

# ╔═╡ 92972793-3445-428e-9ffb-3d52dc62e050
fontsize = 22

# ╔═╡ ad65706c-75ec-4ddb-9d1b-571fbe4998d1
λ = CGLParams{Arb}(3, 1, 1, 0)

# ╔═╡ da71a660-d085-42c6-851a-e7d28f472574
ϵ = Arb("0.1681")

# ╔═╡ 06b5e1d2-e5f2-4e1d-ad61-d086dc37dcd2
ξ₁ = Arb(16)

# ╔═╡ fbbaad41-c446-4c0c-8caa-8137713394b9
md"""
## Existence of a backward self-similar solution

The first step is to prove the existence of a backward self-similar solution, and also compute enclosures of the associated parameters.
"""

# ╔═╡ 6a2cd7c8-2760-42a2-a7d1-5dc89550cdd5
md"""
Our starting point is these rough initial approximations. Note that the approximation for ``\gamma`` is computed automatically by `refine_approximation_fix_epsilon`.
"""

# ╔═╡ 4279b8df-7e4e-457f-b67c-c9637604ca64
μ_approx = Arb(2.2600451)

# ╔═╡ 0e924af2-bb12-4fe2-b0e5-ad4001c9d06d
κ_approx = Arb(0.7912471)

# ╔═╡ 185f4a22-ad1b-4d27-863d-ff68d3c34809
md"""
These approximations are then refined using a Newton-Raphson method.
"""

# ╔═╡ 098d9c62-a9c4-4a3e-a1f7-7d466f136f0f
μ₀, γ₀, κ₀ = CGL2.refine_approximation_fix_epsilon(μ_approx, κ_approx, ϵ, ξ₁, λ)

# ╔═╡ 93669102-c090-4a86-bfb2-9609fbfe924f
md"""
We print 16 digit strings of these approximations, for inclusion in the paper.
"""

# ╔═╡ adc8d719-f492-4633-907e-15b23e28940d
string(μ₀, digits = 16, no_radius = true)

# ╔═╡ 16d83b4f-d023-462a-afa6-4cef99f44196
string(γ₀, digits = 16, no_radius = true)

# ╔═╡ 3e607803-f77c-4467-9967-61bbb4f0cedb
string(κ₀, digits = 16, no_radius = true)

# ╔═╡ 01f5801b-8605-4c1a-b07d-670abfd451e9
md"""
The application of the interval Newton method is handled by `G_solve_fix_epsilon`. Note that this uses a heuristic choice for the box around the given approximation on which we apply the interval Newton operator.
"""

# ╔═╡ baa627cb-3a2c-4a48-89ee-515ceaed7049
μ, γ_real, γ_imag, κ =
    CGL2.G_solve_fix_epsilon(μ₀, real(γ₀), imag(γ₀), κ₀, ϵ, ξ₁, λ, verbose = true)

# ╔═╡ 683c7c63-57df-4bb6-9a3d-96dfbcd85784
γ = Acb(γ_real, γ_imag)

# ╔═╡ bbe48cc0-2e1e-46ec-b9b2-4eccff0c4018
@assert_proof isfinite(μ) && isfinite(γ) && isfinite(κ)

# ╔═╡ 8b16d819-ae9a-407d-a8fc-61bff5048789
md"""
Print formatted enclosures of the zero, for inclusion in the paper.
"""

# ╔═╡ f54a47a9-704c-4536-a38f-e935bad69f01
CGL2.format_interval_precise(μ)

# ╔═╡ 2320539d-6309-4826-b72d-f9b9b5459abf
CGL2.format_interval_precise(γ)

# ╔═╡ c16e015c-d8c6-4684-91c7-09d49c98ccf1
CGL2.format_interval_precise(κ)

# ╔═╡ 3032a818-401b-416a-9d74-1340a98df7db
md"""
Finally, we compute an enclosure of the coefficient ``p_Q`` determining the leading order asymptotic at infinity.

TODO: Rename `p_Q_0` to `p_Q`?
"""

# ╔═╡ f044e3bf-87b5-4d4f-8bcd-b0ec3f816676
p_Q = CGL2.p_Q_0(γ, κ, ϵ, ξ₁, λ)

# ╔═╡ 5888a3a0-7064-4b6b-9428-42fcb61d396a
md"""
Print a formatted enclosure for inclusion in the paper.
"""

# ╔═╡ fd9e8538-381b-44c4-a737-c03f9ad4b590
CGL2.format_interval_precise(p_Q)

# ╔═╡ 48ddc751-4a96-4a04-8e71-040430d203f6
let
    fig = Figure(; fontsize)
    ax = Axis(fig[1, 1], xlabel = L"\xi", ylabel = L"|Q(\xi)|")
    ξs = range(0, 20, 1000)
    Qs = CGL2.Q_zero_float_curve(
        Float64(μ),
        Float64(κ),
        Float64(ϵ),
        ξs[end],
        CGLParams{Float64}(λ),
    ).(
        ξs,
    )
    lines!(ax, ξs, abs.(complex.(getindex.(Qs, 1), getindex.(Qs, 2))), linewidth = 3)
    save_figures && save("figures/Q.pdf", fig)
    fig
end

# ╔═╡ 56193256-54c7-425d-a0dc-56227e009aba
md"""
## Existence of a forward self-similar solution

Next we prove the existance of a forward self-similar solution with the same asymptotic behavior as the above backward self-similar solution. The first step is to determine the parameter ``\gamma_1`` such that the forward solution has the right asymptotic behavior.
"""

# ╔═╡ b4049ea0-ee9f-409c-b2e3-2fdbabe1b8a9
a, b, c = CGL2._abc(κ, ϵ, λ)

# ╔═╡ c2f067f3-4a78-4048-834d-186ec249350e
γ₁ = p_Q / (-c)^-a

# ╔═╡ 40f63eb4-5ece-44e7-835a-b33f97e309bd
@assert_proof isfinite(γ₁)

# ╔═╡ c3ad0d66-d491-42d0-96f4-147350ae6ffa
md"""
Print a formatted enclosures of ``\gamma_1`` for inclusion in the paper.
"""

# ╔═╡ a2feb21f-1cd5-45e2-9f4d-424a23f31788
CGL2.format_interval_precise(γ₁)

# ╔═╡ b50ba676-2872-444e-98a2-1cf97f2dbd3a
md"""
Next we compute a numerical approximation for ``\nu`` and ``\gamma_2``. This is based on a Newton-Raphson method applied to the initial guess ``(\nu, \gamma_2) = (0, 0)``.
"""

# ╔═╡ c5aa3ed9-31a4-477e-9066-b8045b2cb829
ν_approx, γ₂_approx = Acf.(
    CGL2.G_hat_approximate(
        ComplexF64(γ₁),
        Float64(κ),
        Float64(ϵ),
        Float64(ξ₁),
        CGLParams{Float64}(λ),
    ),
)

# ╔═╡ c4fa5393-8708-4643-a95b-9cd398f96364
md"""
We print 16 digit strings of these approximations, for inclusion in the paper.
"""

# ╔═╡ d8d43055-8663-4824-8b92-d6ab45f1cffc
string(ν_approx, digits = 16)

# ╔═╡ 4ca4c796-7cf4-447c-9f1a-2cc775e3eb16
string(γ₂_approx, digits = 16)

# ╔═╡ 6c018849-3827-4a45-bd07-cc8e746a44a4
md"""
The application of the interval Newton method is handled by `G_solve_hat`. Note that this uses a heuristic choice for the box around the given approximation on which we apply the interval Newton operator.
"""

# ╔═╡ 75653ff9-94ca-4c7b-a0c5-363b1de657e0
ν, γ₂ = CGL2.G_hat_solve(ν_approx, γ₁, γ₂_approx, κ, ϵ, ξ₁, λ, verbose = true)

# ╔═╡ 3df5f78f-d780-4db7-bd39-a858e04f650a
@assert_proof isfinite(ν) && isfinite(γ₂)

# ╔═╡ feeea3b3-5de0-42d6-8599-6f19ee3853e3
md"""
Print formatted enclosures of the zero, for inclusion in the paper.
"""

# ╔═╡ db74a6b9-e8ed-4139-af21-5ea7fc8d55a0
CGL2.format_interval_precise(ν)

# ╔═╡ 7bffa0f1-b6ae-4367-945c-a24841a38f1b
CGL2.format_interval_precise(γ₂)

# ╔═╡ 91492771-4441-4e1e-b129-38026188a2d2
let
    fig = Figure(; fontsize)
    ax = Axis(fig[1, 1], xlabel = L"\xi", ylabel = L"|\hat{Q}(\xi)|")
    ξs = range(0, 20, 1000)
    Q_hats = CGL2.Q_hat_zero_float_curve(
        Float64(real(ν)),
        Float64(imag(ν)),
        Float64(κ),
        Float64(ϵ),
        ξs[end],
        CGLParams{Float64}(λ),
    ).(
        ξs,
    )
    lines!(
        ax,
        ξs,
        abs.(complex.(getindex.(Q_hats, 1), getindex.(Q_hats, 2))),
        linewidth = 3,
    )
    save_figures && save("figures/Q-hat.pdf", fig)
    fig
end

# ╔═╡ 3bf117b4-8f2b-4c88-a1f4-a9b599d74af4
md"""
## Existence of an unstable eigenvalue

Finally, we prove the existence of an unstable eigenvalue for the linear operator associated with the above forward self-similar solution.
"""

# ╔═╡ 880df8b3-36ff-4a74-9f9c-428bf9663462
md"""
### Approximate eigenvalue

The first step is to compute an approximate solution using a finite-difference method.
"""

# ╔═╡ f1ea08ed-8969-448f-859b-8aea8ba546de
# Find approximation using finite differences
λsF64 = filter(
    lambda -> imag(lambda) > 0,
    CGL2.linearization_eigenvalues_real_1(
        ComplexF64(ν),
        Float64(κ),
        Float64(ϵ),
        Float64(ξ₁),
        CGLParams{Float64}(λ),
    )[1],
)

# ╔═╡ 21cf1d05-df4a-4853-b1be-66bb66e5b35f
lambda_approx = Acf(λsF64[findmax(real, λsF64)[2]])

# ╔═╡ 71b8a3b1-3227-44b9-8f06-50c9192ac96f
md"""
Print a 16 digit approximation for inclusion in the paper.
"""

# ╔═╡ 7d38f474-0f8a-4bc4-8f33-fbb31c540703
string(lambda_approx, digits = 8)

# ╔═╡ 4ed7e92b-0205-4774-a4a2-42587d385db4
md"""
### Construction of contour

The next step is to take a contour around the approximate solution. For this we take a rectangle with side lengths (approximately) given by
"""

# ╔═╡ 2f6b6452-f2e1-4c07-9511-93f355ba9707
r = Arb("0.015")

# ╔═╡ 0d7d7d0a-ae40-41a4-a41d-074cb617a466
md"""
To plot the square we compute enclosures of the corners.
"""

# ╔═╡ eb8e5ece-2782-4f2d-b519-93a9f336e5d4
corner_bl = lambda_approx + Acb(-r, -r)

# ╔═╡ ea45a2a5-7e79-44fc-a490-ac4889aa3b1f
corner_br = lambda_approx + Acb(r, -r)

# ╔═╡ 0f5d4f04-0aee-4582-a504-868816112aa4
corner_tl = lambda_approx + Acb(-r, r)

# ╔═╡ 933ced1e-a9a6-495b-865e-d9a8a3cf7436
corner_tr = lambda_approx + Acb(r, r)

# ╔═╡ 04b37248-22a5-4fb8-8d0c-a4dbf29a58e4
md"""
We can then plot the square and the approximate eigenvalue.
"""

# ╔═╡ 356dc846-b0fc-4743-997a-73d9b9bf7bcf
let
    fig = Figure(; fontsize)
    ax = Axis(fig[1, 1], xlabel = L"\mathrm{Re}(\lambda)", ylabel = L"\mathrm{Im}(\lambda)")
    xlims!(ax, -0.005, 0.04)
    corner_list = [corner_bl, corner_br, corner_tr, corner_tl, corner_bl]
    lines!(ax, real.(corner_list), imag.(corner_list), linewidth = 3)
    scatter!(
        ax,
        [real(lambda_approx)],
        [imag(lambda_approx)],
        color = :red,
        label = L"\lambda_0",
    )
    axislegend(ax)
    save_figures && save("figures/contour.pdf", fig)
    fig
end

# ╔═╡ 03513e93-3fda-4e9a-b02d-1ed1552811f6
md"""
For printing the enclosure in the paper, the simplest approach is to construct an `Acb` enclosing the square. Due to rounding of the radius, this `Acb` will not be exactly the same as the square but will be slightly larger. For purposes of giving an enclosure this is however fine.
"""

# ╔═╡ f951ad64-9b79-477e-83fb-28fa2bba6e2a
CGL2.format_interval_precise(add_error(Acb(lambda_approx), r), min_digits = 4)

# ╔═╡ 93d9e75b-5861-4653-8786-439a35008d00
md"""
### Computation of winding number

Finally, we compute the winding argument of ``H(C)``.
"""

# ╔═╡ a17f3ef0-9476-4ebd-b52e-657d4fede475
ts_bottom = Arb.(ArbExtras.bisect_interval_recursive(Arb(0), Arb(1), 7))

# ╔═╡ 6b411c4f-1000-4ebb-84a6-c020e9a33b7c
ts_right = Arb.(ArbExtras.bisect_interval_recursive(Arb(1), Arb(2), 7))

# ╔═╡ fd79a36b-af81-427e-aefb-a968bc1cb42a
ts_top = Arb.(ArbExtras.bisect_interval_recursive(Arb(2), Arb(3), 7))

# ╔═╡ ec096027-d898-46e7-abe5-7bcc72229f04
ts_left = Arb.(ArbExtras.bisect_interval_recursive(Arb(3), Arb(4), 7))

# ╔═╡ f72808a0-5019-49eb-9f8d-41d72b970f0c
md"""
We then use these to construct the pieces for the different sides of the square.
"""

# ╔═╡ 3e493c1b-b786-4cd3-b321-0de29ae04bf0
λs_bottom = corner_bl .+ ts_bottom .* (corner_br - corner_bl)

# ╔═╡ 7d7217b1-5186-4d2f-97ab-f2581add502e
λs_right = corner_br .+ (ts_right .- 1) .* (corner_tr - corner_br)

# ╔═╡ c0a0d2fd-9a35-4028-b4e6-dbcf4e777613
λs_top = corner_tr .+ (ts_top .- 2) .* (corner_tl - corner_tr)

# ╔═╡ ca27e948-94c4-44a9-93c8-8b83493eb479
λs_left = corner_tl .+ (ts_left .- 3) .* (corner_bl - corner_tl)

# ╔═╡ 6f7407fb-8d41-451a-a220-80bdcc1fe881
md"""
Compute ``H`` for each piece.
"""

# ╔═╡ 340f6f00-629a-4c52-8f07-baf6568a2957
H_square_bottom = tmap(lambda -> CGL2.H(Acb(lambda), ν, γ₁, γ₂, κ, ϵ, ξ₁, λ), λs_bottom)

# ╔═╡ 060bd72c-1e0f-47cf-bbf1-f03c6e8c3599
H_square_right = tmap(lambda -> CGL2.H(Acb(lambda), ν, γ₁, γ₂, κ, ϵ, ξ₁, λ), λs_right)

# ╔═╡ da7f234b-b681-4956-ba49-7f2ac6b363c2
H_square_top = tmap(lambda -> CGL2.H(Acb(lambda), ν, γ₁, γ₂, κ, ϵ, ξ₁, λ), λs_top)

# ╔═╡ 942b9ffc-f12a-4c3c-9002-ffb69c784fd1
H_square_left = tmap(lambda -> CGL2.H(Acb(lambda), ν, γ₁, γ₂, κ, ϵ, ξ₁, λ), λs_left)

# ╔═╡ 193f2a06-c0d3-45f5-8a2f-ba65d5a41ab3
to_rect(x::Arb, y::Arb) = Rect2d(lbound(x), lbound(y), 2radius(x), 2radius(y))

# ╔═╡ 086a8256-b0a7-41a6-9fff-c99c7f953e45
to_rect(z::Acb) = to_rect(real(z), imag(z))

# ╔═╡ 7cffb3af-9329-4d5a-b848-6aa1f500b0e7
abs_tight(z::Acb) = Arb((abs_lbound(z), abs_ubound(z)))

# ╔═╡ 68b90902-028e-4715-b477-2b38bc0d8e4d
md"""
Next we check the following four conditions:

1. For ``t \in [0, 1]`` (the bottom part) the curve ``H(c(t))`` has positive real part;
2. For ``t \in [1, 2]`` (the right part) the curve ``H(c(t))`` has positive imaginary part;
3. For ``t \in [2, 3]`` (the top part) the curve ``H(c(t))`` has negative real part;
4. For ``t \in [3, 4]`` (the left part) the curve ``H(c(t))`` has negative imaginary part;
"""

# ╔═╡ 5ebbff59-5819-467b-bd72-18b29d49e478
@assert_proof all(Arblib.ispositive, real(H_square_bottom))

# ╔═╡ 44ae8114-ba68-4a4d-8e78-9f462091f6cf
@assert_proof all(Arblib.isnegative, real(H_square_top))

# ╔═╡ 5b6e9013-1951-4fe1-9a9a-618ac4afcd48
@assert_proof all(Arblib.ispositive, imag(H_square_right))

# ╔═╡ 57ea897c-5f5f-4adf-af0c-bf48caf8aaeb
@assert_proof all(Arblib.isnegative, imag(H_square_left))

# ╔═╡ 6af950e3-91a0-4103-aad5-ee45d1d17e14
let
    fig = Figure(; fontsize)
    ax = Axis(
        fig[1, 1],
        xlabel = L"\mathrm{Re}(H(\lambda))",
        ylabel = L"\mathrm{Im}(H(\lambda))",
        xticks = [-2e-20, 0, 2e-20],
        yticks = [-2e-20, 0, 2e-20],
    )
    xlims!(ax, -2e-20, 2e-20)
    ylims!(ax, -2e-20, 2e-20)
    poly!(ax, to_rect.(H_square_top), strokewidth = 1, alpha = 0.5)
    poly!(ax, to_rect.(H_square_bottom), strokewidth = 1, alpha = 0.5)
    poly!(ax, to_rect.(H_square_left), strokewidth = 1, alpha = 0.5)
    poly!(ax, to_rect.(H_square_right), strokewidth = 1, alpha = 0.5)
    #save_figures && save("figures/contour-image.pdf", fig)
    fig
end

# ╔═╡ 9d0a1908-ea5e-4574-b1a6-3a90c79ed00e
let
    fig = Figure(; fontsize)
    ax = Axis(
        fig[1, 1],
        xlabel = L"t",
        ylabel = L"|H(c(t))|",
        yticks = [0, 1e-20, 2e-20, 3e-20],
    )
    ylims!(ax, 0, 3e-20)
    poly!(
        ax,
        to_rect.(ts_bottom, abs_tight.(H_square_bottom)),
        strokewidth = 1,
        alpha = 0.25,
    )
    poly!(ax, to_rect.(ts_right, abs_tight.(H_square_right)), strokewidth = 1, alpha = 0.25)
    poly!(ax, to_rect.(ts_top, abs_tight.(H_square_top)), strokewidth = 1, alpha = 0.25)
    poly!(ax, to_rect.(ts_left, abs_tight.(H_square_left)), strokewidth = 1, alpha = 0.25)
    save_figures && save("figures/H-abs.pdf", fig)
    fig
end

# ╔═╡ f00a20db-3b35-4c0c-8358-e00c98a36385
let
    fig = Figure(; fontsize)
    ax = Axis(
        fig[1, 1],
        xlabel = L"t",
        ylabel = L"\mathrm{Arg}(H(c(t)))",
        yticks = (
            [-π, -π / 2, 0, π / 2, π],
            [L"-\pi", L"-\frac{\pi}{2}", L"0", L"\frac{\pi}{2}", L"\pi"],
        ),
    )
    poly!(ax, to_rect.(ts_bottom, angle.(H_square_bottom)), strokewidth = 1, alpha = 0.25)
    poly!(ax, to_rect.(ts_right, angle.(H_square_right)), strokewidth = 1, alpha = 0.25)
    poly!(ax, to_rect.(ts_top, angle.(H_square_top)), strokewidth = 1, alpha = 0.25)
    poly!(ax, to_rect.(ts_left, angle.(H_square_left)), strokewidth = 1, alpha = 0.25)
    save_figures && save("figures/H-arg-principal.pdf", fig)
    fig
end

# ╔═╡ 72c86cec-d5aa-425e-b82e-f01b677e6dbd
let
    fig = Figure(; fontsize)
    ax = Axis(
        fig[1, 1],
        xlabel = L"t",
        ylabel = L"\mathrm{arg}(H(c(t)))",
        yticks = (
            [-π / 2, 0, π / 2, π, 3π / 2, 2π],
            [L"-\frac{\pi}{2}", L"0", L"\frac{\pi}{2}", L"\pi", L"\frac{3\pi}{2}", L"2\pi"],
        ),
    )
    poly!(ax, to_rect.(ts_bottom, angle.(H_square_bottom)), strokewidth = 1, alpha = 0.25)
    poly!(ax, to_rect.(ts_right, angle.(H_square_right)), strokewidth = 1, alpha = 0.25)
    poly!(
        ax,
        to_rect.(ts_top, angle.(-H_square_top) .+ Arb(π)),
        strokewidth = 1,
        alpha = 0.25,
    )
    poly!(
        ax,
        to_rect.(ts_left, angle.(H_square_left) .+ 2Arb(π)),
        strokewidth = 1,
        alpha = 0.25,
    )
    save_figures && save("figures/H-arg-continuous.pdf", fig)
    fig
end

# ╔═╡ Cell order:
# ╟─8e9f1365-8c0e-41cd-b094-663a3760c0d1
# ╠═e914aa7a-31fe-11f1-b1c4-7dea491a45e5
# ╠═baa97573-bc3e-433b-88f1-a7c0a826d27a
# ╟─abc2a340-ab72-4340-838c-bfff029ddd77
# ╠═92972793-3445-428e-9ffb-3d52dc62e050
# ╠═ad65706c-75ec-4ddb-9d1b-571fbe4998d1
# ╠═da71a660-d085-42c6-851a-e7d28f472574
# ╠═06b5e1d2-e5f2-4e1d-ad61-d086dc37dcd2
# ╟─fbbaad41-c446-4c0c-8caa-8137713394b9
# ╟─6a2cd7c8-2760-42a2-a7d1-5dc89550cdd5
# ╠═4279b8df-7e4e-457f-b67c-c9637604ca64
# ╠═0e924af2-bb12-4fe2-b0e5-ad4001c9d06d
# ╟─185f4a22-ad1b-4d27-863d-ff68d3c34809
# ╠═098d9c62-a9c4-4a3e-a1f7-7d466f136f0f
# ╟─93669102-c090-4a86-bfb2-9609fbfe924f
# ╠═adc8d719-f492-4633-907e-15b23e28940d
# ╠═16d83b4f-d023-462a-afa6-4cef99f44196
# ╠═3e607803-f77c-4467-9967-61bbb4f0cedb
# ╟─01f5801b-8605-4c1a-b07d-670abfd451e9
# ╠═baa627cb-3a2c-4a48-89ee-515ceaed7049
# ╠═683c7c63-57df-4bb6-9a3d-96dfbcd85784
# ╠═bbe48cc0-2e1e-46ec-b9b2-4eccff0c4018
# ╟─8b16d819-ae9a-407d-a8fc-61bff5048789
# ╠═f54a47a9-704c-4536-a38f-e935bad69f01
# ╠═2320539d-6309-4826-b72d-f9b9b5459abf
# ╠═c16e015c-d8c6-4684-91c7-09d49c98ccf1
# ╟─3032a818-401b-416a-9d74-1340a98df7db
# ╠═f044e3bf-87b5-4d4f-8bcd-b0ec3f816676
# ╟─5888a3a0-7064-4b6b-9428-42fcb61d396a
# ╠═fd9e8538-381b-44c4-a737-c03f9ad4b590
# ╟─48ddc751-4a96-4a04-8e71-040430d203f6
# ╟─56193256-54c7-425d-a0dc-56227e009aba
# ╠═b4049ea0-ee9f-409c-b2e3-2fdbabe1b8a9
# ╠═c2f067f3-4a78-4048-834d-186ec249350e
# ╟─40f63eb4-5ece-44e7-835a-b33f97e309bd
# ╟─c3ad0d66-d491-42d0-96f4-147350ae6ffa
# ╠═a2feb21f-1cd5-45e2-9f4d-424a23f31788
# ╟─b50ba676-2872-444e-98a2-1cf97f2dbd3a
# ╠═c5aa3ed9-31a4-477e-9066-b8045b2cb829
# ╟─c4fa5393-8708-4643-a95b-9cd398f96364
# ╠═d8d43055-8663-4824-8b92-d6ab45f1cffc
# ╠═4ca4c796-7cf4-447c-9f1a-2cc775e3eb16
# ╟─6c018849-3827-4a45-bd07-cc8e746a44a4
# ╠═75653ff9-94ca-4c7b-a0c5-363b1de657e0
# ╠═3df5f78f-d780-4db7-bd39-a858e04f650a
# ╟─feeea3b3-5de0-42d6-8599-6f19ee3853e3
# ╠═db74a6b9-e8ed-4139-af21-5ea7fc8d55a0
# ╠═7bffa0f1-b6ae-4367-945c-a24841a38f1b
# ╟─91492771-4441-4e1e-b129-38026188a2d2
# ╟─3bf117b4-8f2b-4c88-a1f4-a9b599d74af4
# ╟─880df8b3-36ff-4a74-9f9c-428bf9663462
# ╠═f1ea08ed-8969-448f-859b-8aea8ba546de
# ╠═21cf1d05-df4a-4853-b1be-66bb66e5b35f
# ╟─71b8a3b1-3227-44b9-8f06-50c9192ac96f
# ╠═7d38f474-0f8a-4bc4-8f33-fbb31c540703
# ╟─4ed7e92b-0205-4774-a4a2-42587d385db4
# ╠═2f6b6452-f2e1-4c07-9511-93f355ba9707
# ╟─0d7d7d0a-ae40-41a4-a41d-074cb617a466
# ╠═eb8e5ece-2782-4f2d-b519-93a9f336e5d4
# ╠═ea45a2a5-7e79-44fc-a490-ac4889aa3b1f
# ╠═0f5d4f04-0aee-4582-a504-868816112aa4
# ╠═933ced1e-a9a6-495b-865e-d9a8a3cf7436
# ╟─04b37248-22a5-4fb8-8d0c-a4dbf29a58e4
# ╟─356dc846-b0fc-4743-997a-73d9b9bf7bcf
# ╟─03513e93-3fda-4e9a-b02d-1ed1552811f6
# ╠═f951ad64-9b79-477e-83fb-28fa2bba6e2a
# ╟─93d9e75b-5861-4653-8786-439a35008d00
# ╠═a17f3ef0-9476-4ebd-b52e-657d4fede475
# ╠═6b411c4f-1000-4ebb-84a6-c020e9a33b7c
# ╠═fd79a36b-af81-427e-aefb-a968bc1cb42a
# ╠═ec096027-d898-46e7-abe5-7bcc72229f04
# ╟─f72808a0-5019-49eb-9f8d-41d72b970f0c
# ╠═3e493c1b-b786-4cd3-b321-0de29ae04bf0
# ╠═7d7217b1-5186-4d2f-97ab-f2581add502e
# ╠═c0a0d2fd-9a35-4028-b4e6-dbcf4e777613
# ╠═ca27e948-94c4-44a9-93c8-8b83493eb479
# ╟─6f7407fb-8d41-451a-a220-80bdcc1fe881
# ╠═340f6f00-629a-4c52-8f07-baf6568a2957
# ╠═060bd72c-1e0f-47cf-bbf1-f03c6e8c3599
# ╠═da7f234b-b681-4956-ba49-7f2ac6b363c2
# ╠═942b9ffc-f12a-4c3c-9002-ffb69c784fd1
# ╠═193f2a06-c0d3-45f5-8a2f-ba65d5a41ab3
# ╠═086a8256-b0a7-41a6-9fff-c99c7f953e45
# ╠═7cffb3af-9329-4d5a-b848-6aa1f500b0e7
# ╟─68b90902-028e-4715-b477-2b38bc0d8e4d
# ╠═5ebbff59-5819-467b-bd72-18b29d49e478
# ╠═44ae8114-ba68-4a4d-8e78-9f462091f6cf
# ╠═5b6e9013-1951-4fe1-9a9a-618ac4afcd48
# ╠═57ea897c-5f5f-4adf-af0c-bf48caf8aaeb
# ╟─6af950e3-91a0-4103-aad5-ee45d1d17e14
# ╟─9d0a1908-ea5e-4574-b1a6-3a90c79ed00e
# ╟─f00a20db-3b35-4c0c-8358-e00c98a36385
# ╟─72c86cec-d5aa-425e-b82e-f01b677e6dbd
