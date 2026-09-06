### A Pluto.jl notebook ###
# v1.0.3

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

# ╔═╡ 1a17e5de-9185-4a57-aa4c-807d64625983
begin
    using Pkg
    Pkg.activate("..", io = devnull)
    using ArbExtras
    using Arblib
    using CairoMakie
    using CGL2
    using LaTeXStrings
    using OhMyThreads
    using PlutoUI

    # Set the bits of precision used for the computations with Arblib
    setprecision(Arb, 128)

    nothing
end

# ╔═╡ dc87e3d2-a59c-11f1-bfc4-030d33b56d0b
md"""
# Branch of unstable eigenvalues
"""

# ╔═╡ 22b3dba2-7578-4cf8-ac1d-7df3f87705d5
TableOfContents()

# ╔═╡ 9fc24b52-545a-4945-b4e6-716bbe7049cd
md"""
Check this box to set the code to save the figures.
- Save figures $(@bind save_figures CheckBox(default = false))
"""

# ╔═╡ 73588240-32de-4e38-abe7-154426af3b76
fontsize = 22

# ╔═╡ 80f2cabf-bb8c-4e45-bc18-b96488019fcb
Λ = CGLParams{Float64}(3, 1, 1, 0)

# ╔═╡ f3606554-1153-4c49-8c0d-de4410aee607
ξ₁ = 15.0

# ╔═╡ 7866a11b-8507-40fe-8909-de990b04cc47
branch = CGL2.CGLBranch.branch_epsilon(CGL2.CGLBranch.sverak_initial(1, 3)...)

# ╔═╡ 14f1af2e-2d7f-417e-9b98-eba2eae8522b
length(branch)

# ╔═╡ eda95376-87db-4afb-8d1e-9eb52854c28d
indices = 90:10:1900

# ╔═╡ 109f3448-d71d-4690-88c7-7cc705e217f8
@time μγκs = map(indices) do i
    μ₀, κ₀, ϵ = branch.μ[i], branch.κ[i], branch.param[i]
    # Refine backward self-similar solution
    μ, γ, κ = CGL2.refine_approximation_fix_epsilon(μ₀, κ₀, ϵ, ξ₁, Λ)
end

# ╔═╡ a4d66ee2-b2b0-4883-8d4e-eaa6f93ecf6f
@time νγ₂s = tmap(indices, μγκs) do i, (μ, γ, κ)
    ϵ = branch.param[i]
    # Compute forward self-similar solution
    p_Q_0 = ComplexF64(CGL2.p_Q_0(Acb(γ), Arb(κ), Arb(ϵ), Arb(ξ₁), CGLParams{Arb}(Λ)))
    a, b, c = CGL2._abc(κ, ϵ, Λ)
    γ₁ = (-1)^a * p_Q_0
    ν, γ₂ = CGL2.G_hat_approximate(γ₁, κ, ϵ, ξ₁, Λ)
end

# ╔═╡ 5088ce66-01cf-4507-999d-128515a09943
@time all_λs = map(indices, μγκs, νγ₂s) do i, (_, _, κ), (ν, γ₂)
    ϵ = branch.param[i]
    # Find approximation using finite differences
    CGL2.eigenvalues_finite_difference(ν, κ, ϵ, ξ₁, Λ, nev = 40)[1]
end

# ╔═╡ f2f3f28e-f41a-4370-a231-e0ac45cbef62
unstable_λs = map(all_λs) do λs
    maximum(real, filter(λ -> imag(λ) > 0, λs))
end

# ╔═╡ 3471e9ed-3a36-4d28-a43b-307a3f6f404c
let
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"\kappa", ylabel = L"Re(\lambda)")
    lines!(ax, branch.κ[indices], unstable_λs, linewidth = 3)
    save_figures && save("figures/kappa-eigenvalue.pdf", fig)
    fig
end

# ╔═╡ 2469c5ed-0ac1-44dc-860f-f9c8a51a4da2
let
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = L"\epsilon", ylabel = L"\kappa")
    xlims!(ax, 0, nothing)
    ylims!(ax, 0, nothing)

    # Full branch
    lines!(
        ax,
        branch.param,
        branch.κ,
        linewidth = 3,
        linestyle = (:dot, :dense),
        color = :grey,
    )

    lines!(
        ax,
        branch.param[indices],
        branch.κ[indices],
        color = ifelse.(unstable_λs .> 0, :red, :blue),
        linewidth = 3,
    )

    save_figures && save("figures/branch.pdf", fig)
    fig
end

# ╔═╡ Cell order:
# ╠═dc87e3d2-a59c-11f1-bfc4-030d33b56d0b
# ╠═1a17e5de-9185-4a57-aa4c-807d64625983
# ╠═22b3dba2-7578-4cf8-ac1d-7df3f87705d5
# ╟─9fc24b52-545a-4945-b4e6-716bbe7049cd
# ╠═73588240-32de-4e38-abe7-154426af3b76
# ╠═80f2cabf-bb8c-4e45-bc18-b96488019fcb
# ╠═f3606554-1153-4c49-8c0d-de4410aee607
# ╠═7866a11b-8507-40fe-8909-de990b04cc47
# ╠═14f1af2e-2d7f-417e-9b98-eba2eae8522b
# ╠═eda95376-87db-4afb-8d1e-9eb52854c28d
# ╠═109f3448-d71d-4690-88c7-7cc705e217f8
# ╠═a4d66ee2-b2b0-4883-8d4e-eaa6f93ecf6f
# ╠═5088ce66-01cf-4507-999d-128515a09943
# ╠═f2f3f28e-f41a-4370-a231-e0ac45cbef62
# ╠═3471e9ed-3a36-4d28-a43b-307a3f6f404c
# ╠═2469c5ed-0ac1-44dc-860f-f9c8a51a4da2
