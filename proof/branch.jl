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
# Eigenvalue along the branch

This notebook produces the numerical approximations of the eigenvalue with largest real part along the first branch of backward self-similar solutions.

This notebook is not part of any proof, it only generates figures showing the expected behavior along the branch.
"""

# ╔═╡ 22b3dba2-7578-4cf8-ac1d-7df3f87705d5
TableOfContents()

# ╔═╡ 9fc24b52-545a-4945-b4e6-716bbe7049cd
md"""
Check this box to set the code to save the figures.
- Save figures $(@bind save_figures CheckBox(default = false))
"""

# ╔═╡ 73588240-32de-4e38-abe7-154426af3b76
fontsize = 25

# ╔═╡ 80f2cabf-bb8c-4e45-bc18-b96488019fcb
Λ = CGLParams{Float64}(3, 1, 1, 0)

# ╔═╡ f3606554-1153-4c49-8c0d-de4410aee607
ξ₁ = 15.0

# ╔═╡ 84a83e41-1759-4091-9ee6-22719ae4dd2c
md"""
Compute numerical approximation of the branch.
"""

# ╔═╡ 7866a11b-8507-40fe-8909-de990b04cc47
branch = CGL2.CGLBranch.branch_epsilon(CGL2.CGLBranch.sverak_initial(1, 3)...)

# ╔═╡ ce60ce14-f869-43f4-bf76-913a5e196c1d
md"""
We don't compute eigenvalues along the entire branch, only for these indices.
"""

# ╔═╡ eda95376-87db-4afb-8d1e-9eb52854c28d
indices = 90:10:1900

# ╔═╡ 463c1a2e-dc8a-4f18-ade1-b1e105b0b5f2
md"""
Find the index that is closest to the $\epsilon$ used in the computer-assisted proof. Since the branch turns there are two points on it for every $\epsilon$. We are look for the one in the top part.
"""

# ╔═╡ 7c53470f-f1f3-4f2e-bc21-2e198afc6c57
ϵ = 0.1681

# ╔═╡ 8b1246af-18b4-4e9f-87c8-132130be95ab
ϵ_index = findmin(eachindex(indices)) do i
    # Ensure that we are in the top part by checking that ϵ is increasing
    if branch.param[indices[i]] < branch.param[indices[i]+1]
        # Compute distance
        abs(branch.param[indices[i]] - ϵ)
    else
        10.0 # Some large value
    end
end[2]

# ╔═╡ bb3ef1b7-bb82-4fe3-ad80-8d58b82ab4a1
md"""
Compute parameters for the backward self-similar solution along branch.
"""

# ╔═╡ 109f3448-d71d-4690-88c7-7cc705e217f8
@time μγκs = map(indices) do i
    μ₀, κ₀, ϵ = branch.μ[i], branch.κ[i], branch.param[i]
    # Refine backward self-similar solution
    μ, γ, κ = CGL2.refine_approximation_fix_epsilon(μ₀, κ₀, ϵ, ξ₁, Λ)
end

# ╔═╡ 77325ef0-2839-4346-af9d-b9334a903689
md"""
Compute parameters for the forward self-similar solution along branch.
"""

# ╔═╡ a4d66ee2-b2b0-4883-8d4e-eaa6f93ecf6f
@time νγ₂s = tmap(indices, μγκs) do i, (μ, γ, κ)
    ϵ = branch.param[i]
    # Compute forward self-similar solution
    p_Q_0 = ComplexF64(CGL2.p_Q_0(Acb(γ), Arb(κ), Arb(ϵ), Arb(ξ₁), CGLParams{Arb}(Λ)))
    a, b, c = CGL2._abc(κ, ϵ, Λ)
    γ₁ = (-1)^a * p_Q_0
    ν, γ₂ = CGL2.G_hat_approximate(γ₁, κ, ϵ, ξ₁, Λ)
end

# ╔═╡ 1ca879d0-9961-4d6e-a1d2-2979086e1fcc
md"""
Compute eigenvalues along branch.
"""

# ╔═╡ 5088ce66-01cf-4507-999d-128515a09943
@time all_λs = map(indices, μγκs, νγ₂s) do i, (_, _, κ), (ν, γ₂)
    ϵ = branch.param[i]
    # Find approximation using finite differences
    CGL2.eigenvalues_finite_difference(ν, κ, ϵ, ξ₁, Λ, nev = 40)[1]
end

# ╔═╡ 3abd8f47-7b3c-4e80-af15-51952088165f
md"""
Extract the eigenvalues with largest real part for each point along the branch.
"""

# ╔═╡ f2f3f28e-f41a-4370-a231-e0ac45cbef62
unstable_λs = map(all_λs) do λs
    maximum(real, filter(λ -> imag(λ) > 0, λs))
end

# ╔═╡ cdefc68a-5e48-46e7-a2ae-ab661dff72d5
md"""
Plot the branch in $\epsilon$-$\kappa$ space. The part where the real part of the eigenvalue with largest real part is positive is drawn in red and where it is negative is drawn in blue. The part of the branch where we do not compute eigenvalues is dotted.
"""

# ╔═╡ 2469c5ed-0ac1-44dc-860f-f9c8a51a4da2
let
    fig = Figure(; fontsize)
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

    scatter!(ax, branch.param[indices[ϵ_index]], branch.κ[indices[ϵ_index]], color = :black)

    save_figures && save("figures/branch.pdf", fig)
    fig
end

# ╔═╡ 1b543107-9f2a-43c7-98bd-84f29d8f17f0
md"""
Plot the real part of the eigenvalue with largest real part as a function of $\kappa$.
"""

# ╔═╡ 3471e9ed-3a36-4d28-a43b-307a3f6f404c
let
    fig = Figure(; fontsize)
    ax = Axis(fig[1, 1], xlabel = L"\kappa", ylabel = L"\mathrm{Re}(\lambda)")
    lines!(ax, branch.κ[indices], unstable_λs, linewidth = 3)
    scatter!(ax, branch.κ[indices[ϵ_index]], unstable_λs[ϵ_index], color = :black)
    save_figures && save("figures/kappa-eigenvalue.pdf", fig)
    fig
end

# ╔═╡ Cell order:
# ╟─dc87e3d2-a59c-11f1-bfc4-030d33b56d0b
# ╠═1a17e5de-9185-4a57-aa4c-807d64625983
# ╠═22b3dba2-7578-4cf8-ac1d-7df3f87705d5
# ╟─9fc24b52-545a-4945-b4e6-716bbe7049cd
# ╠═73588240-32de-4e38-abe7-154426af3b76
# ╠═80f2cabf-bb8c-4e45-bc18-b96488019fcb
# ╠═f3606554-1153-4c49-8c0d-de4410aee607
# ╟─84a83e41-1759-4091-9ee6-22719ae4dd2c
# ╠═7866a11b-8507-40fe-8909-de990b04cc47
# ╟─ce60ce14-f869-43f4-bf76-913a5e196c1d
# ╠═eda95376-87db-4afb-8d1e-9eb52854c28d
# ╟─463c1a2e-dc8a-4f18-ade1-b1e105b0b5f2
# ╠═7c53470f-f1f3-4f2e-bc21-2e198afc6c57
# ╠═8b1246af-18b4-4e9f-87c8-132130be95ab
# ╟─bb3ef1b7-bb82-4fe3-ad80-8d58b82ab4a1
# ╠═109f3448-d71d-4690-88c7-7cc705e217f8
# ╟─77325ef0-2839-4346-af9d-b9334a903689
# ╠═a4d66ee2-b2b0-4883-8d4e-eaa6f93ecf6f
# ╟─1ca879d0-9961-4d6e-a1d2-2979086e1fcc
# ╠═5088ce66-01cf-4507-999d-128515a09943
# ╟─3abd8f47-7b3c-4e80-af15-51952088165f
# ╠═f2f3f28e-f41a-4370-a231-e0ac45cbef62
# ╟─cdefc68a-5e48-46e7-a2ae-ab661dff72d5
# ╟─2469c5ed-0ac1-44dc-860f-f9c8a51a4da2
# ╟─1b543107-9f2a-43c7-98bd-84f29d8f17f0
# ╟─3471e9ed-3a36-4d28-a43b-307a3f6f404c
