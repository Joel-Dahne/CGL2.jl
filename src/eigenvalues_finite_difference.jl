"""
    eigenvalues_finite_difference(ν, κ, ϵ, ξ₁, Λ::CGLParams; n = 2048, nev = 10, return_L = false)

Compute the eigenvalues and eigenvectors of the linearized equation
using a finite difference method.

It uses `n` grid points for the discretization and computes the `nev`
eigenvalues with smallest magnitude. If `return_L` is true it also
returns the matrix given by the discretization as the first argument.

The implementation uses a formulation of the second order real system
given in terms of block matrices.
"""
function eigenvalues_finite_difference(
    ν,
    κ,
    ϵ,
    ξ₁,
    Λ::CGLParams;
    n = 2048,
    nev = 10,
    return_L = false,
)
    (; d, ω, σ, δ) = Λ

    # Grid we discretize the linear operator on
    ξs = range(zero(ξ₁), ξ₁, n + 2)[2:(end-1)]
    Qs = Q_hat_zero_float_curve(real(ν), imag(ν), κ, ϵ, ξ₁, Λ, saveat = ξs).u

    h = step(ξs)

    As = map(Qs, ξs) do Q, ξ
        @SMatrix[ϵ -1; 1 ϵ]
    end

    Bs = map(Qs, ξs) do Q, ξ
        @SMatrix[(ϵ*(d-1)/ξ+κ*ξ) -((d - 1) / ξ); ((d-1)/ξ) (ϵ*(d-1)/ξ+κ*ξ)]
    end

    Cs = map(Qs, ξs) do Q, ξ
        @SMatrix[κ/σ -ω; ω κ/σ]
    end

    J_Ns = map(Qs, ξs) do Q, ξ
        a, b, _, _ = Q
        N₁_a = -(a^2 + b^2)^(σ - 1) * (δ * (1 + 2σ) * a^2 + 2σ * a * b + δ * b^2)
        N₁_b = -(a^2 + b^2)^(σ - 1) * (a^2 + 2δ * σ * a * b + (1 + 2σ) * b^2)
        N₂_a = (a^2 + b^2)^(σ - 1) * ((1 + 2σ) * a^2 - 2δ * σ * a * b + b^2)
        N₂_b = -(a^2 + b^2)^(σ - 1) * (δ * a^2 - 2σ * a * b + δ * (1 + 2σ) * b^2)

        @SMatrix[N₁_a N₁_b; N₂_a N₂_b]
    end

    D1 =
        1 / 2h * Tridiagonal(
            fill(-SMatrix{2,2}(I), n - 1),
            fill(@SMatrix(zeros(2, 2)), n),
            fill(SMatrix{2,2}(I), n - 1),
        )
    D2 =
        1 / h^2 * Tridiagonal(
            fill(SMatrix{2,2}(I), n - 1),
            -2fill(SMatrix{2,2}(I), n),
            fill(SMatrix{2,2}(I), n - 1),
        )

    # Adjust top left corner of matrices for Neumann boundary conditions
    D1[1, 1] = -SMatrix{2,2}(I) / 2h
    D2[1, 1] /= 2

    An = Diagonal(As)
    Bn = Diagonal(Bs)
    Cn = Diagonal(Cs)
    J_Nn = Diagonal(J_Ns)

    L = SparseArrays.sparse(Matrix(BlockArrays.mortar(An * D2 + Bn * D1 + Cn + J_Nn)) / 2κ)

    # Compute eigenvalues
    v0 = ones(size(L, 1)) # Fix v0 to give reproducible results
    Λs, vs = Arpack.eigs(L, which = :SM; nev, v0)

    if return_L
        return L, Λs, vs
    end
    return Λs, vs
end
