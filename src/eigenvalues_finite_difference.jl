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

    A = @SMatrix[ϵ -1; 1 ϵ]
    B₁ = @SMatrix[κ 0; 0 κ]
    B₂ = (d - 1) * A
    C = @SMatrix[κ/σ -ω; ω κ/σ]

    An = Diagonal(fill(A, n))
    B₁n = Diagonal(fill(B₁, n))
    B₂n = Diagonal(fill(B₂, n))
    Cn = Diagonal(fill(C, n))

    Xn = Diagonal([@SMatrix[ξ 0; 0 ξ] for ξ in ξs])

    J_Nn =
        map(Qs) do Q
            a, b, _, _ = Q
            J_N_11 = -(a^2 + b^2)^(σ - 1) * (δ * (1 + 2σ) * a^2 + 2σ * a * b + δ * b^2)
            J_N_12 = -(a^2 + b^2)^(σ - 1) * (a^2 + 2δ * σ * a * b + (1 + 2σ) * b^2)
            J_N_21 = (a^2 + b^2)^(σ - 1) * ((1 + 2σ) * a^2 - 2δ * σ * a * b + b^2)
            J_N_22 = -(a^2 + b^2)^(σ - 1) * (δ * a^2 - 2σ * a * b + δ * (1 + 2σ) * b^2)
            @SMatrix[J_N_11 J_N_12; J_N_21 J_N_22]
        end |> Diagonal

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

    L = SparseArrays.sparse(
        Matrix(BlockArrays.mortar(An * D2 + (B₁n * Xn + B₂n * inv(Xn)) * D1 + Cn + J_Nn)) / 2κ,
    )

    # Compute eigenvalues
    v0 = ones(size(L, 1)) # Fix v0 to give reproducible results
    λs, vs = Arpack.eigs(L, which = :SM; nev, v0)

    if return_L
        return L, λs, vs
    end
    return λs, vs
end
