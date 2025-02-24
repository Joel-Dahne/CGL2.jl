"""
    linearization_eigenvalues_1(ν, κ, ϵ, ξ₁, λ::CGLParams; n = 2048, nev = 10)

Compute the eigenvalues of the linearized equation.

This uses the normalization with `im + ϵ` in front of the second
derivative, which seems to be the "correct" one.

The implementation mimics that from Vladimir's script, though with
some renamed variables.
"""
function linearization_eigenvalues_1(ν, κ, ϵ, ξ₁, λ::CGLParams; n = 2048, nev = 10)
    (; d, ω, σ, δ) = λ

    # Grid we discretize the linear operator on
    ξs = range(zero(ξ₁), ξ₁, n + 2)[2:end-1]
    sol = CGL2.Q_hat_zero_float_curve(real(ν), imag(ν), κ, ϵ, ξ₁, λ, saveat = ξs)
    # Compute complex values at grid points
    Qs = map(y -> complex(y[1], y[2]), sol.u)

    h = step(ξs)
    H1 = 1 / h
    H2 = 1 / h^2

    # A represents the term (1 - im * ϵ) * (q'' + (d - 1) / ξ * q')

    A = SparseArrays.spzeros(ComplexF64, n, n)
    for k = 1:n
        if k > 1
            A[k, k-1] = H2 - (d - 1) / ξs[k] * H1 / 2
        end
        A[k, k] = -2H2
        if k < n
            A[k, k+1] = H2 + (d - 1) / ξs[k] * H1 / 2
        end
    end
    A[1, 1] = -H2 - (d - 1) / ξs[1] * H1 / 2 # Set top left coefficient separately
    A .*= (im + ϵ) # Adjust the coefficient of A

    # B represents the term -im * κ * ξ * q'

    B = SparseArrays.spzeros(ComplexF64, n, n)
    for k = 1:n
        if k > 1
            B[k, k-1] = -ξs[k] * H1 / 2
        end
        if k < n
            B[k, k+1] = ξs[k] * H1 / 2
        end
    end
    B[1, 1] = -ξs[1] * H1 / 2 # Set top left coefficient separately
    B .*= κ # Adjust the coefficient of B

    # C1 represents the "holomorphic part" of the linearization together
    # with ``left-overs" of the constant diagonal linear part of the
    # non-linear term
    # C2 represents the the anti-holomorphic part

    C1 = SparseArrays.spzeros(ComplexF64, n, n)
    C2 = SparseArrays.spzeros(ComplexF64, n, n)
    for k = 1:n
        C1[k, k] = κ / σ + im * ω + (im - δ) * (σ + 1) * (Qs[k] * conj(Qs[k]))^σ
        C2[k, k] = (im - δ) * σ * Qs[k]^2 * (Qs[k] * conj(Qs[k]))^(σ - 1)
    end
    L1 = A + B + C1

    L2 = C2
    M = [real(L1 + L2) imag(-L1 + L2); imag(L1 + L2) real(L1 - L2)]

    # Compute eigenvalues
    return Arpack.eigs(M, which = :SM; nev)
end

"""
    linearization_eigenvalues_2(ν, κ, ϵ, ξ₁, λ::CGLParams; n = 2048, nev = 10)

Compute the eigenvalues of the linearized equation.

This uses the normalization with `im + ϵ` in front of the second
derivative, which seems to be the "correct" one.

The implementation is based on writing all terms as products of
diagonal and tridiagonal matrices.
"""
function linearization_eigenvalues_2(ν, κ, ϵ, ξ₁, λ::CGLParams; n = 2048, nev = 10)
    (; d, ω, σ, δ) = λ

    # Grid we discretize the linear operator on
    ξs = range(zero(ξ₁), ξ₁, n + 2)[2:end-1]
    sol = CGL2.Q_hat_zero_float_curve(real(ν), imag(ν), κ, ϵ, ξ₁, λ, saveat = ξs)
    # Compute complex values at grid points
    Qs = map(y -> complex(y[1], y[2]), sol.u)

    h = step(ξs)

    As = map(Qs, ξs) do Q, ξ
        im + ϵ
    end

    Bs = map(Qs, ξs) do Q, ξ
        (im + ϵ) * (d - 1) / ξ + κ * ξ
    end

    Cs = map(Qs, ξs) do Q, ξ
        κ / σ + im * ω
    end

    N₁s = map(Qs, ξs) do Q, ξ
        (im - δ) * (σ + 1) * (Q * conj(Q))^σ
    end

    N₂s = map(Qs, ξs) do Q, ξ
        (im - δ) * σ * (Q * conj(Q))^(σ - 1) * Q^2
    end

    D1 = 1 / 2h * Tridiagonal(fill(-1.0, n - 1), [-1.0; fill(0.0, n - 1)], fill(1.0, n - 1))
    D2 =
        1 / h^2 * Tridiagonal(fill(1.0, n - 1), [-1.0; fill(-2.0, n - 1)], fill(1.0, n - 1))

    An = Diagonal(As)
    Bn = Diagonal(Bs)
    Cn = Diagonal(Cs)
    N₁n = Diagonal(N₁s)
    N₂n = Diagonal(N₂s)

    L1 = An * D2 + Bn * D1 + Cn + N₁n
    L2 = N₂n

    M = [
        SparseArrays.sparse(real(L1 + L2)) SparseArrays.sparse(imag(-L1 + L2))
        SparseArrays.sparse(imag(L1 + L2)) SparseArrays.sparse(real(L1 - L2))
    ]

    # Compute eigenvalues
    return Arpack.eigs(M, which = :SM; nev)
end

"""
    linearization_eigenvalues_3(ν, κ, ϵ, ξ₁, λ::CGLParams; n = 2048, nev = 10)

Compute the eigenvalues of the linearized equation.

This uses the normalization with `1 - ϵ` in front of the second
derivative, which seems to be the "wrong" one to use.

The implementation mimics [`linearization_eigenvalues_1`](@ref).
"""
function linearization_eigenvalues_3(ν, κ, ϵ, ξ₁, λ::CGLParams; n = 2048, nev = 10)
    (; d, ω, σ, δ) = λ

    # Grid we discretize the linear operator on
    ξs = range(zero(ξ₁), ξ₁, n + 2)[2:end-1]
    sol = CGL2.Q_hat_zero_float_curve(real(ν), imag(ν), κ, ϵ, ξ₁, λ, saveat = ξs)
    # Compute complex values at grid points
    Qs = map(y -> complex(y[1], y[2]), sol.u)

    h = step(ξs)
    H1 = 1 / h
    H2 = 1 / h^2

    # A represents the term (1 - im * ϵ) * (q'' + (d - 1) / ξ * q')

    A = SparseArrays.spzeros(ComplexF64, n, n)
    for k = 1:n
        if k > 1
            A[k, k-1] = H2 - (d - 1) / ξs[k] * H1 / 2
        end
        A[k, k] = -2H2
        if k < n
            A[k, k+1] = H2 + (d - 1) / ξs[k] * H1 / 2
        end
    end
    A .*= (1 - im * ϵ) # Adjust the coefficient of A

    # B represents the term -im * κ * ξ * q'

    B = SparseArrays.spzeros(ComplexF64, n, n)
    for k = 1:n
        if k > 1
            B[k, k-1] = -ξs[k] * H1 / 2
        end
        if k < n
            B[k, k+1] = ξs[k] * H1 / 2
        end
    end
    B .*= -im * κ # Adjust the coefficient of B

    # C1 represents the "holomorphic part" of the linearization together
    # with ``left-overs" of the constant diagonal linear part of the
    # non-linear term
    # C2 represents the the anti-holomorphic part

    C1 = SparseArrays.spzeros(ComplexF64, n, n)
    C2 = SparseArrays.spzeros(ComplexF64, n, n)
    for k = 1:n
        C1[k, k] = -im * κ / σ + ω + (1 + im * δ) * (σ + 1) * (Qs[k] * conj(Qs[k]))^σ
        C2[k, k] = (1 + im * δ) * σ * (Qs[k] * conj(Qs[k]))^(σ - 1) * Qs[k]^2
    end

    L1 = A + B + C1
    L2 = C2
    M = [real(L1 + L2) imag(-L1 + L2); imag(L1 + L2) real(L1 - L2)]

    # Compute eigenvalues
    return Arpack.eigs(M, which = :SM; nev)
end

"""
    linearization_eigenvalues_4(ν, κ, ϵ, ξ₁, λ::CGLParams; n = 2048, nev = 10)

Compute the eigenvalues of the linearized equation.

This uses the normalization with `1 - ϵ` in front of the second
derivative, which seems to be the "wrong" one to use.

The implementation mimics [`linearization_eigenvalues_2`](@ref).
"""
function linearization_eigenvalues_4(ν, κ, ϵ, ξ₁, λ::CGLParams; n = 2048, nev = 10)
    (; d, ω, σ, δ) = λ

    # Grid we discretize the linear operator on
    ξs = range(zero(ξ₁), ξ₁, n + 2)[2:end-1]
    sol = CGL2.Q_hat_zero_float_curve(real(ν), imag(ν), κ, ϵ, ξ₁, λ, saveat = ξs)
    # Compute complex values at grid points
    Qs = map(y -> complex(y[1], y[2]), sol.u)

    h = step(ξs)

    As = map(Qs, ξs) do Q, ξ
        1 - im * ϵ
    end

    Bs = map(Qs, ξs) do Q, ξ
        (1 - im * ϵ) * (d - 1) / ξ - im * κ * ξ
    end

    Cs = map(Qs, ξs) do Q, ξ
        -im * κ / σ + ω
    end

    N₁s = map(Qs, ξs) do Q, ξ
        (1 + im * δ) * (σ + 1) * (Q * conj(Q))^σ
    end

    N₂s = map(Qs, ξs) do Q, ξ
        (1 + im * δ) * σ * (Q * conj(Q))^(σ - 1) * Q^2
    end

    D1 = 1 / 2h * Tridiagonal(fill(-1.0, n - 1), fill(0.0, n), fill(1.0, n - 1))
    D2 = 1 / h^2 * Tridiagonal(fill(1.0, n - 1), fill(-2.0, n), fill(1.0, n - 1))

    An = Diagonal(As)
    Bn = Diagonal(Bs)
    Cn = Diagonal(Cs)
    N₁n = Diagonal(N₁s)
    N₂n = Diagonal(N₂s)

    L1 = An * D2 + Bn * D1 + Cn + N₁n
    L2 = N₂n

    M = [
        SparseArrays.sparse(real(L1 + L2)) SparseArrays.sparse(imag(-L1 + L2))
        SparseArrays.sparse(imag(L1 + L2)) SparseArrays.sparse(real(L1 - L2))
    ]

    # Compute eigenvalues
    return Arpack.eigs(M, which = :SM; nev)
end

"""
    linearization_eigenvalues_real_1(ν, κ, ϵ, ξ₁, λ::CGLParams; n = 2048, nev = 10)

Compute the eigenvalues of the linearized equation.

This uses the normalization with `im + ϵ` in front of the second
derivative, which seems to be the "correct" one.

The implementation uses a formulation of the second order real system
given in terms of block matrices.
"""
function linearization_eigenvalues_real_1(
    ν,
    κ,
    ϵ,
    ξ₁,
    λ::CGLParams;
    n = 2048,
    nev = 10,
    return_L = false,
)
    (; d, ω, σ, δ) = λ

    # Grid we discretize the linear operator on
    ξs = range(zero(ξ₁), ξ₁, n + 2)[2:end-1]
    Qs = CGL2.Q_hat_zero_float_curve(real(ν), imag(ν), κ, ϵ, ξ₁, λ, saveat = ξs).u
    Qs = fill(SVector(0.0, 0.0, 0.0, 0.0), length(ξs))

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

    An = Diagonal(As)
    Bn = Diagonal(Bs)
    Cn = Diagonal(Cs)
    J_Nn = Diagonal(J_Ns)

    L = SparseArrays.sparse(Matrix(BlockArrays.mortar(An * D2 + Bn * D1 + Cn + J_Nn)))
    if return_L
        return L, Arpack.eigs(L, which = :SM; nev)...
    end
    # Compute eigenvalues
    return Arpack.eigs(L, which = :SM; nev)
end

"""
    linearization_eigenvalues_real_2(ν, κ, ϵ, ξ₁, λ::CGLParams; n = 2048, nev = 10)

Compute the eigenvalues of the linearized equation.

This uses the normalization with `1 - ϵ` in front of the second
derivative, which seems to be the "wrong" one to use.

The implementation uses a formulation of the second order real system
given in terms of block matrices.
"""
function linearization_eigenvalues_real_2(ν, κ, ϵ, ξ₁, λ::CGLParams; n = 2048, nev = 10)
    (; d, ω, σ, δ) = λ

    # Grid we discretize the linear operator on
    ξs = range(zero(ξ₁), ξ₁, n + 2)[2:end-1]
    Qs = CGL2.Q_hat_zero_float_curve(real(ν), imag(ν), κ, ϵ, ξ₁, λ, saveat = ξs).u

    h = step(ξs)

    As = map(Qs, ξs) do Q, ξ
        @SMatrix[1 ϵ; -ϵ 1]
    end

    Bs = map(Qs, ξs) do Q, ξ
        @SMatrix[((d-1)/ξ) (ϵ*(d-1)/ξ+κ*ξ); -(ϵ * (d - 1) / ξ + κ * ξ) ((d-1)/ξ)]
    end

    Cs = map(Qs, ξs) do Q, ξ
        @SMatrix[ω κ/σ; -κ/σ ω]
    end

    J_Ns = map(Qs, ξs) do Q, ξ
        a, b, _, _ = Q
        N₁_a = (a^2 + b^2)^(σ - 1) * ((1 + 2σ) * a^2 - 2δ * σ * a * b + b^2)
        N₁_b = -(a^2 + b^2)^(σ - 1) * (δ * a^2 - 2σ * a * b + δ * (1 + 2σ) * b^2)
        N₂_a = (a^2 + b^2)^(σ - 1) * (δ * (1 + 2σ) * a^2 + 2σ * a * b + δ * b^2)
        N₂_b = (a^2 + b^2)^(σ - 1) * (a^2 + 2δ * σ * a * b + (1 + 2σ) * b^2)

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

    An = Diagonal(As)
    Bn = Diagonal(Bs)
    Cn = Diagonal(Cs)
    J_Nn = Diagonal(J_Ns)

    L = SparseArrays.sparse(Matrix(BlockArrays.mortar(An * D2 + Bn * D1 + Cn + J_Nn)))

    # Compute eigenvalues
    return Arpack.eigs(L, which = :SM; nev)
end

"""
    linearization_eigenvalues_real_3(ν, κ, ϵ, ξ₁, λ::CGLParams; n = 2048, nev = 10)

Compute the eigenvalues of the linearized equation.

This uses the normalization with `im + ϵ` in front of the second
derivative, which seems to be the "correct" one.

The implementation uses a formulation of the first order real system
given in terms of block matrices.
"""
function linearization_eigenvalues_real_3(ν, κ, ϵ, ξ₁, λ::CGLParams; n = 2048, nev = 10)
    (; d, ω, σ, δ) = λ

    # Grid we discretize the linear operator on
    ξs = range(zero(ξ₁), ξ₁, n + 2)[2:end-1]
    Qs = CGL2.Q_hat_zero_float_curve(real(ν), imag(ν), κ, ϵ, ξ₁, λ, saveat = ξs).u

    h = step(ξs)

    Bs = map(Qs, ξs) do Q, ξ
        linearization_real_matrix(Q, κ, ϵ, ξ, λ)
    end

    dl = 1 / 2h * fill(-SMatrix{4,4}(I), n - 1)
    du = 1 / 2h * fill(SMatrix{4,4}(I), n - 1)

    L = SparseArrays.sparse(Matrix(BlockArrays.mortar(Tridiagonal(dl, -Bs, du))))

    # Compute eigenvalues
    return Arpack.eigs(L, which = :SM; nev)
end


"""
    linearization_eigenvalues_real_4(ν, κ, ϵ, ξ₁, λ::CGLParams; n = 2048, nev = 10)

Compute the eigenvalues of the linearized equation.

This uses the normalization with `1 - ϵ` in front of the second
derivative, which seems to be the "wrong" one to use.

The implementation uses a formulation of the first order real system
given in terms of block matrices.
"""
function linearization_eigenvalues_real_4(ν, κ, ϵ, ξ₁, λ::CGLParams; n = 2048, nev = 10)
    (; d, ω, σ, δ) = λ

    # Grid we discretize the linear operator on
    ξs = range(zero(ξ₁), ξ₁, n + 2)[2:end-1]
    Qs = CGL2.Q_hat_zero_float_curve(real(ν), imag(ν), κ, ϵ, ξ₁, λ, saveat = ξs).u

    h = step(ξs)

    Bs = map(Qs, ξs) do Q, ξ
        linearization_real_matrix(Q, κ, ϵ, ξ, λ)
    end

    dl = 1 / 2h * fill(-SMatrix{4,4}(I), n - 1)
    du = 1 / 2h * fill(SMatrix{4,4}(I), n - 1)

    L = SparseArrays.sparse(Matrix(BlockArrays.mortar(Tridiagonal(dl, -Bs, du))))

    # Compute eigenvalues
    return Arpack.eigs(L, which = :SM; nev)
end
