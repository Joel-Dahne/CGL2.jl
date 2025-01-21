"""
    linearization_eigenvalues_FEM(ν, γ, κ, ϵ, ξ₁, λ::CGLParams; n = 128)

Compute the eigenvalues of the linearized equation.
"""
function linearization_eigenvalues_FEM(ν, κ, ϵ, ξ₁, λ::CGLParams; n = 128)
    # Grid we discretize the linear operator on
    ξs = range(zero(ξ₁), ξ₁, n + 2)[2:end-1]

    sol = CGL2.Q_hat_zero_float_curve(real(ν), imag(ν), κ, ϵ, ξ₁, λ, saveat = ξs)
    # Compute complex values at grid points
    qs = map(y -> complex(y[1], y[2]), sol.u)

    h = step(ξs)
    H1 = 1 / h
    H2 = 1 / h^2

    # A represents the radial 3D Laplacian

    A = zeros(n, n) # TODO: Should be sparse
    A[1, 1] = -2H2
    A[1, 2] = 2H2
    A[n, n-1] = H2 - H1 / ξs[n]
    A[n, n] = -2H2
    for k = 2:n-1
        A[k, k+1] = H2 + H1 / ξs[k]
        A[k, k] = -2H2
        A[k, k-1] = H2 - H1 / ξs[k]
    end
    A = (im + ϵ) * A # Adjust the coefficient of A

    # B represents the symmetric difference for the term xv' (TODO: which is what?)

    B = zeros(n, n) # TODO: Should be sparse
    B[1, 2] = ξs[1] * H1 / 2
    B[n, n-1] = -ξs[n] * H1 / 2
    for k = 2:n-1
        B[k, k+1] = ξs[k] * H1 / 2
        B[k, k-1] = -ξs[k] * H1 / 2
    end
    B = κ * B # Adjust the coefficient of B

    # C1 represents the "holomorphic part" of the linearization together
    # with ``left-overs" of the constant diagonal linear part of the
    # non-linear term
    # C2 represents the the anti-holomorphic part

    C1 = zeros(ComplexF64, n, n) # TODO: Should be sparse
    C2 = zeros(ComplexF64, n, n) # TODO: Should be sparse
    for k = 1:n
        C1[k, k] = im * λ.ω + κ / λ.σ + im * (λ.σ + 1) * (qs[k] * conj(qs[k]))^λ.σ
        C2[k, k] = im * λ.σ * qs[k]^2 * (qs[k] * conj(qs[k]))^(λ.σ - 1)
    end

    # Rename things (in a slightly confusing way...)

    A = A + B + C1
    B = C2
    M = [real(A + B) imag(-A + B); imag(A + B) real(A - B)]

    # Compute eigenvalues
    return eigvals(M)
end

"""
    linearization_eigenvalues_FEM_2(ν, γ, κ, ϵ, ξ₁, λ::CGLParams; n = 128)

Compute the eigenvalues of the linearized equation.

Using other form of equation.
"""
function linearization_eigenvalues_FEM_2(ν, κ, ϵ, ξ₁, λ::CGLParams; n = 128)
    (; d, ω, σ, δ) = λ

    # Grid we discretize the linear operator on
    ξs = range(zero(ξ₁), ξ₁, n + 2)[2:end-1]

    sol = CGL2.Q_hat_zero_float_curve(real(ν), imag(ν), κ, ϵ, ξ₁, λ, saveat = ξs)
    # Compute complex values at grid points
    Q_hats = map(y -> complex(y[1], y[2]), sol.u)

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
    B .*= -im * κ

    # C1 represents the terms with q
    # C2 represents the term with conj(q)

    C1 = SparseArrays.spzeros(ComplexF64, n, n)
    C2 = SparseArrays.spzeros(ComplexF64, n, n)
    for k = 1:n
        C1[k, k] =
            -im * κ / σ + ω + (1 + im * δ) * (σ + 1) * (Q_hats[k] * conj(Q_hats[k]))^σ
        C2[k, k] = (1 + im * δ) * σ * (Q_hats[k] * conj(Q_hats[k]))^(σ - 1) * Q_hats[k]^2
    end

    L1 = A + B + C1
    L2 = C2
    M = [real(L1 + L2) imag(-L1 + L2); imag(L1 + L2) real(L1 - L2)]

    # Compute eigenvalues
    return Arpack.eigs(M, nev = 10, which = :SM)[1]
end

"""
    linearization_eigenvalues_FEM_3(ν, γ, κ, ϵ, ξ₁, λ::CGLParams; n = 128)

Compute the eigenvalues of the linearized equation.

Using other form of equation.
"""
function linearization_eigenvalues_FEM_3(ν, κ, ϵ, ξ₁, λ::CGLParams; n = 128, nev = 10)
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
    return Arpack.eigs(M, which = :SM; nev)[1]
end

"""
    linearization_eigenvalues_FEM_real_1(ν, γ, κ, ϵ, ξ₁, λ::CGLParams; n = 128)

Compute the eigenvalues of the linearized equation.

This uses the second order equation for the real form of the linear
operator.
"""
function linearization_eigenvalues_FEM_real_1(ν, κ, ϵ, ξ₁, λ::CGLParams; n = 128)
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

    D1 = BlockArrays.mortar(
        1 / 2h * Tridiagonal(
            fill(-SMatrix{2,2}(I), n - 1),
            fill(@SMatrix(zeros(2, 2)), n),
            fill(SMatrix{2,2}(I), n - 1),
        ),
    )
    D2 = BlockArrays.mortar(
        1 / h^2 * Tridiagonal(
            fill(SMatrix{2,2}(I), n - 1),
            -2fill(SMatrix{2,2}(I), n),
            fill(SMatrix{2,2}(I), n - 1),
        ),
    )

    An = BlockArrays.mortar(Diagonal(As))
    Bn = BlockArrays.mortar(Diagonal(Bs))
    Cn = BlockArrays.mortar(Diagonal(Cs))
    J_Nn = BlockArrays.mortar(Diagonal(J_Ns))

    L = An * D2 + Bn * D1 + Cn + J_Nn

    return eigvals(L)
end


"""
    linearization_eigenvalues_FEM_real_2(ν, γ, κ, ϵ, ξ₁, λ::CGLParams; n = 128)

Compute the eigenvalues of the linearized equation.

This uses the first order equation for the real form of the linear
operator.
"""
function linearization_eigenvalues_FEM_real_2(ν, κ, ϵ, ξ₁, λ::CGLParams; n = 128)
    (; d, ω, σ, δ) = λ

    # Grid we discretize the linear operator on
    ξs = range(zero(ξ₁), ξ₁, n + 2)[2:end-1]
    Qs = CGL2.Q_hat_zero_float_curve(real(ν), imag(ν), κ, ϵ, ξ₁, λ, saveat = ξs).u

    h = step(ξs)

    Bs = map(Qs, ξs) do Q, ξ
        linearization_real_matrix(Q, κ, ϵ, ξ, λ)
    end

    dl = fill(-Matrix(I, 4, 4) / 2h, n - 1)
    du = fill(Matrix(I, 4, 4) / 2h, n - 1)

    L = BlockArrays.mortar(Tridiagonal(dl, -Bs, du))

    return eigvals(L)
end
