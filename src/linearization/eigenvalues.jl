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
    for k = 2:n-1
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
    # IMPROVE: Use sparse methods to compute
    return eigvals(Matrix(M))
end
