struct FunctionEnclosures_Y
    Y_12::SMatrix{2,2,Acb}
    Y_12_dξ::SMatrix{2,2,Acb}
    Y_34::SMatrix{2,2,Acb}
    Y_34_dξ::SMatrix{2,2,Acb}
    K_1::SMatrix{2,2,Acb}
    K_2::SMatrix{2,2,Acb}
    J_N::SMatrix{2,2,Acb}

    Y_12_dλ::SMatrix{2,2,Acb}
    Y_12_dλ_dξ::SMatrix{2,2,Acb}
    Y_34_dλ::SMatrix{2,2,Acb}
    Y_34_dλ_dξ::SMatrix{2,2,Acb}
    K_1_dλ::SMatrix{2,2,Acb}
    K_2_dλ::SMatrix{2,2,Acb}

    function FunctionEnclosures_Y(
        lambda::Acb,
        γ₁::Acb,
        γ₂::Acb,
        κ::Arb,
        ϵ::Arb,
        ξ₁::Arb,
        λ::CGLParams{Arb};
        include_dλ::Bool = false,
    )
        Y_12 = hcat(Y_1(ξ₁, lambda, κ, ϵ, λ), Y_2(ξ₁, lambda, κ, ϵ, λ))
        Y_12_dξ = hcat(Y_1_dξ(ξ₁, lambda, κ, ϵ, λ), Y_2_dξ(ξ₁, lambda, κ, ϵ, λ))
        Y_34 = hcat(Y_3(ξ₁, lambda, κ, ϵ, λ), Y_4(ξ₁, lambda, κ, ϵ, λ))
        Y_34_dξ = hcat(Y_3_dξ(ξ₁, lambda, κ, ϵ, λ), Y_4_dξ(ξ₁, lambda, κ, ϵ, λ))

        A = SMatrix{2,2}(ϵ, 1, -1, ϵ)
        K1, K2 = K_1_2(Y_12, Y_12_dξ, Y_34, Y_34_dξ, A)

        # Compute enclosure of forward solution
        Q_hat_ξ₁, _ = Q_hat_infinity(γ₁, γ₂, κ, ϵ, ξ₁, λ)
        JN = J_N(Q_hat_ξ₁, λ)

        # IMPROVE: Only compute these if include_dλ is true
        Y_12_dλ = hcat(Y_1_dλ(ξ₁, lambda, κ, ϵ, λ), Y_2_dλ(ξ₁, lambda, κ, ϵ, λ))
        Y_12_dλ_dξ = hcat(Y_1_dλ_dξ(ξ₁, lambda, κ, ϵ, λ), Y_2_dλ_dξ(ξ₁, lambda, κ, ϵ, λ))
        Y_34_dλ = hcat(Y_3_dλ(ξ₁, lambda, κ, ϵ, λ), Y_4_dλ(ξ₁, lambda, κ, ϵ, λ))
        Y_34_dλ_dξ = hcat(Y_3_dλ_dξ(ξ₁, lambda, κ, ϵ, λ), Y_4_dλ_dξ(ξ₁, lambda, κ, ϵ, λ))

        K1_dλ, K2_dλ = K_1_2_dλ(
            Y_12,
            Y_12_dξ,
            Y_34,
            Y_34_dξ,
            Y_12_dλ,
            Y_12_dλ_dξ,
            Y_34_dλ,
            Y_34_dλ_dξ,
            A,
        )

        F = new(
            Y_12,
            Y_12_dξ,
            Y_34,
            Y_34_dξ,
            K1,
            K2,
            JN,
            Y_12_dλ,
            Y_12_dλ_dξ,
            Y_34_dλ,
            Y_34_dλ_dξ,
            K1_dλ,
            K2_dλ,
        )

        return F
    end
end

struct FunctionEnclosures_Y_new
    E_12::Diagonal{Acb,SVector{2,Acb}}
    E_12_dξ::Diagonal{Acb,SVector{2,Acb}}
    P_12::Diagonal{Acb,SVector{2,Acb}}
    P_12_dξ::Diagonal{Acb,SVector{2,Acb}}
    K_1::SMatrix{2,2,Acb}
    K_2::SMatrix{2,2,Acb}
    J_N::SMatrix{2,2,Acb}

    E_12_dλ::Diagonal{Acb,SVector{2,Acb}}
    E_12_dλ_dξ::Diagonal{Acb,SVector{2,Acb}}
    P_12_dλ::Diagonal{Acb,SVector{2,Acb}}
    P_12_dλ_dξ::Diagonal{Acb,SVector{2,Acb}}
    K_1_dλ::SMatrix{2,2,Acb}
    K_2_dλ::SMatrix{2,2,Acb}

    function FunctionEnclosures_Y_new(
        lambda::Acb,
        γ₁::Acb,
        γ₂::Acb,
        κ::Arb,
        ϵ::Arb,
        ξ₁::Arb,
        λ::CGLParams{Arb};
        include_dλ::Bool = false,
    )
        K1, K2 = K_1_2_new(ξ₁, lambda, κ, ϵ, λ)
        K1_dλ, K2_dλ = K_1_2_dλ_new(ξ₁, lambda, κ, ϵ, λ)

        # Compute enclosure of forward solution
        Q_hat_ξ₁, _ = Q_hat_infinity(γ₁, γ₂, κ, ϵ, ξ₁, λ)
        JN = J_N(Q_hat_ξ₁, λ)

        F = new(
            Diagonal(SVector(E_1(ξ₁, lambda, κ, ϵ, λ), E_2(ξ₁, lambda, κ, ϵ, λ))),
            Diagonal(SVector(E_1_dξ(ξ₁, lambda, κ, ϵ, λ), E_2_dξ(ξ₁, lambda, κ, ϵ, λ))),
            Diagonal(SVector(P_1(ξ₁, lambda, κ, ϵ, λ), P_2(ξ₁, lambda, κ, ϵ, λ))),
            Diagonal(SVector(P_1_dξ(ξ₁, lambda, κ, ϵ, λ), P_2_dξ(ξ₁, lambda, κ, ϵ, λ))),
            K1,
            K2,
            JN,
            Diagonal(SVector(E_1_dλ(ξ₁, lambda, κ, ϵ, λ), E_2_dλ(ξ₁, lambda, κ, ϵ, λ))),
            Diagonal(
                SVector(E_1_dλ_dξ(ξ₁, lambda, κ, ϵ, λ), E_2_dλ_dξ(ξ₁, lambda, κ, ϵ, λ)),
            ),
            Diagonal(SVector(P_1_dλ(ξ₁, lambda, κ, ϵ, λ), P_2_dλ(ξ₁, lambda, κ, ϵ, λ))),
            Diagonal(
                SVector(P_1_dλ_dξ(ξ₁, lambda, κ, ϵ, λ), P_2_dλ_dξ(ξ₁, lambda, κ, ϵ, λ)),
            ),
            K1_dλ,
            K2_dλ,
        )

        return F
    end
end

function J_N(Q_hat_ξ, λ::CGLParams{T}) where {T}
    (; σ) = λ
    a, b = reim(Q_hat_ξ)
    (; M1, M2, M3) = J_N_coeff_matrices(λ)

    return -abs2(Q_hat_ξ)^(σ - 1) * (a^2 * M1 + 2σ * a * b * M2 + b^2 * M3)
end

function K_1_2_new(ξ, lambda, κ, ϵ, λ::CGLParams)
    inv_AM = SMatrix{2,2}(-1 - im * ϵ, -1 + im * ϵ, ϵ - im, ϵ + im) / (2(1 + ϵ^2))
    K_1 = -Diagonal(SVector(J_P_1(ξ, lambda, κ, ϵ, λ), J_P_2(ξ, lambda, κ, ϵ, λ))) * inv_AM
    K_2 = Diagonal(SVector(J_E_1(ξ, lambda, κ, ϵ, λ), J_E_2(ξ, lambda, κ, ϵ, λ))) * inv_AM
    return K_1, K_2
end

function K_1_2_dλ_new(ξ, lambda, κ, ϵ, λ::CGLParams)
    inv_AM = SMatrix{2,2}(-1 - im * ϵ, -1 + im * ϵ, ϵ - im, ϵ + im) / (2(1 + ϵ^2))
    K_1 =
        -Diagonal(SVector(J_P_1_dλ(ξ, lambda, κ, ϵ, λ), J_P_2_dλ(ξ, lambda, κ, ϵ, λ))) *
        inv_AM

    K_2 =
        Diagonal(SVector(J_E_1_dλ(ξ, lambda, κ, ϵ, λ), J_E_2_dλ(ξ, lambda, κ, ϵ, λ))) *
        inv_AM
    return K_1, K_2
end

K_1_2(ξ, lambda, κ, ϵ, λ::CGLParams) = K_1_2(
    hcat(Y_1(ξ, lambda, κ, ϵ, λ), Y_2(ξ, lambda, κ, ϵ, λ)),
    hcat(Y_1_dξ(ξ, lambda, κ, ϵ, λ), Y_2_dξ(ξ, lambda, κ, ϵ, λ)),
    hcat(Y_3(ξ, lambda, κ, ϵ, λ), Y_4(ξ, lambda, κ, ϵ, λ)),
    hcat(Y_3_dξ(ξ, lambda, κ, ϵ, λ), Y_4_dξ(ξ, lambda, κ, ϵ, λ)),
    SMatrix{2,2}(ϵ, 1, -1, ϵ),
)

# TODO: This uses linear algebra with SMatrix. It is probably rigorous
# for 2x2 matrices, but we should make sure to implement a rigorous
# version.
function K_1_2(
    E12::SMatrix{2,2,T},
    E12_dξ::SMatrix{2,2,T},
    P12::SMatrix{2,2,T},
    P12_dξ::SMatrix{2,2,T},
    A::SMatrix{2,2},
) where {T}
    m_K2_inv = A * (P12_dξ - E12_dξ * inv(E12) * P12)

    K1 = inv(E12) * P12 * inv(m_K2_inv)
    K2 = -inv(m_K2_inv)

    return K1, K2
end

K_1_2_dλ(ξ, lambda, κ, ϵ, λ::CGLParams) = K_1_2_dλ(
    hcat(Y_1(ξ, lambda, κ, ϵ, λ), Y_2(ξ, lambda, κ, ϵ, λ)),
    hcat(Y_1_dξ(ξ, lambda, κ, ϵ, λ), Y_2_dξ(ξ, lambda, κ, ϵ, λ)),
    hcat(Y_3(ξ, lambda, κ, ϵ, λ), Y_4(ξ, lambda, κ, ϵ, λ)),
    hcat(Y_3_dξ(ξ, lambda, κ, ϵ, λ), Y_4_dξ(ξ, lambda, κ, ϵ, λ)),
    hcat(Y_1_dλ(ξ, lambda, κ, ϵ, λ), Y_2_dλ(ξ, lambda, κ, ϵ, λ)),
    hcat(Y_1_dλ_dξ(ξ, lambda, κ, ϵ, λ), Y_2_dλ_dξ(ξ, lambda, κ, ϵ, λ)),
    hcat(Y_3_dλ(ξ, lambda, κ, ϵ, λ), Y_4_dλ(ξ, lambda, κ, ϵ, λ)),
    hcat(Y_3_dλ_dξ(ξ, lambda, κ, ϵ, λ), Y_4_dλ_dξ(ξ, lambda, κ, ϵ, λ)),
    SMatrix{2,2}(ϵ, 1, -1, ϵ),
)

function K_1_2_dλ(
    E12::SMatrix{2,2,T},
    E12_dξ::SMatrix{2,2,T},
    P12::SMatrix{2,2,T},
    P12_dξ::SMatrix{2,2,T},
    E12_dλ::SMatrix{2,2,T},
    E12_dλ_dξ::SMatrix{2,2,T},
    P12_dλ::SMatrix{2,2,T},
    P12_dλ_dξ::SMatrix{2,2,T},
    A::SMatrix{2,2},
) where {T}
    # This uses the formula d/dx inv(M) = -inv(M) * M_x * inv(M) for
    # differentiating the inverses.
    m_K2_inv = A * (P12_dξ - E12_dξ * inv(E12) * P12)
    m_K2_inv_dλ =
        A * (
            P12_dλ_dξ - E12_dλ_dξ * inv(E12) * P12 -
            -E12_dξ * inv(E12) * E12_dλ * inv(E12) * P12 - E12_dξ * inv(E12) * P12_dλ
        )

    K1_dλ =
        -inv(E12) * E12_dλ * inv(E12) * P12 * inv(m_K2_inv) +
        inv(E12) * P12_dλ * inv(m_K2_inv) -
        inv(E12) * P12 * inv(m_K2_inv) * m_K2_inv_dλ * inv(m_K2_inv)
    K2_dλ = inv(m_K2_inv) * m_K2_inv_dλ * inv(m_K2_inv)

    return K1_dλ, K2_dλ
end
