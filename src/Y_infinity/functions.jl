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

function J_N(Q_hat_ξ, λ::CGLParams{T}) where {T}
    (; σ) = λ
    a, b = reim(Q_hat_ξ)
    (; M1, M2, M3) = J_N_coeff_matrices(λ)

    return -abs2(Q_hat_ξ)^(σ - 1) * (a^2 * M1 + 2σ * a * b * M2 + b^2 * M3)
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
    Y12::SMatrix{2,2,T},
    Y12_dξ::SMatrix{2,2,T},
    Y34::SMatrix{2,2,T},
    Y34_dξ::SMatrix{2,2,T},
    A::SMatrix{2,2},
) where {T}
    m_K2_inv = A * (Y34_dξ - Y12_dξ * inv(Y12) * Y34)

    K1 = inv(Y12) * Y34 * inv(m_K2_inv)
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
    Y12::SMatrix{2,2,T},
    Y12_dξ::SMatrix{2,2,T},
    Y34::SMatrix{2,2,T},
    Y34_dξ::SMatrix{2,2,T},
    Y12_dλ::SMatrix{2,2,T},
    Y12_dλ_dξ::SMatrix{2,2,T},
    Y34_dλ::SMatrix{2,2,T},
    Y34_dλ_dξ::SMatrix{2,2,T},
    A::SMatrix{2,2},
) where {T}
    # This uses the formula d/dx inv(M) = -inv(M) * M_x * inv(M) for
    # differentiating the inverses.
    m_K2_inv = A * (Y34_dξ - Y12_dξ * inv(Y12) * Y34)
    m_K2_inv_dλ =
        A * (
            Y34_dλ_dξ - Y12_dλ_dξ * inv(Y12) * Y34 -
            -Y12_dξ * inv(Y12) * Y12_dλ * inv(Y12) * Y34 - Y12_dξ * inv(Y12) * Y34_dλ
        )

    K1_dλ =
        -inv(Y12) * Y12_dλ * inv(Y12) * Y34 * inv(m_K2_inv) +
        inv(Y12) * Y34_dλ * inv(m_K2_inv) -
        inv(Y12) * Y34 * inv(m_K2_inv) * m_K2_inv_dλ * inv(m_K2_inv)
    K2_dλ = inv(m_K2_inv) * m_K2_inv_dλ * inv(m_K2_inv)

    return K1_dλ, K2_dλ
end
