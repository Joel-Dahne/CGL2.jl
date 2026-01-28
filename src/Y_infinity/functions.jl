struct FunctionEnclosures_Y
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
        K1, K2 = K_1_2(ξ₁, lambda, κ, ϵ, λ)
        K1_dλ, K2_dλ = K_1_2_dλ(ξ₁, lambda, κ, ϵ, λ)

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

function K_1_2(ξ, lambda, κ, ϵ, λ::CGLParams)
    inv_AM = SMatrix{2,2}(-1 - im * ϵ, -1 + im * ϵ, ϵ - im, ϵ + im) / (2(1 + ϵ^2))
    K_1 = -Diagonal(SVector(J_P_1(ξ, lambda, κ, ϵ, λ), J_P_2(ξ, lambda, κ, ϵ, λ))) * inv_AM
    K_2 = Diagonal(SVector(J_E_1(ξ, lambda, κ, ϵ, λ), J_E_2(ξ, lambda, κ, ϵ, λ))) * inv_AM
    return K_1, K_2
end

function K_1_2_dλ(ξ, lambda, κ, ϵ, λ::CGLParams)
    inv_AM = SMatrix{2,2}(-1 - im * ϵ, -1 + im * ϵ, ϵ - im, ϵ + im) / (2(1 + ϵ^2))
    K_1 =
        -Diagonal(SVector(J_P_1_dλ(ξ, lambda, κ, ϵ, λ), J_P_2_dλ(ξ, lambda, κ, ϵ, λ))) *
        inv_AM

    K_2 =
        Diagonal(SVector(J_E_1_dλ(ξ, lambda, κ, ϵ, λ), J_E_2_dλ(ξ, lambda, κ, ϵ, λ))) *
        inv_AM
    return K_1, K_2
end
