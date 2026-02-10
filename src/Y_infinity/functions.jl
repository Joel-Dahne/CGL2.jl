struct FunctionEnclosures_Y
    E_12::Diagonal{Acb,SVector{2,Acb}}
    E_12_dξ::Diagonal{Acb,SVector{2,Acb}}
    P_12::Diagonal{Acb,SVector{2,Acb}}
    P_12_dξ::Diagonal{Acb,SVector{2,Acb}}
    K_1::Diagonal{Acb,SVector{2,Acb}}
    K_2::Diagonal{Acb,SVector{2,Acb}}
    K_1_dξ::Diagonal{Acb,SVector{2,Acb}}
    K_2_dξ::Diagonal{Acb,SVector{2,Acb}}
    J_N::SMatrix{2,2,Arb}
    J_N_dξ::SMatrix{2,2,Arb}
    I_N::SMatrix{2,2,Acb}
    I_N_dξ::SMatrix{2,2,Acb}

    E_12_dλ::Diagonal{Acb,SVector{2,Acb}}
    E_12_dλ_dξ::Diagonal{Acb,SVector{2,Acb}}
    P_12_dλ::Diagonal{Acb,SVector{2,Acb}}
    P_12_dλ_dξ::Diagonal{Acb,SVector{2,Acb}}
    K_1_dλ::Diagonal{Acb,SVector{2,Acb}}
    K_2_dλ::Diagonal{Acb,SVector{2,Acb}}

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
        K1_dξ, K2_dξ = K_1_2_dξ(ξ₁, lambda, κ, ϵ, λ)
        K1_dλ, K2_dλ = K_1_2_dλ(ξ₁, lambda, κ, ϵ, λ)

        # Compute enclosure of forward solution
        Q_hat_ξ₁, dQ_hat_ξ₁ = Q_hat_infinity(γ₁, γ₂, κ, ϵ, ξ₁, λ)
        JN = J_N(Q_hat_ξ₁, λ)
        JN_dξ = J_N_dξ(Q_hat_ξ₁, dQ_hat_ξ₁, λ)
        IN = I_N(Q_hat_ξ₁, λ)
        IN_dξ = I_N_dξ(Q_hat_ξ₁, dQ_hat_ξ₁, λ)

        F = new(
            Diagonal(SVector(E_1(ξ₁, lambda, κ, ϵ, λ), E_2(ξ₁, lambda, κ, ϵ, λ))),
            Diagonal(SVector(E_1_dξ(ξ₁, lambda, κ, ϵ, λ), E_2_dξ(ξ₁, lambda, κ, ϵ, λ))),
            Diagonal(SVector(P_1(ξ₁, lambda, κ, ϵ, λ), P_2(ξ₁, lambda, κ, ϵ, λ))),
            Diagonal(SVector(P_1_dξ(ξ₁, lambda, κ, ϵ, λ), P_2_dξ(ξ₁, lambda, κ, ϵ, λ))),
            K1,
            K2,
            K1_dξ,
            K2_dξ,
            JN,
            JN_dξ,
            IN,
            IN_dξ,
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

function J_N_dξ(Q_hat_ξ, dQ_hat_ξ, λ::CGLParams{T}) where {T}
    (; σ) = λ
    a, b = reim(Q_hat_ξ)
    a_dξ, b_dξ = reim(dQ_hat_ξ)
    (; M1, M2, M3) = J_N_coeff_matrices(λ)
    # TODO: If we assume that σ is one then the factor
    # abs2(Q_hat_ξ)^(σ - 1) doesn't play a role and the derivative is
    # much simpler. Do we need to care about the general case?
    @assert isone(σ)
    # TODO: Check that this is correct
    return -(2a * a_dξ * M1 + 2σ * (a_dξ * b + a * b_dξ) * M2 + 2b * b_dξ * M3)
end

function I_N(Q_hat_ξ, λ::CGLParams{T}) where {T}
    @assert isone(λ.σ)
    @assert iszero(λ.δ)

    return SMatrix{2,2}(
        2im * abs2(Q_hat_ξ),
        conj(-im * Q_hat_ξ^2),
        -im * Q_hat_ξ^2,
        -2im * abs2(Q_hat_ξ),
    )
end

function I_N_dξ(Q_hat_ξ, dQ_hat_ξ, λ::CGLParams{T}) where {T}
    @assert isone(λ.σ)
    @assert iszero(λ.δ)

    return SMatrix{2,2}(
        4im * (real(Q_hat_ξ) * real(dQ_hat_ξ) + imag(Q_hat_ξ) * imag(dQ_hat_ξ)),
        conj(-2im * Q_hat_ξ * dQ_hat_ξ),
        -2im * Q_hat_ξ * dQ_hat_ξ,
        -4im * (real(Q_hat_ξ) * real(dQ_hat_ξ) + imag(Q_hat_ξ) * imag(dQ_hat_ξ)),
    )
end

function K_1_2(ξ, lambda, κ, ϵ, λ::CGLParams)
    K_1 =
        -Diagonal(
            SVector(
                (ϵ - im) * J_P_1(ξ, lambda, κ, ϵ, λ),
                (ϵ + im) * J_P_2(ξ, lambda, κ, ϵ, λ),
            ),
        ) / (1 + ϵ^2)
    K_2 =
        Diagonal(
            SVector(
                (ϵ - im) * J_E_1(ξ, lambda, κ, ϵ, λ),
                (ϵ + im) * J_E_2(ξ, lambda, κ, ϵ, λ),
            ),
        ) / (1 + ϵ^2)
    return K_1, K_2
end

function K_1_2_dξ(ξ, lambda, κ, ϵ, λ::CGLParams)
    K_1_dξ =
        -Diagonal(
            SVector(
                (ϵ - im) * J_P_1_dξ(ξ, lambda, κ, ϵ, λ),
                (ϵ + im) * J_P_2_dξ(ξ, lambda, κ, ϵ, λ),
            ),
        ) / (1 + ϵ^2)
    K_2_dξ =
        Diagonal(
            SVector(
                (ϵ - im) * J_E_1_dξ(ξ, lambda, κ, ϵ, λ),
                (ϵ + im) * J_E_2_dξ(ξ, lambda, κ, ϵ, λ),
            ),
        ) / (1 + ϵ^2)
    return K_1_dξ, K_2_dξ
end

function K_1_2_dλ(ξ, lambda, κ, ϵ, λ::CGLParams)
    K_1_dλ =
        -Diagonal(
            SVector(
                (ϵ - im) * J_P_1_dλ(ξ, lambda, κ, ϵ, λ),
                (ϵ + im) * J_P_2_dλ(ξ, lambda, κ, ϵ, λ),
            ),
        ) / (1 + ϵ^2)
    K_2_dλ =
        Diagonal(
            SVector(
                (ϵ - im) * J_E_1_dλ(ξ, lambda, κ, ϵ, λ),
                (ϵ + im) * J_E_2_dλ(ξ, lambda, κ, ϵ, λ),
            ),
        ) / (1 + ϵ^2)
    return K_1_dλ, K_2_dλ
end
