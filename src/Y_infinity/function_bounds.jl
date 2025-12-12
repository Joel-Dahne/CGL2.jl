struct FunctionBounds_Y
    J_N::Arb
    Y_12::Arb
    Y_34::Arb
    Y_12_dξ::Arb
    Y_34_dξ::Arb
    K_1::Arb
    K_2::Arb
    Y_12_dλ::Arb
    Y_34_dλ::Arb
    Y_12_dλ_dξ::Arb
    Y_34_dλ_dξ::Arb
    K_1_dλ::Arb
    K_2_dλ::Arb

    function FunctionBounds_Y(
        lambda::Acb,
        γ₁::Acb,
        γ₂::Acb,
        κ::Arb,
        ϵ::Arb,
        ξ₁::Arb,
        λ::CGLParams{Arb};
        include_dλ::Bool = false,
    )
        C = new(
            C_J_N(lambda, γ₁, γ₂, κ, ϵ, ξ₁, λ),
            C_Y_12(lambda, κ, ϵ, ξ₁, λ),
            C_Y_34(lambda, κ, ϵ, ξ₁, λ),
            C_Y_12_dξ(lambda, κ, ϵ, ξ₁, λ),
            C_Y_34_dξ(lambda, κ, ϵ, ξ₁, λ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
        )

        C.K_1[] = C_K_1(lambda, κ, ϵ, ξ₁, λ, C)
        C.K_2[] = C_K_2(lambda, κ, ϵ, ξ₁, λ, C)

        if include_dλ
            C.Y_12_dλ[] = C_Y_12_dλ(lambda, κ, ϵ, ξ₁, λ, C)
            C.Y_34_dλ[] = C_Y_34_dλ(lambda, κ, ϵ, ξ₁, λ, C)
            C.Y_12_dλ_dξ[] = C_Y_12_dλ_dξ(lambda, κ, ϵ, ξ₁, λ, C)
            C.Y_34_dλ_dξ[] = C_Y_34_dλ_dξ(lambda, κ, ϵ, ξ₁, λ, C)
            C.K_1_dλ[] = C_K_1_dλ(lambda, κ, ϵ, ξ₁, λ, C)
            C.K_2_dλ[] = C_K_2_dλ(lambda, κ, ϵ, ξ₁, λ, C)
        end

        return C
    end
end

# TODO: Decide if we want to split the bounds in _1 and _2 versions or
# not.
struct FunctionBounds_Y_new
    J_N::Arb
    E_1::Arb
    E_2::Arb
    E_12::Arb
    P_1::Arb
    P_2::Arb
    P_12::Arb
    E_1_dξ::Arb
    E_2_dξ::Arb
    E_12_dξ::Arb
    P_1_dξ::Arb
    P_2_dξ::Arb
    P_12_dξ::Arb
    K_1::Arb
    K_2::Arb
    E_1_dλ::Arb
    E_2_dλ::Arb
    E_12_dλ::Arb
    P_1_dλ::Arb
    P_2_dλ::Arb
    P_12_dλ::Arb
    E_1_dλ_dξ::Arb
    E_2_dλ_dξ::Arb
    E_12_dλ_dξ::Arb
    P_1_dλ_dξ::Arb
    P_2_dλ_dξ::Arb
    P_12_dλ_dξ::Arb
    K_1_dλ::Arb
    K_2_dλ::Arb

    # TODO: Update the bounds
    function FunctionBounds_Y_new(
        lambda::Acb,
        γ₁::Acb,
        γ₂::Acb,
        κ::Arb,
        ϵ::Arb,
        ξ₁::Arb,
        λ::CGLParams{Arb};
        include_dλ::Bool = false,
    )
        a, b, c = _abc(κ, ϵ, λ)
        CU = UBounds(a - lambda / 2κ, b, -c, ξ₁, include_da = include_dλ)
        CU_conj = UBounds(conj(a) - lambda / 2κ, b, -conj(c), ξ₁, include_da = include_dλ)

        C = new(
            C_J_N(lambda, γ₁, γ₂, κ, ϵ, ξ₁, λ),
            C_E_1(lambda, κ, ϵ, ξ₁, λ, CU),
            C_E_2(lambda, κ, ϵ, ξ₁, λ, CU_conj),
            indeterminate(κ),
            C_P_1(lambda, κ, ϵ, ξ₁, λ, CU),
            C_P_2(lambda, κ, ϵ, ξ₁, λ, CU_conj),
            indeterminate(κ),
            C_E_1_dξ(lambda, κ, ϵ, ξ₁, λ, CU),
            C_E_2_dξ(lambda, κ, ϵ, ξ₁, λ, CU_conj),
            indeterminate(κ),
            C_P_1_dξ(lambda, κ, ϵ, ξ₁, λ, CU),
            C_P_2_dξ(lambda, κ, ϵ, ξ₁, λ, CU_conj),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
            indeterminate(κ),
        )

        C.E_12[] = max(C.E_1, C.E_2)
        C.P_12[] = max(C.P_1, C.P_2)
        C.E_12_dξ[] = max(C.E_1_dξ, C.E_2_dξ)
        C.P_12_dξ[] = max(C.P_1_dξ, C.P_2_dξ)
        C.K_1[] = C_K_1(lambda, κ, ϵ, ξ₁, λ, C)
        C.K_2[] = C_K_2(lambda, κ, ϵ, ξ₁, λ, C)

        if include_dλ
            C.E_1_dλ[] = C_E_1_dλ(lambda, κ, ϵ, ξ₁, λ, CU)
            C.E_2_dλ[] = C_E_2_dλ(lambda, κ, ϵ, ξ₁, λ, CU_conj)
            C.E_12_dλ[] = max(C.E_1_dλ, C.E_2_dλ)
            C.P_1_dλ[] = C_P_1_dλ(lambda, κ, ϵ, ξ₁, λ, CU)
            C.P_2_dλ[] = C_P_2_dλ(lambda, κ, ϵ, ξ₁, λ, CU_conj)
            C.P_12_dλ[] = max(C.P_1_dλ, C.P_2_dλ)
            C.E_1_dλ_dξ[] = C_E_1_dλ_dξ(lambda, κ, ϵ, ξ₁, λ, CU)
            C.E_2_dλ_dξ[] = C_E_2_dλ_dξ(lambda, κ, ϵ, ξ₁, λ, CU_conj)
            C.E_12_dλ_dξ[] = max(C.E_1_dλ_dξ, C.E_2_dλ_dξ)
            C.P_1_dλ_dξ[] = C_P_1_dλ_dξ(lambda, κ, ϵ, ξ₁, λ, CU)
            C.P_2_dλ_dξ[] = C_P_2_dλ_dξ(lambda, κ, ϵ, ξ₁, λ, CU_conj)
            C.P_12_dλ_dξ[] = max(C.P_1_dλ_dξ, C.P_2_dλ_dξ)
            C.K_1_dλ[] = C_K_1_dλ(lambda, κ, ϵ, ξ₁, λ, C)
            C.K_2_dλ[] = C_K_2_dλ(lambda, κ, ϵ, ξ₁, λ, C)
        end

        return C
    end
end

function C_J_N(lambda::Acb, γ₁::Acb, γ₂::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    (; σ, δ) = λ

    # FIXME
    Q_hat, _ = Q_hat_infinity(γ₁, γ₂, κ, ϵ, ξ₁, λ)
    a, b = reim(Q_hat)
    C_ab = 1.1max(abs(a), abs(b)) / ξ₁^(-1 / σ)

    return 2^(σ - 1) *
           C_ab^2σ *
           (
               1 +
               abs(δ) * (1 + 2σ) +
               2σ * (1 + abs(δ)) +
               max(1 + 2σ * abs(δ), 1 + abs(δ) * (1 + 2σ))
           )
end

# FIXME
function C_Y_12(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    return 1.5max(norm_inf(Y_1(ξ₁, lambda, κ, ϵ, λ)), norm_inf(Y_2(ξ₁, lambda, κ, ϵ, λ))) /
           (exp(real_a12(κ, ϵ) * ξ₁^2) * ξ₁^-real_s12(lambda, κ, λ))
end

# FIXME
function C_Y_34(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    return 1.5max(norm_inf(Y_3(ξ₁, lambda, κ, ϵ, λ)), norm_inf(Y_4(ξ₁, lambda, κ, ϵ, λ))) /
           ξ₁^-real_s34(lambda, κ, λ)
end

# FIXME
function C_Y_12_dξ(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    return 1.5max(
        norm_inf(Y_1_dξ(ξ₁, lambda, κ, ϵ, λ)),
        norm_inf(Y_2_dξ(ξ₁, lambda, κ, ϵ, λ)),
    ) / (exp(real_a12(κ, ϵ) * ξ₁^2) * ξ₁^(-real_s12(lambda, κ, λ) + 1))
end

# FIXME
function C_Y_34_dξ(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb})
    return 1.5max(
        norm_inf(Y_3_dξ(ξ₁, lambda, κ, ϵ, λ)),
        norm_inf(Y_4_dξ(ξ₁, lambda, κ, ϵ, λ)),
    ) / ξ₁^(-real_s34(lambda, κ, λ) - 1)
end

# FIXME
function C_Y_12_dλ(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    return 1.5max(
        norm_inf(Y_1_dλ(ξ₁, lambda, κ, ϵ, λ)),
        norm_inf(Y_2_dλ(ξ₁, lambda, κ, ϵ, λ)),
    ) / (exp(real_a12(κ, ϵ) * ξ₁^2) * log(ξ₁) * ξ₁^-real_s12(lambda, κ, λ))
end

function C_Y_34_dλ(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    return 1.5max(
        norm_inf(Y_3_dλ(ξ₁, lambda, κ, ϵ, λ)),
        norm_inf(Y_4_dλ(ξ₁, lambda, κ, ϵ, λ)),
    ) / (log(ξ₁) * ξ₁^-real_s34(lambda, κ, λ))
end

# FIXME
function C_Y_12_dλ_dξ(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    return 1.5max(
        norm_inf(Y_1_dλ_dξ(ξ₁, lambda, κ, ϵ, λ)),
        norm_inf(Y_2_dλ_dξ(ξ₁, lambda, κ, ϵ, λ)),
    ) / (exp(real_a12(κ, ϵ) * ξ₁^2) * log(ξ₁) * ξ₁^(-real_s12(lambda, κ, λ) + 1))
end

# FIXME
function C_Y_34_dλ_dξ(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    return 1.5max(
        norm_inf(Y_3_dλ_dξ(ξ₁, lambda, κ, ϵ, λ)),
        norm_inf(Y_4_dλ_dξ(ξ₁, lambda, κ, ϵ, λ)),
    ) / (log(ξ₁) * ξ₁^(-real_s34(lambda, κ, λ) - 1))
end

# TODO: Update this
function C_E_1(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU::UBounds)
    (; d, σ) = λ
    a, b, c = _abc(κ, ϵ, λ)

    return CU.U_bma_b * abs(c^((a - lambda / 2κ) - b))
end

function C_E_2(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU_conj::UBounds)
    (; d, σ) = λ
    a, b, c = _abc(κ, ϵ, λ)

    return CU_conj.U_bma_b * abs(conj(c)^((conj(a) - lambda / 2κ) - b))
end

# PROVE: Verify this
function C_P_1(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU::UBounds)
    (; d, σ) = λ
    a, b, c = _abc(κ, ϵ, λ)

    return CU.U_a_b * abs((-c)^(lambda / 2κ - a))
end

# PROVE: Verify this
function C_P_2(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU_conj::UBounds)
    (; d, σ) = λ
    a, b, c = _abc(κ, ϵ, λ)

    return CU_conj.U_a_b * abs(conj(-c)^(lambda / 2κ - conj(a)))
end

# PROVE: Verify this
function C_E_1_dξ(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU::UBounds)
    a, b, c = _abc(κ, ϵ, λ)

    C1 = abs(c^((a - lambda / 2κ) - b)) * CU.U_bma_b
    C2 = abs(c^((a - lambda / 2κ) - b - 1)) * CU.U_dz_bma_b

    return abs(2c) * C1 + abs(2c) * C2 * ξ₁^-2
end

function C_E_2_dξ(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU_conj::UBounds)
    a, b, c = _abc(κ, ϵ, λ)

    C1 = abs(conj(c)^((conj(a) - lambda / 2κ) - b)) * CU_conj.U_bma_b
    C2 = abs(conj(c)^((conj(a) - lambda / 2κ) - b - 1)) * CU_conj.U_dz_bma_b

    return abs(2c) * C1 + abs(2c) * C2 * ξ₁^-2
end

# PROVE: Verify this
function C_P_1_dξ(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU::UBounds)
    a, b, c = _abc(κ, ϵ, λ)

    return abs(2(-c)^(lambda / 2κ - a)) * CU.U_dz_a_b
end

# PROVE: Verify this
function C_P_2_dξ(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU_conj::UBounds)
    a, b, c = _abc(κ, ϵ, λ)

    return abs(2conj(-c)^(lambda / 2κ - conj(a))) * CU_conj.U_dz_a_b
end

function C_E_1_dλ(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU::UBounds)
    a, b, c = _abc(κ, ϵ, λ)

    return CU.U_da_bma_b * abs(c^((a - lambda / 2κ) - b)) / 2κ * (2 + abs(log(c)) / log(ξ₁))
end

function C_E_2_dλ(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU_conj::UBounds)
    a, b, c = _abc(κ, ϵ, λ)

    return CU_conj.U_da_bma_b * abs(conj(c)^((conj(a) - lambda / 2κ) - b)) / 2κ *
           (2 + abs(log(conj(c))) / log(ξ₁))
end

function C_P_1_dλ(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU::UBounds)
    a, b, c = _abc(κ, ϵ, λ)

    return CU.U_da_a_b * abs((-c)^-(a - lambda / 2κ)) / 2κ * (2 + abs(log(-c)) / log(ξ₁))
end

function C_P_2_dλ(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU_conj::UBounds)
    a, b, c = _abc(κ, ϵ, λ)

    return CU_conj.U_da_a_b * abs(conj(-c)^-(conj(a) - lambda / 2κ)) / 2κ *
           (2 + abs(log(conj(-c))) / log(ξ₁))
end

# PROVE: Verify this
function C_E_1_dλ_dξ(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU::UBounds)
    a, b, c = _abc(κ, ϵ, λ)

    #C1 = abs((-c)^((a - lambda / 2κ) - b)) * CU.U_da_bma_b
    #C2 = abs((-c)^((a - lambda / 2κ) - b - 1)) * CU.U_dzda_bma_b

    return Arb(1) # FIXME
    #return (abs(2c) * C1 + abs(2c) * C2 * ξ₁^-2) / 2κ
end

function C_E_2_dλ_dξ(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    CU_conj::UBounds,
)
    a, b, c = _abc(κ, ϵ, λ)

    #C1 = abs((-c)^((a - lambda / 2κ) - b)) * CU.U_da_bma_b
    #C2 = abs((-c)^((a - lambda / 2κ) - b - 1)) * CU.U_dzda_bma_b

    return Arb(1) # FIXME
    #return (abs(2c) * C1 + abs(2c) * C2 * ξ₁^-2) / 2κ
end

# FIXME
function C_P_1_dλ_dξ(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, CU::UBounds)
    a, b, c = _abc(κ, ϵ, λ)

    return Arb(1) # FIXME
    #return abs(2c^-(a - lambda / 2κ)) * CU.U_dz_a_b
end

# FIXME
function C_P_2_dλ_dξ(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    CU_conj::UBounds,
)
    a, b, c = _abc(κ, ϵ, λ)

    return Arb(1) # FIXME
    #return abs(2c^-(a - lambda / 2κ)) * CU.U_dz_a_b
end

# TODO: Remove
function C_K_1(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, C::FunctionBounds_Y)
    K_1, _ = K_1_2(ξ₁, lambda, κ, ϵ, λ)
    return 1.5norm_inf(K_1) /
           (exp(-real_a12(κ, ϵ) * ξ₁^2) * ξ₁^(real_s12(lambda, κ, λ) - 1))
end

# FIXME
function C_K_1(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y_new,
)
    K_1, _ = K_1_2_new(ξ₁, lambda, κ, ϵ, λ)
    return 1.5norm_inf(K_1) /
           (exp(-real_a12(κ, ϵ) * ξ₁^2) * ξ₁^(real_s12(lambda, κ, λ) - 1))
end

# FIXME
function C_K_2(lambda::Acb, κ::Arb, ϵ::Arb, ξ₁::Arb, λ::CGLParams{Arb}, C::FunctionBounds_Y)
    _, K_2 = K_1_2(ξ₁, lambda, κ, ϵ, λ)
    return 1.5norm_inf(K_2) / ξ₁^(real_s34(lambda, κ, λ) - 1)
end

# FIXME
function C_K_2(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y_new,
)
    _, K_2 = K_1_2_new(ξ₁, lambda, κ, ϵ, λ)
    return 1.5norm_inf(K_2) / ξ₁^(real_s34(lambda, κ, λ) - 1)
end

# FIXME
function C_K_1_dλ(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    K_1_dλ, _ = K_1_2_dλ(ξ₁, lambda, κ, ϵ, λ)
    return 1.5norm_inf(K_1_dλ) /
           (exp(-real_a12(κ, ϵ) * ξ₁^2) * log(ξ₁) * ξ₁^(real_s12(lambda, κ, λ) - 1))
end

# FIXME
function C_K_1_dλ(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y_new,
)
    K_1_dλ, _ = K_1_2_dλ_new(ξ₁, lambda, κ, ϵ, λ)
    return 1.5norm_inf(K_1_dλ) /
           (exp(-real_a12(κ, ϵ) * ξ₁^2) * log(ξ₁) * ξ₁^(real_s12(lambda, κ, λ) - 1))
end

# FIXME
function C_K_2_dλ(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y,
)
    _, K_2_dλ = K_1_2_dλ(ξ₁, lambda, κ, ϵ, λ)
    return 1.5norm_inf(K_2_dλ) / (log(ξ₁) * ξ₁^(real_s34(lambda, κ, λ) - 1))
end

# FIXME
function C_K_2_dλ(
    lambda::Acb,
    κ::Arb,
    ϵ::Arb,
    ξ₁::Arb,
    λ::CGLParams{Arb},
    C::FunctionBounds_Y_new,
)
    _, K_2_dλ = K_1_2_dλ_new(ξ₁, lambda, κ, ϵ, λ)
    return 1.5norm_inf(K_2_dλ) / (log(ξ₁) * ξ₁^(real_s34(lambda, κ, λ) - 1))
end
