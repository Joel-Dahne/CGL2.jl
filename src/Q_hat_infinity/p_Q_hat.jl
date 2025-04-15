function p_Q_hat(γ₁, κ, ϵ, λ::CGLParams{T}) where {T}
    a, b, c = CGL2._abc(κ, ϵ, λ)
    return γ₁ * (-c)^-a
end
