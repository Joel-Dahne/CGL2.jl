using ArbExtras
using Arblib
using CGL2
using LinearAlgebra
using SpecialFunctions
using StaticArrays

using Test
using FiniteDifferences

import CGL2: P, P_dξ, P_dξ_dξ, P_dξ_dξ_dξ
import CGL2: P_dκ, P_dξ_dκ, P_dξ_dξ_dκ
import CGL2: P_dϵ, P_dξ_dϵ, P_dξ_dξ_dϵ
import CGL2: E, E_dξ, E_dξ_dξ, E_dξ_dξ_dξ
import CGL2: E_dκ, E_dξ_dκ
import CGL2: E_dϵ, E_dξ_dϵ
import CGL2: W
import CGL2: J_E, J_E_dξ, J_E_dξ_dξ, J_E_dκ, J_E_dϵ
import CGL2: J_P, J_P_dξ, J_P_dξ_dξ, J_P_dκ, J_P_dϵ
import CGL2: D, D_dξ, D_dξ_dξ, H, H_dξ, H_dξ_dξ

@testset "CGL2" verbose = true begin
    include("arb.jl")
    include("verify_and_refine_root.jl")

    include("U.jl")
    include("U_expansion.jl")

    @testset "Q" begin
        include("Q_zero/equation.jl")
        include("Q_zero/Q.jl")

        include("Q_infinity/parameters.jl")
        include("Q_infinity/functions.jl")
        include("Q_infinity/function_bounds.jl")
        include("Q_infinity/Q.jl")
        include("Q_infinity/p_Q_0.jl")

        include("refine_approximation.jl")
        include("G.jl")
        include("G_solve.jl")
    end

    @testset "Q_hat" begin
        include("Q_hat_zero/Q_hat.jl")

        include("Q_hat_infinity/functions.jl")
        include("Q_hat_infinity/I_E_hat.jl")
    end

    @testset "Y" begin
        include("Y_zero/Y.jl")

        include("Y_infinity/functions.jl")
        include("Y_infinity/function_bounds.jl")
        include("Y_infinity/Y.jl")
    end
end
