using ArbExtras
using Arblib
using CGL2
using LinearAlgebra
using SpecialFunctions
using StaticArrays

using Test
using FiniteDifferences

@testset "CGL2" verbose = true begin
    include("arb.jl")
    include("verify_and_refine_root.jl")
    include("special-functions.jl")

    include("U.jl")
    include("U_expansion.jl")

    @testset "Q" begin
        include("Q_zero/equation.jl")
        include("Q_zero/Q.jl")

        include("Q_infinity/parameters.jl")
        include("Q_infinity/functions.jl")
        include("Q_infinity/function_bounds.jl")
        include("Q_infinity/Q.jl")

        include("refine_approximation.jl")
        include("G.jl")
        include("G_solve.jl")
    end

    @testset "Q_hat" begin
        include("Q_hat_zero/Q_hat.jl")
    end

    @testset "Y" begin
        include("Y_zero/Y.jl")

        #include("Y_infinity/parameters.jl") # TODO: Update these tests
        include("Y_infinity/functions.jl")
        include("Y_infinity/function_bounds.jl")
        include("Y_infinity/Y.jl")
    end
end
