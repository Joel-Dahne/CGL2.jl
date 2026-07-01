module CGL2

using ArbExtras
using Arblib
using LinearAlgebra
using NonlinearSolve
using OhMyThreads: tmap
using OrdinaryDiffEqVerner
using SpecialFunctions
using StaticArrays

import Arpack
import BlockArrays
import DiffEqBase
import ForwardDiff
import IntervalArithmetic: Interval, interval, nai, inf, sup
import ProgressLogging: @progress, @withprogress, @logprogress
import SciMLLogging
import SparseArrays

include("CGLBranch/CGLBranch.jl")

include("assert_proof.jl")

include("arb.jl")
include("helper.jl")
include("det.jl")

include("verify_and_refine_root.jl")

include("U.jl")
include("U_expansion.jl")

include("CGLParams.jl")

include("Q_zero/equation.jl")
include("Q_zero/Q.jl")
include("Q_zero/Q_float.jl")
include("Q_zero/Q_taylor.jl")
include("Q_zero/Q_capd.jl")

include("Q_infinity/parameters.jl")
include("Q_infinity/functions.jl")
include("Q_infinity/function_bounds.jl")
include("Q_infinity/I_bounds.jl")
include("Q_infinity/norm_bounds_constants.jl")
include("Q_infinity/norm_bounds.jl")
include("Q_infinity/I_P.jl")
include("Q_infinity/Q.jl")
include("Q_infinity/p_Q_0.jl")

include("refine_approximation.jl")
include("G.jl")
include("G_solve.jl")

include("Q_hat_zero/Q_hat.jl")
include("Q_hat_zero/Q_hat_float.jl")
include("Q_hat_zero/Q_hat_taylor.jl")
include("Q_hat_zero/Q_hat_capd.jl")

include("Q_hat_infinity/functions.jl")
include("Q_hat_infinity/function_bounds.jl")
include("Q_hat_infinity/norm_bounds.jl")
include("Q_hat_infinity/Q_hat_bounds.jl")
include("Q_hat_infinity/I_E_hat.jl")
include("Q_hat_infinity/Q_hat.jl")

include("G_hat.jl")
include("G_hat_solve.jl")

include("eigenvalues_finite_difference.jl")

include("Y_zero/Y.jl")
include("Y_zero/Y_float.jl")
include("Y_zero/Y_taylor.jl")
include("Y_zero/Y_capd.jl")
include("Y_zero/equation.jl")

include("Y_infinity/functions.jl")
include("Y_infinity/function_bounds.jl")
include("Y_infinity/norm_bounds.jl")
include("Y_infinity/I_K_2.jl")
include("Y_infinity/Y.jl")

include("H.jl")

end # module CGL2
