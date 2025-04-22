module CGL2

using Arblib
using ArbExtras
using CSV
using DataFrames
using LinearAlgebra
using NonlinearSolve
using OhMyThreads: tmap, tforeach
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqVerner
using SpecialFunctions
using StaticArrays

import Dates
import Distributed
import ForwardDiff
import IntervalArithmetic
import IntervalArithmetic:
    BareInterval, Interval, bareinterval, interval, inf, sup, isempty_interval, nai
import ProgressLogging: @progress, @withprogress, @logprogress
import SparseArrays
import BlockArrays
import Arpack

include("CGLBranch/CGLBranch.jl")

include("arb.jl")
include("interval.jl")
include("helper.jl")
include("special-functions.jl")

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
include("Q_infinity/norm_bounds_constants.jl")
include("Q_infinity/norm_bounds.jl")
include("Q_infinity/I_bounds.jl")
include("Q_infinity/I.jl")
include("Q_infinity/check_existence.jl")
include("Q_infinity/Q.jl")
include("Q_infinity/verify_monotonicity.jl")
include("Q_infinity/p_Q.jl")

include("refine_approximation.jl")
include("G.jl")
include("G_solve.jl")
include("count_critical_points.jl")

include("branch/branch_points.jl")
include("branch/branch_existence.jl")
include("branch/branch_segment_existence.jl")
include("branch/branch_continuation.jl")
include("branch/branch_critical_points.jl")
include("branch/data_handling.jl")
include("branch/proof_witness.jl")

include("orchestration/run_branch_points.jl")
include("orchestration/run_branch_points_verification.jl")
include("orchestration/run_branch_existence.jl")
include("orchestration/run_branch_continuation.jl")
include("orchestration/run_branch_critical_points.jl")

include("Q_hat_zero/equation.jl")
include("Q_hat_zero/Q_hat.jl")
include("Q_hat_zero/Q_hat_float.jl")
include("Q_hat_zero/Q_hat_taylor.jl")
include("Q_hat_zero/Q_hat_capd.jl")

include("Q_hat_infinity/functions.jl")
include("Q_hat_infinity/function_bounds.jl")
include("Q_hat_infinity/check_existence.jl")
include("Q_hat_infinity/Q_hat.jl")
include("Q_hat_infinity/p_Q_hat.jl")

include("G_hat.jl")
include("G_hat_solve.jl")

include("linearization/linearization.jl")
include("linearization/eigenvalues.jl")

include("Y_zero/Y.jl")
include("Y_zero/Y_float.jl")
include("Y_zero/equation.jl")

include("Y_infinity/parameters.jl")
include("Y_infinity/functions.jl")
include("Y_infinity/function_bounds.jl")
#include(Y_infinity/norm_bounds_constants.jl")
include("Y_infinity/norm_bounds.jl")
#include(Y_infinity/I_bounds.jl")
#include(Y_infinity/I.jl")
include("Y_infinity/check_existence.jl")
include("Y_infinity/Y.jl")
#include(Y_infinity/verify_monotonicity.jl")

include("H.jl")
include("H_solve.jl")

include("eigenvalue_solve.jl")

include("branch_eigenvalue/branch_eigenvalue_approximation.jl")
include("branch_eigenvalue/branch_eigenvalue_points.jl")

end # module CGL2
