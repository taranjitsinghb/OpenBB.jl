# @Author: Massimo De Mauri <massimo>
# @Date:   2019-03-04T03:18:32+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: GUROBI_types.jl
# @Last modified by:   massimo
# @Last modified time: 2019-06-12T22:34:03+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

const GUROBImaxint = 2000000000


mutable struct GUROBIsettings <:AbstractSettings
    # algorithm and IO
    Method::Int                    #1->simplex, 2->barrier
    OutputFlag::Bool               #Enables or disables screen output
    DisplayInterval::Int           #Frequency at which log lines are printed
    # termination
    IterationLimit::Float64        # max number of iterations allowed
    TimeLimit::Float64             # max time allowed
    # tolerances
    FeasibilityTol::Float64 	   #Primal feasibility tolerance
    OptimalityTol::Float64 	       #Dual feasibility tolerance
    PSDTol::Float64 	           #Positive semi-definite tolerance
    # simplex algorithm
    MarkowitzTol::Float64	       #Threshold pivoting tolerance
    InfUnbdInfo::Int 	           #Generate additional info for infeasible/unbounded models
    NormAdjust::Int  	           #Simplex pricing norm
    ObjScale::Float64 	           #Objective scaling
    PerturbValue::Float64 	       #Simplex perturbation magnitude
    Quad::Int 	                   #Quad precision computation in simplex
    ScaleFlag::Int 	               #Model scaling
    Sifting::Int 	               #Sifting within dual simplex
    SiftMethod::Int 	           #LP method used to solve sifting sub-problems
    SimplexPricing::Int 	       #Simplex variable pricing strategy
    # barrier algorithm
    BarIterLimit::Int 	           #Barrier iteration limit
    BarConvTol::Float64 	       #Barrier convergence tolerance
    BarQCPConvTol::Float64 	       #Barrier QCP convergence tolerance
    BarCorrectors::Int 	           #Central correction limit
    BarHomogeneous::Int 	       #Barrier homogeneous algorithm
    BarOrder::Int 	               #Barrier ordering algorithm
    Crossover::Int  	           #Barrier crossover strategy
    CrossoverBasis::Int 	       #Crossover initial basis construction strategy
    QCPDual::Int 	               #Compute dual variables for QCP models
    # thead count
    Threads::Int                   # number of threads to use
end


function GUROBIsettings(;   Method::Int=2,
                            OutputFlag::Bool=false,
                            DisplayInterval::Int=10,
                            IterationLimit::Float64=Inf,
                            TimeLimit::Float64=Inf,
                            FeasibilityTol::Float64=1e-6,
                            OptimalityTol::Float64=1e-6,
                            PSDTol::Float64=1e-6,
                            MarkowitzTol::Float64=0.0078125,
                            InfUnbdInfo::Int=0,
                            NormAdjust::Int=-1,
                            ObjScale::Float64=0.0,
                            PerturbValue::Float64=0.0002,
                            Quad::Int=-1,
                            ScaleFlag::Int=-1,
                            Sifting::Int=-1,
                            SiftMethod::Int=-1,
                            SimplexPricing::Int=-1,
                            BarIterLimit::Int=1000,
                            BarConvTol::Float64=1e-8,
                            BarQCPConvTol::Float64=1e-6,
                            BarCorrectors::Int=-1,
                            BarHomogeneous::Int=-1,
                            BarOrder::Int=-1,
                            Crossover::Int=-1,
                            CrossoverBasis::Int=0,
                            QCPDual::Int=0,
                            Threads::Int=1)::GUROBIsettings

    return GUROBIsettings(Method, OutputFlag, DisplayInterval, IterationLimit, TimeLimit,
                          FeasibilityTol, OptimalityTol, PSDTol, MarkowitzTol,
                          InfUnbdInfo, NormAdjust, ObjScale, PerturbValue, Quad,
                          ScaleFlag, Sifting, SiftMethod, SimplexPricing,
                          BarIterLimit, BarConvTol, BarQCPConvTol, BarCorrectors,
                          BarHomogeneous, BarOrder, Crossover, CrossoverBasis,
                          QCPDual, Threads)
end
# structure used for storing data for OSQP solver
mutable struct GUROBIworkspace <: AbstractWorkspace
    # objective
    Q::SparseMatrixCSC{Float64}
    L::Array{Float64,1}

    # constraints
    A::SparseMatrixCSC{Float64}
    cnsLoBs::Array{Float64,1}
    cnsUpBs::Array{Float64,1}

    # variables
    varLoBs::Array{Float64,1}
    varUpBs::Array{Float64,1}

    # workspace
    environment::Gurobi.Env
    settings::GUROBIsettings
end


# getters and setters
function get_primalTolerance(workspace::GUROBIworkspace)::Float64
    return workspace.settings.FeasibilityTol
end
