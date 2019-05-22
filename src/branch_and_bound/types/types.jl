# @Author: Massimo De Mauri <massimodemauri>
# @Date:   2019-02-06T16:19:47+01:00
# @Email:  massimo.demauri@gmail.com
# @Project: OpenBB
# @Filename: types.jl
# @Last modified by:   massimo
# @Last modified time: 2019-05-22T10:41:10+02:00
# @License: apache 2.0
# @Copyright: {{copyright}}


######### types definition ##########
# define abstract and null types
abstract type AbstractSettings end; struct NullSettings  <: AbstractSettings  end
abstract type AbstractWorkspace end; struct NullWorkspace <: AbstractWorkspace end
abstract type AbstractBBnode end; struct NullBBnode <: AbstractBBnode end

# this is the set of data that distinguishes a node from another
mutable struct BBnode <: AbstractBBnode
    branchLoBs::Dict{Int,Float64}
    branchUpBs::Dict{Int,Float64}
    pseudoCosts::Array{Float64,1}
    primal::Array{Float64,1}
    bndDual::Array{Float64,1}
    cnsDual::Array{Float64,1}
    avgFrac::Float64
    objVal::Float64
    reliable::Bool
end

# this node is used to arrest the branch and bound process
struct KillerNode <:AbstractBBnode
    count::Int
end


mutable struct BBstatus
    numSolutions::Int
    objLoB::Float64
    objUpB::Float64
    absoluteGap::Float64
    relativeGap::Float64
    totalTime::Float64
    numRelaxationsSolved::Int64
    description::String
end

function BBstatus(; objLoB::Float64=-Inf,objUpB::Float64=Inf,
                    absoluteGap::Float64= Inf,relativeGap::Float64=Inf,
                    totalTime::Float64=0.0,numRelaxationsSolved::Int=0,numSolutions::Int=0,
                    description::String="new")::BBstatus
    return BBstatus(numSolutions,objLoB,objUpB,absoluteGap,relativeGap,totalTime,numRelaxationsSolved,description)
end



mutable struct BBsettings <: AbstractSettings
    # execution modifiers
    verbose::Bool                           # print info during execution
    iterationInfoFreq::Int                  # frequency of iteration info print
    numProcesses::Int                       # max number of process to launch
    stopAfterSolution::Bool                 # stop workers after the solution for the current problem has been found
    dynamicMode::Bool                       # store in memory suboptimal solutions and nodes to allow later updates
    useSosConstraints::Bool                 # consider the sos constraints
    # problem bounds
    integerTolerance::Float64               # integer tolerance
    primalTolerance::Float64                # constraint violation tolerance
    objectiveCutoff::Float64                # look only for solutions that are better than the provided upper bound
    # priority rules
    expansion_priority_rule::Function       # ordering of the nodes in the activeQueue
    branching_priority_rule::Function       # ordering of the discrete variables for branching
    unreliable_subps_priority::Int          # activeQueue insertion priority for unreliable nodes (-1->low, 0->normal, 1->high)
    # stopping criteria
    custom_stopping_rule::Function          # user-defined stopping criterion
    timeLimit::Float64                      # max running time
    numSolutionsLimit::Int                  # stop after the given number of integral solutions have been found
    absoluteGapTolerance::Float64                # stop if the absolute optimality gap is below the given level
    relativeGapTolerance::Float64                # stop if the relative optimality gap is below the given level
    # heuristic search
    roundingHeuristicsThreshold::Float64    # use the simple rounding heuristics whenever the average fractionality is under the threshold
end



function BBsettings(;verbose::Bool=false,
                     iterationInfoFreq::Int = 10,
                     numProcesses::Int=0,
                     stopAfterSolution::Bool = true,
                     dynamicMode::Bool=false,
                     useSosConstraints::Bool=true,
                     integerTolerance::Float64=1e-4,
                     primalTolerance::Float64=1e-4,
                     objectiveCutoff::Float64=Inf,
                     expansion_priority_rule::Function=lower_mixed,
                     branching_priority_rule::Function=pseudo_cost,
                     unreliable_subps_priority::Int=0,
                     custom_stopping_rule::Function=x->false,
                     timeLimit::Float64=Inf,
                     numSolutionsLimit::Int=0,
                     absoluteGapTolerance::Float64=1e-4,
                     relativeGapTolerance::Float64=1e-6,
                     roundingHeuristicsThreshold::Float64 = -1.
                     )::BBsettings


    return BBsettings(verbose,iterationInfoFreq,numProcesses,stopAfterSolution,dynamicMode,
                      useSosConstraints,integerTolerance,primalTolerance,objectiveCutoff,
                      expansion_priority_rule,branching_priority_rule,unreliable_subps_priority,
                      custom_stopping_rule,timeLimit,numSolutionsLimit,absoluteGapTolerance,relativeGapTolerance,roundingHeuristicsThreshold)
end


# this structure collects all the informations needed to start and resume a B&B
# run
struct BBworkspace{T<:AbstractWorkspace} <: AbstractWorkspace
    # problem description
    subsolverWS::T
    dscIndices::Array{Int64,1}
    sos1Groups::Array{Int64,1}
    sosConstraints::LinearCns
    # branch and bound status
    activeQueue::Array{BBnode,1}
    solutionPool::Array{BBnode,1}
    unactivePool::Array{BBnode,1}
    status::BBstatus
    # multiprocessing communication
    inputChannel::Union{RemoteChannel{Channel{AbstractBBnode}},Nothing} # receive BBnodes from other workers
    outputChannel::Union{RemoteChannel{Channel{AbstractBBnode}},Nothing} # send BBnodes to other workers
    globalInfo::Union{SharedArray{Float64,1},Nothing} # info on the global status of the process
    # user settings
    settings::BBsettings
end


struct SubSolution
    primal::Array{Float64,1}
    bndDual::Array{Float64,1}
    cnsDual::Array{Float64,1}
    objVal::Float64
    status::Int
    time::Float64
end
