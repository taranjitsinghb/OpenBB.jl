# @Author: Massimo De Mauri <massimo>
# @Date:   2019-05-23T11:28:34+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: BBsettings.jl
# @Last modified by:   massimo
# @Last modified time: 2019-06-03T19:33:50+02:00
# @License: apache 2.0
# @Copyright: {{copyright}}


# define abstract and null types
abstract type AbstractSettings end; struct NullSettings  <: AbstractSettings  end


mutable struct BBsettings <: AbstractSettings
    # execution modifiers
    verbose::Bool                           # print info during execution
    statusInfoPeriod::Float64                 # frequency of status info print
    numProcesses::Int                       # max number of process to launch
    stopAfterSolution::Bool                 # stop workers after the solution for the current problem has been found
    dynamicMode::Bool                       # store in memory suboptimal solutions and nodes to allow later updates
    # problem bounds
    integerTolerance::Float64               # integer tolerance
    primalTolerance::Float64                # constraint violation tolerance
    objectiveCutoff::Float64                # look only for solutions that are better than the provided upper bound
    # priority rules
    expansionPriorityRule::Function       # ordering of the nodes in the activeQueue
    branchingPriorityRule::Function       # ordering of the discrete variables for branching
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
                     statusInfoPeriod::Float64 = 1.,
                     numProcesses::Int=0,
                     stopAfterSolution::Bool = true,
                     dynamicMode::Bool=false,
                     integerTolerance::Float64=1e-4,
                     primalTolerance::Float64=1e-4,
                     objectiveCutoff::Float64=Inf,
                     expansionPriorityRule::Function=lower_mixed,
                     branchingPriorityRule::Function=most_fractional,
                     unreliable_subps_priority::Int=0,
                     custom_stopping_rule::Function=x->false,
                     timeLimit::Float64=Inf,
                     numSolutionsLimit::Int=0,
                     absoluteGapTolerance::Float64=1e-4,
                     relativeGapTolerance::Float64=1e-6,
                     roundingHeuristicsThreshold::Float64 = -1.
                     )::BBsettings


    return BBsettings(verbose,statusInfoPeriod,numProcesses,stopAfterSolution,dynamicMode,
                      integerTolerance,primalTolerance,objectiveCutoff,
                      expansionPriorityRule,branchingPriorityRule,unreliable_subps_priority,
                      custom_stopping_rule,timeLimit,numSolutionsLimit,absoluteGapTolerance,relativeGapTolerance,roundingHeuristicsThreshold)
end
