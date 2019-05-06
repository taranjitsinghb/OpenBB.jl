
# @Author: Massimo De Mauri <massimo>
# @Date:   2019-02-25T17:07:48+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: inspect.jl
# @Last modified by:   massimo
# @Last modified time: 2019-04-30T19:28:00+02:00
# @License: apache 2.0
# @Copyright: {{copyright}}

# function to print BB status info to screen
function print_status(workspace::BBworkspace)::Nothing
    println("time: ",round(workspace.status.totalTime,digits = 2),
            " | best obj.: ", workspace.status.objUpB,
            " | best possible: ", workspace.status.objLoB,
            " | abs. gap.: ",workspace.status.absoluteGap,
            " | rel. gap.: ",round(workspace.status.relativeGap, digits = 2))
    return
end


###### getters and setters


# finds the best subproblem in the array
function best_subproblem(subproblem_pool::Array{BBsubproblem,1})::Union{BBsubproblem,Nothing}
    if length(subproblem_pool) != 0
        for i in length(subproblem_pool):-1:1
            if subproblem_pool[i].reliable
                return subproblem_pool[i]
            end
        end
    end
    return
end

# returns the best solution
function get_best_solution(workspace::BBworkspace)::Union{BBsubproblem,Nothing}

    subproblem = nothing
    if length(workspace.solutionPool) != 0
        subproblem = best_subproblem(workspace.solutionPool)
    end

    return subproblem
end

# returns the best subproblem
function get_best_subproblem(workspace::BBworkspace)::Union{BBsubproblem,Nothing}

    subproblem = nothing
    if length(workspace.solutionPool) != 0
        subproblem = best_subproblem(workspace.solutionPool)
    end

    if best_subproblem == nothing
        subproblem = best_subproblem(workspace.activeQueue)
    end

    if best_subproblem == nothing
        subproblem = best_subproblem(workspace.unactivePool)
    end

    return subproblem
end


# returns the number of variables
function get_numVariables(workspace::BBworkspace)::Int
    return get_numVariables(workspace.subsolverWS)
end

# returns the number of constraints
function get_numConstraints(workspace::BBworkspace)::Int
    return get_numConstraints(workspace.subsolverWS)
end

# this function returns the sparsity pattern of the constraint set
function get_constraints_sparsity(workspace::BBworkspace)::Tuple{Array{Int,1},Array{Int,1}}
    return get_constraints_sparsity(workspace.subsolverWS)
end

# this function returns the sparsity pattern a constraint in the constraint set
function get_constraint_sparsity(workspace::BBworkspace,index::Int)::Array{Int,1}
    return get_constraint_sparsity(workspace.subsolverWS,index)
end


# this function returns the sparsity pattern of the objective function
function get_objective_sparsity(workspace::BBworkspace)::Tuple{Array{Int,1},Array{Int,1}}
    return get_objective_sparsity(workspace.subsolverWS)
end


#
function get_variableBounds(workspace::BBworkspace)::Tuple{Array{Float64,1},Array{Float64,1}}
    return get_variableBounds(workspace.subsolverWS)
end

#
function get_constraintBounds(workspace::BBworkspace)::Tuple{Array{Float64,1},Array{Float64,1}}
    return get_constraintsBounds(workspace.subsolverWS)
end
