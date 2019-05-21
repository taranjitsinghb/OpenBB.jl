
# @Author: Massimo De Mauri <massimo>
# @Date:   2019-02-25T17:07:48+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: inspect.jl
# @Last modified by:   massimo
# @Last modified time: 2019-05-21T16:33:25+02:00
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


# finds the best node in the array
function find_best_node(node_pool::Array{BBnode,1})::Union{BBnode,Nothing}
    if length(node_pool) != 0
        for i in length(node_pool):-1:1
            if node_pool[i].reliable
                return node_pool[i]
            end
        end
    end
    return
end

# returns the best solution
function get_best_solution(workspace::BBworkspace;localOnly::Bool=false)::Union{BBnode,Nothing}

    if length(workspace.solutionPool) != 0
        solution = workspace.solutionPool[end]
    else
        solution = nothing
    end

    if localOnly || workspace.globalInfo == nothing
        return solution
    elseif solution != nothing && solution.objVal == workspace.status.objLoB
        return solution
    else
        for p in 2:workspace.settings.numProcesses
            node = remotecall_fetch(Main.eval,p,:(OpenBB.get_best_solution(workspace,localOnly=true)))
            if node != nothing && node.objVal == workspace.globalInfo[1]
                solution = node
            end
        end
    end
    return solution
end

# returns the best node
function get_best_node(workspace::BBworkspace;localOnly::Bool=false)::Union{BBnode,Nothing}

    # define dummy best node
    bestNode = nothing

    @sync if !localOnly && workspace.globalInfo != nothing

        nodes = Array{BBnode,1}(undef,workspace.settings.numProcesses)
        # call the local version of the function on the remote workers
        for p in 2:workspace.settings.numProces
            @async nodes[p] = remotecall_fetch(Main.eval,p,:(OpenBB.get_best_node(workspace,localOnly=true)))
        end

        # call the local version of the function on the current process
        nodes[1] = get_best_node(workspace,localOnly=true)

        # choose the best of the returned nodes
        for node in nodes
            if bestNode == nothing ||
               (node.avgFrac==0 && node.avgFrac >0) ||
               workspace.settings.expansion_priority_rule(node,bestNode,workspace.status)
               bestNode = node
           end
       end
   else
        # take the last solution generated
        if length(workspace.solutionPool) > 0
            bestNode = workspace.solutionPool[end]
        end
        # take the first problem in the active queue
        if bestNode == nothing && length(activeQueue) > 0
            bestNode = activeQueue[end]
        end
        # check the unactivePool for better nodes
        if bestNode == nothing
            for node in workspace.unactivePool
                if bestNode == nothing || workspace.settings.expansion_priority_rule(bestNode,node,workspace.status)
                    bestNode = node
                end
            end
        end
    end

    return deepcopy(bestNode)

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
