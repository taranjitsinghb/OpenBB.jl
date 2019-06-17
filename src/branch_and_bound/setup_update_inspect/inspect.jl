
# @Author: Massimo De Mauri <massimo>
# @Date:   2019-02-25T17:07:48+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: inspect.jl
# @Last modified by:   massimo
# @Last modified time: 2019-06-08T13:46:02+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

# function to print BB status info to screen
function print_status(workspace::BBworkspace{T1,T2})::Nothing where T1<:AbstractWorkspace where T2<:AbstractSharedMemory

    if workspace.sharedMemory isa NullSharedMemory

        globalObjectiveUpB = workspace.status.objUpB
        globalObjectiveLoB = workspace.status.objLoB
        globalAbsGap = workspace.status.absoluteGap
        globalRelGap = workspace.status.relativeGap

    else
        globalObjectiveUpB = workspace.sharedMemory.objectiveBounds[end]
        globalObjectiveLoB = minimum(workspace.sharedMemory.objectiveBounds[1:end-1])

        # recompute optimality gaps
        if globalObjectiveUpB == Inf || globalObjectiveLoB == -Inf
            globalAbsGap = globalRelGap =  Inf
        else
            globalAbsGap = globalObjectiveUpB - globalObjectiveLoB
            globalRelGap = workspace.status.absoluteGap/abs(1e-10 + globalObjectiveUpB)
        end

    end

    println("time: ",round(workspace.status.totalTime,digits = 2),
            " | best obj.: ", workspace.status.objUpB,
            " | best possible: ", globalObjectiveLoB,
            " | abs. gap.: ", globalAbsGap,
            " | rel. gap.: ",round(globalRelGap, digits = 2))

    return
end


# returns the best solution
function get_best_solution(workspace::BBworkspace{T1,T2};localOnly::Bool=false)::AbstractBBnode where T1<:AbstractWorkspace where T2<:AbstractSharedMemory

    solution = NullBBnode()

    if length(workspace.solutionPool) != 0
        solution = workspace.solutionPool[end]
    end

    # check the other workers if needed/required
    if !(localOnly || workspace.sharedMemory isa NullSharedMemory) &&
       (solution isa NullBBnode || solution.objective > workspace.status.objLoB)
        for p in 2:workspace.settings.numProcesses
            node = remotecall_fetch(Main.eval,p,:(OpenBB.get_best_solution(workspace,localOnly=true)))
            if !(node isa NullBBnode) && node.objective == workspace.sharedMemory.objectiveBounds[end]
                solution = node
                break
            end
        end
    end
    return solution
end

# returns the best node
function get_best_node(workspace::BBworkspace{T1,T2};localOnly::Bool=false)::AbstractBBnode  where T1<:AbstractWorkspace where T2<:AbstractSharedMemory

    # define dummy best node
    bestNode = NullBBnode()

    if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)

        nodes = Array{AbstractBBnode,1}(undef,workspace.settings.numProcesses)

        # call the local version of the function on the remote workers
        @sync for p in 2:workspace.settings.numProcesses
            @async nodes[p] = remotecall_fetch(Main.eval,p,:(OpenBB.get_best_node(workspace,localOnly=true)))
        end

        # call the local version of the function on the current process
        nodes[1] = get_best_node(workspace,localOnly=true)

        # choose the best of the returned nodes
        for node in nodes
            if bestNode isa NullBBnode || # there is no best node yet
               (node.avgAbsFrac==0 && bestNode.avgAbsFrac > 0) || # the new node is a solution while the best so far isn't
               expansion_priority_rule(workspace.settings.expansionPriorityRule,node,bestNode,workspace.status) # the new node is better than the best so far

                # set the new node as the best one
                bestNode = node
           end
       end
   else
        # take the last solution generated if any
        if length(workspace.solutionPool) > 0
            bestNode = workspace.solutionPool[end]
        end
        # otherwise, take the first problem in the active queue
        if bestNode isa NullBBnode && length(workspace.activeQueue) > 0
            bestNode = activeQueue[end]
        end
        # otherwise, check the unactivePool for better nodes
        if bestNode isa NullBBnode && length(workspace.unactivePool) > 0
            for node in workspace.unactivePool
                if bestNode isa NullBBnode || expansion_priority_rule(workspace.settings.expansionPriorityRule,bestNode,node,workspace.status)
                    bestNode = node
                end
            end
        end
    end

    return deepcopy(bestNode)
end


# returns the number of variables
function get_numVariables(workspace::BBworkspace{T1,T2})::Int where T1<:AbstractWorkspace where T2<:AbstractSharedMemory
    return get_numVariables(workspace.subsolverWS)
end

# returns the number of discrete variables
function get_numDiscreteVariables(workspace::BBworkspace{T1,T2})::Int where T1<:AbstractWorkspace where T2<:AbstractSharedMemory
    return length(workspace.dscIndices)
end

# returns the number of constraints
function get_numConstraints(workspace::BBworkspace{T1,T2})::Int where T1<:AbstractWorkspace where T2<:AbstractSharedMemory
    return get_numConstraints(workspace.subsolverWS)
end

# this function returns the sparsity pattern of the constraint set
function get_constraints_sparsity(workspace::BBworkspace{T1,T2})::Tuple{Array{Int,1},Array{Int,1}} where T1<:AbstractWorkspace where T2<:AbstractSharedMemory
    return get_constraints_sparsity(workspace.subsolverWS)
end

# this function returns the sparsity pattern a constraint in the constraint set
function get_constraint_sparsity(workspace::BBworkspace{T1,T2},index::Int)::Array{Int,1}  where T1<:AbstractWorkspace where T2<:AbstractSharedMemory
    return get_constraint_sparsity(workspace.subsolverWS,index)
end

# this function returns the sparsity pattern of the objective function
function get_objective_sparsity(workspace::BBworkspace{T1,T2})::Tuple{Array{Int,1},Array{Int,1}}  where T1<:AbstractWorkspace where T2<:AbstractSharedMemory
    return get_objective_sparsity(workspace.subsolverWS)
end

#
function get_variableBounds(workspace::BBworkspace{T1,T2})::Tuple{Array{Float64,1},Array{Float64,1}}  where T1<:AbstractWorkspace where T2<:AbstractSharedMemory
    return get_variableBounds(workspace.subsolverWS)
end

#
function get_constraintBounds(workspace::BBworkspace{T1,T2})::Tuple{Array{Float64,1},Array{Float64,1}}  where T1<:AbstractWorkspace where T2<:AbstractSharedMemory
    return get_constraintsBounds(workspace.subsolverWS)
end
