# @Author: Massimo De Mauri <massimo>
# @Date:   2019-02-19T12:14:43+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: update_problem.jl
# @Last modified by:   massimo
# @Last modified time: 2019-07-05T17:21:20+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

#
function append_constraints!(workspace::BBworkspace{T1,T2},
                             A::Union{Array{Float64,2},SparseMatrixCSC{Float64}},
                             cnsLoBs::Array{Float64,1},
                             cnsUpBs::Array{Float64,1};
                             suppressWarnings::Bool=false,
                             suppressUpdate::Bool=false,
                             localOnly::Bool=false)::Nothing where T1<:AbstractWorkspace where T2<:AbstractSharedMemory

    return insert_constraints!(workspace,A,cnsLoBs,cnsUpBs,
                               get_numConstraints(workspace)+1,
                               suppressWarnings=suppressWarnings,
                               suppressUpdate=suppressUpdate,
                               localOnly=localOnly)
end

#
function insert_constraints!(workspace::BBworkspace{T1,T2},
                             A::Union{Array{Float64,2},SparseMatrixCSC{Float64}},
                             cnsLoBs::Array{Float64,1},
                             cnsUpBs::Array{Float64,1},
                             index::Int;
                             suppressWarnings::Bool=false,
                             suppressUpdate::Bool=false,
                             localOnly::Bool=false)::Nothing where T1<:AbstractWorkspace where T2<:AbstractSharedMemory

    # check if it is possible to make changes
    if !suppressWarnings && workspace.status.description != "new" && !workspace.settings.dynamicMode && myid() == 1
        @warn "In order to correctly manipulate the problem formulation, OpenBB must be run in dynamic mode"
    end

    # call the same function on the other workers
    @sync if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)
        # call the local version of the function on the remote workers
        for p in 2:workspace.settings.numProcesses
            @async remotecall_fetch(Main.eval,p,:(OpenBB.insert_constraints!(workspace,$A,$cnsLoBs,$cnsUpBs,$index,
                                                                             suppressWarnings=$suppressWarnings,
                                                                             suppressUpdate=true,localOnly=true)))
        end

        # call the local version of the function on the current process
        insert_constraints!(workspace,A,cnsLoBs,cnsUpBs,index,
                            suppressWarnings=suppressWarnings,
                            suppressUpdate=true,localOnly=true)

    else

        # propagate the changes to the subsolver
        insert_constraints!(workspace.subsolverWS,A,cnsLoBs,cnsUpBs,index,suppressUpdate=true)

        # update all the problems in the activeQueue
        if length(workspace.activeQueue)>0
            for i in 1:length(workspace.activeQueue)
                splice!(workspace.activeQueue[i].cnsDual,index:index-1,zeros(size(A,1)))
            end
        end
        if length(workspace.unactivePool)>0
            for i in  1:length(workspace.unactivePool)
                splice!(workspace.unactivePool[i].cnsDual,index:index-1,zeros(size(A,1)))
            end
        end
        if length(workspace.solutionPool)>0
            for i in  1:length(workspace.solutionPool)
                splice!(workspace.solutionPool[i].cnsDual,index:index-1,zeros(size(A,1)))
            end
        end

        # adapt the workspace to the changes
        if !suppressUpdate
            update!(workspace,localOnly=true)
        end
    end

    return
end


#
function remove_constraints!(workspace::BBworkspace{T1,T2},indices::Array{Int,1};
                             suppressWarnings::Bool=false,
                             suppressUpdate::Bool=true,
                             localOnly::Bool=true)::Nothing where T1<:AbstractWorkspace where T2<:AbstractSharedMemory

    ## check if it is possible to make changes
    if !suppressWarnings && workspace.status.description != "new" && myid() == 1
        @warn "Removing constraints will invalidate the results found in the last solution process"
    end

    @sync if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)
        # call the local version of the function on the remote workers
        for p in 2:workspace.settings.numProcesses
            @async remotecall_fetch(Main.eval,p,:(OpenBB.remove_constraints!(workspace,$indices,
                                                                             suppressWarnings=$suppressWarnings,
                                                                             suppressUpdate=true,localOnly=true)))
        end
        # call the local version of the function on the main process
        remove_constraints!(workspace,indices,
                            suppressWarnings=suppressWarnings,
                            suppressUpdate=true,localOnly=true)
    else
        # propagate the changes to the nodes solver
        remove_constraints!(workspace.subsolverWS,indices,suppressUpdate=true)

        # update all the problems in the activeQueue
        for i in 1:length(workspace.activeQueue)
            deleteat!(workspace.activeQueue[i].cnsDual,indices)
        end
        # update all the problems in the solutionPool
        for i in 1:length(workspace.solutionPool)
            deleteat!(workspace.solutionPool[i].cnsDual,indices)
        end
        # update all the problems in the activeQueue
        for i in 1:length(workspace.unactivePool)
            deleteat!(workspace.unactivePool[i].cnsDual,indices)
        end

        # update the workspace
        if !suppressUpdate
            update!(workspace,localOnly=true)
        end
    end

    return
end


#
function permute_constraints!(workspace::BBworkspace{T1,T2},permutation::Array{Int,1};
                              suppressWarnings::Bool=false,
                              suppressUpdate::Bool=true,
                              localOnly::Bool=false)::Nothing where T1<:AbstractWorkspace where T2<:AbstractSharedMemory

    @sync if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)
        # call the local version of the function on the remote workers
        for p in 2:workspace.settings.numProcesses
            @async remotecall_fetch(Main.eval,p,:(OpenBB.permute_constraints!(workspace,$permutation,
                                                                              suppressWarnings=$suppressWarnings,
                                                                              suppressUpdate=true,localOnly=true)))
        end

        # call the local version of the function on the main process
        permute_constraints!(workspace,permutation,
                             suppressWarnings=suppressWarnings,
                             suppressUpdate=true,localOnly=true)

    else

        # propagate the change to the subsolver
        permute_constraints!(workspace.subsolverWS,permutation,suppressUpdate=true)

        # permute the optimality results
        if length(workspace.activeQueue)>0
            for i in 1:length(workspace.activeQueue)
                permute!(workspace.activeQueue[i].cnsDual,permutation)
            end
        end
        if length(workspace.unactivePool)>0
            for i in 1:length(workspace.unactivePool)
                permute!(workspace.unactivePool[i].cnsDual,permutation)
            end
        end
        if length(workspace.solutionPool)>0
            for i in 1:length(workspace.solutionPool)
                permute!(workspace.solutionPool[i].cnsDual,permutation)
            end
        end

        # update the subsolver workspace only
        if !suppressUpdate
            update!(workspace.subsolverWS,localOnly=true)
        end
    end

    return
end


#
function update_bounds!(workspace::BBworkspace{T1,T2};
                        cnsLoBs::Array{Float64,1}=Float64[],
                        cnsUpBs::Array{Float64,1}=Float64[],
                        varLoBs::Array{Float64,1}=Float64[],
                        varUpBs::Array{Float64,1}=Float64[],
                        suppressWarnings::Bool=false,
                        suppressUpdate::Bool=false,
                        localOnly::Bool=false)::Nothing where T1<:AbstractWorkspace where T2<:AbstractSharedMemory

    # check the correctness of the input
    @assert length(cnsLoBs)==0 || length(cnsLoBs)==length(workspace.subsolverWS.cnsLoBs)
    @assert length(cnsUpBs)==0 || length(cnsUpBs)==length(workspace.subsolverWS.cnsUpBs)
    @assert length(varLoBs)==0 || length(varLoBs)==length(workspace.subsolverWS.varLoBs)
    @assert length(varUpBs)==0 || length(varUpBs)==length(workspace.subsolverWS.varUpBs)


    # check if it is possible to make changes
    if !suppressWarnings && workspace.status.description != "new" && !workspace.settings.dynamicMode && myid() == 1
        @warn "In order to correctly manipulate the problem formulation, OpenBB must be run in dynamic mode"

        # check if a bounds relaxation was requested
        if (length(cnsLoBs) > 0 && any(@. cnsLoBs < workspace.subsolverWS.cnsLoBs)) ||
           (length(cnsUpBs) > 0 && any(@. cnsUpBs > workspace.subsolverWS.cnsUpBs)) ||
           (length(varLoBs) > 0 && any(@. varLoBs < workspace.subsolverWS.varLoBs)) ||
           (length(varUpBs) > 0 && any(@. varUpBs > workspace.subsolverWS.varLoBs))

            @warn "It is not possible to update the BB tree for a bound relaxation, a restart might be necessary"
        end
    end

    # ensure consistency in the inputs
    if length(varLoBs) == 0
        varLoBs = copy(workspace.subsolverWS.varLoBs)
    end
    if length(varUpBs) == 0
        varUpBs = copy(workspace.subsolverWS.varUpBs)
    end
    if length(cnsLoBs) == 0
        cnsLoBs = copy(workspace.subsolverWS.cnsLoBs)
    end
    if length(cnsUpBs) == 0
        cnsUpBs = copy(workspace.subsolverWS.cnsUpBs)
    end


    @sync if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)
        # call the local version of the function on the remote workers
        for p in 2:workspace.settings.numProcesses
            @async remotecall_fetch(Main.eval,p,:(OpenBB.update_bounds!(workspace,
                                                                        cnsLoBs=$cnsLoBs,cnsUpBs=$cnsUpBs,
                                                                        varLoBs=$varLoBs,varUpBs=$varUpBs,
                                                                        suppressWarnings=$suppressWarnings,
                                                                        suppressUpdate=true,localOnly=true)))
        end

        # call the local version of the function on the main process
        update_bounds!(workspace,
                       cnsLoBs=cnsLoBs,cnsUpBs=cnsUpBs,
                       varLoBs=varLoBs,varUpBs=varUpBs,
                       suppressWarnings=suppressWarnings,
                       suppressUpdate=true,localOnly=true)
    else

        # propagate the changes to the nodes solver
        update_bounds!(workspace.subsolverWS,cnsLoBs,cnsUpBs,varLoBs,varUpBs,suppressUpdate=true)

        # adapt the workspace to the changes
        if !suppressUpdate
            update!(workspace,localOnly=true)
        end
    end

    return
end

# ...
function set_objective!(workspace::BBworkspace{T1,T2},newObjective::T3;
                        suppressWarnings::Bool=false,
                        suppressUpdate::Bool=false,
                        localOnly::Bool=false)::Nothing where T1<:AbstractWorkspace where T2<:AbstractSharedMemory where T3 <: AbstractObjective

    # check if it is possible to make changes
    if !suppressWarnings
        if workspace.status.description != "new" && !workspace.settings.dynamicMode && myid() == 1
            @warn "In order to correctly manipulate the problem formulation, OpenBB must be run in dynamic mode"
        else
            # warn the user about the potential dangers
            @warn "changing the objective after some iterations is potentially destructive, please be sure that the new objective is always greater or equal to the old one"
        end
    end

    @sync if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)
        # call the local version of the function on the remote workers
        for p in 2:workspace.settings.numProcesses
            @async remotecall_fetch(Main.eval,p,:(OpenBB.set_objective!(workspace,$newObjective,
                                                                         suppressWarnings=true,
                                                                         suppressUpdate=true,localOnly=true)))
        end

        # call the local version of the function on the main process
        set_objective!(workspace,newObjective,
                        suppressWarnings=true,
                        suppressUpdate=true,localOnly=true)
    else

        # propagate the change to the subsolver
        set_objective!(workspace.subsolverWS,newObjective,suppressUpdate=true)

        # update the workspace
        if !suppressUpdate
            update!(workspace,localOnly=true)
        end
    end

    return
end

#...
function set_constraintSet!(workspace::BBworkspace{T1,T2},newConstraintSet::T3;
                           suppressWarnings::Bool=false,
                           suppressUpdate::Bool=false,
                           localOnly::Bool=false)::Nothing where T1<:AbstractWorkspace where T2<:AbstractSharedMemory where T3 <: AbstractConstraintSet

    # check if it is possible to make changes
    if !suppressWarnings
        if workspace.status.description != "new" && !workspace.settings.dynamicMode && myid() == 1
            @warn "In order to correctly manipulate the problem formulation, OpenBB must be run in dynamic mode"
        else
            # warn the user about the potential dangers
            @warn "changing the constraintSet after some iterations is potentially destructive, please be sure that the new constraint set is more restrictive than the old"
        end
    end

    @sync if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)
        # call the local version of the function on the remote workers
        for p in 2:workspace.settings.numProcesses
            @async remotecall_fetch(Main.eval,p,:(OpenBB.set_constraintSet!(workspace,$newConstraintSet,
                                                                            suppressWarnings=true,
                                                                            suppressUpdate=true,localOnly=true)))
        end

        # call the local version of the function on the main process
        set_constraintSet!(workspace,newConstraintSet,
                           suppressWarnings=true,
                           suppressUpdate=true,localOnly=true)
    else

        # reset the constraints dual
        if length(workspace.activeQueue)>0
            for i in 1:length(workspace.activeQueue)
                workspace.activeQueue[i].cnsDual = zeros(get_numConstraints(newConstraintSet))
            end
        end
        if length(workspace.unactivePool)>0
            for i in 1:length(workspace.unactivePool)
                workspace.unactivePool[i].cnsDual = zeros(get_numConstraints(newConstraintSet))
            end
        end
        if length(workspace.solutionPool)>0
            for i in 1:length(workspace.solutionPool)
                workspace.solutionPool[i].cnsDual = zeros(get_numConstraints(newConstraintSet))
            end
        end

        # propagate the change to the subsolver
        set_constraintSet!(workspace.subsolverWS,newConstraintSet,suppressUpdate=true)

        # update the workspace
        if !suppressUpdate
            update!(workspace,localOnly=true)
        end
    end

    return
end

#
function append_problem!(workspace::BBworkspace{T1,T2},problem::Problem;
                         suppressWarnings::Bool=false,
                         suppressUpdate::Bool=false,
                         localOnly::Bool=false)::Nothing where T1<:AbstractWorkspace where T2<:AbstractSharedMemory

    # check if it is possible to make changes
    if !suppressWarnings && workspace.status.description != "new" && !workspace.settings.dynamicMode && myid() == 1
        @warn "In order to correctly manipulate the problem formulation, OpenBB must be run in dynamic mode"
    end

    @sync if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)
        # call the local version of the function on the remote workers
        for p in 2:workspace.settings.numProcesses
            @async remotecall_fetch(Main.eval,p,:(OpenBB.append_problem!(workspace,$problem,
                                                                         suppressWarnings=$suppressWarnings,
                                                                         suppressUpdate=true,localOnly=true)))
        end
        # call the local version of the function on the main process
        append_problem!(workspace,problem,
                        suppressWarnings=suppressWarnings,
                        suppressUpdate=true,localOnly=true)
    else
        # collect info on the problem
        nVars1 = length(workspace.subsolverWS.varLoBs)
        nVars2 = length(problem.varSet.loBs)
        nCnss1 = length(workspace.subsolverWS.cnsLoBs)
        nCnss2 = length(problem.cnsSet.loBs)
        # propagate the changes to the nodes solver
        reliableObjLoBs = append_problem!(workspace.subsolverWS,problem,suppressUpdate=true)
        # update all the sub-problems
        newPrimal = problem.varSet.vals
        for i in  1:length(workspace.activeQueue)
            # mark all the problems as unreliable if necessary
            workspace.activeQueue[i].reliable = reliableObjLoBs && workspace.activeQueue[i].reliable
            # extend branching bounds
            append!(workspace.activeQueue[i].branchLoBs,-Infs(length(problem.varSet.dscIndices)))
            append!(workspace.activeQueue[i].branchUpBs, Infs(length(problem.varSet.dscIndices)))
            # extend primal and dual optimization results
            append!(workspace.activeQueue[i].primal,copy(newPrimal))
            append!(workspace.activeQueue[i].bndDual,zeros(nVars2))
            append!(workspace.activeQueue[i].cnsDual,zeros(nCnss2))
        end
        for i in  1:length(workspace.solutionPool)
            # mark all the problems as unreliable if necessary
            workspace.solutionPool[i].reliable = reliableObjLoBs && workspace.solutionPool[i].reliable
            # extend branching bounds
            append!(workspace.solutionPool[i].branchLoBs,-Infs(length(problem.varSet.dscIndices)))
            append!(workspace.solutionPool[i].branchUpBs, Infs(length(problem.varSet.dscIndices)))
            # extend primal and dual optimization results
            append!(workspace.solutionPool[i].primal,copy(newPrimal))
            append!(workspace.solutionPool[i].bndDual,zeros(nVars2))
            append!(workspace.solutionPool[i].cnsDual,zeros(nCnss2))
        end
        for i in  1:length(workspace.unactivePool)
            # mark all the problems as unreliable if necessary
            workspace.unactivePool[i].reliable = reliableObjLoBs && workspace.unactivePool[i].reliable
            # extend branching bounds
            append!(workspace.unactivePool[i].branchLoBs,-Infs(length(problem.varSet.dscIndices)))
            append!(workspace.unactivePool[i].branchUpBs, Infs(length(problem.varSet.dscIndices)))
            # extend primal and dual optimization results
            append!(workspace.unactivePool[i].primal,copy(newPrimal))
            append!(workspace.unactivePool[i].bndDual,zeros(nVars2))
            append!(workspace.unactivePool[i].cnsDual,zeros(nCnss2))
        end
        # set objective lower bound to -Inf if the requested update invalidates the already computed lower bounds
        if !reliableObjLoBs
            workspace.status.objLoB = -Inf
        end
        # update discrete variables
        append!(workspace.dscIndices,@. problem.varSet.dscIndices + nVars1)
        append!(workspace.sos1Groups,@. problem.varSet.sos1Groups + maximum(workspace.sos1Groups) + 1)
        workspace.pseudoCosts = (vcat(workspace.pseudoCosts[1],problem.varSet.pseudoCosts),
                                 vcat(workspace.pseudoCosts[2],repeat([0 0],length(problem.varSet.dscIndices),1)))

         # adapt the workspace to the changes
         if !suppressUpdate
             update!(workspace,localOnly=true)
         end
    end

    return
end


#
function integralize_variables!(workspace::BBworkspace{T1,T2},newDscIndices::Array{Int,1};
                                newSos1Groups::Array{Int,1}=Int[],
                                suppressWarnings::Bool=false,
                                suppressUpdate::Bool=false,
                                localOnly::Bool=false)::Nothing where T1<:AbstractWorkspace where T2<:AbstractSharedMemory

    # check the correctness of the input
    @assert length(intersect(newDscIndices,workspace.dscIndices)) == 0


    @sync if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)
        # call the local version of the function on the remote workers
        for p in 2:workspace.settings.numProcesses
            @async remotecall_fetch(Main.eval,p,:(OpenBB.integralize_variables!(workspace,$newDscIndices,
                                                                                newSos1Groups=$newSos1Groups,
                                                                                suppressWarnings=$suppressUpdate,
                                                                                suppressUpdate=true,localOnly=true)))
        end
        # call the local version of the function on the main process
        integralize_variables!(workspace,newDscIndices,
                               newSos1Groups=newSos1Groups,
                               suppressWarnings=suppressUpdate,
                               suppressUpdate=true,localOnly=true)
    else
        # check if it is possible to make changes
        if !suppressWarnings && workspace.status.description != "new" && !workspace.settings.dynamicMode && myid() == 1
            @warn "In order to correctly manipulate the problem formulation, OpenBB must be run in dynamic mode"
        end
        # check correctness of the inputs
        if length(newSos1Groups) == 0
            newSos1Groups = repeat([-1],length(newDscIndices))
        elseif length(newSos1Groups) != length(newDscIndices)
            @error "newSos1Groups should either be empty or have the same dimension than newDscIndices"
        end

        # add the new discrete variables
        append!(workspace.dscIndices,copy(newDscIndices))
        append!(workspace.sos1Groups,copy(newSos1Groups))
        workspace.pseudoCosts = (vcat(workspace.pseudoCosts[1],repeat(sum(workspace.pseudoCosts[1],dims=1)/size(workspace.pseudoCosts[1],1),length(newDscIndices),1)),
                                 vcat(workspace.pseudoCosts[2],repeat([0 0],length(newDscIndices),1)))

        tmpPerm = sortperm(workspace.dscIndices)
        permute!(workspace.dscIndices,tmpPerm)
        permute!(workspace.sos1Groups,tmpPerm)
        workspace.pseudoCosts = (workspace.pseudoCosts[1][tmpPerm,:],workspace.pseudoCosts[2][tmpPerm,:])

        # propagate the change to the nodes
        for k in 1:length(workspace.activeQueue)
            append!(workspace.activeQueue[k].branchLoBs,-Inf)
            append!(workspace.activeQueue[k].branchUpBs, Inf)
            permute!(workspace.activeQueue[k].branchLoBs,tmpPerm)
            permute!(workspace.activeQueue[k].branchUpBs,tmpPerm)
        end
        for k in 1:length(workspace.unactivePool)
            append!(workspace.unactivePool[k].branchLoBs,-Inf)
            append!(workspace.unactivePool[k].branchUpBs, Inf)
            permute!(workspace.unactivePool[k].branchLoBs,tmpPerm)
            permute!(workspace.unactivePool[k].branchUpBs,tmpPerm)
        end
        for k in 1:length(workspace.solutionPool)
            append!(workspace.solutionPool[k].branchLoBs,-Inf)
            append!(workspace.solutionPool[k].branchUpBs, Inf)
            permute!(workspace.solutionPool[k].branchLoBs,tmpPerm)
            permute!(workspace.solutionPool[k].branchUpBs,tmpPerm)
        end

        # adapt the workspace to the changes
        if !suppressUpdate
            update!(workspace,localOnly=true)
        end
    end

    return
end
