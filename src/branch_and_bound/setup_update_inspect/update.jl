# @Author: Massimo De Mauri <massimo>
# @Date:   2019-02-19T12:14:43+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: update.jl
# @Last modified by:   massimo
# @Last modified time: 2019-04-30T17:54:16+02:00
# @License: apache 2.0
# @Copyright: {{copyright}}


# return the workspace to the initial state
function reset!(workspace::BBworkspace)::Nothing

    deleteat!(workspace.activeQueue,1:length(workspace.activeQueue))
    deleteat!(workspace.solutionPool,1:length(workspace.solutionPool))
    deleteat!(workspace.unactivePool,1:length(workspace.unactivePool))
    push!(workspace.activeQueue,BBnode(Dict{Int,Float64}(),Dict{Int,Float64}(),
                                             zeros(get_numVariables(workspace)),
                                             zeros(get_numVariables(workspace)),
                                             zeros(get_numConstraints(workspace)),
                                             1.0,-Inf,false))

    # reset the status
    defaultStatus = BBstatus()
    for field in fieldnames(BBstatus)
        setfield!(workspace.status,field,getfield(defaultStatus,field))
    end

    return
end

#
function update!(workspace::BBworkspace)::Nothing

    # update the subsolver workspace
    update!(workspace.subsolverWS)

    # reset the explored sub-problems
    reset_explored_nodes!(workspace)

    return
end


#
function reset_explored_nodes!(workspace::BBworkspace)::Nothing


    # adapt the workspace to the changes
    append!(workspace.activeQueue,
            sort(workspace.solutionPool,
                 lt=(l,r)->workspace.settings.expansion_priority_rule(l,r,workspace.status),
                 rev=true))

    append!(workspace.activeQueue,
            sort(workspace.unactivePool,
                 lt=(l,r)->workspace.settings.expansion_priority_rule(l,r,workspace.status),
                 rev=true))

    deleteat!(workspace.solutionPool,1:length(workspace.solutionPool))
    deleteat!(workspace.unactivePool,1:length(workspace.unactivePool))

    workspace.status.objUpB = Inf
    workspace.status.absoluteGap = Inf
    workspace.status.relativeGap = Inf
    workspace.status.description = "interrupted"

    sort!(workspace.activeQueue,lt=(l,r)->workspace.settings.expansion_priority_rule(l,r,workspace.status),rev=true,alg=MergeSort)

    return
end


#
function append_constraints!(workspace::BBworkspace,
                             A::Union{Array{Float64,2},SparseMatrixCSC{Float64}},
                             cnsLoBs::Array{Float64,1},
                             cnsUpBs::Array{Float64,1};
                             suppressWarnings::Bool=false,
                             suppressUpdate::Bool=false)::Nothing

    return insert_constraints!(workspace,A,cnsLoBs,cnsUpBs,
                               get_numConstraints(workspace.subsolverWS)+1,
                               suppressWarnings=suppressWarnings,
                               suppressUpdate=suppressWarnings)
end

#
function insert_constraints!(workspace::BBworkspace,
                             A::Union{Array{Float64,2},SparseMatrixCSC{Float64}},
                             cnsLoBs::Array{Float64,1},
                             cnsUpBs::Array{Float64,1},
                             index::Int;
                             suppressWarnings::Bool=false,
                             suppressUpdate::Bool=false)::Nothing

    # check if it is possible to make changes
    if !suppressWarnings && workspace.status.description != "new" && !workspace.settings.dynamicMode
        @warn "In order to correctly manipulate the problem formulation, OpenBB must be run in dynamic mode"
    end

    # propagate the changes to the nodes solver
    insert_constraints!(workspace.subsolverWS,A,cnsLoBs,cnsUpBs,index,suppressUpdate=true)



    # update all the problems in the activeQueue
    if length(workspace.activeQueue)>0
        for i in  1:length(workspace.activeQueue)
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
        update!(workspace)
    end


    return
end


#
function remove_constraints!(workspace::BBworkspace,indices::Array{Int,1};
                             suppressWarnings::Bool=false,suppressUpdate::Bool=true)::Nothing

    ## check if it is possible to make changes
    if !suppressWarnings && workspace.status.description != "new"
        @warn "Removing constraints will invalidate the results found in the last solution process"
    end

    # propagate the changes to the nodes solver
    remove_constraints!(workspace.subsolverWS,indices,suppressUpdate=suppressUpdate)


    # adapt the workspace to the changes
    reset_explored_nodes!(workspace)

    # update all the problems in the activeQueue
    for i in 1:length(workspace.activeQueue)
        deleteat!(workspace.activeQueue[i].cnsDual,indices)
    end


    return
end


#
function permute_constraints!(workspace::BBworkspace,indices::Array{Int,1};
                              suppressWarnings::Bool=false,suppressUpdate::Bool=true)::Nothing

    # propagate the change to the subsolver
    permute_constraints!(workspace.subsolverWS,indices,suppressUpdate=suppressUpdate)

    # permute the optimality results
    if length(workspace.activeQueue)>0
        for i in 1:length(workspace.activeQueue)
            permute!(workspace.activeQueue[i].cnsDual,indices)
        end
    end
    if length(workspace.unactivePool)>0
        for i in 1:length(workspace.unactivePool)
            permute!(workspace.unactivePool[i].cnsDual,indices)
        end
    end
    if length(workspace.solutionPool)>0
        for i in 1:length(workspace.solutionPool)
            permute!(workspace.solutionPool[i].cnsDual,indices)
        end
    end


    # update the workspace
    if !suppressUpdate
        update!(workspace.subsolverWS)
    end


    return
end


#
function update_bounds!(workspace::BBworkspace;
                            cnsLoBs::Array{Float64,1}=Float64[],
                            cnsUpBs::Array{Float64,1}=Float64[],
                            varLoBs::Array{Float64,1}=Float64[],
                            varUpBs::Array{Float64,1}=Float64[],
                            suppressWarnings::Bool=false,
                            suppressUpdate::Bool=false
                            )::Nothing

    # check the correctness of the input
    @assert length(cnsLoBs)==0 || length(cnsLoBs)==length(workspace.subsolverWS.cnsLoBs)
    @assert length(cnsUpBs)==0 || length(cnsUpBs)==length(workspace.subsolverWS.cnsUpBs)
    @assert length(varLoBs)==0 || length(varLoBs)==length(workspace.subsolverWS.varLoBs)
    @assert length(varUpBs)==0 || length(varUpBs)==length(workspace.subsolverWS.varUpBs)


    # check if it is possible to make changes
    if !suppressWarnings && workspace.status.description != "new" && !workspace.settings.dynamicMode
        @warn "In order to correctly manipulate the problem formulation, OpenBB must be run in dynamic mode"

        # check if a bounds relaxation was requested
        if (length(cnsLoBs) > 0 && any(@. cnsLoBs < workspace.subsolverWS.cnsLoBs)) ||
           (length(cnsUpBs) > 0 && any(@. cnsUpBs > workspace.subsolverWS.cnsUpBs)) ||
           (length(varLoBs) > 0 && any(@. varLoBs < workspace.subsolverWS.varLoBs)) ||
           (length(varUpBs) > 0 && any(@. varUpBs > workspace.subsolverWS.varLoBs))

            @warn "It is not possible to update the BB tree for a bounds relaxation, a restart might be necessary"
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

    # propagate the changes to the nodes solver
    update_bounds!(workspace.subsolverWS,cnsLoBs,cnsUpBs,varLoBs,varUpBs,suppressUpdate=suppressUpdate)

    # adapt the workspace to the changes
    if !suppressUpdate
        reset_explored_nodes!(workspace)
    end

    return
end

#
function append_problem!(workspace::BBworkspace,problem::Problem;
                         suppressWarnings::Bool=false,
                         suppressUpdate::Bool=false)::Nothing

    # check if it is possible to make changes
    if !suppressWarnings && workspace.status.description != "new" && !workspace.settings.dynamicMode
        @warn "In order to correctly manipulate the problem formulation, OpenBB must be run in dynamic mode"
    end

    # collect info on the problem
    nVars1 = length(workspace.subsolverWS.varLoBs)
    nVars2 = length(problem.varSet.loBs)
    nCnss1 = length(workspace.subsolverWS.cnsLoBs)
    nCnss2 = length(problem.cnsSet.loBs)

    # propagate the changes to the nodes solver
    reliableObjLoBs = append_problem!(workspace.subsolverWS,problem,suppressUpdate=suppressUpdate)



    # update all the sub-problems
    newPrimal = problem.varSet.val
    newPseudoCosts = problem.varSet.pseudoCosts
    for i in  1:length(workspace.activeQueue)
        # mark all the problems as unreliable if necessary
        workspace.activeQueue[i].reliable = reliableObjLoBs && workspace.activeQueue[i].reliable
        # extend primal and dual optimization results
        append!(workspace.activeQueue[i].primal,copy(newPrimal))
        append!(workspace.activeQueue[i].pseudoCosts,copy(newPseudoCosts))
        append!(workspace.activeQueue[i].bndDual,zeros(nVars2))
        append!(workspace.activeQueue[i].cnsDual,zeros(nCnss2))
    end
    for i in  1:length(workspace.solutionPool)
        # mark all the problems as unreliable if necessary
        workspace.solutionPool[i].reliable = reliableObjLoBs && workspace.solutionPool[i].reliable
        # extend primal and dual optimization results
        append!(workspace.solutionPool[i].primal,copy(newPrimal))
        append!(workspace.solutionPool[i].pseudoCosts,copy(newPseudoCosts))
        append!(workspace.solutionPool[i].bndDual,zeros(nVars2))
        append!(workspace.solutionPool[i].cnsDual,zeros(nCnss2))
    end
    for i in  1:length(workspace.unactivePool)
        # mark all the problems as unreliable if necessary
        workspace.unactivePool[i].reliable = reliableObjLoBs && workspace.unactivePool[i].reliable
        # extend primal and dual optimization results
        append!(workspace.unactivePool[i].primal,copy(newPrimal))
        append!(workspace.unactivePool[i].pseudoCosts,copy(newPseudoCosts))
        append!(workspace.unactivePool[i].bndDual,zeros(nVars2))
        append!(workspace.unactivePool[i].cnsDual,zeros(nCnss2))
    end




    # set objective lower bound to -Inf if the requested update
    # invalidates the already computed lower bounds
    if !reliableObjLoBs
        workspace.status.objLoB = -Inf
    end

    # update discrete variables
    append!(workspace.dscIndices,problem.varSet.dscIndices .+ nVars1)
    append!(workspace.sos1Groups,problem.varSet.sos1Groups)

    # adapt the workspace to the chanes
    if !suppressUpdate
        reset_explored_nodes!(workspace)
    end

    return
end


#
function integralize_variables!(workspace::BBworkspace,newDscIndices::Array{Int,1};
                                newSos1Groups=Int[],newPseudoCosts=Float64,
                                suppressWarnings::Bool=false,suppressUpdate::Bool=false)::Nothing

    # check if it is possible to make changes
    if !suppressWarnings && workspace.status.description != "new" && !workspace.settings.dynamicMode
        @warn "In order to correctly manipulate the problem formulation, OpenBB must be run in dynamic mode"
    end

    # check correctness of the inputs
    if length(newSos1Groups) == 0
        newSos1Groups = repeat([-1],length(newDscIndices))
    elseif length(newSos1Groups) != length(newDscIndices)
        @error "newSos1Groups should either be empty or have the same dimension than newDscIndices"
    end

    if length(newPseudoCosts) == 0
        newPseudoCosts = ones(length(newDscIndices))
    elseif length(newPseudoCosts) != length(newDscIndices)
        @error "newPseudoCosts should either be empty or have the same dimension than newDscIndices"
    end



    # add the new discrete variables
    append!(workspace.dscIndices,copy(newDscIndices))
    append!(workspace.sos1Groups,copy(newSos1Groups))
    tmpPerm = sortperm!(workspace.dscIndices)
    permute!(workspace.sos1Groups,tmpPerm)

    # update the nodes pseudoCosts
    for i in 1:length(workspace.activeQueue)
        append!(workspace.activeQueue[i].pseudoCosts,newPseudoCosts)
        permute!(workspace.activeQueue[i].pseudoCosts,tmpPerm)
    end
    for i in 1:length(workspace.solutionPool)
        append!(workspace.solutionPool[i].pseudoCosts,newPseudoCosts)
        permute!(workspace.solutionPool[i].pseudoCosts,tmpPerm)
    end
    for i in 1:length(workspace.unactivePool)
        append!(workspace.unactivePool[i].pseudoCosts,newPseudoCosts)
        permute!(workspace.unactivePool[i].pseudoCosts,tmpPerm)
    end
    


    # adapt the workspace to the changes
    if !suppressUpdate
        reset_explored_nodes!(workspace)
    end

    return
end



#
function update_settings!(workspace::BBworkspace;bb_settings::AbstractSettings=NullSettings(),ss_settings::AbstractSettings=NullSettings())::Nothing


    @error "Not implemented yet"

    # check if it is possible to make changes
    @assert workspace.settings.dynamicMode == true || workspace.status.description == "new"


    # update the BB settings
    if !(bb_settings isa NullSettings)
        if workspace.status.totalTime > 0.

            # adapt the workspace to the new settings


            # pass from dynamic mode to static mode
            if workspace.dynamicMode == true && bb_settings.dynamicMode == false
                for k in length(workspace.solutionPool):-1:1
                    if workspace.solutionPool[k].reliable
                        deleteat!(workspace.solutionPool,k+1:length(workspace.solutionPool))
                        deleteat!(workspace.solutionPool,1:k-1)
                        break
                    elseif k == 1
                        deleteat!(workspace.solutionPool,1:length(workspace.solutionPool)-1)
                    end
                end
                deleteat!(workspace.unactivePool,1:length(workspace.unactivePool))
            end
        end

        for field in fieldnames(BBsettings)
            setfield!(workspace.settings,field,getfield(bb_settings,field))
        end
    end

    # update the subsolver settings
    if !(ss_settings isa NullSettings)
        update_settings!(workspace.subsolverWS,ss_settings,bb_primalTolerance=workspace.settings.primalTolerance,bb_timeLimit=workspace.settings.timeLimit)
    end

    return
end
