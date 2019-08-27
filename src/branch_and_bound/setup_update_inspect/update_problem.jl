# @Author: Massimo De Mauri <massimo>
# @Date:   2019-02-19T12:14:43+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: update_problem.jl
# @Last modified by:   massimo
# @Last modified time: 2019-08-27T15:31:36+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

#
function append_constraints!(workspace::BBworkspace{T1,T2},
    						 constraintSet::T3;
                             suppressWarnings::Bool=false,
                             suppressUpdate::Bool=false,
                             localOnly::Bool=false)::Nothing where T1<:AbstractWorkspace where T2<:AbstractSharedMemory where T3 <: AbstractConstraintSet

    return insert_constraints!(workspace,constraintSet,
                               get_numConstraints(workspace)+1,
                               suppressWarnings=suppressWarnings,
                               suppressUpdate=suppressUpdate,
                               localOnly=localOnly)
end

#
function insert_constraints!(workspace::BBworkspace{T1,T2},
                             constraintSet::T3,index::Int;
                             suppressWarnings::Bool=false,
                             suppressUpdate::Bool=false,
                             localOnly::Bool=false)::Nothing where T1<:AbstractWorkspace where T2<:AbstractSharedMemory where T3 <: AbstractConstraintSet

    # call the same function on the other workers
    @sync if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)
        # call the local version of the function on the remote workers
        for p in 2:workspace.settings.numProcesses
            @async remotecall_fetch(Main.eval,p,:(OpenBB.insert_constraints!(workspace,$constraintSet,$index,
                                                                             suppressWarnings=$suppressWarnings,
                                                                             suppressUpdate=$suppressUpdate,localOnly=true)))
        end

        # call the local version of the function on the current process
        insert_constraints!(workspace,constraintSet,index,
                            suppressWarnings=suppressWarnings,
                            suppressUpdate=suppressUpdate,localOnly=true)

		# reset the information shared among the processes
		reset_global_info!(workspace)

    else

		# check if it is possible to make changes
		if myid() == 1 && !suppressWarnings && workspace.status.description != "new" && !workspace.settings.interactiveMode && myid() == 1
			@warn "In order to correctly manipulate the problem formulation, OpenBB must be run in interactive mode"
		end


        # propagate the changes to the subsolver
        insert_constraints!(workspace.subsolverWS,constraintSet,index,suppressUpdate=true)

        # update all the problems in the activeQueue
		numConstraints = get_size(constraintSet)
		newBounds = get_bounds(constraintSet)
        if length(workspace.activeQueue)>0
            for i in 1:length(workspace.activeQueue)
                splice!(workspace.activeQueue[i].cnsDual,index:index-1,zeros(numConstraints))
                splice!(workspace.activeQueue[i].cnsLoBs,index:index-1,newBounds[1])
				splice!(workspace.activeQueue[i].cnsUpBs,index:index-1,newBounds[2])
            end
        end
        if length(workspace.unactivePool)>0
            for i in  1:length(workspace.unactivePool)
                splice!(workspace.unactivePool[i].cnsDual,index:index-1,zeros(numConstraints))
				splice!(workspace.unactivePool[i].cnsLoBs,index:index-1,newBounds[1])
				splice!(workspace.unactivePool[i].cnsUpBs,index:index-1,newBounds[2])
            end
        end
        if length(workspace.solutionPool)>0
            for i in  1:length(workspace.solutionPool)
                splice!(workspace.solutionPool[i].cnsDual,index:index-1,zeros(numConstraints))
				splice!(workspace.solutionPool[i].cnsLoBs,index:index-1,newBounds[1])
				splice!(workspace.solutionPool[i].cnsUpBs,index:index-1,newBounds[2])
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


    @sync if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)
        # call the local version of the function on the remote workers
        for p in 2:workspace.settings.numProcesses
            @async remotecall_fetch(Main.eval,p,:(OpenBB.remove_constraints!(workspace,$indices,
                                                                             suppressWarnings=true,
                                                                             suppressUpdate=$suppressUpdate,
																			 localOnly=true)))
        end

        # call the local version of the function on the main process
        remove_constraints!(workspace,indices,
                            suppressWarnings=true,
                            suppressUpdate=suppressUpdate,
							localOnly=true)


		# reset the information shared among the processes
		reset_global_info!(workspace)

    else

		## check if it is possible to make changes
		if myid() == 1 && !suppressWarnings && workspace.status.description != "new"
			@warn "Removing constraints after some iterations is potentially destructive, I hope you know what you are doing."
		end

        # propagate the changes to the nodes solver
        remove_constraints!(workspace.subsolverWS,indices,suppressUpdate=true)

        # update all the problems in the activeQueue
        for i in 1:length(workspace.activeQueue)
            deleteat!(workspace.activeQueue[i].cnsDual,indices)
			deleteat!(workspace.activeQueue[i].cnsLoBs,indices)
			deleteat!(workspace.activeQueue[i].cnsUpBs,indices)
        end
        # update all the problems in the solutionPool
        for i in 1:length(workspace.solutionPool)
            deleteat!(workspace.solutionPool[i].cnsDual,indices)
			deleteat!(workspace.activeQueue[i].cnsLoBs,indices)
			deleteat!(workspace.activeQueue[i].cnsUpBs,indices)
        end
        # update all the problems in the activeQueue
        for i in 1:length(workspace.unactivePool)
            deleteat!(workspace.unactivePool[i].cnsDual,indices)
			deleteat!(workspace.activeQueue[i].cnsLoBs,indices)
			deleteat!(workspace.activeQueue[i].cnsUpBs,indices)
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
                                                                              suppressWarnings=true,
                                                                              suppressUpdate=$suppressUpdate,
																			  localOnly=true)))
        end

        # call the local version of the function on the main process
        permute_constraints!(workspace,permutation,
                             suppressWarnings=true,
                             suppressUpdate=suppressUpdate,
							 localOnly=true)

 		# reset the information shared among the processes
 		reset_global_info!(workspace)

    else

        # propagate the change to the subsolver
        permute_constraints!(workspace.subsolverWS,permutation,suppressUpdate=true)

        # permute the optimality results
        if length(workspace.activeQueue)>0
            for i in 1:length(workspace.activeQueue)
                permute!(workspace.activeQueue[i].cnsDual,permutation)
				permute!(workspace.activeQueue[i].cnsLoBs,permutation)
				permute!(workspace.activeQueue[i].cnsUpBs,permutation)
            end
        end
        if length(workspace.unactivePool)>0
            for i in 1:length(workspace.unactivePool)
                permute!(workspace.unactivePool[i].cnsDual,permutation)
				permute!(workspace.unactivePool[i].cnsLoBs,permutation)
				permute!(workspace.unactivePool[i].cnsUpBs,permutation)
            end
        end
        if length(workspace.solutionPool)>0
            for i in 1:length(workspace.solutionPool)
                permute!(workspace.solutionPool[i].cnsDual,permutation)
				permute!(workspace.solutionPool[i].cnsLoBs,permutation)
				permute!(workspace.solutionPool[i].cnsUpBs,permutation)
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
						varLoBs::Array{Float64,1}=Array{Float64,1}(),
						varUpBs::Array{Float64,1}=Array{Float64,1}(),
                        cnsLoBs::Array{Float64,1}=Array{Float64,1}(),
                        cnsUpBs::Array{Float64,1}=Array{Float64,1}(),
                        suppressWarnings::Bool=false,
                        suppressUpdate::Bool=false,
                        localOnly::Bool=false)::Nothing where T1<:AbstractWorkspace where T2<:AbstractSharedMemory


    # ensure the correctness of the input
	@assert length(varLoBs)==length(workspace.subsolverWS.varLoBs) || length(varLoBs)==0
	@assert length(varUpBs)==length(workspace.subsolverWS.varUpBs) || length(varUpBs)==0
	@assert length(cnsLoBs)==length(workspace.subsolverWS.cnsLoBs) || length(cnsLoBs)==0
	@assert length(cnsUpBs)==length(workspace.subsolverWS.cnsUpBs) || length(cnsUpBs)==0


    @sync if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)
        # call the local version of the function on the remote workers
        for p in 2:workspace.settings.numProcesses
            @async remotecall_fetch(Main.eval,p,:(OpenBB.update_bounds!(workspace,
                                                                        cnsLoBs=$cnsLoBs,cnsUpBs=$cnsUpBs,
                                                                        varLoBs=$varLoBs,varUpBs=$varUpBs,
                                                                        suppressWarnings=true,
                                                                        suppressUpdate=$suppressUpdate,
																		localOnly=true)))
        end

        # call the local version of the function on the main process
        update_bounds!(workspace,
                       cnsLoBs=cnsLoBs,cnsUpBs=cnsUpBs,
                       varLoBs=varLoBs,varUpBs=varUpBs,
                       suppressWarnings=true,
                       suppressUpdate=suppressUpdate,
					   localOnly=true)


		# reset the information shared among the processes
		reset_global_info!(workspace)
    else

		# check if it is possible to make changes
		if myid() == 1 && !suppressWarnings && workspace.status.description != "new"
			if !workspace.settings.interactiveMode
				@warn "In order to correctly manipulate the problem formulation, OpenBB must be run in interactive mode"
			# check if a bounds relaxation was requested
			elseif (length(cnsLoBs) > 0 && any(@. cnsLoBs < workspace.subsolverWS.cnsLoBs)) ||
			   	   (length(cnsUpBs) > 0 && any(@. cnsUpBs > workspace.subsolverWS.cnsUpBs)) ||
			   	   (length(varLoBs) > 0 && any(@. varLoBs < workspace.subsolverWS.varLoBs)) ||
			       (length(varUpBs) > 0 && any(@. varUpBs > workspace.subsolverWS.varLoBs))

				@warn "Relaxing the variable bounds after some iterations is potentially destructive, please make sure that the new bound set is more restrictive than the old one."
			end
		end


        # update the bounds into the BB tree
        if length(workspace.activeQueue)>0
            for i in 1:length(workspace.activeQueue)
				if length(varLoBs)>0 @. workspace.activeQueue[i].varLoBs = max(workspace.activeQueue[i].varLoBs,varLoBs) end
				if length(varUpBs)>0 @. workspace.activeQueue[i].varUpBs = min(workspace.activeQueue[i].varUpBs,varUpBs) end
				if length(cnsLoBs)>0 @. workspace.activeQueue[i].cnsLoBs = max(workspace.activeQueue[i].cnsLoBs,cnsLoBs) end
				if length(cnsUpBs)>0 @. workspace.activeQueue[i].cnsUpBs = min(workspace.activeQueue[i].cnsUpBs,cnsUpBs) end
            end
        end
        if length(workspace.unactivePool)>0
            for i in 1:length(workspace.unactivePool)
				if length(varLoBs)>0 @. workspace.unactivePool[i].varLoBs = max(workspace.unactivePool[i].varLoBs,varLoBs) end
				if length(varUpBs)>0 @. workspace.unactivePool[i].varUpBs = min(workspace.unactivePool[i].varUpBs,varUpBs) end
				if length(cnsLoBs)>0 @. workspace.unactivePool[i].cnsLoBs = max(workspace.unactivePool[i].cnsLoBs,cnsLoBs) end
				if length(cnsUpBs)>0 @. workspace.unactivePool[i].cnsUpBs = min(workspace.unactivePool[i].cnsUpBs,cnsUpBs) end
            end
        end
        if length(workspace.solutionPool)>0
            for i in 1:length(workspace.solutionPool)
				if length(varLoBs)>0 @. workspace.solutionPool[i].varLoBs = max(workspace.solutionPool[i].varLoBs,varLoBs) end
				if length(varUpBs)>0 @. workspace.solutionPool[i].varUpBs = min(workspace.solutionPool[i].varUpBs,varUpBs) end
				if length(cnsLoBs)>0 @. workspace.solutionPool[i].cnsLoBs = max(workspace.solutionPool[i].cnsLoBs,cnsLoBs) end
				if length(cnsUpBs)>0 @. workspace.solutionPool[i].cnsUpBs = min(workspace.solutionPool[i].cnsUpBs,cnsUpBs) end
            end
        end

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

    @sync if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)
        # call the local version of the function on the remote workers
        for p in 2:workspace.settings.numProcesses
            @async remotecall_fetch(Main.eval,p,:(OpenBB.set_objective!(workspace,$newObjective,
                                                                         suppressWarnings=true,
                                                                         suppressUpdate=$suppressUpdate,
																		 localOnly=true)))
        end

        # call the local version of the function on the main process
        set_objective!(workspace,newObjective,
                        suppressWarnings=true,
                        suppressUpdate=suppressUpdate,
						localOnly=true)

		# reset the information shared among the processes
		reset_global_info!(workspace)
    else


		if myid() == 1 && !suppressWarnings && workspace.status.description != "new"
			if !workspace.settings.interactiveMode
				@warn "In order to correctly manipulate the problem formulation, OpenBB must be run in interactive mode"
			else
				# warn the user about the potential dangers
				@warn "changing the objective after some iterations is potentially destructive, please be sure that the new objective is always greater or equal to the old one"
			end
		end

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


    @sync if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)
        # call the local version of the function on the remote workers
        for p in 2:workspace.settings.numProcesses
            @async remotecall_fetch(Main.eval,p,:(OpenBB.set_constraintSet!(workspace,$newConstraintSet,
                                                                            suppressWarnings=true,
                                                                            suppressUpdate=$suppressUpdate,
																			localOnly=true)))
        end

        # call the local version of the function on the main process
        set_constraintSet!(workspace,newConstraintSet,
                           suppressWarnings=true,
						   suppressUpdate=suppressUpdate,
						   localOnly=true)

		# reset the information shared among the processes
		reset_global_info!(workspace)

    else
		# check if it is possible to make changes
		if myid() == 1 && !suppressWarnings && workspace.status.description != "new"
			if !workspace.settings.interactiveMode
				@warn "In order to correctly manipulate the problem formulation, OpenBB must be run in interactive mode"
			else
				# warn the user about the potential dangers
				@warn "Changing the constraintSet after some iterations is potentially destructive, please be sure that the new constraint set is more restrictive than the old one"
			end
		end

        # reset the constraints dual and bounds
		newBounds = get_bounds(newConstraintSet)
        if length(workspace.activeQueue)>0
            for i in 1:length(workspace.activeQueue)
                workspace.activeQueue[i].cnsDual = zeros(get_numConstraints(newConstraintSet))
				workspace.activeQueue[i].cnsLoBs = copy(newBounds[1])
				workspace.activeQueue[i].cnsUpBs = copy(newBounds[2])
            end
        end
        if length(workspace.unactivePool)>0
            for i in 1:length(workspace.unactivePool)
                workspace.unactivePool[i].cnsDual = zeros(get_numConstraints(newConstraintSet))
				workspace.unactivePool[i].cnsLoBs = copy(newBounds[1])
				workspace.unactivePool[i].cnsUpBs = copy(newBounds[2])
            end
        end
        if length(workspace.solutionPool)>0
            for i in 1:length(workspace.solutionPool)
                workspace.solutionPool[i].cnsDual = zeros(get_numConstraints(newConstraintSet))
				workspace.solutionPool[i].cnsLoBs = copy(newBounds[1])
				workspace.solutionPool[i].cnsUpBs = copy(newBounds[2])
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


    @sync if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)
        # call the local version of the function on the remote workers
        for p in 2:workspace.settings.numProcesses
            @async remotecall_fetch(Main.eval,p,:(OpenBB.append_problem!(workspace,$problem,
                                                                         suppressWarnings=$suppressWarnings,
                                                                         suppressUpdate=$suppressUpdate,localOnly=true)))
        end
        # call the local version of the function on the main process
        append_problem!(workspace,problem,
                        suppressWarnings=suppressWarnings,
                        suppressUpdate=suppressUpdate,localOnly=true)

		# adapt the shared memory to the new problem
		update_sharedMemory!(workspace)

		# reset the information shared among the processes
		reset_global_info!(workspace)
    else

		# check if it is possible to make changes
		if myid() == 1 && !suppressWarnings && workspace.status.description != "new" && !workspace.settings.interactiveMode && myid() == 1
			@warn "In order to correctly manipulate the problem formulation, OpenBB must be run in interactive mode"
		end

        # collect info on the problem
        nVars1 = length(workspace.subsolverWS.varLoBs)
        nVars2 = length(problem.varSet.loBs)
        nCnss1 = length(workspace.subsolverWS.cnsLoBs)
        nCnss2 = length(problem.cnsSet.loBs)
        # propagate the changes to the subsolver
        append_problem!(workspace.subsolverWS,problem,suppressUpdate=true)

        # update all the sub-problems
        newPrimal = problem.varSet.vals
		newVarBounds = get_bounds(problem.varSet)
		newCnsBounds = get_bounds(problem.cnsSet)
        for i in  1:length(workspace.activeQueue)
            # extend branching bounds
            append!(workspace.activeQueue[i].varLoBs,newVarBounds[1])
            append!(workspace.activeQueue[i].varUpBs,newVarBounds[2])
			append!(workspace.activeQueue[i].cnsLoBs,newCnsBounds[1])
            append!(workspace.activeQueue[i].cnsUpBs,newCnsBounds[2])
            # extend primal and dual optimization results
            append!(workspace.activeQueue[i].primal,newPrimal)
            append!(workspace.activeQueue[i].bndDual,zeros(nVars2))
            append!(workspace.activeQueue[i].cnsDual,zeros(nCnss2))
        end
        for i in  1:length(workspace.solutionPool)
            # extend branching bounds
			append!(workspace.solutionPool[i].varLoBs,newVarBounds[1])
            append!(workspace.solutionPool[i].varUpBs,newVarBounds[2])
			append!(workspace.solutionPool[i].cnsLoBs,newCnsBounds[1])
            append!(workspace.solutionPool[i].cnsUpBs,newCnsBounds[2])
            # extend primal and dual optimization results
            append!(workspace.solutionPool[i].primal,newPrimal)
            append!(workspace.solutionPool[i].bndDual,zeros(nVars2))
            append!(workspace.solutionPool[i].cnsDual,zeros(nCnss2))
        end
        for i in  1:length(workspace.unactivePool)
            # extend branching bounds
			append!(workspace.unactivePool[i].varLoBs,newVarBounds[1])
            append!(workspace.unactivePool[i].varUpBs,newVarBounds[2])
			append!(workspace.unactivePool[i].cnsLoBs,newCnsBounds[1])
            append!(workspace.unactivePool[i].cnsUpBs,newCnsBounds[2])
            # extend primal and dual optimization results
            append!(workspace.unactivePool[i].primal,newPrimal)
            append!(workspace.unactivePool[i].bndDual,zeros(nVars2))
            append!(workspace.unactivePool[i].cnsDual,zeros(nCnss2))
        end

        # update discrete variables
        append!(workspace.dscIndices,@. problem.varSet.dscIndices + nVars1)
        append!(workspace.sos1Groups,@. problem.varSet.sos1Groups + maximum(workspace.sos1Groups) + 1)
        workspace.pseudoCosts = (vcat(workspace.pseudoCosts[1],problem.varSet.pseudoCosts[1]),
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
                                                                                suppressWarnings=true,
                                                                                suppressUpdate=$suppressUpdate,
																				localOnly=true)))
        end
        # call the local version of the function on the main process
        integralize_variables!(workspace,newDscIndices,
                               newSos1Groups=newSos1Groups,
                               suppressWarnings=true,
                               suppressUpdate=suppressUpdate,
							   localOnly=true)

	    # adapt the shared memory to the new problem
		update_sharedMemory!(workspace)

		# reset the information shared among the processes
		reset_global_info!(workspace)

    else
        # check if it is possible to make changes
        if myid() == 1 && !suppressWarnings && workspace.status.description != "new" && !workspace.settings.interactiveMode && myid() == 1
            @warn "In order to correctly manipulate the problem formulation, OpenBB must be run in interactive mode"
        end
        # check correctness of the inputs
        if length(newSos1Groups) == 0
            newSos1Groups = repeat([0],length(newDscIndices))
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
			@. workspace.activeQueue[k].varLoBs =  ceil(workspace.activeQueue[k].varLoBs-workspace.settings.primalTolerance)
			@. workspace.activeQueue[k].varUpBs = floor(workspace.activeQueue[k].varUpBs+workspace.settings.primalTolerance)
        end
        for k in 1:length(workspace.unactivePool)
			@. workspace.unactivePool[k].varLoBs =  ceil(workspace.unactivePool[k].varLoBs-workspace.settings.primalTolerance)
			@. workspace.unactivePool[k].varUpBs = floor(workspace.unactivePool[k].varUpBs+workspace.settings.primalTolerance)
        end
        for k in 1:length(workspace.solutionPool)
			@. workspace.solutionPool[k].varLoBs =  ceil(workspace.solutionPool[k].varLoBs-workspace.settings.primalTolerance)
			@. workspace.solutionPool[k].varUpBs = floor(workspace.solutionPool[k].varUpBs+workspace.settings.primalTolerance)
        end

        # adapt the workspace to the changes
        if !suppressUpdate
            update!(workspace,localOnly=true)
        end
    end

    return
end
