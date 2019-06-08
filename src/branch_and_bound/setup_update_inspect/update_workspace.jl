# @Author: Massimo De Mauri <massimo>
# @Date:   2019-05-20T10:04:17+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: update_nodes.jl
# @Last modified by:   massimo
# @Last modified time: 2019-06-06T14:22:53+02:00
# @License: apache 2.0
# @Copyright: {{copyright}}


#
function reset_explored_nodes!(workspace::BBworkspace{T1,T2};localOnly::Bool=false)::Nothing where T1<:AbstractWorkspace where T2<:AbstractSharedMemory

    @sync if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)

        # call the local version of the function on the remote workers
        for p in 2:workspace.settings.numProcesses
            @async remotecall_fetch(Main.eval,p,:(OpenBB.reset_explored_nodes!(workspace,localOnly=true)))
        end

        # call the local version of the function on the main process
        reset_explored_nodes!(workspace,localOnly=true)

        # reset the global info
        workspace.sharedMemory.objectiveBounds[end] = Inf
        @. workspace.sharedMemory.stats = 0
		@. workspace.sharedMemory.arrestable = false

    else


        # adapt the node pools to the changes (the solutions go on top of all)
        append!(workspace.activeQueue,workspace.unactivePool)
		append!(workspace.activeQueue,workspace.solutionPool)
        sort!(workspace.activeQueue,alg=MergeSort,rev=true,
              lt=(l,r)->workspace.settings.expansionPriorityRule(l,r,workspace.status))

        deleteat!(workspace.solutionPool,1:length(workspace.solutionPool))
        deleteat!(workspace.unactivePool,1:length(workspace.unactivePool))

        # adapt the status to the changes
        workspace.status.objUpB = Inf
        workspace.status.absoluteGap = Inf
        workspace.status.relativeGap = Inf
        workspace.status.description = "interrupted"

    end

    return
end



#
function update_sharedMemory(workspace::BBworkspace{T1,T2})::Nothing where T1<:AbstractWorkspace where T2<:AbstractSharedMemory

	if workspace.sharedMemory isa BBsharedMemory{BBnodeChannel}
		numVars = get_numVariables(workspace)
		numCnss = get_numConstraints(workspace)
		numDscVars = get_numDiscreteVariables(workspace)

		# construct new communication Channels
		communicationChannels = Array{BBnodeChannel,1}(undef,workspace.settings.numProcesses)
		for k in 1:workspace.settings.numProcesses
			communicationChannels[k] = BBnodeChannel(flat_size(numVars,numDscVars,numCnss))
		end
		@sync for k in 2:workspace.settings.numProcesses
			@async if k < workspace.settings.numProcesses
				remotecall_fetch(Main.eval,k,:(workspace.sharedMemory.inputChannel = $(communicationChannels[k]);
											   workspace.sharedMemory.outputChannel = $(communicationChannels[k+1]);
											   nothing))
			else
				remotecall_fetch(Main.eval,k,:(workspace.sharedMemory.inputChannel = $(communicationChannels[k]);
											   workspace.sharedMemory.outputChannel = $(communicationChannels[1]);
											   nothing))
			end
		end
	end
	return
end


#
function update!(workspace::BBworkspace{T1,T2};localOnly::Bool=false)::Nothing where T1<:AbstractWorkspace where T2<:AbstractSharedMemory

    @sync if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)

        # call the local version of the function on the remote workers
        for p in 2:workspace.settings.numProcesses
            @async remotecall_fetch(Main.eval,p,:(OpenBB.update!(workspace,localOnly=true)))
        end

        # call the local version of the function on the current process
        update!(workspace,localOnly=true)

		# adapt the shared memory to the new problem
		update_sharedMemory(workspace)

        # reset the global info
        workspace.sharedMemory.objectiveBounds[end] = Inf
        @. workspace.sharedMemory.stats = 0
		@. workspace.sharedMemory.arrestable = false


        # update the communication channels (if needed)
		if workspace.sharedMemory.inputChannel isa BBnodeChannel

			workersList = workers()
	        communicationChannels = Array{BBnodeChannel,1}(undef,workspace.settings.numProcesses)
			for k in 1:workspace.settings.numProcesses
				communicationChannels[k] = BBnodeChannel(flat_size(get_numVariables(workspace),
																   get_numDiscreteVariables(workspace),
																   get_numConstraints(workspace)))
			end
			workspace.sharedMemory.inputChannel = communicationChannels[1]
			workspace.sharedMemory.outputChannel = communicationChannels[2]
			@sync for k in 2:workspace.settings.numProcesses
				if k < workspace.settings.numProcesses
					@async remotecall_fetch(Main.eval,workersList[k-1],:(workspace.sharedMemory.inputChannel = $(communicationChannels[k]);
																	     workspace.sharedMemory.outputChannel = $(communicationChannels[k+1])))
				else
					@async remotecall_fetch(Main.eval,workersList[k-1],:(workspace.sharedMemory.inputChannel = $(communicationChannels[k]);
																	     workspace.sharedMemory.outputChannel = $(communicationChannels[1])))
				end
			end
		end


    else

        # update the subsolver workspace
        update!(workspace.subsolverWS)

        # reset the explored sub-problems
        reset_explored_nodes!(workspace,localOnly=true)
    end

    return
end


# return the workspace to the initial state
function reset!(workspace::BBworkspace{T1,T2};localOnly::Bool=false)::Nothing where T1<:AbstractWorkspace where T2<:AbstractSharedMemory

    @sync if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)

        # remove all nodes in the remote workers
        for p in 2:workspace.settings.numProcesses
            @async remotecall_fetch(Main.eval,p,:(OpenBB.clear!(workspace,localOnly=true)))
        end

        # call the local version of the function on the main process
        reset!(workspace,localOnly=true)

        # reset the global info
        @. workspace.sharedMemory.objectiveBounds[1:end-1] = -Inf
        workspace.sharedMemory.objectiveBounds[end] = Inf
        @. workspace.sharedMemory.stats = 0
		@. workspace.sharedMemory.arrestable = false

    else
        # eliminate all the generated nodes and reinsert the root of the BB tree
        clear!(workspace,localOnly=true)
        push!(workspace.activeQueue,BBnode(Dict{Int,Float64}(),Dict{Int,Float64}(),
                                             zeros(get_numDiscreteVariables(workspace)),
                                             zeros(get_numVariables(workspace)),
											 zeros(get_numVariables(workspace)),
                                             zeros(get_numConstraints(workspace)),
                                             1.0,-Inf,false))
    end

    return
end


# eliminates all the generated nodes from the workspace
function clear!(workspace::BBworkspace{T1,T2};localOnly::Bool=false)::Nothing where T1<:AbstractWorkspace where T2<:AbstractSharedMemory

    @sync if !localOnly && !(workspace.sharedMemory isa NullSharedMemory)

        # call function on the remote workers
        for p in 2:workspace.settings.numProcesses
            @async remotecall_fetch(Main.eval,p,:(OpenBB.clear!(workspace,localOnly=true)))
        end

        # call function on the main process
        clear!(workspace,localOnly=true)

        # update the global info
        @. workspace.sharedMemory.objectiveBounds[1:end-1] = Inf
        workspace.sharedMemory.objectiveBounds[end] = Inf
        @. workspace.sharedMemory.stats = 0
		@. workspace.sharedMemory.arrestable = false

    else
        deleteat!(workspace.activeQueue, 1:length(workspace.activeQueue))
        deleteat!(workspace.solutionPool,1:length(workspace.solutionPool))
        deleteat!(workspace.unactivePool,1:length(workspace.unactivePool))
        # reset the status
        defaultStatus = BBstatus()
        for field in fieldnames(BBstatus)
            setfield!(workspace.status,field,getfield(defaultStatus,field))
        end
		workspace.status.objLoB = -Inf
    end

    return
end
