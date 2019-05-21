# @Author: Massimo De Mauri <massimo>
# @Date:   2019-05-20T10:04:17+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: update_nodes.jl
# @Last modified by:   massimo
# @Last modified time: 2019-05-21T14:15:28+02:00
# @License: apache 2.0
# @Copyright: {{copyright}}


function update!(workspace::BBworkspace;localOnly::Bool=false)::Nothing

    @sync if !localOnly && workspace.globalInfo != nothing
        # call the local version of the function on the remote workers
        for p in 2:workspace.settings.numProcesses
            @async remotecall_fetch(Main.eval,p,:(OpenBB.update!(workspace,localOnly=true)))
        end

        # call the local version of the function on the current process
        update!(workspace,localOnly=true)

    else
        # update the subsolver workspace
        update!(workspace.subsolverWS)

        # reset the explored sub-problems
        reset_explored_nodes!(workspace,localOnly=true)

        # reset the global info
        if workspace.globalInfo != nothing
            workspace.globalInfo[1] = Inf
            workspace.globalInfo[2] = 0.
            workspace.globalInfo[3] = 0.
        end
    end

    return
end


# eliminates all the generated nodes from the workspace
function clear!(workspace;localOnly::Bool=localOnly)::Nothing

    @sync if !localOnly && workspace.globalInfo != nothing
        # call function on the remote workers
        for p in 2:workspace.settings.numProcesses
            @async remotecall_fetch(Main.eval,p,:(OpenBB.clear!(workspace,localOnly=true)))
        end

        # call function on the main process
        clear!(workspace,localOnly=true)

        # reset the global info
        if workspace.globalInfo != nothing
            workspace.globalInfo[1] = Inf
            workspace.globalInfo[2] = 0.
            workspace.globalInfo[3] = 0.
        end
    else
        deleteat!(workspace.activeQueue, 1:length(workspace.activeQueue))
        deleteat!(workspace.solutionPool,1:length(workspace.solutionPool))
        deleteat!(workspace.unactivePool,1:length(workspace.unactivePool))
        # reset the status
        defaultStatus = BBstatus()
        for field in fieldnames(BBstatus)
            setfield!(workspace.status,field,getfield(defaultStatus,field))
        end
    end

    return
end

# return the workspace to the initial state
function reset!(workspace::BBworkspace;localOnly::Bool=false)::Nothing

    @sync if !localOnly && workspace.globalInfo != nothing
        # remove all nodes in the remote workers
        for p in 2:workspace.settings.numProcesses
            @async remotecall_fetch(Main.eval,p,:(OpenBB.clear!(workspace,localOnly=true)))
        end

        # call the local version of the function on the main process
        reset!(workspace,localOnly=true)

        # reset the global info
        if workspace.globalInfo != nothing
            workspace.globalInfo[1] = Inf
            workspace.globalInfo[2] = 0.
            workspace.globalInfo[3] = 0.
        end
    else
        # eliminate all the generated nodes and reinsert the root of the BB tree
        clear!(workspace,localOnly=true)
        push!(workspace.activeQueue,BBnode(Dict{Int,Float64}(),Dict{Int,Float64}(),
                                             zeros(get_numVariables(workspace)),
                                             zeros(get_numVariables(workspace)),
                                             zeros(get_numConstraints(workspace)),
                                             1.0,-Inf,false))
    end

    return
end

# function reset!(;localOnly::Bool=false)::Nothing
#     return reset!(workspace,localOnly=localOnly)
# end


#
function reset_explored_nodes!(workspace::BBworkspace;localOnly::Bool=false)::Nothing

    @sync if !localOnly && workspace.globalInfo != nothing
        # call the local version of the function on the remote workers
        for p in 2:workspace.settings.numProcesses
            @async remotecall_fetch(Main.eval,p,:(OpenBB.reset_explored_nodes!(workspace,localOnly=true)))
        end
        # call the local version of the function on the main process
        reset_explored_nodes!(workspace,localOnly=true)

        # reset the global info
        if workspace.globalInfo != nothing
            workspace.globalInfo[1] = Inf
            workspace.globalInfo[2] = 0.
            workspace.globalInfo[3] = 0.
        end

    else

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
    end

    return
end
