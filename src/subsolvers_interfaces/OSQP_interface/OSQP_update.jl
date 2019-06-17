# @Author: Massimo De Mauri <massimo>
# @Date:   2019-02-19T20:26:37+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: OSQP_interface_update.jl
# @Last modified by:   massimo
# @Last modified time: 2019-06-17T15:36:37+02:00
# @License: apache 2.0
# @Copyright: {{copyright}}


#
function update!(workspace::OSQPworkspace)::Nothing

    # reformat the settings for OSQP
    settings_dict = Dict{Symbol,Any}()
    for field in fieldnames(OSQPsettings)
        settings_dict[field] = getfield(workspace.settings,field)
    end

    # setup OSQP with the new constraint set
    OSQP.setup!(workspace.model;
                    P=sparse(workspace.Q),
                    q=workspace.L,
                    A=vcat(speye(size(workspace.A,2)),sparse(workspace.A)),
                    l=vcat(workspace.varLoBs,workspace.cnsLoBs),
                    u=vcat(workspace.varUpBs,workspace.cnsUpBs),
                    settings_dict...)
    return
end


#
function insert_constraints!(workspace::OSQPworkspace,
                            A::Union{Array{Float64,2},SparseMatrixCSC{Float64}},
                            cnsLoBs::Array{Float64,1},
                            cnsUpBs::Array{Float64,1},
                            index::Int;
                            suppressUpdate::Bool=false)::Nothing


    # check correctness of the input
    @assert size(A,2) == size(workspace.A,2)
    @assert size(A,1) == length(cnsLoBs)
    @assert size(A,1) == length(cnsUpBs)

    # perform the insertion
    workspace.A = vcat(workspace.A[1:index-1,:],A,workspace.A[index:end,:])
    splice!(workspace.cnsLoBs,index:index-1,cnsLoBs)
    splice!(workspace.cnsUpBs,index:index-1,cnsUpBs)

    # update the osqp workspace
    if !suppressUpdate
        update!(workspace)
    end

end

#
function remove_constraints!(workspace::OSQPworkspace,indices::Array{Int,1};
                            suppressUpdate::Bool=false)::Nothing



    # perform the removal
    indicesCnssToKeep = filter(x->!(x in indices),collect(1:size(workspace.A,1)))
    workspace.A = workspace.A[indicesCnssToKeep,:]
    deleteat!(workspace.cnsLoBs,indices)
    deleteat!(workspace.cnsUpBs,indices)

    # update the GUROBI workspace
    if !suppressUpdate
        update!(workspace)
    end

    return
end

#
function permute_constraints!(workspace::OSQPworkspace,indices::Array{Int,1};suppressUpdate::Bool=false)::Nothing

    # sort the constraints in memory
    workspace.A = workspace.A[indices,:]
    permute!(workspace.cnsLoBs,indices)
    permute!(workspace.cnsUpBs,indices)

    # update the osqp workspace
    if !suppressUpdate
        update!(workspace)
    end

    return
end


#
function update_bounds!(workspace::OSQPworkspace,
                        cnsLoBs::Array{Float64,1},cnsUpBs::Array{Float64,1},
                        varLoBs::Array{Float64,1},varUpBs::Array{Float64,1};
                        suppressUpdate::Bool=false)::Nothing


    @. workspace.cnsLoBs = cnsLoBs
    @. workspace.cnsUpBs = cnsUpBs
    @. workspace.varLoBs = varLoBs
    @. workspace.varUpBs = varUpBs

    # propagate the changes to the OSQP solver
    if !suppressUpdate
        OSQP.update!(workspace.model;l=vcat(workspace.varLoBs,workspace.cnsLoBs),u=vcat(workspace.varUpBs,workspace.cnsLoBs))
    end

    return
end


#
function append_problem!(workspace::OSQPworkspace,
                         problem::Problem{LinearObjective,LinearConstraintSet{T}};
                         suppressUpdate::Bool=false)::Bool where T

    # test the future validity of the already computed lower bounds
    testWorkspace = setup(problem,workspace.settings)
    testSolution = solve!(testWorkspace)
    if testSolution[1] < -workspace.settings.eps_prim_inf
        reliableObjLoBs = false
    else
        reliableObjLoBs = true
    end


    append!(workspace.L, problem.objFun.L)


    workspace.A = vcat(hcat( workspace.A,                                         zeros(size(workspace.A,1),size(problem.cnsSet.A,2))    ),
                       hcat( zeros(size(problem.cnsSet.A,1),size(workspace.A,2)), problem.cnsSet.A                                       ))
    append!(workspace.cnsLoBs, problem.cnsSet.loBs)
    append!(workspace.cnsUpBs, problem.cnsSet.upBs)

    append!(workspace.varLoBs, problem.varSet.loBs)
    append!(workspace.varUpBs, problem.varSet.upBs)

    # update the osqp workspace
    if !suppressUpdate
        update!(workspace)
    end

    return reliableObjLoBs
end

#
function append_problem!(workspace::OSQPworkspace,
                         problem::Problem{QuadraticObjective{T1},LinearConstraintSet{T2}};
                         suppressUpdate::Bool=false)::Bool where T1 where T2

    # test the future validity of the already computed lower bounds
    testWorkspace = setup(problem,workspace.settings)
    testSolution = solve!(testWorkspace)
    if testSolution[1] < -workspace.settings.eps_prim_inf
        reliableObjLoBs = false
    else
        reliableObjLoBs = true
    end


    # update the node data
    workspace.Q = vcat(hcat( workspace.Q,                                         zeros(size(workspace.Q,1),size(problem.objFun.Q,2))    ),
                       hcat( zeros(size(problem.objFun.Q,1),size(workspace.Q,2)), problem.objFun.Q                                       ))

    append!(workspace.L, problem.objFun.L)


    workspace.A = vcat(hcat( workspace.A,                                         zeros(size(workspace.A,1),size(problem.cnsSet.A,2))    ),
                       hcat( zeros(size(problem.cnsSet.A,1),size(workspace.A,2)), problem.cnsSet.A                                       ))
    append!(workspace.cnsLoBs, problem.cnsSet.loBs)
    append!(workspace.cnsUpBs, problem.cnsSet.upBs)

    append!(workspace.varLoBs, problem.varSet.loBs)
    append!(workspace.varUpBs, problem.varSet.upBs)

    # update the osqp workspace
    if !suppressUpdate
        update!(workspace)
    end

    return reliableObjLoBs
end


# this function updates the settings in the osqp model
function update_settings!(workspace::OSQPworkspace,settings::OSQPsettings;bb_primalTolerance::Float64=Inf,bb_timeLimit=Inf)

    # reformat the settings for OSQP
    settings_dict = Dict{Symbol,Any}()
    for field in fieldnames(OSQPsettings)
        settings_dict[field] = getfield(settings,field)
    end

    # overwrite the osqp setting depending on the branch and bound settings
    settings_dict[:eps_prim_inf] = min(settings_dict[:eps_prim_inf],bb_primalTolerance)

    if bb_timeLimit < Inf
        if settings_dict[:timeLimit] == 0.
            settings_dict[:timeLimit] = bb_timeLimit
        else
            settings_dict[:timeLimit] = min(settings_dict[:timeLimit],bb_timeLimit)
        end
    end


    #TODO: some settings can not be changed without the need for a new setup phase
    OSQP.setup!(workspace.model;P=workspace.Q,
                                q=workspace.L,
                                A=vcat(speye(length(workspace.varLoBs)),sparse(workspace.A)),
                                l=vcat(workspace.varLoBs,workspace.cnsLoBs),
                                u=vcat(workspace.varUpBs,workspace.cnsUpBs),
                                settings_dict...)
    return
end
