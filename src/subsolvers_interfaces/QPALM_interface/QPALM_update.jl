# @Author: Massimo De Mauri <massimo>
# @Date:   2019-02-19T20:26:37+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: QPALM_interface_update.jl
# @Last modified by:   massimo
# @Last modified time: 2019-08-22T17:58:29+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


#
function update!(workspace::QPALMworkspace)::Nothing

    # reformat the settings for QPALM
    settings_dict = Dict{Symbol,Any}()
    for field in fieldnames(QPALMsettings)
        settings_dict[field] = getfield(workspace.settings,field)
    end

    # setup QPALM with the new constraint set
    QPALM.setup!(workspace.model;
                    Q=sparse(workspace.Q),
                    q=workspace.L,
                    A=vcat(speye(size(workspace.A,2)),sparse(workspace.A)),
                    bmin=vcat(workspace.varLoBs,workspace.cnsLoBs),
                    bmax=vcat(workspace.varUpBs,workspace.cnsUpBs),
                    settings_dict...)
    return
end


#
function insert_constraints!(workspace::QPALMworkspace,
                            constraintSet::T,index::Int;
                            suppressUpdate::Bool=false)::Nothing where T <: Union{NullConstraintSet,LinearConstraintSet}


    # if a null workspace was given as input nothing has to be done
    if constraintSet isa NullConstraintSet
        return
    end

    # check correctness of the input
    @assert size(constraintSet.A,2) == size(workspace.A,2)
    @assert size(constraintSet.A,1) == length(constraintSet.loBs)
    @assert size(constraintSet.A,1) == length(constraintSet.upBs)

    # perform the insertion
    workspace.A = vcat(workspace.A[1:index-1,:],constraintSet.A,workspace.A[index:end,:])
    splice!(workspace.cnsLoBs,index:index-1,constraintSet.loBs)
    splice!(workspace.cnsUpBs,index:index-1,constraintSet.upBs)

    # update the osqp workspace
    if !suppressUpdate
        update!(workspace)
    end

end

#
function remove_constraints!(workspace::QPALMworkspace,indices::Array{Int,1};
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
function permute_constraints!(workspace::QPALMworkspace,indices::Array{Int,1};suppressUpdate::Bool=false)::Nothing

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
function update_bounds!(workspace::QPALMworkspace,
                        cnsLoBs::Array{Float64,1},cnsUpBs::Array{Float64,1},
                        varLoBs::Array{Float64,1},varUpBs::Array{Float64,1};
                        suppressUpdate::Bool=false)::Nothing


    @. workspace.cnsLoBs = cnsLoBs
    @. workspace.cnsUpBs = cnsUpBs
    @. workspace.varLoBs = varLoBs
    @. workspace.varUpBs = varUpBs

    # propagate the changes to the QPALM solver
    if !suppressUpdate
        QPALM.update!(workspace.model;
                      bmin=vcat(workspace.varLoBs,workspace.cnsLoBs),
                      bmax=vcat(workspace.varUpBs,workspace.cnsLoBs))
    end

    return
end


# ...
function set_objective!(workspace::QPALMworkspace,newObjective::T;suppressUpdate::Bool=false)::Nothing where T <: AbstractObjective

    if newObjective isa NullObjective
        @. workspace.Q = 0.
        @. workspace.L = 0.
    elseif newObjective isa LinearObjective
        @. workspace.Q = 0.
        @. workspace.L = newObjective.L
    elseif newObjective isa QuadraticObjective
        @. workspace.Q = newObjective.Q
        @. workspace.L = newObjective.L
    else
        @error "QPALM cannot deal with the given objective function"
    end

    # update the qpalm workspace
    if !suppressUpdate
        update!(workspace)
    end

    return
end


# ...
function set_constraintSet!(workspace::QPALMworkspace,newConstraintSet::T;suppressUpdate::Bool=false)::Nothing where T <: AbstractConstraintSet

    if newConstraintSet isa NullConstraintSet
        numVariables = get_numVariables(workspace)
        @. workspace.A = zeros(0,numVariables)
        @. workspace.cnsLoBs = Float64[]
        @. workspace.cnsUpBs = Float64[]
    elseif newConstraintSet isa LinearConstraintSet
        @assert get_numVariables(newObjectiveTerm,2) == get_numVariables(workspace)
        @. workspace.A = newConstraintSet.A
        @. workspace.cnsLoBs = newConstraintSet.upBs
        @. workspace.cnsUpBs = newConstraintSet.loBs
    else
        @error "QPALM cannot deal with the given constraint set"
    end

    # update the gurobi workspace
    if !suppressUpdate
        update!(workspace)
    end

    return
end



#
function append_problem!(workspace::QPALMworkspace,
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
function append_problem!(workspace::QPALMworkspace,
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

    # update the qpalm workspace
    if !suppressUpdate
        update!(workspace)
    end

    return reliableObjLoBs
end


# this function updates the settings in the osqp model
function update_settings!(workspace::QPALMworkspace,settings::QPALMsettings;bb_primalTolerance::Float64=Inf,bb_timeLimit=Inf)

    # reformat the settings for QPALM
    settings_dict = Dict{Symbol,Any}()
    for field in fieldnames(QPALMsettings)
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
    QPALM.setup!(workspace.model;Q=workspace.Q,
                                q=workspace.L,
                                A=vcat(speye(length(workspace.varLoBs)),sparse(workspace.A)),
                                bmin=vcat(workspace.varLoBs,workspace.cnsLoBs),
                                bmax=vcat(workspace.varUpBs,workspace.cnsUpBs),
                                settings_dict...)
    return
end
