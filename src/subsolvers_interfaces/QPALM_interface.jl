# @Author: Massimo De Mauri <massimo>
# @Date:   2019-02-19T20:29:29+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: QPALM_interface.jl
# @Last modified by:   massimo
# @Last modified time: 2019-09-27T17:17:52+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

using QPALM

## Settings ##########################################################

# structure used to hold the settings for QPALM
mutable struct QPALMsettings <: AbstractSettings
    max_iter::Int                       # Maximum number of interations 0 < max_iter
    eps_abs::Float64                    # Absolute tolerance 	0 <= eps_abs
    eps_rel::Float64                    # Relative tolerance 	0 <= eps_rel
    eps_abs_in::Float64                 # Intermediate absolute convergence tolerance 	0 <= eps_abs_in
    eps_rel_in::Float64                 # Intermediate relative convergence tolerance 	0 <= eps_rel_in
    rho::Float64                        # Tolerance scaling factor step 	0 < rho < 1
    eps_prim_inf::Float64               # Primal infeasibility tolerance    0 <= eps_prim_inf
    eps_dual_inf::Float64               # Dual infeasibility tolerance      0 <= eps_dual_inf
    theta::Float64                      # Penalty update criterion          theta <= 1
    delta::Float64                      # Penalty update factor             1 < delta
    tau_init::Float64                   # Initial stepsize in backtracking  0 < tau_init
    proximal::Bool                      # Use proximal method of multipliers  True/False
    gamma_init::Float64                 # Injtial proximal penalty          0 < gamma_init
    gamma_upd::Float64                  # Proximal penalty update factor    1 <= gamma_upd
    gamma_max::Float64                  # Proximal penalty parameter cap    gamma_init < gamma_max
    scaling::Int                        # Number of scaling iterations 	0 (disabled) or 0 < scaling (integer)
    nonconvex::Bool                     # Indicate if QP is convex          True/False
    verbose::Bool                       # Print output 	True/False
    warm_start::Bool                    # Perform warm starting 	True/False
end


function QPALMsettings(; max_iter::Int=10000,
                         eps_abs::Float64=1e-8,
                         eps_rel::Float64=1e-8,
                         eps_abs_in::Float64=1.0,
                         eps_rel_in::Float64=1.0,
                         eps_prim_inf::Float64=1e-8,
                         eps_dual_inf::Float64=1e-8,
                         rho::Float64=0.1,
                         theta::Float64=0.25,
                         delta::Float64=10.0,
                         tau_init::Float64=1.0,
                         proximal::Bool=true,
                         gamma_init::Float64=1e6,
                         gamma_upd::Float64=10.0,
                         gamma_max::Float64=1e8,
                         scaling::Int=10,
                         nonconvex::Bool=false,
                         verbose::Bool=false,
                         warm_start::Bool=true)::QPALMsettings



    return QPALMsettings(max_iter,eps_abs,eps_rel,eps_abs_in,eps_rel_in,rho,eps_prim_inf,
                        eps_dual_inf,theta,delta,tau_init,proximal,gamma_init,
                        gamma_upd,gamma_max,scaling,nonconvex,verbose,warm_start)
end




## Workspace ##########################################################
# structure used for storing data for QPALM solver
mutable struct QPALMworkspace <: AbstractWorkspace
    # objective
    problem::Problem
    # memory
    model::QPALM.Model
    settings::QPALMsettings
    # outdated flag
    outdated::Bool
end

## Setup & Update ##########################################################
# this function creates an QPALM.Model representing the given CvxQproblem
function setup(problem::Problem,settings::QPALMsettings;bb_primalTolerance::Float64=Inf,bb_timeLimit=Inf)::QPALMworkspace

    # check the problem
    @assert problem.objFun isa NullObjective || problem.objFun isa LinearObjective || problem.objFun isa QuadraticObjective
    @assert problem.cnsSet isa NullConstraintSet || problem.cnsSet isa LinearConstraintSet


    # overwrite the qpalm settings depending on the branch and bound settings
    settings.eps_prim_inf = min(settings.eps_prim_inf,bb_primalTolerance*1e-1)
    if bb_timeLimit < Inf
        if settings.timeLimit == 0.
            settings.timeLimit = bb_timeLimit
        else
            settings.timeLimit = min(settings.timeLimit,bb_timeLimit)
        end
    end

    # reformat the settings
    settings_dict = Dict{Symbol,Any}()
    for field in fieldnames(QPALMsettings)
        settings_dict[field] = getfield(settings,field)
    end

    # ensure type consistency
    objFun = QuadraticObjective{SparseMatrixCSC{Float64,Int64},Array{Float64,1}}(problem.objFun)
    cnsSet = LinearConstraintSet{SparseMatrixCSC{Float64,Int64}}(problem.cnsSet)


    # create the QPALMworkspace
    model = QPALM.Model()
    if length(problem.varSet.loBs) > 0
        QPALM.setup!(model;Q=sparse(objFun.Q),q=objFun.L,
                           A=vcat(speye(get_size(problem.varSet)),sparse(cnsSet.A)),
                           bmin=vcat(problem.varSet.loBs,cnsSet.loBs),
                           bmax=vcat(problem.varSet.upBs,cnsSet.upBs),
                           settings_dict...)
    end

    return QPALMworkspace(problem,model,settings,false)
end

# it marks the workspace as outdated
function make_outdated!(workspace::QPALMworkspace)::Nothing
    workspace.outdated = true
    return
end

#
function update!(workspace::QPALMworkspace)::Nothing

    # reformat the settings
    settings_dict = Dict{Symbol,Any}()
    for field in fieldnames(QPALMsettings)
        settings_dict[field] = getfield(workspace.settings,field)
    end

    # ensure type consistency
    objFun = QuadraticObjective{SparseMatrixCSC{Float64,Int64},Array{Float64,1}}(workspace.problem.objFun)
    cnsSet = LinearConstraintSet{SparseMatrixCSC{Float64,Int64}}(workspace.problem.cnsSet)

    # setup QPALM for the new problem
    QPALM.setup!(workspace.model;Q=sparse(objFun.Q),q=objFun.L,
                 A=vcat(speye(get_size(workspace.problem.varSet)),sparse(cnsSet.A)),
                 bmin=vcat(workspace.problem.varSet.loBs,cnsSet.loBs),
                 bmax=vcat(workspace.problem.varSet.upBs,cnsSet.upBs),
                 settings_dict...)

    # mark the workspace as up to date
    workspace.outdated = false
    return
end

## Solve ##########################################################
function solve!(node::BBnode,workspace::QPALMworkspace)::Tuple{Int8,Float64}

    # update the problem formulation if needed
    if workspace.outdated
        update!(workspace)
    end

    # collect info on the problem
    numVars = get_size(workspace.problem.varSet)

    # update bounds in the the qpalm model
    QPALM.update!(workspace.model;bmin=vcat(node.varLoBs,node.cnsLoBs),bmax=vcat(node.varUpBs,node.cnsUpBs))

    # set hotstart info
    if length(node.primal) > 0 && length(node.bndDual) > 0 && length(node.cnsDual) > 0
        QPALM.warm_start!(workspace.model; x_warm_start=node.primal, y_warm_start=vcat(node.bndDual,node.cnsDual))
    end

    # solve problem
    sol = QPALM.solve!(workspace.model)

    # output sol info
    if  sol.info.status_val == 1
        status = 0 # "solved"
        @. node.primal = sol.x
        @. node.bndDual = sol.y[1:numVars]
        @. node.cnsDual = sol.y[numVars+1:end]
        objFun = QuadraticObjective{SparseMatrixCSC{Float64,Int64},Array{Float64,1}}(workspace.problem.objFun)
        node.objVal = 1/2 * transpose(node.primal) * objFun.Q * node.primal + transpose(objFun.L) * node.primal
        node.objGap = max(workspace.settings.eps_abs,
                          workspace.settings.eps_rel*abs(node.objVal))
    elseif sol.info.status_val == -3
        status = 1 # "infeasible"
        @. node.primal = @. min(max(sol.x,node.varLoBs),node.varUpBs)
        @. node.bndDual = sol.y[1:numVars]
        @. node.cnsDual = sol.y[numVars+1:end]
        node.objVal = Inf
        node.objGap = 0.0
    elseif sol.info.status_val in [2,3,4,-6,-2]
        status = 2 # "unreliable"
        @. node.primal = min(max(sol.x,node.varLoBs),node.varUpBs)
        @. node.bndDual = sol.y[1:numVars]
        @. node.cnsDual = sol.y[numVars+1:end]
        objFun = QuadraticObjective{SparseMatrixCSC{Float64,Int64},Array{Float64,1}}(workspace.problem.objFun)
        newObjVal = 1/2 * transpose(node.primal) * objFun.Q * node.primal + transpose(objFun.L) * node.primal
        if newObjVal >= node.ObjVal - node.objGap
            node.objGap = newObjVal - node.objVal + node.objGap #TODO: recopute the gap if possible
            node.objVal = newObjVal
        else
            node.objGap = Inf #TODO: recopute the gap if possible
            @warn "Inaccuracy in node sol, status: "*string(sol.info.status)*" (code: "*string(status)*")"
        end
        @warn "Inaccuracy in node sol, message: "*string(sol.info.status)*" (code: "*string(sol.info.status_val)*")"
    elseif sol.info.status_val in [-7,-10]
        status = 3 # "error"
        @error "Subsover error, status: "*string(sol.info.status)*" (code: "*string(sol.info.status_val)*")"
    else
        @error "Subsolver unknown status: "*string(sol.info.status)*"("*string(sol.info.status_val)*")"
    end

    return (status, sol.info.run_time)
end
