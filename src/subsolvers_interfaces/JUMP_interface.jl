using JuMP
using Mosek

mutable struct JUMP_SDPsettings <: AbstractSettings
    max_iter::Int                       # Maximum number of iterations 	0 < max_iter (integer)
    eps_abs::Float64                    # Absolute tolerance 	0 <= eps_abs
    eps_rel::Float64                    # Relative tolerance 	0 <= eps_rel
    eps_prim_inf::Float64               # Primal infeasibility tolerance 	0 <= eps_prim_inf
    eps_dual_inf::Float64               # Dual infeasibility tolerance 	0 <= eps_dual_inf
    verbose::Bool                       # Print output 	True/False
end

function JUMP_SDPsettings(; max_iter::Int=10000,
                        eps_abs::Float64=1e-6,
                        eps_rel::Float64=1e-6,
                        eps_prim_inf::Float64=1e-4,
                        eps_dual_inf::Float64=1e-4,
                        verbose::Bool=false)::JUMP_SDPsettings

    return JUMP_SDPsettings(max_iter,eps_abs,eps_rel,eps_prim_inf,eps_dual_inf,verbose)
end

mutable struct JUMPworkspace <: AbstractWorkspace
    # objective
    problem::Problem
    # memory
    model::JuMP.Model
    settings::JUMP_SDPsettings
    # outdated flag
    outdated::Bool
end


## Setup & Update ##########################################################
# this function creates an QPALM.Model representing the given CvxQproblem
function setup(problem::Problem,settings::JUMP_SDPsettings;bb_primalTolerance::Float64=Inf,bb_timeLimit=Inf)::JUMPworkspace

    # check the problem
    @assert problem.objFun isa NullObjective || problem.objFun isa LinearObjective || problem.objFun isa QuadraticObjective
    @assert problem.cnsSet isa NullConstraintSet || problem.cnsSet isa LinearConstraintSet || problem.cnsSet isa SDPConstraintSet


    # overwrite the qpalm settings depending on the branch and bound settings
    # settings.eps_prim_inf = min(settings.eps_prim_inf,bb_primalTolerance*1e-1)
    # if bb_timeLimit < Inf
    #     if settings.timeLimit == 0.
    #         settings.timeLimit = bb_timeLimit
    #     else
    #         settings.timeLimit = min(settings.timeLimit,bb_timeLimit)
    #     end
    # end

    # reformat the settings
    settings_dict = Dict{Symbol,Any}()
    for field in fieldnames(JUMP_SDPsettings)
        settings_dict[field] = getfield(settings,field)
    end

    # ensure type consistency
    objFun = LinearObjective{SparseMatrixCSC{Float64,Int64},Array{Float64,1}}(problem.objFun)
    cnsSet = SDPConstraintSet{Symmetric{Float64,Array{Float64,2}}}(problem.cnsSet)


    # create the QPALMworkspace
    # model = QPALM.Model()
    # if length(problem.varSet.loBs) > 0
    #     QPALM.setup!(model;Q=sparse(objFun.Q),q=objFun.L,
    #                        A=vcat(speye(get_size(problem.varSet)),sparse(cnsSet.A)),
    #                        bmin=vcat(problem.varSet.loBs,cnsSet.loBs),
    #                        bmax=vcat(problem.varSet.upBs,cnsSet.upBs),
    #                        settings_dict...)
    # end

    return JUMPworkspace(problem,model,settings,false)
end

# it marks the workspace as outdated
function make_outdated!(workspace::JUMPworkspace)::Nothing
    workspace.outdated = true
    return
end

#
function update!(workspace::JUMPworkspace)::Nothing

    # reformat the settings
    settings_dict = Dict{Symbol,Any}()
    for field in fieldnames(JUMPsettings)
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
