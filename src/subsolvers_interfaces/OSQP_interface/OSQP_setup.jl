# @Author: Massimo De Mauri <massimo>
# @Date:   2019-02-19T22:28:53+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: OSQP_setup.jl
# @Last modified by:   massimo
# @Last modified time: 2019-09-02T12:21:29+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# this function creates an OSQP.Model representing the given CvxQproblem
function setup(problem::Problem,settings::OSQPsettings;bb_primalTolerance::Float64=Inf,bb_timeLimit=Inf)::OSQPworkspace

    # collect info on the problem
    numVars = get_numVariables(problem)

    # use the workspace info to update the OSQP settings
    settings.eps_prim_inf = min(settings.eps_prim_inf,bb_primalTolerance*1e-1)
    if bb_timeLimit < Inf
        if settings.timeLimit == 0.
            settings.timeLimit = bb_timeLimit
        else
            settings.timeLimit = min(settings.timeLimit,bb_timeLimit)
        end
    end

    # check the objective function
    if problem.objFun isa NullObjective
        Q = spzeros(nVars,nVars)
        L = zeros(nVars)
    elseif problem.objFun isa LinearObjective
        Q = spzeros(nVars,nVars)
        L = problem.objFun.L
    elseif problem.objFun isa QuadraticObjective
        Q = dropzeros(sparse(problem.objFun.Q))
        L = problem.objFun.L
    else
        @error "OSQP cannot deal with the given objective function"
    end

    # check the constraint set
    if problem.cnsSet isa NullConstraintSet
        println(problem.cnsSet)
        A = spzeros(0,length(problem.varSet.loBs))
        cnsLoBs = Float64[]
        cnsUpBs = Float64[]
    elseif problem.cnsSet isa LinearConstraintSet
        A = dropzeros(sparse(problem.cnsSet.A))
        cnsLoBs = copy(problem.cnsSet.loBs)
        cnsUpBs = copy(problem.cnsSet.upBs)
    else
        @error "OSQP cannot deal with the given constraint set"
    end

    # reformat the settings for OSQP
    settings_dict = Dict{Symbol,Any}()
    for field in fieldnames(OSQPsettings)
        settings_dict[field] = getfield(settings,field)
    end

    # create the subsolver OSQPworkspace
    model = OSQP.Model()
    if length(problem.varSet.loBs) > 0
        OSQP.setup!(model;P=Q,q=L,
                          A=vcat(speye(length(problem.varSet.loBs)),sparse(A)),
                          l=vcat(problem.varSet.loBs,cnsLoBs),
                          u=vcat(problem.varSet.upBs,cnsUpBs),
                          settings_dict...)
    end

    return OSQPworkspace(copy(Q),copy(L),
                     copy(A),copy(cnsLoBs),copy(cnsUpBs),
                     copy(problem.varSet.loBs),copy(problem.varSet.upBs),
                     model,settings)
end
