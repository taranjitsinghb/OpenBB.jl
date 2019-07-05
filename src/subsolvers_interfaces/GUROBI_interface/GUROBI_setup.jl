# @Author: Massimo De Mauri <massimo>
# @Date:   2019-03-03T19:26:23+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: Gurobi_setup.jl
# @Last modified by:   massimo
# @Last modified time: 2019-07-01T18:45:35+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


function setup(problem::Problem,settings::GUROBIsettings;bb_primalTolerance::Float64=Inf,bb_timeLimit=Inf)::GUROBIworkspace

    # reformat the settings for GUROBI
    settings_dict = Dict{Symbol,Any}()
    for field in fieldnames(GUROBIsettings)
        settings_dict[field] = getfield(settings,field)
    end

    # overwrite the osqp setting depending on the branch and bound settings
    settings_dict[:FeasibilityTol] = min(settings_dict[:FeasibilityTol],bb_primalTolerance)

    if bb_timeLimit < Inf
        if settings_dict[:TimeLimit] == 0.
            settings_dict[:TimeLimit] = bb_timeLimit
        else
            settings_dict[:TimeLimit] = min(settings_dict[:TimeLimit],bb_timeLimit)
        end
    end


    nVars = length(problem.varSet.loBs)


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
        @error "GUROBI cannot deal with the given objective function"
    end

    # check the constraint set
    if problem.cnsSet isa NullConstraintSet
        # adapt the constraint set to accomodate for variables bounds
        A = spzeros((0,length(problem.varSet.loBs)))
        cnsLoBs = Float64[]
        cnsUpBs = Float64[]

    elseif problem.cnsSet isa LinearConstraintSet
        # adapt the constraint set to accomodate for variables bounds
        A = dropzeros(sparse(problem.cnsSet.A))
        cnsLoBs = copy(problem.cnsSet.loBs)
        cnsUpBs = copy(problem.cnsSet.upBs)
    else
        @error "GUROBI cannot deal with the given constraint set"
    end


    # create the subsolver OSQPworkspace
    env = Gurobi.Env()
    Gurobi.setparams!(env;settings_dict...)

    return GUROBIworkspace(copy(Q),copy(L),
                     copy(A),copy(cnsLoBs),copy(cnsUpBs),
                     copy(problem.varSet.loBs),copy(problem.varSet.upBs),
                     env,settings)
end
