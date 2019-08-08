# this function creates an JUMP_MOSEK.Model representing the given JUMPProblem

function setup(problem::Problem,settings::JUMP_MOSEKsettings;bb_primalTolerance::Float64=Inf,bb_timeLimit=Inf)::MOSEKworkspace


    settings_dict = Dict{Symbol,Any}()
    for field in fieldnames(JUMP_MOSEKsettings)
        settings_dict[field] = getfield(settings,field)
    end


    settings_dict[:eps_prim_inf] = min(settings_dict[:eps_prim_inf],bb_primalTolerance)

    if bb_timeLimit < Inf
        if settings_dict[:timeLimit] == 0.
            settings_dict[:timeLimit] = bb_timeLimit
        else
            settings_dict[:timeLimit] = min(settings_dict[:timeLimit],bb_timeLimit)
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
    else
        @error "MOSEKSDP cannot deal with the given objective function"
    end

    # check the constraint set
    if problem.cnsSet isa NullConstraintSet
        A = spzeros(0,length(problem.varSet.loBs))
        cnsLoBs = Float64[]
        cnsUpBs = Float64[]

    elseif problem.cnsSet isa LinearConstraintSet
        A = dropzeros(sparse(problem.cnsSet.A))
        cnsLoBs = copy(problem.cnsSet.loBs)
        cnsUpBs = copy(problem.cnsSet.upBs)
    else
        @error "MOSEKSDP cannot deal with the given constraint set"
    end


    # create the subsolver MOSEKworkspace


end
