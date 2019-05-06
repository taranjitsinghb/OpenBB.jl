# @Author: Massimo De Mauri <massimodemauri>
# @Date:   2019-02-06T16:22:41+01:00
# @Email:  massimo.demauri@gmail.com
# @Project: OpenBB
# @Filename: setup.jl
# @Last modified by:   massimo
# @Last modified time: 2019-04-30T17:10:07+02:00
# @License: apache 2.0
# @Copyright: {{copyright}}


function setup(problem::Problem, bb_settings::BBsettings=BBsettings(), ss_settings::AbstractSettings=NullSettings())::BBworkspace

    # check correctness of the inputs
    @assert bb_settings.maxProcesses>0

    # load default settings
    if ss_settings isa NullSettings
        if problem isa Problem{LinearObj,LinearCns} # linear Problem
            ss_settings = OSQPsettings()
        elseif problem isa Problem{QuadraticObj,LinearCns} # quadratic Problem
            ss_settings = OSQPsettings()
        else
            @error "Type of the problem not understood"
        end
    end

    # create the continuous relaxation for the problem
    subsolverWS = setup(problem,ss_settings,bb_primalTolerance=bb_settings.primalTolerance,bb_timeLimit=bb_settings.timeLimit)    # check if at least one process is allowed_settings)

    # collect some data for the BBworkspace
    nVars = length(problem.varSet.loBs)

    if problem.cnsSet isa NullConstraintSet || problem isa NullProblem
        Ncnss = 0
        sosConstraints = LinearCns(sparse(Int[],Int[],Float64[]),Float64[],Float64[],Int[])
    else
        Ncnss = length(problem.cnsSet.loBs)
        dscIndices = problem.varSet.dscIndices
        sosIndices = problem.cnsSet.sosIndices
        sosConstraints = LinearCns( sparse(problem.cnsSet.A[sosIndices,dscIndices]),
                                    problem.cnsSet.loBs[sosIndices],
                                    problem.cnsSet.upBs[sosIndices],Int[])
    end



    # construct the BBworkspace
    return BBworkspace( subsolverWS,
                        problem.varSet.dscIndices,
                        problem.varSet.sos1Groups,
                        sosConstraints,
                        [BBsubproblem(Dict{Int,Float64}(),Dict{Int,Float64}(),problem.varSet.pseudoCosts,problem.varSet.val,zeros(nVars),zeros(Ncnss),1.,NaN,false)],
                        Array{BBsubproblem,1}(),
                        Array{BBsubproblem,1}(),
                        BBstatus(),
                        bb_settings)

end

function setup(problem::NullProblem,bb_settings::BBsettings, ss_settings::AbstractSettings=NullSettings())::BBworkspace
    return setup(Problem(NullObjectiveFun(),NullConstraintSet(),EmptyVarSet()),bb_settings,ss_settings)
end

function setup(bb_settings::BBsettings, ss_settings::AbstractSettings=NullSettings())::BBworkspace
    return setup(Problem(NullObjectiveFun(),NullConstraintSet(),EmptyVarSet()),bb_settings,ss_settings)
end
