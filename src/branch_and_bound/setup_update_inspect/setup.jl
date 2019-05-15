# @Author: Massimo De Mauri <massimodemauri>
# @Date:   2019-02-06T16:22:41+01:00
# @Email:  massimo.demauri@gmail.com
# @Project: OpenBB
# @Filename: setup.jl
# @Last modified by:   massimo
# @Last modified time: 2019-05-15T10:04:16+02:00
# @License: apache 2.0
# @Copyright: {{copyright}}


function setup(problem::Problem, bb_settings::BBsettings=BBsettings(), ss_settings::AbstractSettings=NullSettings())::BBworkspace


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


    # check correctness of the inputs
    @assert bb_settings.numProcesses>=0

	if bb_settings.numProcesses == 0
	# set the processes to use to the number of cpus in the machine
    	bb_settings.numProcesses = length(Sys.cpu_info())
	else
		bb_settings.numProcesses = min(length(Sys.cpu_info()),bb_settings.numProcesses)
	end

	if nprocs() < bb_settings.numProcesses
	# add new processes if there aren't enough
		addprocs(bb_settings.numProcesses - nprocs())
	end



	if bb_settings.numProcesses > 1

		# send load OpenBB in the workers global scope
		workersList = workers()
		@sync for k in 1:length(workersList)
			# load OpenBB in the workers
			@async remotecall_fetch(Main.eval,workersList[k],:(using OpenBB))
		end


		# construct the communication channels
		communicationChannels = Array{RemoteChannel,1}(undef,bb_settings.numProcesses+1)
		for k in 1:bb_settings.numProcesses+1
			communicationChannels[k] = RemoteChannel(()->Channel{AbstractBBnode}(5))
		end

		# create the remote workspaces
		expressions = Array{Expr,1}(undef,length(workersList))
		globalInfo = SharedArray{Float64,1}([Inf,0.,0.])
		for k in 1:length(workersList)
			# create an empty workspace in each of the workers global scope
			expressions[k]  = :(workspace = OpenBB.BBworkspace(OpenBB.setup($problem,$ss_settings,
																		bb_primalTolerance=$(bb_settings.primalTolerance),
																		bb_timeLimit=$(bb_settings.timeLimit)
																		),
		                        						   $(problem.varSet.dscIndices),$(problem.varSet.sos1Groups),$sosConstraints,
														   Array{OpenBB.BBnode,1}(),Array{OpenBB.BBnode,1}(),Array{OpenBB.BBnode,1}(),OpenBB.BBstatus(),
														   $(communicationChannels[k]),$(communicationChannels[k+1]),$globalInfo,
														   $bb_settings);
							nothing)
		end
		@sync for k in 1:length(workersList)
			@async remotecall_fetch(Main.eval,workersList[k],expressions[k])
		end



		# construct the master BBworkspace
		workspace = BBworkspace(setup(problem,ss_settings,
									  bb_primalTolerance=bb_settings.primalTolerance,
									  bb_timeLimit=bb_settings.timeLimit),
								problem.varSet.dscIndices,problem.varSet.sos1Groups,sosConstraints,
								[BBnode(Dict{Int,Float64}(),Dict{Int,Float64}(),
									   problem.varSet.pseudoCosts,problem.varSet.val,
									   zeros(nVars),zeros(Ncnss),1.,NaN,false)],
								Array{BBnode,1}(),Array{BBnode,1}(),BBstatus(),
								communicationChannels[end],communicationChannels[1],globalInfo,
								bb_settings)

	else
	# only one process: no communication channels needed

		# construct the master BBworkspace
		workspace = BBworkspace(setup(problem,ss_settings,
									  bb_primalTolerance=bb_settings.primalTolerance,
									  bb_timeLimit=bb_settings.timeLimit),
								problem.varSet.dscIndices,problem.varSet.sos1Groups,sosConstraints,
								[BBnode(Dict{Int,Float64}(),Dict{Int,Float64}(),
									   problem.varSet.pseudoCosts,problem.varSet.val,
									   zeros(nVars),zeros(Ncnss),1.,NaN,false)],
								Array{BBnode,1}(),Array{BBnode,1}(),BBstatus(),
								nothing,nothing,nothing,bb_settings)

	end





    return workspace

end

function setup(problem::NullProblem,bb_settings::BBsettings, ss_settings::AbstractSettings=NullSettings())::BBworkspace
    return setup(Problem(NullObjectiveFun(),NullConstraintSet(),EmptyVarSet()),bb_settings,ss_settings)
end

function setup(bb_settings::BBsettings, ss_settings::AbstractSettings=NullSettings())::BBworkspace
    return setup(Problem(NullObjectiveFun(),NullConstraintSet(),EmptyVarSet()),bb_settings,ss_settings)
end
