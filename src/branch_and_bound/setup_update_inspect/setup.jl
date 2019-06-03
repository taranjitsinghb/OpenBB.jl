# @Author: Massimo De Mauri <massimodemauri>
# @Date:   2019-02-06T16:22:41+01:00
# @Email:  massimo.demauri@gmail.com
# @Project: OpenBB
# @Filename: setup.jl
# @Last modified by:   massimo
# @Last modified time: 2019-06-03T16:09:08+02:00
# @License: apache 2.0
# @Copyright: {{copyright}}


function setup(problem::Problem, bb_settings::BBsettings=BBsettings(), ss_settings::AbstractSettings=NullSettings())::BBworkspace


    # load default settings
    if ss_settings isa NullSettings
        if problem isa Problem{LinearObjective,LinearConstraintSet} # linear Problem
            ss_settings = OSQPsettings()
        elseif problem isa Problem{QuadraticObjective,LinearConstraintSet} # quadratic Problem
            ss_settings = OSQPsettings()
        else
            @error "Type of the problem not understood"
        end
    end

    # collect some data for the BBworkspace
    numVars = get_numVariables(problem)
	numCnss = get_numConstraints(problem)
	numDscVars = get_numDiscreteVariables(problem)

    if !(problem.cnsSet isa NullConstraintSet) && !(problem isa NullProblem)
        dscIndices = problem.varSet.dscIndices
        sosIndices = problem.cnsSet.sosIndices
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
		@everywhere Main.eval(:(using OpenBB))
		workersList = workers()[1:bb_settings.numProcesses-1]


		# construct the communication channels
		communicationChannels = Array{BBnodeChannel,1}(undef,bb_settings.numProcesses)
		for k in 1:bb_settings.numProcesses
			communicationChannels[k] = BBnodeChannel(flat_size(numVars,numDscVars,numCnss))
		end
		# communicationChannels = Array{RemoteChannel,1}(undef,bb_settings.numProcesses)
		# @sync for k in 1:bb_settings.numProcesses
		# 	@async communicationChannels[k] = RemoteChannel(()->Channel{AbstractBBnode}(2),k)
		# end

		# construct shared Memory
		objectiveBounds = SharedArray{Float64,1}(vcat([-Inf],repeat([Inf],bb_settings.numProcesses)))
		stats = SharedArray{Int,1}([0,0,0])

		# construct the master BBworkspace
		workspace = BBworkspace(setup(problem,ss_settings,
									  bb_primalTolerance=bb_settings.primalTolerance,
									  bb_timeLimit=bb_settings.timeLimit),
								problem.varSet.dscIndices,problem.varSet.sos1Groups,problem.varSet.pseudoCosts,
								[BBnode(-Infs(numDscVars),Infs(numDscVars),problem.varSet.val,
									   zeros(numVars),zeros(numCnss),1.,NaN,false)],
								Array{BBnode,1}(),Array{BBnode,1}(),
								BBstatus(),BBsharedMemory(communicationChannels[1],communicationChannels[2],objectiveBounds,stats),bb_settings)


		# construct the remote workspaces
		expressions = Array{Expr,1}(undef,length(workersList))
		for k in 2:workspace.settings.numProcesses
			if k < workspace.settings.numProcesses
				sharedMemory = BBsharedMemory(communicationChannels[k],communicationChannels[k+1],objectiveBounds,stats)
			else
				sharedMemory = BBsharedMemory(communicationChannels[k],communicationChannels[1],objectiveBounds,stats)
			end
			expressions[k-1] = :(workspace = OpenBB.BBworkspace(OpenBB.setup($problem,$ss_settings,
																		bb_primalTolerance=$(bb_settings.primalTolerance),
																		bb_timeLimit=$(bb_settings.timeLimit)
																		),
														   $(problem.varSet.dscIndices),$(problem.varSet.sos1Groups),$(problem.varSet.pseudoCosts),
														   Array{OpenBB.BBnode,1}(),Array{OpenBB.BBnode,1}(),Array{OpenBB.BBnode,1}(),
														   OpenBB.BBstatus(objLoB=Inf,description="empty"),$sharedMemory,$bb_settings))
	    end
		@sync for k in 1:length(workersList)
			@async remotecall_fetch(Main.eval,workersList[k],expressions[k])
		end

	else # only one process: no communication channels needed

		# construct the master BBworkspace
		workspace = BBworkspace(setup(problem,ss_settings,
									  bb_primalTolerance=bb_settings.primalTolerance,
									  bb_timeLimit=bb_settings.timeLimit),
								problem.varSet.dscIndices,problem.varSet.sos1Groups,problem.varSet.pseudoCosts,
								[BBnode(-Infs(numDscVars),Infs(numDscVars),problem.varSet.val,
									   zeros(numVars),zeros(numCnss),1.,NaN,false)],
								Array{BBnode,1}(),Array{BBnode,1}(),
								BBstatus(),NullSharedMemory(),bb_settings)

	end

    return workspace

end

function setup(problem::NullProblem,bb_settings::BBsettings, ss_settings::AbstractSettings=NullSettings())::BBworkspace
    return setup(Problem(NullObjectiveFunction(),NullConstraintSet(),EmptyVarSet()),bb_settings,ss_settings)
end

function setup(bb_settings::BBsettings, ss_settings::AbstractSettings=NullSettings())::BBworkspace
    return setup(Problem(NullObjectiveFunction(),NullConstraintSet(),EmptyVarSet()),bb_settings,ss_settings)
end
