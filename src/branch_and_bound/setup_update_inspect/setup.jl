# @Author: Massimo De Mauri <massimodemauri>
# @Date:   2019-02-06T16:22:41+01:00
# @Email:  massimo.demauri@gmail.com
# @Project: OpenBB
# @Filename: setup.jl
# @Last modified by:   massimo
# @Last modified time: 2019-09-25T23:15:29+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


function setup(problem::Problem, bbSettings::BBsettings=BBsettings(), ssSettings::AbstractSettings=NullSettings())::BBworkspace


    # load default settings
    if ssSettings isa NullSettings
        if problem isa Problem{LinearObjective,LinearConstraintSet} # linear Problem
            ssSettings = OSQPsettings()
        elseif problem isa Problem{QuadraticObjective,LinearConstraintSet} # quadratic Problem
            ssSettings = OSQPsettings()
        else
            @error "Type of the problem not understood"
        end
    end

    # collect some data for the BBworkspace
    numVars = get_numVariables(problem)
	numCnss = get_numConstraints(problem)
	numDscVars = get_numDiscreteVariables(problem)

	localProblem = deepcopy(problem)

	if bbSettings.numProcesses == 1	# single-process setup

		# construct the master BBworkspace
		workspace = BBworkspace(localProblem,
							    setup(localProblem,ssSettings,
									  bb_primalTolerance=bbSettings.primalTolerance,
									  bb_timeLimit=bbSettings.timeLimit),
								NullSharedMemory(),
								Array{BBnode,1}(),Array{BBnode,1}(),Array{BBnode,1}(),
								BBstatus(),bbSettings,false)

		# build the root node and solve it
		push!(workspace.activeQueue,BBroot(workspace))
		solve_node!(workspace.activeQueue[1],workspace)
		workspace.status.objLoB = workspace.activeQueue[1].objVal - workspace.activeQueue[1].objGap

		# initialize the pseudo costs
		initialize_pseudoCosts!(workspace.settings.pseudoCostsInitialization,workspace.problem.varSet.pseudoCosts,workspace.activeQueue[1])

	else # multi-process setup

		# add new processes if there aren't enough
		if nprocs() < bbSettings.numProcesses
			addprocs(bbSettings.numProcesses - nprocs())
		end

		# send load OpenBB in the workers global scope
		@everywhere Main.eval(:(using OpenBB))
		workersList = workers()[1:bbSettings.numProcesses-1]

		# construct the communication channels
		communicationChannels = Array{BBnodeChannel,1}(undef,bbSettings.numProcesses)
		for k in 1:bbSettings.numProcesses
			communicationChannels[k] = BBnodeChannel(flat_size(numVars,numCnss))
		end

		# construct shared Memory
		objectiveBounds = SharedArray{Float64,1}(vcat([-Inf],repeat([Inf],bbSettings.numProcesses)))
		stats = SharedArray{Int,1}([0])
		arrestable = SharedArray{Bool,1}(repeat([false],bbSettings.numProcesses))

		# construct the master BBworkspace
		workspace = BBworkspace(localProblem,
								setup(localProblem,ssSettings,
									  bb_primalTolerance=bbSettings.primalTolerance,
									  bb_timeLimit=bbSettings.timeLimit),
							    BBsharedMemory(communicationChannels[1],communicationChannels[2],objectiveBounds,stats,arrestable),
								Array{BBnode,1}(),Array{BBnode,1}(),Array{BBnode,1}(),
								BBstatus(),bbSettings,false)

		# construct the remote workspaces
		expressions = Array{Expr,1}(undef,length(workersList))
		for k in 2:workspace.settings.numProcesses
			if k < workspace.settings.numProcesses
				sharedMemory = BBsharedMemory(communicationChannels[k],communicationChannels[k+1],objectiveBounds,stats,arrestable)
			else
				sharedMemory = BBsharedMemory(communicationChannels[k],communicationChannels[1],objectiveBounds,stats,arrestable)
			end
			expressions[k-1] = :(workspace = OpenBB.BBworkspace($localProblem,
																OpenBB.setup($localProblem,$ssSettings,
																			 bb_primalTolerance=$(bbSettings.primalTolerance),
																			 bb_timeLimit=$(bbSettings.timeLimit)),
																$sharedMemory,
															    Array{OpenBB.BBnode,1}(),Array{OpenBB.BBnode,1}(),Array{OpenBB.BBnode,1}(),
															    OpenBB.BBstatus(objLoB=Inf,description="empty"),$bbSettings,false))
	    end

		@sync for k in 1:length(workersList)
			@async remotecall_fetch(Main.eval,workersList[k],expressions[k])
		end

		# build the root node and solve it
		push!(workspace.activeQueue,BBroot(workspace))
		solve!(workspace.activeQueue[1],workspace.subsolverWS)
		workspace.status.objLoB = workspace.activeQueue[1].objVal - workspace.activeQueue[1].objGap

		# initialize the pseudo costs in the master process
		initialize_pseudoCosts!(workspace.settings.pseudoCostsInitialization,workspace.problem.varSet.pseudoCosts,workspace.activeQueue[1])

		# initialize the pseudoCosts in the remote workers
		@sync for k in 2:workspace.settings.numProcesses
			@async remotecall_fetch(Main.eval,k,:(workspace.problem.varSet.pseudoCosts[1] .= $(workspace.problem.varSet.pseudoCosts[1]);workspace.problem.varSet.pseudoCosts[2] .= $(workspace.problem.varSet.pseudoCosts[2])))
		end
	end

    return workspace

end

function setup(problem::NullProblem,bbSettings::BBsettings, ssSettings::AbstractSettings=NullSettings())::BBworkspace
    return setup(Problem(NullObjective(),NullConstraintSet(),EmptyVarSet()),bbSettings,ssSettings)
end

function setup(bbSettings::BBsettings, ssSettings::AbstractSettings=NullSettings())::BBworkspace
    return setup(Problem(NullObjective(),NullConstraintSet(),EmptyVarSet()),bbSettings,ssSettings)
end
