# @Author: Massimo De Mauri <massimodemauri>
# @Date:   2019-02-06T16:22:41+01:00
# @Email:  massimo.demauri@gmail.com
# @Project: OpenBB
# @Filename: setup.jl
# @Last modified by:   massimo
# @Last modified time: 2019-08-14T11:53:16+02:00
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

    if !(problem isa NullProblem)
        dscIndices = problem.varSet.dscIndices
	end


	# add new processes if there aren't enough
	if nprocs() < bbSettings.numProcesses
		addprocs(bbSettings.numProcesses - nprocs())
	end



	if bbSettings.numProcesses > 1

		# send load OpenBB in the workers global scope
		@everywhere Main.eval(:(using OpenBB))
		workersList = workers()[1:bbSettings.numProcesses-1]

		# build the root node
		rootNode = BBroot(numVars,numCnss)

		# construct the communication channels
		communicationChannels = Array{BBnodeChannel,1}(undef,bbSettings.numProcesses)
		for k in 1:bbSettings.numProcesses
			communicationChannels[k] = BBnodeChannel(flat_size(numVars,numCnss))
		end
		# communicationChannels = Array{RemoteChannel,1}(undef,bbSettings.numProcesses)
		# @sync for k in 1:bbSettings.numProcesses
		# 	@async communicationChannels[k] = RemoteChannel(()->Channel{AbstractBBnode}(2),k)
		# end

		# construct shared Memory
		objectiveBounds = SharedArray{Float64,1}(vcat([-Inf],repeat([Inf],bbSettings.numProcesses)))
		stats = SharedArray{Int,1}([0])
		arrestable = SharedArray{Bool,1}(repeat([false],bbSettings.numProcesses))

		# construct the master BBworkspace
		workspace = BBworkspace(setup(problem,ssSettings,
									  bb_primalTolerance=bbSettings.primalTolerance,
									  bb_timeLimit=bbSettings.timeLimit),
								problem.varSet.dscIndices,problem.varSet.sos1Groups,
								(problem.varSet.pseudoCosts,Array{Int,2}(undef,size(problem.varSet.pseudoCosts))),
								[rootNode],Array{BBnode,1}(),Array{BBnode,1}(),
								BBstatus(),BBsharedMemory(communicationChannels[1],communicationChannels[2],objectiveBounds,stats,arrestable),bbSettings)


		# construct the remote workspaces
		expressions = Array{Expr,1}(undef,length(workersList))
		for k in 2:workspace.settings.numProcesses
			if k < workspace.settings.numProcesses
				sharedMemory = BBsharedMemory(communicationChannels[k],communicationChannels[k+1],objectiveBounds,stats,arrestable)
			else
				sharedMemory = BBsharedMemory(communicationChannels[k],communicationChannels[1],objectiveBounds,stats,arrestable)
			end
			expressions[k-1] = :(workspace = OpenBB.BBworkspace(OpenBB.setup($problem,$ssSettings,
																		bb_primalTolerance=$(bbSettings.primalTolerance),
																		bb_timeLimit=$(bbSettings.timeLimit)
																		),
														   copy($(problem.varSet.dscIndices)),copy($(problem.varSet.sos1Groups)),
														   ($problem.varSet.pseudoCosts,Array{Int,2}(undef,size($problem.varSet.pseudoCosts))),
														   Array{OpenBB.BBnode,1}(),Array{OpenBB.BBnode,1}(),Array{OpenBB.BBnode,1}(),
														   OpenBB.BBstatus(objLoB=Inf,description="empty"),$sharedMemory,deepcopy($bbSettings)))
	    end
		@sync for k in 1:length(workersList)
			@async remotecall_fetch(Main.eval,workersList[k],expressions[k])
		end

	else # only one process: no communication channels needed

		# build the root node
		rootNode = BBroot(numVars,numCnss)

		# construct the master BBworkspace
		workspace = BBworkspace(setup(problem,ssSettings,
									  bb_primalTolerance=bbSettings.primalTolerance,
									  bb_timeLimit=bbSettings.timeLimit),
								problem.varSet.dscIndices,problem.varSet.sos1Groups,
								(problem.varSet.pseudoCosts,Array{Int,2}(undef,size(problem.varSet.pseudoCosts))),
								[rootNode],Array{BBnode,1}(),Array{BBnode,1}(),
								BBstatus(),NullSharedMemory(),bbSettings)
	end

    return workspace

end

function setup(problem::NullProblem,bbSettings::BBsettings, ssSettings::AbstractSettings=NullSettings())::BBworkspace
    return setup(Problem(NullObjective(),NullConstraintSet(),EmptyVarSet()),bbSettings,ssSettings)
end

function setup(bbSettings::BBsettings, ssSettings::AbstractSettings=NullSettings())::BBworkspace
    return setup(Problem(NullObjective(),NullConstraintSet(),EmptyVarSet()),bbSettings,ssSettings)
end
