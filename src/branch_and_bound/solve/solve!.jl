# @Author: Massimo De Mauri <massimodemauri>
# @Date:   2019-02-06T18:10:22+01:00
# @Email:  massimo.demauri@gmail.com
# @Project: OpenBB
# @Filename: solve!.jl
# @Last modified by:   massimo
# @Last modified time: 2019-06-10T19:10:33+02:00
# @License: apache 2.0
# @Copyright: {{copyright}}


# include referred code
include("./branch_and_solve!.jl")
include("./run!.jl")


# This is the main function called to solve a branch and bound problem
function solve!(workspace::BBworkspace)::Nothing

	# solve the root node
	@sync if workspace.status.description == "new"
		solve!(workspace.activeQueue[1],workspace)
		workspace.status.objLoB = workspace.activeQueue[1].objective

		# initialize the pseudo costs
		initialize_pseudoCosts!(workspace.settings.pseudoCostsInitialization,workspace.pseudoCosts,workspace.activeQueue[1])
		if workspace.settings.numProcesses > 1
			# initialize the pseudoCosts in the remote workers
			for k in 2:workspace.settings.numProcesses
				@async remotecall_fetch(Main.eval,k,:(@. workspace.pseudoCosts = $(workspace.pseudoCosts)))
			end
		end
	end


	if workspace.settings.numProcesses > 1
		# start the remote branch and bound processes
		for k in 2:workspace.settings.numProcesses
			@async remotecall_fetch(Main.eval,k,:(OpenBB.run!(workspace)))
		end
	end

	# start the local BB process
	run!(workspace)

	@sync return

end
