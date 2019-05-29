# @Author: Massimo De Mauri <massimodemauri>
# @Date:   2019-02-06T18:10:22+01:00
# @Email:  massimo.demauri@gmail.com
# @Project: OpenBB
# @Filename: solve!.jl
# @Last modified by:   massimo
# @Last modified time: 2019-05-21T19:37:59+02:00
# @License: apache 2.0
# @Copyright: {{copyright}}


# include referred code
include("./solve_and_branch!.jl")
include("./run!.jl")


# This is the main function called to solve a branch and bound problem
function solve!(workspace::BBworkspace)::Nothing

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
