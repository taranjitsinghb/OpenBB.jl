# @Author: Massimo De Mauri <massimodemauri>
# @Date:   2019-02-06T18:10:22+01:00
# @Email:  massimo.demauri@gmail.com
# @Project: OpenBB
# @Filename: solve!.jl
# @Last modified by:   massimo
# @Last modified time: 2019-02-19T14:44:42+01:00
# @License: apache 2.0
# @Copyright: {{copyright}}


# include referred code
include("./solve_and_branch!.jl")
include("./run!.jl")


# This is the main function called to solve a branch and bound problem
function solve!(workspace::BBworkspace)
    if workspace.settings.maxProcesses == 1
        return run!(workspace)
    else
        @error "not_implemented yet"
    end
end
