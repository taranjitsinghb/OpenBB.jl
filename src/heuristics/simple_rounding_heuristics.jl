# @Author: Massimo De Mauri <massimo>
# @Date:   2019-03-12T15:41:14+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: simple_rounding_heuristics.jl
# @Last modified by:   massimo
# @Last modified time: 2019-04-30T17:07:55+02:00
# @License: apache 2.0
# @Copyright: {{copyright}}


function simple_rounding_heuristics(subproblem::BBsubproblem, workspace::BBworkspace)::BBsubproblem

    # round the primal info and fix the discrete variables
    primal = copy(subproblem.primal)
    newBranchLoBs = copy(subproblem.branchLoBs)
    newBranchUpBs = copy(subproblem.branchUpBs)

    for ind in workspace.dscIndices
        newBranchLoBs[ind] = newBranchUpBs[ind] = primal[ind] = Int(round(subproblem.primal[ind]))

    end

    # return the resulting subproblem
    BBsubproblem(newBranchLoBs,newBranchUpBs,copy(subproblem.pseudoCosts),
                 primal,copy(subproblem.bndDual),copy(subproblem.cnsDual),
                 0,subproblem.objVal,true)

end
