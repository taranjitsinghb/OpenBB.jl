# @Author: Massimo De Mauri <massimo>
# @Date:   2019-03-12T15:41:14+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: simple_rounding_heuristics.jl
# @Last modified by:   massimo
# @Last modified time: 2019-04-30T17:07:55+02:00
# @License: apache 2.0
# @Copyright: {{copyright}}


function simple_rounding_heuristics(node::BBnode, workspace::BBworkspace)::BBnode

    # round the primal info and fix the discrete variables
    primal = copy(node.primal)
    newBranchLoBs = copy(node.branchLoBs)
    newBranchUpBs = copy(node.branchUpBs)

    for ind in workspace.dscIndices
        newBranchLoBs[ind] = newBranchUpBs[ind] = primal[ind] = Int(round(node.primal[ind]))

    end

    # return the resulting node
    BBnode(newBranchLoBs,newBranchUpBs,copy(node.pseudoCosts),
                 primal,copy(node.bndDual),copy(node.cnsDual),
                 0,node.objVal,true)

end
