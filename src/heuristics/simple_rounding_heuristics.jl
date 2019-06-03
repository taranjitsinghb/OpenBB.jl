# @Author: Massimo De Mauri <massimo>
# @Date:   2019-03-12T15:41:14+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: simple_rounding_heuristics.jl
# @Last modified by:   massimo
# @Last modified time: 2019-06-03T14:36:26+02:00
# @License: apache 2.0
# @Copyright: {{copyright}}


function simple_rounding_heuristics(node::BBnode, workspace::BBworkspace)::BBnode

    # round the primal info and fix the discrete variables
    primal = copy(node.primal)
    newBranchLoBs = copy(node.branchLoBs)
    newBranchUpBs = copy(node.branchUpBs)

    @. primal[workspace.dscIndices] = round(node.primal[workspace.dscIndices])
    @. newBranchLoBs = newBranchUpBs = Int(primal[workspace.dscIndices])

    # return the resulting node
    BBnode(newBranchLoBs,newBranchUpBs,
           primal,copy(node.bndDual),copy(node.cnsDual),
           0,node.objVal,true)

end
