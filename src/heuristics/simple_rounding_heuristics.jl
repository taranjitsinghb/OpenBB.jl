# @Author: Massimo De Mauri <massimo>
# @Date:   2019-03-12T15:41:14+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: simple_rounding_heuristics.jl
# @Last modified by:   massimo
# @Last modified time: 2019-09-02T15:12:06+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


function simple_rounding_heuristics(node::BBnode, workspace::BBworkspace)::BBnode

    # round the primal info and fix the discrete variables
    newPrimal = copy(node.primal)
    newLoBs = -Infs(length(node.primal))
    newUpBs =  Infs(length(node.primal))
    @. newPrimal[workspace.problem.varSet.dscIndices]  =
       newLoBs[workspace.problem.varSet.dscIndices]    =
       newUpBs[workspace.problem.varSet.dscIndices]    = round(newPrimal[workspace.problem.varSet.dscIndices])

    # return the resulting node
    BBnode(newPrimal[workspace.problem.varSet.dscIndices],newPrimal[workspace.problem.varSet.dscIndices],
           newPrimal,copy(node.bndDual),copy(node.cnsDual),
           0.0,node.objVal,node.pseudoObjective,false)

end
