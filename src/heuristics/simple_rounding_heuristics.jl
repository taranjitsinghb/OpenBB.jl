# @Author: Massimo De Mauri <massimo>
# @Date:   2019-03-12T15:41:14+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: simple_rounding_heuristics.jl
# @Last modified by:   massimo
# @Last modified time: 2019-07-10T19:02:53+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


function simple_rounding_heuristics(node::BBnode, workspace::BBworkspace)::BBnode

    # round the primal info and fix the discrete variables
    newPrimal = copy(node.primal)
    @. newPrimal[workspace.dscIndices] = round(newPrimal[workspace.dscIndices])

    # return the resulting node
    BBnode(newPrimal[workspace.dscIndices],newPrimal[workspace.dscIndices],
           newPrimal,copy(node.bndDual),copy(node.cnsDual),
           0.0,node.objective,node.pseudoObjective,false)

end
