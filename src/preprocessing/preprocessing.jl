# @Author: Massimo De Mauri <massimo>
# @Date:   2019-03-11T12:27:34+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: preprocessing.jl
# @Last modified by:   massimo
# @Last modified time: 2019-04-04T00:15:38+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


include("./linear_bounds_propagation.jl")

function preprocess!(workspace::BBworkspace)::Bool
    globalVarLoBs, globalVarUpBs = OpenBB.get_variableBounds(workspace)
    globalCnsLoBs, globalCnsUpBs = OpenBB.get_constraintBounds(workspace)

    return OpenBB.bounds_propagation!(workspace.subsolverWS.A,
                                     globalCnsLoBs,
                                     globalCnsUpBs,
                                     globalVarLoBs,
                                     globalVarUpBs,
                                     workspace.dscIndices)
end

function preprocess!(node::BBnode, updatedVars::Array{Int64,1})::Bool
    try
        OpenBB.variable_bounds_propagation!(
            updatedVars, node
        )
        return true
    catch
        return false
    end
end
