# @Author: Wim Van Roy
# @Date:   2019-03-11T12:27:34+01:00
# @Email:  wim.vanroy@kuleuven.be
# @Filename: preprocessing.jl
# @Last modified by:   massimo
# @Last modified time: 2019-09-02T13:44:21+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

struct InfeasibleError <: Exception
   msg::String
end

include("./linear_bounds_propagation.jl")


function preprocess!(node::BBnode, workspace::BBworkspace, updatedVars::Array{Int64,1};
                     withBoundsPropagation::Bool=true)::Bool

   feasible = true
   if withBoundsPropagation
       feasible, updatedVars = OpenBB.bounds_propagation!(
            node, workspace.subsolverWS.A,
            workspace.dscIndices,
            updatedVars
         )
   end

   return feasible
end
