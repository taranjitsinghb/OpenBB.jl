# @Author: Wim Van Roy
# @Date:   2019-03-11T12:27:34+01:00
# @Email:  wim.van.roy@atlascopco.com
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

   if withBoundsPropagation
      try
         OpenBB.bounds_propagation!(
            node, workspace.subsolverWS.A,
            workspace.problem.varSet.dscIndices,
            updatedVars
         )
      catch e
         if isa(e, InfeasibleError)
            return false
         else
            rethrow(e)
         end
      end
   end

   return true
end
