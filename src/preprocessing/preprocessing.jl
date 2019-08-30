# @Author: Massimo De Mauri <massimo>
# @Date:   2019-03-11T12:27:34+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: preprocessing.jl
# @Last modified by:   massimo
# @Last modified time: 2019-04-04T00:15:38+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

struct InfeasibleError <: Exception
   msg::String
end

include("./linear_bounds_propagation.jl")


function preprocess!(node::BBnode, workspace::BBworkspace, updatedVars::Array{Int64,1})::Bool
   try
      OpenBB.bounds_propagation!(
         node, workspace.subsolverWS.A,
         workspace.dscIndices,
         updatedVars
      )
      return true
   catch e
      if isa(e, InfeasibleError)
         return false
      else
         rethrow(e)
      end
   end
end
