# @Author: Massimo De Mauri <massimo>
# @Date:   2019-05-20T10:04:17+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: update_settings.jl
# @Last modified by:   massimo
# @Last modified time: 2019-06-03T18:04:56+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


#TODO: implement update_settings!
# function update_settings!(workspace::BBworkspace{T1,T2};bb_settings::AbstractSettings=NullSettings(),ss_settings::AbstractSettings=NullSettings())::Nothing
#          where T1<:AbstractWorkspace where T2<:AbstractSharedMemory
#
#
#     @error "Not implemented yet"
#
#     # check if it is possible to make changes
#     @assert workspace.settings.dynamicMode == true || workspace.status.description == "new"
#
#
#     # update the BB settings
#     if !(bb_settings isa NullSettings)
#         if workspace.status.totalTime > 0.
#
#             # adapt the workspace to the new settings
#
#
#             # pass from dynamic mode to static mode
#             if workspace.dynamicMode == true && bb_settings.dynamicMode == false
#                 for k in length(workspace.solutionPool):-1:1
#                     if workspace.solutionPool[k].reliable
#                         deleteat!(workspace.solutionPool,k+1:length(workspace.solutionPool))
#                         deleteat!(workspace.solutionPool,1:k-1)
#                         break
#                     elseif k == 1
#                         deleteat!(workspace.solutionPool,1:length(workspace.solutionPool)-1)
#                     end
#                 end
#                 deleteat!(workspace.unactivePool,1:length(workspace.unactivePool))
#             end
#         end
#
#         for field in fieldnames(BBsettings)
#             setfield!(workspace.settings,field,getfield(bb_settings,field))
#         end
#     end
#
#     # update the subsolver settings
#     if !(ss_settings isa NullSettings)
#         update_settings!(workspace.subsolverWS,ss_settings,bb_primalTolerance=workspace.settings.primalTolerance,bb_timeLimit=workspace.settings.timeLimit)
#     end
#
#     return
# end
