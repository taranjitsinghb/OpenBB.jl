# @Author: Massimo De Mauri <massimo>
# @Date:   2019-05-23T11:28:46+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: BBstatus.jl
# @Last modified by:   massimo
# @Last modified time: 2019-07-03T18:01:30+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


mutable struct BBstatus
    objLoB::Float64
    objUpB::Float64
    absoluteGap::Float64
    relativeGap::Float64
    totalTime::Float64
    waitingTime::Float64
    numSolutions::Int
    numRelaxationsSolved::Int64
    description::String
end

function BBstatus(; objLoB::Float64=-Inf,objUpB::Float64=Inf,
                    absoluteGap::Float64= Inf,relativeGap::Float64=Inf,
                    totalTime::Float64=0.0,waitingTime::Float64=0.0,
                    numSolutions::Int=0,numRelaxationsSolved::Int=0,
                    description::String="new")::BBstatus
    return BBstatus(objLoB,objUpB,absoluteGap,relativeGap,totalTime,waitingTime,numSolutions,numRelaxationsSolved,description)
end


function BBstatus(status::BBstatus)::BBstatus
    return BBstatus(status.objLoB,status.objUpB,
                    status.absoluteGap,status.relativeGap,
                    status.totalTime,status.waitingTime,
                    status.numSolutions,status.numRelaxationsSolved,
                    status.description)
end
