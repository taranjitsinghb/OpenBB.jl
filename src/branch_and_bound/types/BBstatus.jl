# @Author: Massimo De Mauri <massimo>
# @Date:   2019-05-23T11:28:46+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: BBstatus.jl
# @Last modified by:   massimo
# @Last modified time: 2019-06-17T14:17:59+02:00
# @License: apache 2.0
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
                    waitingTime::Float64=0.0,totalTime::Float64=0.0,
                    numSolutions::Int=0,numRelaxationsSolved::Int=0,
                    description::String="new")::BBstatus
    return BBstatus(objLoB,objUpB,absoluteGap,relativeGap,waitingTime,totalTime,numSolutions,numRelaxationsSolved,description)
end
