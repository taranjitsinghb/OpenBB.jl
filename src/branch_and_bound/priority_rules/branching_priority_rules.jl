# @Author: Massimo De Mauri <massimodemauri>
# @Date:   2019-02-06T16:19:41+01:00
# @Email:  massimo.demauri@gmail.com
# @Project: OpenBB
# @Filename: branching_priority_functions.jl
# @Last modified by:   massimo
# @Last modified time: 2019-06-17T12:51:00+02:00
# @License: apache 2.0
# @Copyright: {{copyright}}


# wrapper for branching priority rules
function branching_priority_rule(functionTuple::Tuple,dscValues::Array{Float64,1},pseudoCosts::Tuple{Array{Float64,2},Array{Int,2}},primalTolerance::Float64)::Tuple{Int,Array{Float64,1}}
    @assert functionTuple[1] isa Function
    return functionTuple[1](dscValues,pseudoCosts,primalTolerance,functionTuple[2:end]...)
end

# actual priority rules
function pseudoIncrements_mean(dscValues::Array{Float64,1},pseudoCosts::Tuple{Array{Float64,2},Array{Int,2}},primalTolerance::Float64)::Tuple{Int,Array{Float64,1}}
    return pseudoIncrements_mean(dscValues,pseudoCosts,primalTolerance,1,1)
end

function pseudoIncrements_mean(dscValues::Array{Float64,1},pseudoCosts::Tuple{Array{Float64,2},Array{Int,2}},primalTolerance::Float64,
                          ratio::Rational)::Tuple{Int,Array{Float64,1}}
    return pseudoIncrements_mean(dscValues,pseudoCosts,primalTolerance,numerator(ratio),denominator(ratio)-numerator(ratio))
end

function pseudoIncrements_mean(dscValues::Array{Float64,1},pseudoCosts::Tuple{Array{Float64,2},Array{Int,2}},primalTolerance::Float64,
                          gainOfMin::Int64,gainOfMax::Int64)::Tuple{Int,Array{Float64,1}}

    scores = zeros(length(dscValues))
    bestIndex = 1
    for k in 1: length(dscValues)
        deltaMinus = dscValues[k] - floor(dscValues[k]+primalTolerance)
        deltaMinus = pseudoCosts[1][k,1]*deltaMinus*(deltaMinus>primalTolerance)
        deltaPlus  = ceil(dscValues[k]-primalTolerance) - dscValues[k]
        deltaPlus  = pseudoCosts[1][k,2]*deltaPlus*(deltaPlus>primalTolerance)

        if deltaMinus <= deltaPlus
            scores[k] = deltaMinus*gainOfMin + deltaPlus*gainOfMax
        else
            scores[k] = deltaPlus*gainOfMin + deltaMinus*gainOfMax
        end

        if scores[k] > scores[bestIndex]
            bestIndex = k
        end
    end
    return bestIndex, scores
end

function pseudoIncrements_geomean(dscValues::Array{Float64,1},pseudoCosts::Tuple{Array{Float64,2},Array{Int,2}},primalTolerance::Float64)::Tuple{Int,Array{Float64,1}}
    return pseudoIncrements_geomean(dscValues,pseudoCosts,primalTolerance,1,1)
end

function pseudoIncrements_geomean(dscValues::Array{Float64,1},pseudoCosts::Tuple{Array{Float64,2},Array{Int,2}},primalTolerance::Float64,
                                      ratio::Rational)::Tuple{Int,Array{Float64,1}}
    return pseudoIncrements_geomean(dscValues,pseudoCosts,primalTolerance,numerator(ratio),denominator(ratio)-numerator(ratio))
end

function pseudoIncrements_geomean(dscValues::Array{Float64,1},pseudoCosts::Tuple{Array{Float64,2},Array{Int,2}},primalTolerance::Float64,
                             gainOfMin::Int,gainOfMax::Int)::Tuple{Int,Array{Float64,1}}

    scores = zeros(length(dscValues))
    bestIndex = 1
    for k in 1:length(dscValues)
        deltaMinus = dscValues[k] - floor(dscValues[k])
        deltaMinus = pseudoCosts[1][k,1]*deltaMinus*(deltaMinus>primalTolerance)
        deltaPlus  = ceil(dscValues[k]) - dscValues[k]
        deltaPlus  = pseudoCosts[1][k,2]*deltaPlus*(deltaPlus>primalTolerance)

        if deltaMinus <= deltaPlus
            scores[k] = deltaMinus^gainOfMin*deltaPlus^gainOfMax
        else
            scores[k] = deltaPlus^gainOfMin*deltaMinus^gainOfMax
        end

        if scores[k] > scores[bestIndex]
            bestIndex = k
        end
    end
    return bestIndex, scores
end



function absolute_fractionality(dscValues::Array{Float64,1},pseudoCosts::Tuple{Array{Float64,2},Array{Int,2}},primalTolerance::Float64)::Tuple{Int,Array{Float64,1}}

    scores = zeros(length(dscValues))
    bestIndex = 1
    for k in 1:length(dscValues)
        scores[k] = abs(dscValues[k] - round(dscValues[k]))
        scores[k] = scores[k]*(scores[k]>primalTolerance)

        if scores[k] > scores[bestIndex]
            bestIndex = k
        end
    end

    return bestIndex, scores
end



function order_of_appearance(dscValues::Array{Float64,1},pseudoCosts::Tuple{Array{Float64,2},Array{Int,2}},primalTolerance::Float64;reverse::Bool=false)::Tuple{Int,Array{Float64,1}}

    scores = zeros(length(dscValues))
    if reverse
        for k in length(dscValues):-1:1
            if abs(dscValues[k] - round(dscValues[k])) > primalTolerance
                scores[k] = Inf
                return k, scores
            end
        end
    else
        for k in 1:length(dscValues)
            if abs(dscValues[k] - round(dscValues[k])) > primalTolerance
                scores[k] = Inf
                return k, scores
            end
        end
    end

    return 0, scores
end
