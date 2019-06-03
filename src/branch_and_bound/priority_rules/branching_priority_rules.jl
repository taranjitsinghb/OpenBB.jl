# @Author: Massimo De Mauri <massimodemauri>
# @Date:   2019-02-06T16:19:41+01:00
# @Email:  massimo.demauri@gmail.com
# @Project: OpenBB
# @Filename: branching_priority_functions.jl
# @Last modified by:   massimo
# @Last modified time: 2019-06-03T19:12:52+02:00
# @License: apache 2.0
# @Copyright: {{copyright}}


# priority functions for variables to branch on

function pseudo_cost(fractionality::Array{Float64,1},pseudoCosts::Array{Float64,2})::Int
    (_,index) = findmax(@. (1e-12 + fractionality)*pseudoCosts)
    return index
end

function most_fractional(fractionality::Array{Float64,1},pseudoCosts::Array{Float64,2})::Int
    (_,index) = findmax(fractionality)
    if fractionality[index] == 0
        return 0
    else
        return index
    end
end

function order_of_appearance(fractionality::Array{Float64,1},pseudoCosts::Array{Float64,2})::Int
    for i in 1:length(fractionality)
        if fractionality[i] > 0
            return i
        end
    end
    return 0
end

function reverse_order_of_appearance(fractionality::Array{Float64,1},pseudoCosts::Array{Float64,2})::Int
    for i in length(fractionality):-1:1
        if fractionality[i] > 0
            return i
        end
    end
    return 0
end
