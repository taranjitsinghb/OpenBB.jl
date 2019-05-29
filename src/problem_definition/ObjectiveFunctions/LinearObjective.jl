# @Author: Massimo De Mauri <massimo>
# @Date:   2019-05-23T18:33:51+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: LinearObjective.jl
# @Last modified by:   massimo
# @Last modified time: 2019-05-23T18:54:39+02:00
# @License: apache 2.0
# @Copyright: {{copyright}}


# linear objective
mutable struct LinearObjective <: AbstractObjectiveFunction
    L::Array{Float64,1}
end

function LinearObjective(;L::Array{Float64,1})::LinearObjective
    return LinearObjective(L)
end
