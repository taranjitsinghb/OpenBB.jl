# @Author: Massimo De Mauri <massimo>
# @Date:   2019-02-26T15:03:06+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: ConstraintSets.jl
# @Last modified by:   massimo
# @Last modified time: 2019-08-27T12:11:20+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# abstract and null types
abstract type AbstractConstraintSet end
struct NullConstraintSet <: AbstractConstraintSet end


# linear constraints set
mutable struct LinearConstraintSet{T<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}} <: AbstractConstraintSet
    A::T
    loBs::Array{Float64,1}
    upBs::Array{Float64,1}
end

include("./LinearConstraintSet.jl")
