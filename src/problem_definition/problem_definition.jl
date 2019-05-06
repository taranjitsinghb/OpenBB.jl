# @Author: Massimo De Mauri <massimodemauri>
# @Date:   2019-02-06T16:19:50+01:00
# @Email:  massimo.demauri@gmail.com
# @Project: OpenBB
# @Filename: problem_definitions.jl
# @Last modified by:   massimo
# @Last modified time: 2019-03-14T16:05:49+01:00
# @License: apache 2.0
# @Copyright: {{copyright}}


# load the possible componets definitions
include("./ObjectiveFunctions.jl")
include("./ConstraintSets.jl")
include("./VariableSet.jl")


# abstract and null types
abstract type AbstractProblem end
struct NullProblem <: AbstractProblem end

# concrete problem
struct Problem{T1<:AbstractObjectiveFun,T2<:AbstractConstraintSet}
    objFun::T1
    cnsSet::T2
    varSet::VariableSet
end

function Problem(;objFun::AbstractObjectiveFun,cnsSet::AbstractConstraintSet,varSet::VariableSet)::Problem
    return Problem(objFun,cnsSet,varSet)
end
