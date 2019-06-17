# @Author: Massimo De Mauri <massimodemauri>
# @Date:   2019-02-06T16:19:50+01:00
# @Email:  massimo.demauri@gmail.com
# @Project: OpenBB
# @Filename: problem_definitions.jl
# @Last modified by:   massimo
# @Last modified time: 2019-05-23T19:00:11+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# load the possible componets definitions
include("./ObjectiveFunctions/ObjectiveFunctions.jl")
include("./ConstraintSets/ConstraintSets.jl")
include("./VariableSet.jl")


# abstract and null types
abstract type AbstractProblem end
struct NullProblem <: AbstractProblem end

# concrete problem
struct Problem{T1<:AbstractObjectiveFunction,T2<:AbstractConstraintSet}
    objFun::T1
    cnsSet::T2
    varSet::VariableSet
end

function Problem(;objFun::AbstractObjectiveFunction,cnsSet::AbstractConstraintSet,varSet::VariableSet)::Problem
    return Problem(objFun,cnsSet,varSet)
end

function get_numVariables(problem::Problem)::Int
    return get_numVariables(problem.varSet)
end

function get_numDiscreteVariables(problem::Problem)::Int
    return get_numDiscreteVariables(problem.varSet)
end

function get_numConstraints(problem::Problem)::Int
    return get_numConstraints(problem.cnsSet)
end
