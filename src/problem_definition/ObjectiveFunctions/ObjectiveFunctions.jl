# @Author: Massimo De Mauri <massimo>
# @Date:   2019-02-26T15:03:06+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: ObjectiveFunctions.jl
# @Last modified by:   massimo
# @Last modified time: 2019-05-23T19:01:19+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# abstract and null types
abstract type AbstractObjectiveFunction end

struct NullObjectiveFunction <: AbstractObjectiveFunction end

include("./LinearObjective.jl")
include("./QuadraticObjective.jl")
