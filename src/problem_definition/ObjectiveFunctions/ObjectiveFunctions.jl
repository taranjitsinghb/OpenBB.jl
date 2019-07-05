# @Author: Massimo De Mauri <massimo>
# @Date:   2019-02-26T15:03:06+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: ObjectiveFunctions.jl
# @Last modified by:   massimo
# @Last modified time: 2019-07-02T15:21:36+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# abstract and null types
abstract type AbstractObjective end

struct NullObjective <: AbstractObjective end

# ...
function get_numVariables(objectiveFun::NullObjective)::Int
    return 0
end

include("./LinearObjective.jl")
include("./QuadraticObjective.jl")
